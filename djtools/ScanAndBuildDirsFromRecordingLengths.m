function [results, scaninfo] = ScanAndBuildDirsFromRecordingLengths(varargin)
%ScanAndBuildDirsFromRecordingLengths - scan a rig root directory (e.g.
%/Volumes/Projects/5XFAD/Rig3Phys/) for every recording_lengths.json
%under it -- whether it lives directly in a Bonsai (outer) folder or one
%level down in that Bonsai folder's Ephys (inner) folder, since the
%split function has been observed writing it to either location --
%sanity-check the group of sessions each one describes, and (optionally)
%call BuildDirsFromRecordingLengths on it to write dirs.mat/Bdirs.mat.
%
%Because the grouped folders behind a recording_lengths.json are picked
%by hand in a GUI during kilosorting, and that's subject to human error,
%each group is sanity-checked before anything is written:
%   - every member folder name must parse as "YYYY-MM-DD_H-MM-SS_mouse-ID"
%   - all members must have the SAME mouse ID
%   - all members must fall on the SAME calendar day (recordings are
%     assumed to never span midnight)
%   - the members must be CHRONOLOGICALLY CONTIGUOUS among all the
%     Bonsai folders under root -- i.e. no other, non-grouped session's
%     folder may fall between the earliest and latest member timestamps
%A group that fails a sanity check is reported (with the reason) and
%skipped -- nothing is ever written for it, dry run or not. A group
%whose members span more than 'maxhours' (default 6) is not failed, but
%prints/records a note, since this is meant to be typical rather than
%required.
%
%All console output is also echoed (via diary) to a log file,
%ScanAndBuildDirsFromRecordingLengths.log, in root. Unlike the
%delete-and-restart-each-run logs used elsewhere in djtools (e.g.
%ProcessTC_LFP2.m), this log is APPENDED to on every run -- since this
%function is meant to be re-run periodically as new data/groups show up
%-- with a timestamped header on each run recording whether it ran in
%dry-run or write mode.
%
%Because it's meant to be re-run periodically, a passing group whose
%dirs.mat (in every member's Ephys directory) and Bdirs.mat (in every
%member's Bonsai directory) already exist AND already contain exactly
%the dirs/Bdirs/DataRoot that recording_lengths.json currently describes
%is left alone -- reported as "already up to date" rather than
%rewritten. This applies in both dry-run and write mode (dry-run just
%reports it instead of writing regardless); only a group that's new,
%changed, or only partially written gets an actual write in write mode.
%
%usage:
%   ScanAndBuildDirsFromRecordingLengths()
%       dry run (write=false) over the default root
%   ScanAndBuildDirsFromRecordingLengths(root)
%   ScanAndBuildDirsFromRecordingLengths(root, 'write', true)
%       actually writes dirs.mat/Bdirs.mat for every group that passes
%       the sanity checks
%   ScanAndBuildDirsFromRecordingLengths('write', true)
%       same, using the default root
%   ScanAndBuildDirsFromRecordingLengths(..., 'maxhours', 8)
%   ScanAndBuildDirsFromRecordingLengths(..., 'mouse', 3622)
%   ScanAndBuildDirsFromRecordingLengths(..., 'mouse', {'3622','3555'})
%       only look for recording_lengths.json in these mice's session folders
%       (number, string, or a list); default [] = all mice. The contiguity
%       check still uses every folder under root, so it is unchanged.
%       change the (soft, non-failing) time-span warning threshold
%
%returns:
%   results  - one entry per recording_lengths.json found (including
%       parse failures and duplicates), with every fact that gets
%       printed to the console/log also captured as a field, so results
%       is a structured, isomorphic record of the run:
%         hostfolder   - the folder recording_lengths.json was found in
%         jsonpath     - full path to that recording_lengths.json
%         members      - cell array of grouped Bonsai folder names
%         ephysdirs    - cell array of the corresponding Ephys folder
%                        names (same order as members)
%         dataroot     - DataRoot as returned by BuildDirsFromRecordingLengths
%         mouseid      - shared mouse ID, if the group passed that check
%         date         - shared 'yyyy-mm-dd' date, if the group passed
%         spanhours    - hours between earliest and latest member (NaN
%                        if the group failed before this could be computed)
%         note         - the "spans more than maxhours" note text, or ''
%         passed       - true/false
%         reason       - failure reason, or '' if passed
%         written      - true if, as of the end of this call (with
%                        'write',true), dirs.mat/Bdirs.mat are correctly
%                        in place for this group -- whether because this
%                        run actually wrote them, they were already up
%                        to date, or (for a duplicate) an earlier entry
%                        in this run already covered it. Always false on
%                        a dry run.
%         alreadyuptodate - true if dirs.mat/Bdirs.mat already existed
%                        and already matched recording_lengths.json
%                        BEFORE this call -- i.e. nothing needed to be
%                        written for this group, in either dry-run or
%                        write mode
%         duplicate    - true if this recording_lengths.json describes
%                        the same group as an earlier one in this run
%         duplicateof  - hostfolder of the first occurrence, if duplicate
%   scaninfo - scalar struct with the scan-level facts that aren't tied
%       to any one group: root, write, maxhours, mouse, logfile, runtime,
%       njsonfiles, nunparsedfolders, unparsedfolders
%
%-mike/claude 9.2.26 ('mouse' option 10.13.26)

DEFAULT_ROOT='/Volumes/Projects/5XFAD/Rig3Phys/';

%if called with only name-value pairs and no root (e.g.
%ScanAndBuildDirsFromRecordingLengths('write', true)), inputParser would
%otherwise try to consume 'write'/'maxhours' positionally as the
%(optional) root value. Detect that case and insert the default root
%explicitly.
knownparams={'write', 'maxhours', 'mouse'};
if ~isempty(varargin) && (ischar(varargin{1})||isstring(varargin{1})) && any(strcmpi(varargin{1}, knownparams))
    varargin=[{DEFAULT_ROOT}, varargin];
end

p=inputParser;
addOptional(p, 'root', DEFAULT_ROOT, @(x) ischar(x)||isstring(x));
addParameter(p, 'write', false, @(x) islogical(x)||isnumeric(x));
addParameter(p, 'maxhours', 6, @isnumeric);
addParameter(p, 'mouse', [], @(x) isempty(x)||isnumeric(x)||ischar(x)||isstring(x)||iscell(x));
parse(p, varargin{:});
root=char(p.Results.root);
dowrite=logical(p.Results.write);
maxhours=p.Results.maxhours;
mice=p.Results.mouse; if ~isempty(mice), mice=strtrim(string(mice(:)')); end   %[] = all mice

if ~isfolder(root)
    error('ScanAndBuildDirsFromRecordingLengths:noroot', 'root directory not found: %s', root);
end

%echo everything printed below to a log file in root, appending across
%runs (diary appends by default when the file already exists -- no
%delete() here, unlike the fresh-each-run logs elsewhere in djtools).
logfilename='ScanAndBuildDirsFromRecordingLengths.log';
logfullpath=fullfile(root, logfilename);
diary(logfullpath);
logCleanupObj=onCleanup(@() diary('off')); %guarantees diary is turned off even if this function errors out

runtime=datestr(now);
if dowrite
    modestr='WRITE (dirs.mat/Bdirs.mat will be created for passing groups)';
else
    modestr='DRY RUN (write=false, nothing will be written)';
end
fprintf('\n\n===== %s run started: %s =====\n', mfilename, runtime);
fprintf('root: %s\n', root);
fprintf('mode: %s\n', modestr);
fprintf('maxhours (time-span warning threshold): %g\n', maxhours);
if isempty(mice), fprintf('mice: all\n'); else, fprintf('mice: %s\n', strjoin(mice, ', ')); end

%Build the master list of every immediate subdirectory of root, parsed
%into a timestamp (as a datenum) + mouseID where possible. This is what
%the chronological-contiguity check below is run against, so it has to
%cover ALL sessions at this rig, not just the ones already grouped into
%a recording_lengths.json.
d=dir(root);
d=d([d.isdir]);
d=d(~ismember({d.name}, {'.','..'}));
nall=numel(d);

allnames=cell(1,nall);
alldnum=nan(1,nall);
unparsedfolders={};
for i=1:nall
    [dnum, ~, ~, ok]=parse_bonsai_folder_name(d(i).name);
    allnames{i}=d(i).name;
    if ok
        alldnum(i)=dnum;
    else
        unparsedfolders{end+1}=d(i).name; %#ok<AGROW>
    end
end
nunparsed=numel(unparsedfolders);
if nunparsed>0
    fprintf('\n(note: %d of %d folders under %s did not match the expected "YYYY-MM-DD_H-MM-SS_mouse-ID" naming and were excluded from the contiguity check)\n', nunparsed, nall, root);
end
parsedmask=~isnan(alldnum);
parsednames=allnames(parsedmask);
parseddnum=alldnum(parsedmask);

%find every recording_lengths.json under root: either directly in a
%Bonsai (outer) folder, or one level down inside one of its immediate
%subfolders (typically the Ephys (inner) folder -- the split function
%has been observed writing recording_lengths.json to either location).
%hostfolders entries are paths relative to root, so a Bonsai-level hit
%is just the Bonsai folder name, and an Ephys-level hit is
%"<Bonsai folder>/<Ephys folder>".
%with 'mouse', only those mice's session folders are searched (the master list above stays complete)
scanidx=1:nall;
if ~isempty(mice)
    folderids=regexp(allnames, '_mouse-(.+)$', 'tokens', 'once'); folderids(cellfun(@isempty, folderids))={{''}};
    folderids=string(cellfun(@(c) c{1}, folderids, 'UniformOutput', false));
    scanidx=find(ismember(folderids, mice));
    nf=arrayfun(@(m) sum(folderids==m), mice);
    fprintf('%d session folders for %d mice', numel(scanidx), numel(mice));
    if any(nf==0), fprintf(' (no folders for mouse %s)', strjoin(mice(nf==0), ', ')); end
    fprintf('\n');
end
hostfolders={};
foundlocation={}; %parallel to hostfolders, for the summary counts below
%in-place progress counter; the diary is paused so the \b updates don't go into the log
diary('off');
nscan=numel(scanidx);
fprintf('scanning %d folders for recording_lengths.json:', nscan); nbytes=0;
for k=1:nscan
    i=scanidx(k);
    fprintf(repmat('\b', 1, nbytes)); nbytes=fprintf(' %d/%d', k, nscan);
    bonsaidir=fullfile(root, d(i).name);

    if isfile(fullfile(bonsaidir, 'recording_lengths.json'))
        hostfolders{end+1}=d(i).name; %#ok<AGROW>
        foundlocation{end+1}='Bonsai folder'; %#ok<AGROW>
    end

    sd=dir(bonsaidir);
    sd=sd([sd.isdir]);
    sd=sd(~ismember({sd.name}, {'.','..'}));
    for j=1:numel(sd)
        if isfile(fullfile(bonsaidir, sd(j).name, 'recording_lengths.json'))
            hostfolders{end+1}=fullfile(d(i).name, sd(j).name); %#ok<AGROW>
            foundlocation{end+1}='Ephys folder'; %#ok<AGROW>
        end
    end
end
fprintf('\n'); diary(logfullpath);   %resume the log
nfoundbonsai=sum(strcmp(foundlocation, 'Bonsai folder'));
nfoundephys=sum(strcmp(foundlocation, 'Ephys folder'));
fprintf('\nfound %d recording_lengths.json file(s) under %s (%d in a Bonsai folder, %d in an Ephys subfolder)\n', ...
    numel(hostfolders), root, nfoundbonsai, nfoundephys);

results=struct('hostfolder',{},'jsonpath',{},'members',{},'ephysdirs',{},'dataroot',{}, ...
    'mouseid',{},'date',{},'spanhours',{},'note',{},'passed',{},'reason',{},'written',{}, ...
    'alreadyuptodate',{},'duplicate',{},'duplicateof',{});
sigmap=containers.Map(); %signature (sorted members joined) -> index into results, for de-duplication

for hi=1:numel(hostfolders)
    hostfolder=hostfolders{hi};
    jsonpath=fullfile(root, hostfolder, 'recording_lengths.json');
    fprintf('\n--- %s ---', fullfile(hostfolder, 'recording_lengths.json'));

    try
        [dirs_i, Bdirs_i, DataRoot_i] = BuildDirsFromRecordingLengths(jsonpath); %read-only here -- never writes
    catch err
        fprintf('\n  FAILED to parse: %s\n', err.message);
        results(end+1)=make_entry(hostfolder, jsonpath, {}, {}, '', '', '', NaN, '', false, err.message, false, false, false, ''); %#ok<AGROW>
        continue
    end

    members=cellfun(@local_basename, Bdirs_i, 'UniformOutput', false);
    signature=strjoin(sort(members), '|');

    if isKey(sigmap, signature)
        orig=results(sigmap(signature));
        fprintf('\n  (same group already handled via %s -- skipping duplicate)\n', orig.hostfolder);
        entry=orig;
        entry.hostfolder=hostfolder;
        entry.jsonpath=jsonpath;
        entry.duplicate=true;
        entry.duplicateof=orig.hostfolder;
        results(end+1)=entry; %#ok<AGROW>
        continue
    end

    fprintf('\n  members (n=%d):', numel(members));
    for i=1:numel(members)
        fprintf('\n    %s', members{i});
    end

    [passed, reason, mouseid, thedate, spanhours, note] = sanity_check_group(members, parsednames, parseddnum, maxhours);

    if passed
        fprintf('\n  sanity check: OK (mouseID %s, %s)', mouseid, thedate);
    else
        fprintf('\n  sanity check: FAILED -- %s', reason);
    end
    if ~isempty(note)
        fprintf('\n  NOTE: %s', note);
    end

    written=false;
    alreadyuptodate=false;
    if passed
        isuptodate=already_uptodate(dirs_i, Bdirs_i, DataRoot_i);
        if dowrite
            if isuptodate
                fprintf('\n  already up to date -- dirs.mat/Bdirs.mat already match recording_lengths.json, nothing written');
                alreadyuptodate=true;
            else
                BuildDirsFromRecordingLengths(jsonpath, 'write', true);
                fprintf('\n  wrote dirs.mat/Bdirs.mat for %d directories', numel(members));
            end
            written=true;
        else
            if isuptodate
                fprintf('\n  [DRY RUN] already up to date -- dirs.mat/Bdirs.mat already match recording_lengths.json');
                alreadyuptodate=true;
            else
                fprintf('\n  [DRY RUN] would write dirs.mat/Bdirs.mat for %d directories', numel(members));
            end
        end
    else
        fprintf('\n  SKIPPING (sanity check failed, nothing written)');
    end

    results(end+1)=make_entry(hostfolder, jsonpath, members, dirs_i, DataRoot_i, mouseid, thedate, spanhours, note, passed, reason, written, alreadyuptodate, false, ''); %#ok<AGROW>
    sigmap(signature)=numel(results);
end

if isempty(results)
    npassed=0; nfailed=0; nwritten=0; nuptodate=0;
else
    npassed=sum([results.passed]);
    nfailed=numel(results)-npassed;
    nwritten=sum([results.written]);
    nuptodate=sum([results.alreadyuptodate]);
end
fprintf('\n\n=== summary: %d group(s) found, %d passed, %d failed', numel(results), npassed, nfailed);
if dowrite
    fprintf(', %d written (%d newly written, %d already up to date)', nwritten, nwritten-nuptodate, nuptodate);
else
    fprintf(' (dry run -- nothing written; %d already up to date; call with ''write'', true to actually write)', nuptodate);
end
fprintf(' ===\n');

scaninfo=struct('root',root, 'write',dowrite, 'maxhours',maxhours, 'mouse',{cellstr(string(mice))}, 'logfile',logfullpath, ...
    'runtime',runtime, 'njsonfiles',numel(hostfolders), 'nunparsedfolders',nunparsed, ...
    'unparsedfolders',{unparsedfolders});

%summarize failures to console
fprintf('\nlist of failures (passed==0) with results idx, mouseid, reason')
for i=1:length(results); 
    if results(i).passed==0 
        fprintf('\n%d n=%d mouse %s %s', i,length(results(i).members),results(i).mouseid, results(i).reason) 
    end 
end

end

function entry=make_entry(hostfolder, jsonpath, members, ephysdirs, dataroot, mouseid, thedate, spanhours, note, passed, reason, written, alreadyuptodate, duplicate, duplicateof)
entry=struct('hostfolder',hostfolder, 'jsonpath',jsonpath, 'members',{members}, 'ephysdirs',{ephysdirs}, ...
    'dataroot',dataroot, 'mouseid',mouseid, 'date',thedate, 'spanhours',spanhours, 'note',note, ...
    'passed',passed, 'reason',reason, 'written',written, 'alreadyuptodate',alreadyuptodate, ...
    'duplicate',duplicate, 'duplicateof',duplicateof);
end

function [dnum, datestr_only, mouseid, ok] = parse_bonsai_folder_name(name)
%parse a Bonsai-style folder name "YYYY-MM-DD_H-MM-SS_mouse-ID" into a
%datenum (dnum), a plain 'yyyy-mm-dd' date string, and the mouse ID.
%Uses datenum/regexp only (no datetime/duration classes) for maximum
%compatibility across MATLAB versions.
pat='^(?<y>\d{4})-(?<mo>\d{2})-(?<d>\d{2})_(?<h>\d{1,2})-(?<mi>\d{1,2})-(?<s>\d{1,2})_mouse-(?<mouseid>.+)$';
tok=regexp(name, pat, 'names', 'once');
if isempty(tok)
    dnum=NaN; datestr_only=''; mouseid=''; ok=false;
    return
end
y=str2double(tok.y); mo=str2double(tok.mo); dy=str2double(tok.d);
h=str2double(tok.h); mi=str2double(tok.mi); s=str2double(tok.s);
dnum=datenum(y,mo,dy,h,mi,s); %#ok<DATNM>
datestr_only=sprintf('%04d-%02d-%02d', y, mo, dy);
mouseid=tok.mouseid;
ok=true;
end

function [passed, reason, mouseid, thedate, spanhours, note] = sanity_check_group(members, allnames, alldnum, maxhours)
passed=true;
reason='';
mouseid='';
thedate='';
spanhours=NaN;
note='';

n=numel(members);
dnums=nan(1,n);
mouseids=cell(1,n);
dates=cell(1,n);
for i=1:n
    [dnum,dstr,mid,ok]=parse_bonsai_folder_name(members{i});
    if ~ok
        passed=false;
        reason=sprintf('could not parse folder name "%s" (expected "YYYY-MM-DD_H-MM-SS_mouse-ID")', members{i});
        return
    end
    dnums(i)=dnum;
    mouseids{i}=mid;
    dates{i}=dstr;
end

%same mouse ID
umouse=unique(mouseids);
if numel(umouse)>1
    passed=false;
    reason=sprintf('grouped folders have different mouse IDs: %s', strjoin(umouse, ', '));
    return
end
mouseid=umouse{1};

%same calendar day (recordings are assumed to never span midnight)
udate=unique(dates);
if numel(udate)>1
    passed=false;
    reason=sprintf('grouped folders span more than one calendar day: %s', strjoin(udate, ', '));
    return
end
thedate=udate{1};

%time span sanity -- soft warning only, not a failure
tmin=min(dnums);
tmax=max(dnums);
spanhours=(tmax-tmin)*24;
if spanhours>maxhours
    [~,imin]=min(dnums); [~,imax]=max(dnums);
    note=sprintf('this group spans %.1f hours (%s .. %s) -- longer than the typical few hours, double check this grouping', ...
        spanhours, members{imin}, members{imax});
end

%chronological contiguity: no non-member folder's timestamp may fall
%between this group's earliest and latest member (inclusive)
inwindow = alldnum>=tmin & alldnum<=tmax;
windownames = allnames(inwindow);
intervening = setdiff(windownames, members);
if ~isempty(intervening)
    passed=false;
    reason=sprintf('%d non-grouped folder(s) found chronologically between the earliest and latest grouped session: %s', numel(intervening), strjoin(intervening, ', '));
    return
end

end

function b=local_basename(p)
%last path component of p, tolerant of either separator -- avoids
%fileparts' extension-splitting behavior, which folder names here don't
%need or want.
p=regexprep(p, '[/\\]+$', ''); %strip trailing separator, if any
parts=regexp(p, '[/\\]', 'split');
b=parts{end};
end

function tf=already_uptodate(dirs_i, Bdirs_i, DataRoot_i)
%true only if, for EVERY grouped session, both dirs.mat (in its Ephys
%directory) and Bdirs.mat (in its Bonsai directory) already exist AND
%already contain exactly this dirs/Bdirs/DataRoot triple -- i.e. calling
%BuildDirsFromRecordingLengths(jsonpath,'write',true) right now would be
%a pure no-op. A missing directory, a missing .mat file, or any mismatch
%(including a different session order) makes this false, which is the
%safe direction -- it only ever causes an extra (harmless, identical)
%rewrite, never a skipped one.
tf=true;
n=numel(dirs_i);
for i=1:n
    ephysdir=fullfile(Bdirs_i{i}, dirs_i{i});
    bonsaidir=Bdirs_i{i};

    if ~mat_matches(fullfile(ephysdir, 'dirs.mat'), dirs_i, Bdirs_i, DataRoot_i)
        tf=false;
        return
    end
    if ~mat_matches(fullfile(bonsaidir, 'Bdirs.mat'), dirs_i, Bdirs_i, DataRoot_i)
        tf=false;
        return
    end
end
end

function tf=mat_matches(matpath, dirs_i, Bdirs_i, DataRoot_i)
%true if matpath exists, loads cleanly, and its dirs/Bdirs/DataRoot are
%identical (isequal, ignoring case -- macOS paths are case-insensitive, e.g.
%5xfad/rig3phys vs 5XFAD/Rig3Phys) to the ones just parsed from recording_lengths.json.
tf=false;
if ~isfile(matpath)
    return
end
try
    s=load(matpath, 'dirs', 'Bdirs', 'DataRoot');
catch
    return %unreadable/corrupt .mat -- treat as not up to date
end
if ~isfield(s,'dirs') || ~isfield(s,'Bdirs') || ~isfield(s,'DataRoot')
    return
end
lc=@(x) lower(string(x));
tf = isequal(lc(s.dirs), lc(dirs_i)) && isequal(lc(s.Bdirs), lc(Bdirs_i)) && isequal(lc(s.DataRoot), lc(DataRoot_i));
end
