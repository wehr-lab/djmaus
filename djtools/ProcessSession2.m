function PS = ProcessSession2(BonsaiPath, varargin)
%ProcessSession2 - pre-process one session recorded in the new Open Ephys format
%(function version of the ProcessSession script). Run on a Bonsai folder.
%
%Steps (each skipped if its output already exists, unless forced):
%   1. locate EphysPath / EphysPath_KS (from OEinfo-*.mat, else by searching)
%   2. Open Ephys session object: load session-<EphysPath>.mat, or build it with
%      Session() (open-ephys-matlab-tools), extract stim / soundcard-trigger / laser /
%      running-wheel traces and timestamps, and save it (without continuous data)
%   3. messages, TTL, soundcard-trigger count check; DLC check
%   4. save OEinfo-<BonsaiFolder>.mat
%   5. SortedUnits via ProcessSpikes -- re-run automatically when the existing file is
%      STALE (see below)
%   6. Behavior file via ProcessCams; 7. Assimilation.mat via AssimilateSignals
%
%Stale SortedUnits ('reprocessSpikes','auto', the default) -- ProcessSpikes is re-run if
%   - there is no SortedUnits_*.mat, or more than one
%   - it is older than this session's spike_times.csv or the penetration's master
%     kilosort4/cluster_info.tsv (i.e. it came from an earlier sort or curation)
%   - this session uses spike_times.csv but SortedUnits lacks the 'shank' field, i.e. it
%     predates the 10.2.26 ProcessSpikes fix that stores the true probe channel
%   Old files are renamed *.stale-<timestamp> (not deleted) and restored if ProcessSpikes fails.
%
%usage:
%   PS = ProcessSession2()                     %Bonsai folder = pwd
%   PS = ProcessSession2(BonsaiPath)
%   PS = ProcessSession2(BonsaiPath, 'reprocessSession', true, 'reprocessSpikes', true, ...
%                        'reprocessCams', true, 'reprocessAssimilation', true)
%   ('reprocessSpikes' accepts 'auto' | true | false)
%
%returns PS, a struct whose fields have the same names as the variables the old script
%left in the workspace: BonsaiPath, BonsaiFolder, LocalDataRoot, DataRoot, EphysPath,
%EphysPath_KS, Bdirs, dirs, sessionfilename, session, stimtrace, soundcardtrigger,
%lasertrace, runningwheeltrace, timestamps, num_channels, bit_volts, OEsamplerate,
%OEversion, messages, TTL, numsoundcardtriggers, nummessages, DLC_exists,
%OEinfofilename, SortedUnitsFile, behaviorfilename, assimilationfilename, plus
%PS.log (what each step did; also appended to ProcessSession2.log in the Bonsai folder,
%PS.logfile) and PS.errors (non-fatal step failures).
%Fatal problems (no EphysPath, unsupported channel configuration, wrong sample rate) error.
%Supports the current data format only: new Open Ephys hierarchy, 128-ch probe (142 or 136
%channels). Older formats (64-ch NeuroNexus, tetrodes, Rig2 behavior, old kilosort layouts)
%are handled by ProcessSession_legacy_script.m.
%The working directory is restored on return.
%-mike/claude 10.2.26 (from ProcessSession script, mw 04.27.2024)

p=inputParser;
addOptional(p, 'BonsaiPath', pwd, @(x) ischar(x)||isstring(x));
addParameter(p, 'reprocessSession', false, @(x) islogical(x)||isnumeric(x));
addParameter(p, 'reprocessSpikes', 'auto', @(x) islogical(x)||isnumeric(x)||strcmpi(x,'auto'));
addParameter(p, 'reprocessCams', false, @(x) islogical(x)||isnumeric(x));
addParameter(p, 'reprocessAssimilation', false, @(x) islogical(x)||isnumeric(x));
if nargin<1, BonsaiPath=pwd; end
parse(p, BonsaiPath, varargin{:}); P=p.Results;

startdir=pwd; restoredir=onCleanup(@() cd(startdir));
PS=struct(); PS.log={}; PS.errors={};

%% paths
BonsaiPath=char(P.BonsaiPath);
if ~isfolder(BonsaiPath), error('ProcessSession2: not a folder: %s', BonsaiPath); end
cd(BonsaiPath); BonsaiPath=pwd; cd(startdir);    %absolute, no trailing separator
[LocalDataRoot, BonsaiFolder]=fileparts(BonsaiPath);
PS.BonsaiPath=BonsaiPath; PS.BonsaiFolder=BonsaiFolder; PS.LocalDataRoot=LocalDataRoot;
%log file in the Bonsai folder, appended on every run (each line written as it happens, so the
%log survives a fatal error part-way through). Only ProcessSession2's own messages go here --
%no diary, so callers' diaries (e.g. ProcessAllSessions.log) are left alone.
PS.logfile=fullfile(BonsaiPath, 'ProcessSession2.log');
PS=note(PS, '===== ProcessSession2 run started %s =====', char(datetime('now')));
PS=note(PS, 'options: reprocessSession=%d reprocessSpikes=%s reprocessCams=%d reprocessAssimilation=%d', ...
    P.reprocessSession, string(P.reprocessSpikes), P.reprocessCams, P.reprocessAssimilation);
if isempty(dir(fullfile(BonsaiPath,'Sky_*.csv'))) || isempty(dir(fullfile(BonsaiPath,'TTL_*.csv')))
    PS=note(PS, 'WARNING: no Sky_* or TTL_* csv files -- this may not be a Bonsai folder');
end

%% 1. EphysPath, EphysPath_KS
OEinfofilename=fullfile(BonsaiPath, sprintf('OEinfo-%s.mat', BonsaiFolder));
EphysPath=''; EphysPath_KS='';
if isfile(OEinfofilename)
    oe=load(OEinfofilename);
    if isfield(oe,'EphysPath'), EphysPath=oe.EphysPath; end
    if isfield(oe,'EphysPath_KS'), EphysPath_KS=oe.EphysPath_KS; end
    PS=note(PS, 'loaded %s', OEinfofilename);
end
if isempty(EphysPath)   %identify the Ephys folder by its notebook file (one level down first, then anywhere)
    ed=[dir(fullfile(BonsaiPath,'*','notebook.mat')); dir(fullfile(BonsaiPath,'*','notebook.json'))];
    fol=unique({ed.folder});
    if isempty(fol), error('ProcessSession2: could not find EphysPath (no notebook file) in %s', BonsaiPath); end
    if numel(fol)>1, error('ProcessSession2: more than one candidate EphysPath in %s:\n  %s', BonsaiPath, strjoin(fol, '\n  ')); end
    [~,EphysPath]=fileparts(fol{1});
    PS=note(PS, 'found EphysPath %s', EphysPath);
end
ephysdir=fullfile(BonsaiPath, EphysPath);
PS.EphysPath=EphysPath;

%Bdirs/dirs (session grouping), if present
[PS.Bdirs, PS.dirs]=deal({}); PS.DataRoot=LocalDataRoot;
if isfile(fullfile(BonsaiPath,'Bdirs.mat'))
    b=load(fullfile(BonsaiPath,'Bdirs.mat'));
    if isfield(b,'Bdirs'), PS.Bdirs=b.Bdirs; end
    if isfield(b,'dirs'), PS.dirs=b.dirs; end
    if isfield(b,'DataRoot'), PS.DataRoot=b.DataRoot; end
end

%% 2. Open Ephys session object
sessionfilename=fullfile(ephysdir, ['session-' EphysPath '.mat']);
PS.sessionfilename=sessionfilename;
tracevars={'stimtrace','soundcardtrigger','lasertrace','runningwheeltrace','timestamps','num_channels','bit_volts'};
if isfile(sessionfilename) && ~P.reprocessSession
    S=load(sessionfilename);
    PS=note(PS, 'loaded session object %s', sessionfilename);
    PS.session=S.session;
    for v=tracevars   %variables missing from older session files are returned EMPTY (never left over)
        if isfield(S, v{1}), PS.(v{1})=S.(v{1}); else, PS.(v{1})=[]; end
    end
    if isempty(PS.runningwheeltrace), PS=note(PS, 'session file has no runningwheeltrace (older file; use ''reprocessSession'',true to rebuild)'); end
else
    PS=note(PS, 'building Open Ephys session object (takes a few minutes)...');
    session=Session(ephysdir);   %open-ephys-matlab-tools
    rec=session.recordNodes{1}.recordings{1};
    info=rec.info.continuous(1);   %(1) = the 30 kHz data; Rig2 behavior data has a 2nd memory-usage stream
    num_channels=info.num_channels;
    bit_volts=arrayfun(@(c) c.bit_volts, info.channels(1:num_channels));
    if info.sample_rate~=30000, error('ProcessSession2: sample rate %g, expected 30000 -- check which continuous node is used', info.sample_rate); end
    keys=rec.continuous.keys(); cont=rec.continuous(keys{1});
    timestamps=cont.timestamps;
    c=aux_channels(num_channels);
    stimtrace=cont.samples(c.stim,:);              %read only the channels we need (not all samples)
    soundcardtrigger=cont.samples(c.sct,:);
    lasertrace=cont.samples(c.laser,:);
    if isempty(c.wheel), runningwheeltrace=[]; else, runningwheeltrace=cont.samples(c.wheel,:); end
    PS=note(PS, '%s running wheel trace', ternary(isempty(runningwheeltrace), 'NO', 'extracted'));
    %drop bulky continuous/spike data before saving the session object
    for n=1:numel(session.recordNodes)
        for r=1:numel(session.recordNodes{n}.recordings)
            try session.recordNodes{n}.recordings{r}.continuous=[]; catch, end
            try session.recordNodes{n}.recordings{r}.spikes=[]; catch, end
        end
    end
    save(sessionfilename, 'session', 'stimtrace','soundcardtrigger','lasertrace','runningwheeltrace','timestamps','num_channels','bit_volts');
    PS=note(PS, 'saved %s', sessionfilename);
    PS.session=session;
    [PS.stimtrace, PS.soundcardtrigger, PS.lasertrace, PS.runningwheeltrace, PS.timestamps, PS.num_channels, PS.bit_volts]= ...
        deal(stimtrace, soundcardtrigger, lasertrace, runningwheeltrace, timestamps, num_channels, bit_volts);
end

%% 3. messages, TTL, DLC
rec=PS.session.recordNodes{1}.recordings{1};
PS.OEsamplerate=rec.info.continuous(1).sample_rate;
PS.OEversion=rec.info.GUIVersion;
PS.messages=rec.messages('MessageCenter');
PS.TTL=rec.ttlEvents(erase(rec.info.continuous(1).folder_name, '/'));
PS.numsoundcardtriggers=height(PS.TTL)/2;   %TTL has rising and falling edges
PS.nummessages=height(PS.messages);
if PS.numsoundcardtriggers~=PS.nummessages
    PS=note(PS, 'WARNING: numsoundcardtriggers (%g) ~= nummessages (%d)', PS.numsoundcardtriggers, PS.nummessages);
end
PS.DLC_exists=~isempty(dir(fullfile(BonsaiPath,'*.pickle'))) && ~isempty(dir(fullfile(BonsaiPath,'Sky_m*0_el_filtered.csv'))) ...
    && isempty(dir(fullfile(BonsaiPath,'Sky_mouse*_preycap*_el_filtered.csv')));
if ~PS.DLC_exists, PS=note(PS, 'no (usable) DLC tracks for this session'); end

%kilosort folder: only the penetration's master session has one (phy.log); batch-sorted
%sessions just have spike_times.csv in the Ephys folder
if isempty(EphysPath_KS)
    ksd=[dir(fullfile(BonsaiPath,'*','phy.log')); dir(fullfile(BonsaiPath,'*','*','phy.log'))];
    if isempty(ksd) && ~isfile(fullfile(ephysdir,'spike_times.csv'))
        ksd=dir(fullfile(BonsaiPath,'**','phy.log'));   %last resort: deep search
    end
    if numel(ksd)==1
        EphysPath_KS=ksd.folder; PS=note(PS, 'found kilosort folder %s', EphysPath_KS);
    elseif numel(ksd)>1
        error('ProcessSession2: more than one kilosort folder (phy.log) in %s', BonsaiPath);
    elseif isfile(fullfile(ephysdir,'spike_times.csv'))
        PS=note(PS, 'no kilosort folder here; using batch-sorted spike_times.csv');
    else
        PS=note(PS, 'WARNING: no kilosort data found -- ProcessSpikes will fail');
    end
end
PS.EphysPath_KS=EphysPath_KS;

%% 4. OEinfo
PS.OEinfofilename=OEinfofilename;
OEsamplerate=PS.OEsamplerate; OEversion=PS.OEversion; messages=PS.messages; TTL=PS.TTL; numsoundcardtriggers=PS.numsoundcardtriggers; %#ok<NASGU>
save(OEinfofilename, 'EphysPath', 'EphysPath_KS', 'BonsaiPath', 'BonsaiFolder', 'OEversion', 'messages', 'OEsamplerate', 'numsoundcardtriggers', 'TTL');

%% 5. SortedUnits
PS.SortedUnitsFile='';
sud=dir(fullfile(ephysdir, 'SortedUnits_*.mat'));
[dorun, why]=need_spikes(P.reprocessSpikes, sud, ephysdir, PS, LocalDataRoot);
if ~dorun
    PS.SortedUnitsFile=fullfile(sud.folder, sud.name);
    PS=note(PS, 'using existing %s', sud.name);
else
    PS=note(PS, 'running ProcessSpikes: %s', why);
    stamp=char(datetime('now','Format','yyyyMMdd-HHmmss'));
    moved={};
    for k=1:numel(sud)   %set old files aside (restored if ProcessSpikes fails)
        old=fullfile(sud(k).folder, sud(k).name); movefile(old, [old '.stale-' stamp]); moved{end+1}=old; %#ok<AGROW>
    end
    try
        cd(ephysdir);
        f=ProcessSpikes(BonsaiPath, EphysPath_KS, EphysPath, LocalDataRoot);
        PS.SortedUnitsFile=char(f);
        if ~isempty(moved), PS=note(PS, 'previous SortedUnits kept as *.stale-%s', stamp); end
    catch err
        PS.errors{end+1}=['ProcessSpikes failed: ' err.message];
        PS=note(PS, 'ERROR: ProcessSpikes failed: %s', err.message);
        for k=1:numel(moved), if ~isfile(moved{k}), movefile([moved{k} '.stale-' stamp], moved{k}); end, end
        if ~isempty(moved), PS=note(PS, 'restored previous SortedUnits file(s); SortedUnitsFile left EMPTY so callers fail loudly'); end
    end
    cd(BonsaiPath);
end

%% 6. Behavior (ProcessCams)
PS.behaviorfilename='';
bd=dir(fullfile(BonsaiPath,'Behavior*.mat'));
if isempty(bd) || P.reprocessCams
    try
        cd(BonsaiPath); ProcessCams;
        bd=dir(fullfile(BonsaiPath,'Behavior*.mat'));
        PS=note(PS, 'ran ProcessCams');
    catch err
        PS.errors{end+1}=['ProcessCams failed: ' err.message]; PS=note(PS, 'ERROR: ProcessCams failed: %s', err.message);
    end
end
if ~isempty(bd), PS.behaviorfilename=fullfile(bd(1).folder, bd(1).name); end

%% 7. Assimilation
PS.assimilationfilename='';
af=fullfile(BonsaiPath,'Assimilation.mat');
if ~isfile(af) || P.reprocessAssimilation
    try
        cd(BonsaiPath);
        [vids,units,chans]=AssimilateSignals; %#ok<ASGLU>
        save(af, 'vids','units','chans');
        PS=note(PS, 'ran AssimilateSignals');
    catch err
        PS.errors{end+1}=['AssimilateSignals failed: ' err.message]; PS=note(PS, 'AssimilateSignals failed (possibly no camera data): %s', err.message);
    end
end
if isfile(af), PS.assimilationfilename=af; end

PS=note(PS, 'ProcessSession2 done%s', ternary(isempty(PS.errors), '', sprintf(' with %d error(s)', numel(PS.errors))));
end

%% ---------------------------------------------------------------------
function [dorun, why]=need_spikes(mode, sud, ephysdir, PS, LocalDataRoot)
%decide whether ProcessSpikes must (re)run
if islogical(mode) || isnumeric(mode)
    dorun=logical(mode); why=ternary(dorun, 'forced', 'skipped (reprocessSpikes=false)');
    if ~dorun && numel(sud)~=1, dorun=true; why=sprintf('%d SortedUnits files found', numel(sud)); end
    return
end
dorun=true;
if isempty(sud), why='no SortedUnits file'; return, end
if numel(sud)>1, why=sprintf('%d SortedUnits files found', numel(sud)); return, end
t=sud.datenum;
stf=dir(fullfile(ephysdir,'spike_times.csv'));
if ~isempty(stf) && stf.datenum>t, why='SortedUnits is older than spike_times.csv'; return, end
ci=master_cluster_info(PS, LocalDataRoot);
if ~isempty(ci) && ci.datenum>t, why='SortedUnits is older than the master cluster_info.tsv'; return, end
if ~isempty(stf)
    s=load(fullfile(sud.folder, sud.name), 'SortedUnits');
    if ~isfield(s.SortedUnits, 'shank'), why='SortedUnits predates the probe-channel fix (no shank field)'; return, end
    if ~isempty(s.SortedUnits) && all(isnan([s.SortedUnits.channel])), why='SortedUnits has no probe channels (all NaN -- channel lookup failed last time)'; return, end
end
dorun=false; why='';
end

function ci=master_cluster_info(PS, LocalDataRoot)
%dir() entry of the penetration's master kilosort4/cluster_info.tsv, or [] if not found
ci=[];
cands={};
if ~isempty(PS.EphysPath_KS), cands{end+1}=fullfile(PS.EphysPath_KS, 'cluster_info.tsv'); end
if ~isempty(PS.Bdirs)
    b1=regexp(regexprep(PS.Bdirs{1}, '[/\\]+$', ''), '[^/\\]+$', 'match', 'once');
    cands{end+1}=fullfile(LocalDataRoot, b1, 'kilosort4', 'cluster_info.tsv');
end
for k=1:numel(cands), if isfile(cands{k}), ci=dir(cands{k}); return, end, end
end

function c=aux_channels(n)
%auxiliary (ADC) channel indices by total channel count
switch n
    case 142, c=struct('stim',135, 'sct',136, 'laser',137, 'wheel',138);    %128-ch Diagnostic Biochips probe
    case 136, c=struct('stim',129, 'sct',130, 'laser',131, 'wheel',{[]});   %128-ch probe, AUX channels not recorded
    otherwise
        error(['ProcessSession2: unsupported channel configuration (%d channels). Only the 128-ch probe ' ...
            '(142 or 136 channels) is supported; for older formats (64-ch NeuroNexus, tetrodes, Rig2 ' ...
            'behavior) use ProcessSession_legacy_script.m'], n);
end
end

function out=ternary(c, a, b)
if c, out=a; else, out=b; end
end

function PS=note(PS, varargin)
%print, keep in PS.log, and append to PS.logfile (if set)
msg=sprintf(varargin{:}); fprintf('\n%s', msg); PS.log{end+1}=msg;
if isfield(PS, 'logfile')
    fid=fopen(PS.logfile, 'a');
    if fid>0, fprintf(fid, '%s\n', msg); fclose(fid); end
end
end
