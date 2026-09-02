function [dirs, Bdirs, DataRoot] = BuildDirsFromRecordingLengths(varargin)
%BuildDirsFromRecordingLengths - reconstruct the historical dirs.mat /
%Bdirs.mat variables (dirs, Bdirs, DataRoot) from a recording_lengths.json
%file written by the new split function, for backwards compatibility with
%code (ProcessTC_LFP2, ProcessGPIAS_PSTH2, PlotTC_LFP2, PlotGPIAS_PSTH2,
%etc.) that expects to find these variables in dirs.mat/Bdirs.mat.
%
%Background: each recording session has a Bonsai folder (outer, written
%by Bonsai) containing an Ephys folder (inner, written by OpenEphys).
%Several sessions recorded from the same site (typically 6, one per
%stimulus) are kilosorted together, then split back into their
%respective directories. dirs.mat/Bdirs.mat historically tracked this
%grouping -- dirs holds the Ephys folder name for each grouped session,
%Bdirs holds the corresponding Bonsai folder path, and DataRoot is the
%shared parent directory of the Bonsai folders (useful because the same
%NAS hierarchy looks different from different computers -- see
%macifypath.m / FixDataRoot.m). This grouping is what lets you, e.g.,
%find depths.mat from the CSD session and apply it to the other sessions
%in the group. dirs.mat and Bdirs.mat both contain all three variables
%(dirs, Bdirs, DataRoot), but historically live in different places for
%each grouped session: dirs.mat goes in the (inner) Ephys directory, and
%Bdirs.mat goes in the (outer) Bonsai directory.
%
%The new split function no longer writes dirs.mat/Bdirs.mat directly;
%instead it writes recording_lengths.json, a flat map of
%    "<full path to continuous.dat>": {"length_samples": <n>, "folder": "<Bonsai folder path>"}
%with one entry per grouped session. This function parses that file and
%reconstructs dirs/Bdirs/DataRoot from it.
%
%usage:
%   [dirs, Bdirs, DataRoot] = BuildDirsFromRecordingLengths()
%       looks for recording_lengths.json in the current directory
%   [dirs, Bdirs, DataRoot] = BuildDirsFromRecordingLengths(jsonpath)
%       jsonpath can be the .json file itself, or a directory containing
%       recording_lengths.json
%   BuildDirsFromRecordingLengths(..., 'write', true)
%       also writes, for every grouped session i (each containing dirs,
%       Bdirs, and DataRoot):
%         dirs.mat  -> fullfile(Bdirs{i}, dirs{i})  (the Ephys directory)
%         Bdirs.mat -> Bdirs{i}                     (the Bonsai directory)
%
%example:
%   BuildDirsFromRecordingLengths('/Volumes/Projects/5XFAD/Rig3Phys/2025-11-19_6-22-38_mouse-4091', 'write', true)
%
%-mike/claude 9.2.26

%if called with only name-value pairs and no jsonpath (e.g.
%BuildDirsFromRecordingLengths('write', true)), inputParser would
%otherwise try to consume 'write' positionally as the (optional)
%jsonpath value, then choke on `true` where it expects a parameter
%name. Detect that case and insert the default jsonpath explicitly.
knownparams={'write'};
if ~isempty(varargin) && (ischar(varargin{1})||isstring(varargin{1})) && any(strcmpi(varargin{1}, knownparams))
    varargin=[{pwd}, varargin];
end

p=inputParser;
addOptional(p, 'jsonpath', pwd, @(x) ischar(x)||isstring(x));
addParameter(p, 'write', false, @(x) islogical(x)||isnumeric(x));
parse(p, varargin{:});
jsonpath=char(p.Results.jsonpath);
dowrite=logical(p.Results.write);

if isfolder(jsonpath)
    jsonpath=fullfile(jsonpath, 'recording_lengths.json');
end
if ~isfile(jsonpath)
    error('BuildDirsFromRecordingLengths:notfound', 'could not find recording_lengths.json at %s', jsonpath);
end

txt=fileread(jsonpath);

%Parse the flat {key: {length_samples, folder}} structure directly with
%regexp rather than jsondecode. jsondecode would mangle the original key
%strings (which contain '/', '\', and spaces) into valid MATLAB
%fieldnames, and we need the exact key text to recover the Ephys folder
%name -- it isn't stored as its own field in the json.
pat='"((?:[^"\\]|\\.)*)"\s*:\s*\{\s*"length_samples"\s*:\s*(-?\d+)\s*,\s*"folder"\s*:\s*"((?:[^"\\]|\\.)*)"\s*\}';
tok=regexp(txt, pat, 'tokens');
if isempty(tok)
    error('BuildDirsFromRecordingLengths:parse', 'could not parse any entries out of %s -- unexpected format?', jsonpath);
end

n=numel(tok);
dirs=cell(1,n);
Bdirs=cell(1,n);
for i=1:n
    key=unescape_json_string(tok{i}{1});
    folder=unescape_json_string(tok{i}{3});

    %the Ephys folder name is the first path component of key
    %immediately following the Bonsai folder (=folder)
    if ~startsWith(key, folder)
        error('BuildDirsFromRecordingLengths:mismatch', 'entry %d: key does not start with its own folder:\n key:    %s\n folder: %s', i, key, folder);
    end
    rest=key(length(folder)+1:end);
    rest=regexprep(rest, '^[/\\]+', ''); %strip leading separator(s)
    parts=regexp(rest, '[/\\]+', 'split');
    if isempty(parts) || isempty(parts{1})
        error('BuildDirsFromRecordingLengths:mismatch', 'entry %d: could not find an Ephys folder name after the Bonsai folder in key:\n%s', i, key);
    end
    dirs{i}=parts{1};
    Bdirs{i}=resolve_local_path(folder);
end

%DataRoot = the shared parent of all the (resolved) Bonsai folders
parentdirs=cell(1,n);
for i=1:n
    parentdirs{i}=fileparts(Bdirs{i});
end
DataRoot=parentdirs{1};
if ~all(strcmp(parentdirs, DataRoot))
    warning('BuildDirsFromRecordingLengths:datarootmismatch', ...
        'not all grouped folders in %s share the same parent directory -- using the first one (%s). Found:\n%s', ...
        jsonpath, DataRoot, strjoin(unique(parentdirs), '\n'));
end

fprintf('\nfound %d grouped recording(s) in %s:', n, jsonpath);
for i=1:n
    fprintf('\n  dirs{%d}  = %s', i, dirs{i});
    fprintf('\n  Bdirs{%d} = %s', i, Bdirs{i});
end
fprintf('\nDataRoot = %s\n', DataRoot);

if dowrite
    for i=1:n
        ephysdir=fullfile(Bdirs{i}, dirs{i});
        bonsaidir=Bdirs{i};

        if ~isfolder(ephysdir)
            warning('BuildDirsFromRecordingLengths:nodir', 'Ephys directory does not exist, skipping dirs.mat: %s', ephysdir);
        else
            save(fullfile(ephysdir, 'dirs.mat'), 'dirs', 'Bdirs', 'DataRoot');
            fprintf('\nwrote dirs.mat to %s', ephysdir);
        end

        if ~isfolder(bonsaidir)
            warning('BuildDirsFromRecordingLengths:nodir', 'Bonsai directory does not exist, skipping Bdirs.mat: %s', bonsaidir);
        else
            save(fullfile(bonsaidir, 'Bdirs.mat'), 'dirs', 'Bdirs', 'DataRoot');
            fprintf('\nwrote Bdirs.mat to %s', bonsaidir);
        end
    end
    fprintf('\n');
end

end

function s=unescape_json_string(s)
%minimal JSON string unescaping for the characters that actually appear
%in file paths (backslash, quote); extend if needed.
s=strrep(s, '\\', '\'); %JSON-escaped backslash -> literal backslash
s=strrep(s, '\"', '"'); %JSON-escaped quote -> literal quote (shouldn't occur in a path, but just in case)
end

function localpath=resolve_local_path(rawpath)
%best-effort conversion of a raw (as-collected) NAS path to this local
%machine's view of it. Defers to the existing macifypath.m convention
%used throughout djtools/djmaus; adds a case-insensitive fallback for
%the wehr-nas share since recording_lengths.json capitalizes it
%('Projects') differently than macifypath.m's hardcoded pattern
%('projects') -- strrep in macifypath.m is case-sensitive, so without
%this fallback the conversion silently would not fire on a Mac.
localpath=macifypath(rawpath);
if ismac
    localpath=regexprep(localpath, '//wehr-nas\.uoregon\.edu/projects/', '/Volumes/Projects/', 'ignorecase');
end
end
