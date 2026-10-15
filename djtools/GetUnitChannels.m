function U = GetUnitChannels(BonsaiPath, EphysPath, varargin)
%GetUnitChannels - probe channel, shank and depth of every Kilosort cluster in a session.
%
%All sessions of a penetration are sorted together; the master Kilosort output is
%in Bdirs{1}/kilosort4 (Bdirs from the session's dirs.mat). Its cluster_info.tsv
%(written by Phy) gives, per cluster:
%   ch    = acquisition (Open Ephys) channel, 0-based -- a value of channel_map.npy
%   depth = y position on the probe (um above the tip)
%   sh    = shank (1 or 2)
%depths.mat (and sink_chans) index channels 1..128 in DEPTH order: 1..64 = shank 1
%from the tip up, 65..128 = shank 2. That index is the row of ch in channel_map.npy
%(1-based), which equals round(y/pitch)+1 + 64*(shank-1) for this probe.
%
%usage:
%   U = GetUnitChannels(BonsaiPath, EphysPath)         %EphysPath: folder name or full path
%   U = GetUnitChannels(..., 'ksdir', path)            %explicit kilosort4 folder
%returns table U, one row per cluster in cluster_info.tsv:
%   cluster, group (Phy curation: good/mua/noise), KSLabel, ch_OE (0-based acquisition
%   channel), chanIdx (1..128 depth-ordered; index into depths.mat), shank, yPos (um above
%   tip), corrected_depth, angle_corrected_depth (from depths.mat in the Ephys folder, NaN
%   if absent/empty), amp, fr, n_spikes, ContamPct
%-mike/claude 10.2.26

p=inputParser;
addParameter(p, 'ksdir', '', @(x) ischar(x)||isstring(x));
addParameter(p, 'pitch', 20, @isnumeric);
addParameter(p, 'chansPerShank', 64, @isnumeric);
parse(p, varargin{:}); P=p.Results;
BonsaiPath=char(BonsaiPath); EphysPath=char(EphysPath);
ed=EphysPath; if ~isfolder(ed), ed=fullfile(BonsaiPath, EphysPath); end

%% locate the master kilosort4 folder
ks=char(P.ksdir);
if isempty(ks)
    d=load(fullfile(ed, 'dirs.mat'), 'Bdirs');
    b1=regexp(regexprep(d.Bdirs{1}, '[/\\]+$', ''), '[^/\\]+$', 'match', 'once');
    ks=fullfile(fileparts(BonsaiPath), b1, 'kilosort4');   %Bdirs{1} relative to the local data root
    if ~isfolder(ks) && isfolder(fullfile(BonsaiPath, 'kilosort4')), ks=fullfile(BonsaiPath, 'kilosort4'); end
end
assert(isfile(fullfile(ks, 'cluster_info.tsv')), 'no cluster_info.tsv in %s', ks);

%% clusters
T=readtable(fullfile(ks, 'cluster_info.tsv'), 'FileType','text', 'Delimiter','\t', 'TextType','string');
U=table(T.cluster_id, T.group, T.KSLabel, T.ch, 'VariableNames', {'cluster','group','KSLabel','ch_OE'});
U.shank=double(T.sh); U.yPos=double(T.depth);

%% depth-ordered channel index (1..128)
cmf=fullfile(ks, 'channel_map.npy'); cpf=fullfile(ks, 'channel_positions.npy');
if exist('readNPY', 'file') && isfile(cmf) && isfile(cpf)
    cm=double(readNPY(cmf)); cp=double(readNPY(cpf));
    [tf, row]=ismember(U.ch_OE, cm(:));
    assert(all(tf), 'cluster_info ch values not found in channel_map.npy');
    U.chanIdx=row;                                         %1-based row = depth-ordered index
    if any(abs(cp(row,2)-U.yPos)>1e-6), warning('channel_positions y does not match cluster_info depth for %d clusters', sum(abs(cp(row,2)-U.yPos)>1e-6)); end
else   %no readNPY: same mapping from position
    U.chanIdx=round(U.yPos/P.pitch)+1+P.chansPerShank*(U.shank-1);
end

%% cortical depth from depths.mat
[U.corrected_depth, U.angle_corrected_depth]=deal(nan(height(U),1));
df=fullfile(ed, 'depths.mat');
if isfile(df)
    D=load(df);
    if isfield(D,'corrected_depth') && ~isempty(D.corrected_depth), x=D.corrected_depth(:); U.corrected_depth=x(U.chanIdx); end   %(:) first: depths.mat stores rows
    if isfield(D,'angle_corrected_depth') && ~isempty(D.angle_corrected_depth), x=D.angle_corrected_depth(:); U.angle_corrected_depth=x(U.chanIdx); end
end
U.amp=T.amp; U.fr=T.fr; U.n_spikes=T.n_spikes; U.ContamPct=T.ContamPct;
U.Properties.UserData.ksdir=ks;
end
