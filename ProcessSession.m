% ProcessOpenEphysSession
%
% pre-process data recorded in new open ephys data format using Session function from open-ephys-matlab-tools
% also processes and aligns camera and behavioral data
% run this from a Bonsai folder
%
%   should extend to work with dirs/bdirs if those are part of the workflow
%
% in the future we should expand this to scan and processes all data directories contained in the current
% directory
%mw 04.27.2024
%
ForceReprocess=0;

% check to make sure we're in the bonsai folder 
BonsaiPath= pwd;
cd ..
LocalDataRoot=pwd; %parent directory of BonsaiPath
cd(BonsaiPath)

dsky=dir('Sky_mouse-*');
dttl=dir('TTL_mouse-*');
if isempty(dsky) | isempty(dttl)
    error('I don''t think this is a Bonsai folder because I can''t find any Sky_mouse=* or TTL_mouse-* files')
end

%get EphysPath
EphysPath_exists=1;
% % ed=dir('2024-*'); %this will obviously only work for 2024 data
% % %we should add deeper checking for EphysPath
% % if length(ed)==1 & ed.isdir
% %     EphysPath=ed.name;
% % else
% %     error ('cannot find ephys dir')
% %     EphysPath_exists=0;
% % end
%Here's a better way: look for notebook.mat file to identify EphysPath
cd(BonsaiPath)
ed=dir('**/notebook.mat');
if isempty(d)
    error('could not find EphysPath in this Bonsai directory')
    EphysPath_exists=0;
    EphysPath=[];
elseif length(ed)==1
    [~,EphysPath,~]=fileparts(ed.folder);
else error('more than one candidate EphysPath folder found')
end
fprintf('\nfound EphysPath %s\n', EphysPath)

%load OpenEphys session information
sessionfilename=['session-',EphysPath];
cd(EphysPath)
try
    load(sessionfilename)
catch
    fprintf('\ndid not find saved Open Ephys session object for this directory ...')
    fprintf('\nloading Open Ephys data (this will take a couple minutes) ...\n')
    session = Session(pwd); %open-ephys-matlab-tools function
    %it is slow, takes 3 minutes. so we save for fast loading in the
    %future.
    % manually delete any continuous data fields prior to saving a copy of
    % the Session to avoid creating a massive session file
    session.recordNodes{1}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.spikes=[];
    fprintf('\nsaving Open Ephys session object in OpenEphys folder...')
    save(sessionfilename, 'session')
    fprintf(' done. ')
end
OEversion = session.recordNodes{1}.recordings{1}.info.GUIVersion;
messages=session.recordNodes{1}.recordings{1}.messages('MessageCenter');
try
    %this ttl map key works for Rig4
    TTL=session.recordNodes{1}.recordings{1}.ttlEvents('Acquisition_Board-100.Rhythm Data-A');
catch
    %this ttl map key works for Rig3
    TTL=session.recordNodes{1}.recordings{1}.ttlEvents('OE_FPGA_Acquisition_Board-100.Rhythm Data');
    if isempty(TTL)
        error('empty TTL, what is wrong?')
    end
    %if this fails we might need to add another try...catch block
end
% TTL has both the rising and falling times
% messages includes GetRecordingPath
% so height(TTL)==2*(height(messages)-1)
numsoundcardtriggers=height(TTL)/2;
nummessages=height(messages);
if numsoundcardtriggers ~= nummessages
    warning('numsoundcardtriggers ~= nummessages')
end
cd(BonsaiPath) % return to where we started

%check for DLC
DLC_exists=1;
cd(BonsaiPath)
if length(dir('*.pickle'))==0
    warning('no DLC tracks for this session') %analyze video step
    DLC_exists=0;
end
if length(dir('Sky_m*0_el_filtered.csv'))==0
    warning('no DLC tracks for this session') %tracks output to csv step
    DLC_exists=0;
end
if length(dir('Sky_mouse*_preycap*_el_filtered.csv'))>0
    %change this to check for any DLC model you want to exclude
    warning('wrong DLC model used for this session') %tracks output to csv step
    DLC_exists=0;
end
        
%get path to kilosorted data 
cd(BonsaiPath)
d=dir('**/phy.log');
if isempty(d)
    warning('could not find kilosort data')
    EphysPath_KS=[];
elseif length(d)==1
    EphysPath_KS=d.folder;
else error('more than one candidate KS folder found')
end
fprintf('\nfound EphysPath_KS kilosort folder: \n%s\n', EphysPath_KS)

% dirName=split(session.recordNodes{1}.recordings{1}.directory, filesep);
% EphysPath_KS=fullfile(filesep, dirName{[1:end-3, end-4]}); %absolute path

%get SortedUnits from ProcessSpikes
try
    cd(EphysPath_KS) %abs path
sudir=dir('./SortedUnits_*.mat'); %check if SortedUnits file output from ProcessSpikes already exists
if isempty(sudir)
%    [SortedUnitsFile] = ProcessSpikes(EphysPath_KS,BonsaiPath, LocalDataRoot); %Process the spikes
    % [SortedUnitsFile] = ProcessSpikes(EphysPath_KS, LocalDataRoot); %Process the spikes
    [SortedUnitsFile] = ProcessSpikes(BonsaiPath, EphysPath_KS, EphysPath, LocalDataRoot) %new way for new OE hierarchy
end
catch
    warning('ProcessSpikes failed')
    if isempty(EphysPath_KS) warning ('ProcessSpikes failed. Probably because there is no kilosort data.'); end
end

% these are the things I will want to add to Sky
% I'm saving them here for convenience, because sometimes we want to re-add them if we regenerate Sky in the future
[~,BonsaiFolder,~]=fileparts(BonsaiPath); %BonsaiFolder is just the timestamp-mouseID identifier of the Bonsai directory (excluding the abolute path)
OEinfofilename=sprintf('OEinfo-%s', BonsaiFolder);
OEsamplerate=session.recordNodes{1}.recordings{1}.info.continuous.sample_rate;
cd(BonsaiPath)
save(OEinfofilename, 'EphysPath', 'EphysPath_KS', 'EphysPath_KS', 'BonsaiPath', 'BonsaiFolder', 'OEversion', 'messages', 'OEsamplerate', 'numsoundcardtriggers', 'TTL');

behaviorfilename = dir('Behavior*.mat');
if isempty(behaviorfilename) | ForceReprocess
    fprintf('\nno Behavior file found, calling ProcessCams')
    behaviorfile=ProcessCams;
    fprintf('done\n')
else
    fprintf('\nBehavior file already exists, skipping ProcessCams')
end

cd(LocalDataRoot)
cd(BonsaiPath)
assimilationfilename = dir('Assimilation.mat');
if isempty(assimilationfilename)  | ForceReprocess
    [vids,units,chans] = AssimilateSignals;
    save Assimilation vids units chans
end
    
fprintf('\nProcessSession done\n')




