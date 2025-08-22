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
% This really should be converted into a function because as a script it's
% brittle and could be carrying over the wrong info.
%this would require changing all calls to ProcessSession to return outputs
% function [Bdirs, dirs, BonsaiFolder, OEsamplerate, ...
%     BonsaiPath, OEversion, soundcardtrigger, SortedUnitsFile, ...
%     lasertrace, stimtrace, TTL, messages, EphysPath, assimilationfilename, num_channels, timestamps, ...
%     EphysPath_KS, behaviorfilename, nummessages, bit_volts, numsoundcardtriggers, LocalDataRoot]=ProcessSession
%
% these are the variables we potentially might want to return
% Bdirs                 OEinfofilename        dirs                  session
% BonsaiFolder          OEsamplerate          dsky                  sessionfilename
% BonsaiPath            OEversion             dttl                  soundcardtrigger
% DLC_exists            SortedUnitsFile       lasertrace            stimtrace
% DataRoot              TTL                   messages              sudir
% EphysPath             assimilationfilename  num_channels          timestamps
% EphysPath_KS          behaviorfilename      nummessages
% ForceReprocess        bit_volts             numsoundcardtriggers
% LocalDataRoot         d                     oe

clear EphysPath

ForceReprocess=0;

% check to make sure we're in the bonsai folder
BonsaiPath= pwd;
cd ..
LocalDataRoot=pwd; %parent directory of BonsaiPath
cd(BonsaiPath)

dsky=dir('Sky_*.csv');
dttl=dir('TTL_*.csv');
if isempty(dsky) | isempty(dttl)
    warning('I don''t think this is a Bonsai folder because I can''t find any Sky_mouse=* or TTL_mouse-* files')
end

%get EphysPath
cd(BonsaiPath)
[~,BonsaiFolder,~]=fileparts(BonsaiPath); %BonsaiFolder is just the timestamp-mouseID identifier of the Bonsai directory (excluding the abolute path)
OEinfofilename=sprintf('OEinfo-%s.mat', BonsaiFolder);
if exist(OEinfofilename)==2
    oe=load(OEinfofilename);
    fprintf('\nfound and loaded %s\n', OEinfofilename)
    EphysPath=oe.EphysPath;
    EphysPath_KS=oe.EphysPath_KS;
else
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
    ed=dir('**/notebook.mat');
    if isempty(ed)
        error('could not find EphysPath in this Bonsai directory')
        EphysPath_exists=0;
        EphysPath=[];
    elseif length(ed)==1
        [~,EphysPath,~]=fileparts(ed.folder);
    else error('more than one candidate EphysPath folder found')
    end
    fprintf('\nfound EphysPath %s\n', EphysPath)
end

%load OpenEphys session information
sessionfilename=['session-',EphysPath];
cd(EphysPath)
try
    load(sessionfilename)
catch
    fprintf('\ndid not find saved Open Ephys session object for this directory ...')
    fprintf('\nloading Open Ephys data (this will take a couple minutes) ...\n')
    session = Session(pwd); %open-ephys-matlab-tools function

    num_channels=session.recordNodes{1}.recordings{1}.info.continuous.num_channels;
    for ch=1:num_channels
        %bit_volts(ch)=session.recordNodes{1}.recordings{1}.info.continuous.channels(ch).bit_volts;
        bit_volts(ch)=session.recordNodes{1}.recordings{1}.info.continuous(1).channels(ch).bit_volts; %changed because gap detection behavior data on rig2 with new OE has 2 continuous nodes. (1) should be the data at 30kHz, (2) is some memory usage stream
    end
    % double-check samplerate to make sure we're getting the right continuous node 
    OEsamplerate=session.recordNodes{1}.recordings{1}.info.continuous(1).sample_rate;
    if OEsamplerate~= 30000, error('wrong samplerate, check which OE continuous node we''re extracting from'), end

keys=session.recordNodes{1}.recordings{1}.continuous.keys();
    key=keys{1};
    timestamps=session.recordNodes{1}.recordings{1}.continuous(key).timestamps;
    samples=session.recordNodes{1}.recordings{1}.continuous(key).samples(:,:);
    if num_channels==78 %for 64 channel neuronexus probe
        stimtracech=71;
        soundcardtriggerch=72;
        lasertracech=73;
    elseif num_channels==142 %for 128 channel diagnostic biochips probe
        stimtracech=135;
        soundcardtriggerch=136;
        lasertracech=137;
    elseif num_channels==136 %for 128 channel diagnostic biochips probe where we forgot to record the 6 AUX channels
        stimtracech=129;
        soundcardtriggerch=130;
        lasertracech=131;
    elseif num_channels==43 %32-ch tetrode config
        stimtracech=36;
        soundcardtriggerch=37;
        lasertracech=38;
    elseif num_channels==11 %Rig2 gap detection behavior config
        stimtracech=4;
        soundcardtriggerch=5;
        lasertracech=6; %I think, but haven't confirmed
        %piezo data from mice 1-4 are on chans 8,9,10,11. Chans 1-3 are empty
    end
    stimtrace=session.recordNodes{1}.recordings{1}.continuous(key).samples(stimtracech,:);
    soundcardtrigger=session.recordNodes{1}.recordings{1}.continuous(key).samples(soundcardtriggerch,:);
    lasertrace=session.recordNodes{1}.recordings{1}.continuous(key).samples(lasertracech,:);


    %it is slow, takes 3 minutes. so we save for fast loading in the
    %future.
    % manually delete any continuous data fields prior to saving a copy of
    % the Session to avoid creating a massive session file
    session.recordNodes{1}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.spikes=[];
    fprintf('\nsaving Open Ephys session object in OpenEphys folder...')
    save(sessionfilename, 'session', 'stimtrace','soundcardtrigger','lasertrace','timestamps','num_channels', 'bit_volts')
    %save('samples', 'samples')
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
if ~exist('EphysPath_KS') | isempty(EphysPath_KS)
    cd(BonsaiPath)
    fprintf('\nsearching for kilosort folder...')
    tic
    ksd=dir('**/phy.log');
    if isempty(ksd)
        warning('could not find kilosort data')
        EphysPath_KS=[];
        fprintf('\nProcess Spikes will fail because there is no kilosort data')

    elseif length(ksd)==1
        fprintf('\tfound it. (took %.0fs to find)', toc)
        EphysPath_KS=ksd.folder;
        fprintf('\nfound EphysPath_KS kilosort folder: \n%s\n', EphysPath_KS)
    else error('more than one candidate KS folder found')
    end
end

% dirName=split(session.recordNodes{1}.recordings{1}.directory, filesep);
% EphysPath_KS=fullfile(filesep, dirName{[1:end-3, end-4]}); %absolute path

%check for dirs and bdirs files, if needed create them by calling
%SettingYourStage. Alternatively, if you're only kilosorting individual
%sessions (so there's no need to group sessions that were kilosorted
%together) then you can just use the current directory
if 0 %group sessions that were kilosorted together
    SettingYourStage;
else %assume we're only kilosorting individual sessions
    Bdirs{1}=BonsaiPath;
    dirs{1}=EphysPath;
    DataRoot=LocalDataRoot;
    cd(BonsaiPath)
    save('Bdirs.mat','dirs','Bdirs','DataRoot');
    cd(EphysPath);
    save('dirs.mat', 'Bdirs', 'dirs','DataRoot');
    cd(BonsaiPath)
end

%get SortedUnits from ProcessSpikes
try
    cd(BonsaiPath)
    cd(EphysPath)
    sudir=dir('./SortedUnits_*.mat'); %check if SortedUnits file output from ProcessSpikes already exists

    if isempty(sudir)
        %    [SortedUnitsFile] = ProcessSpikes(EphysPath_KS,BonsaiPath, LocalDataRoot); %Process the spikes
        % [SortedUnitsFile] = ProcessSpikes(EphysPath_KS, LocalDataRoot); %Process the spikes
        SortedUnitsFile = ProcessSpikes(BonsaiPath, EphysPath_KS, EphysPath, LocalDataRoot) %new way for new OE hierarchy
    else
        fprintf('\nfound SortedUnits file, skipping ProcessSpikes')
        SortedUnitsFile = fullfile(sudir.folder,sudir.name);
    end

catch
    lasterr
    warning('ProcessSpikes failed')
    if isempty(EphysPath_KS) warning('ProcessSpikes failed. Probably because there is no kilosort data.'); end

end

% these are the things I will want to add to Sky
% I'm saving them here for convenience, because sometimes we want to re-add them if we regenerate Sky in the future
[~,BonsaiFolder,~]=fileparts(BonsaiPath); %BonsaiFolder is just the timestamp-mouseID identifier of the Bonsai directory (excluding the abolute path)
OEinfofilename=sprintf('OEinfo-%s', BonsaiFolder);
OEsamplerate=session.recordNodes{1}.recordings{1}.info.continuous(1).sample_rate;
cd(BonsaiPath)
save(OEinfofilename, 'EphysPath', 'EphysPath_KS', 'EphysPath_KS', 'BonsaiPath', 'BonsaiFolder', 'OEversion', 'messages', 'OEsamplerate', 'numsoundcardtriggers', 'TTL');

try
    behaviorfilename = dir('Behavior*.mat');
    if isempty(behaviorfilename) | ForceReprocess
        fprintf('\nno Behavior file found, calling ProcessCams')
        behaviorfile=ProcessCams;
        fprintf('done\n')
    else
        fprintf('\nBehavior file already exists, skipping ProcessCams')
    end
catch
    warning('\n ProcessCams failed, probably because there is no camera data')
    lasterr
end


cd(LocalDataRoot)
cd(BonsaiPath)
try
    assimilationfilename = dir('Assimilation.mat');
    if isempty(assimilationfilename)  | ForceReprocess
        [vids,units,chans] = AssimilateSignals;
        save Assimilation vids units chans
    end
catch
    warning('\n AssimilateSignals failed, possibly because there is no camera data')
end

fprintf('\nProcessSession done\n')




