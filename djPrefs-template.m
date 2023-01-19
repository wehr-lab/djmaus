function djPrefs
%Local preferences
%specifics for the local machine are at top
%preferences for individual users at the bottom
global pref SP

% machine-wide default prefs:
pref.users={'lab','mike-home','mike-kombi','mike-linux'};
try
    pref.soundcarddeviceID=GetAsioLynxDevice; %note: using the asio device for both input and output
catch
    pref.soundcarddeviceID=42;
end
pref.num_soundcard_outputchannels=4;
pref.Soundchannel1=1; %the sound
pref.Soundchannel2=nan; %no binaural
pref.Soundcardtriggerchannel=3;
pref.Laserchannel=4;
pref.Shockchannel=2;
pref.CampulseOff=nan;
pref.CampulseOn=nan;
pref.SCT_digital_line_in=1;
pref.reqlatencyclass=0; %Rig1: saw dropouts with 0, fewer 1, fewer with 2, still some with 3
pref.suggestedLatency=.01;
pref.SoundFs=192000;
pref.maxSPL=80;
pref.allmouseIDs='';
pref.root=fileparts(which(mfilename));
pref.windowpos=[2853 861  420  643]; %djmaus GUI position
pref.local =0; %1 for local communication (djmaus and open-ephys on same
%computer), 0 for remote (djmaus and open-ephys on different computers)
pref.stimuli=fullfile(pref.root, 'stimuli');
pref.mesoscope_mode=0; %0=use intan/openephys and save metadata to OE directory, 1=use arduino SCT and write metadata locally, no intan/openephys
if pref.local %same computer
    switch computer
        case 'MACI64'
            pref.zmqurl='tcp://localhost:5556'; %seems to work for mac
            pref.mexpath='mac';
        case 'GLNXA64'
            pref.zmqurl='tcp://127.0.0.1:5556'; %seems to work for linux
            pref.mexpath='unix';
        case 'PCWIN64'
            pref.zmqurl='tcp://127.0.0.1:5556'; %seems to work for windows
            pref.mexpath='windows';
    end
    pref.datapath=fullfile(pref.root, 'Data\lab'); %open-ephys data is acquired on this computer
    pref.remotedatapath=pref.datapath; %these are the same thing for local acquisition
        pref.datahost='\';
else %different computer
    switch computer
        case 'MACI64'
            pref.mexpath='mac';
        case 'GLNXA64'
            pref.mexpath='unix';
        case 'PCWIN64'
            pref.mexpath='windows';
    end
    %specific zmq url for the open-ephys computer
    pref.zmqurl='tcp://184.171.85.38:5556';
    pref.datahost='\\wehrrig1b';
    pref.remotedatapath='d:\lab\djmaus\Data\lab'; %what the open-ephys datapath looks like on the other computer
    pref.datapath='\\wehrrig1b\d\lab\djmaus\Data\lab'; %what the open-ephys datapath looks like to this computer
end

% individual user prefs:
if ~isfield(SP, 'user')
    SP.user='lab'; %default user
end

try
    switch SP.user

        case 'mike-home'
            %pref.stimuli='/home/mike/djmaus/stimuli/mike';
            pref.datapath='/Users/mikewehr/Documents/Data';
        case 'mike-linux'
            pref.stimuli='/home/mike/djmaus/stimuli/mike';
            pref.datapath='/home/Data';
        case 'mike-kombi'
            pref.stimuli='/Users/mikewehr/Documents/Analysis/djmaus-master/stimuli/mike';
            pref.datapath='/Users/mikewehr/Documents/Data2';
            pref.soundcarddeviceID=3;
    end
end



%the following sections are generated automatically, please don't edit below here








%saved mouseIDs for lab
switch SP.user
case  'lab'
	pref.allmouseIDs={'','Fred','PV3','pv3','test'};
end




















