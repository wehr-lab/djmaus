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
pref.reqlatencyclass=2;
pref.SoundFs=192000;
pref.maxSPL=80;
pref.allmouseIDs='';
pref.root=fileparts(which(mfilename));
pref.stimuli=fullfile(pref.root, 'stimuli');
pref.stimuli=fullfile(pref.root, 'stimuli');
pref.datapath=fullfile(pref.root, 'Data');
pref.OEpath='D:\lab\plugin-GUI\Builds\VisualStudio2013\x64\Release64\bin\open-ephys.exe &'; %make sure to put & at the end

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




















