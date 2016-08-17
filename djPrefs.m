function djPrefs
%Local preferences
%specifics for the local machine are at top
%preferences for individual users at the bottom
global pref SP

% machine-wide default prefs:

pref.users={'lab', 'mike-home', 'mike-kombi', 'mike-linux'};
pref.soundcarddeviceID=3;
pref.num_soundcard_outputchannels=2;
pref.reqlatencyclass=2;
pref.SoundFs=44100;
pref.maxSPL=80;
pref.allmouseIDs='';
pref.root=fileparts(which(mfilename));
pref.stimuli=fullfile(pref.root, 'stimuli');
pref.stimuli=fullfile(pref.root, 'stimuli');
pref.datapath=fullfile(pref.root, 'Data');

%hack for exper compatibility:
pref.home='/Users/mikewehr/Documents/Analysis/';

% individual user prefs:

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




%the following sections are generated automatically, please don't edit below here







