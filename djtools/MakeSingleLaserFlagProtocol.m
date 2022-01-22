function MakeSingleLaserFlagProtocol
%simple function to creat a protocol with a single laser flag
%the idea is that you use djmaus front panel to specifiy laser params, and
%the play button to control timing
n=1;

stimuli(n).type='silentsound';
stimuli(n).param.loop_flg=0;
stimuli(n).param.duration=100;
stimuli(n).param.ramp=0;
stimuli(n).param.next=1000;
stimuli(n).protocol_name='SingleLaserFlag';
stimuli(n).stimulus_description=GetParamStr(stimuli(n));
stimuli(n).protocol_description='single laser ON flag for front-panel controlled laser stim';
stimuli(n).PlottingFunction='none';
stimuli(n).version='djmaus';
stimuli(n).param.laser=1;

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli) %where stimulus protocols are saved
warning off MATLAB:MKDIR:DirectoryExists
mkdir('Prey Capture Laser protocols')
cd('Prey Capture Laser protocols')
filename='SingleLaserFlag';
save(filename, 'stimuli')

fprintf('\nwrote file %s \n in directory %s', filename, pwd)
