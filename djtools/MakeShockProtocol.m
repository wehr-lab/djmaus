function [filename,path]=MakeShockProtocol






name= sprintf('Shock1');
description=sprintf('Shock1');
filename=sprintf('Shock1');


nn=0;
for rep=1:10
    
    
    nn=nn+1;
    stimuli(nn).type='whitenoise'; %use nn because stimuli(1) is name/description
    stimuli(nn).param.amplitude=80;
    stimuli(nn).param.duration=25;
    stimuli(nn).param.laser=0;
    stimuli(nn).param.ramp=0;
    stimuli(nn).param.next=1000;
    
    stimuli(nn).param.Shock=1; %whether to deliver a shock on this stimulus
    stimuli(nn).param.Shockstart=-25; %ms relative to sound onset
    stimuli(nn).param.Shockpulsewidth=1; %ms
    stimuli(nn).param.Shocknumpulses=20; % for shock trains
    stimuli(nn).param.Shockisi=5; %ms, for shock trains
    %note: the soundcard channel to deliver shock TTL is specified in djPrefs as pref.Shockchannel
    
    
    stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
    stimuli(nn).protocol_name=name;
    stimuli(nn).protocol_description=description;
    stimuli(nn).PlottingFunction='none';
    stimuli(nn).version='djmaus';
    
end


global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('Shock protocols')
cd ('Shock protocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);

