% script to make an ArchPV laser protocol for a reviewer at Neuron
% Nov 1 2016
% he suggests turning on the laser as the response is starting, so we can see if it takes some time to kick in
% I will just manually create this protocol, a one-off
% white noise, with variable laser onset

%0:2:50 relative to tone onset

amp=60;
dur=25;

name='ArchPVRev1';
description='custom protocol for reviewer#1 on Allie''s paper';
nn=0;

for rep=1:50
    
    nn=nn+1;
    %one laser OFF trial
    stimuli(nn).type='whitenoise'; %use nn because stimuli(1) is name/description
    stimuli(nn).param.amplitude=amp;
    stimuli(nn).param.duration=dur;
    stimuli(nn).param.laser=0;
    stimuli(nn).param.ramp=1;
    stimuli(nn).param.next=850;
    stimuli(nn).param.VarLaser=0;
    stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
    stimuli(nn).protocol_name=name;
    stimuli(nn).protocol_description=description;
    stimuli(nn).version='djmaus';
    
    for laserstart=-50:2:50
        nn=nn+1;
        stimuli(nn).type='whitenoise'; %use nn because stimuli(1) is name/description
        stimuli(nn).param.amplitude=amp;
        stimuli(nn).param.duration=dur;
        stimuli(nn).param.laser=1;
        stimuli(nn).param.ramp=1;
        stimuli(nn).param.next=850;
        stimuli(nn).param.VarLaser=1;
        stimuli(nn).param.VarLaserstart=laserstart; %relative to tone onset
        stimuli(nn).param.VarLaserpulsewidth=150; %for example
        stimuli(nn).param.VarLasernumpulses=1;
        stimuli(nn).param.VarLaserisi=0; %not used
        stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
        stimuli(nn).protocol_name=name;
        stimuli(nn).protocol_description=description;
        stimuli(nn).version='djmaus';
        
    end
end

stimuli=stimuli(randperm(nn));
global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli) %where stimulus protocols are saved
warning off MATLAB:MKDIR:DirectoryExists
mkdir('ArchPVRevProtocols')
cd('ArchPVRevProtocols')
save(name, 'stimuli')