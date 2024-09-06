function [filename,path]=MakeClicktraindjProtocol(amplitude, trainduration, icis, clickduration, next, ramp, interleave_laser, nrepeats)
%usage:  MakeClicktraindjProtocol(amplitude, trainduration, icis, clickduration, next, ramp, interleave_laser, nrepeats)
%
%
% modified from MakeWNTrainProtocol2, and ported to djmaus
% we have a single amplitude and a fixed train duration
% (with variable number of clicks per train)
% creates an djmaus protocol file for a train of white noise bursts (click train)
%
% inputs:
%   amplitude:  amplitude in dB SPL (requires system to be calibrated)
%   trainduration   -   duration of train, in ms
%   icis            -   inter-click interval, i.e. interval between the
%                       start of previous click and start of the next click
%                       (use an array for multiple ICIs)
%   clickduration   -   duration of an individual click (ms)
%   next            -   inter-click-train-interval, i.e. when the next
%                       click train should follow the previous one (ms)
%   ramp            -   rising/falling edge of individual clicks (can be zero)
%   interleave_laser-   0 or 1 to duplicate all stimuli and interleave laser and non-laser trials in random order
%   nrepeats: number of repetitions (different pseudorandom orders)
%
% outputs:
%       - creates a suitably named stimulus protocol in djprefs.stimuli/ClicktrainProtocols
%
%
%example call:
% amplitude=80; trainduration=10e3; icis=[1 2 4 8 16 32 64 128 256];
% clickduration=.5; next=3000; ramp=0; interleave_laser=1; nrepeats=20;
% MakeClicktraindjProtocol(amplitude, trainduration, icis, clickduration, next, ramp, interleave_laser, nrepeats)



numicis=length(icis);
icistring=sprintf('%d-', icis);icistring=icistring(1:end-1);
if interleave_laser
    interleave_laserstr='-IL';
else
    interleave_laserstr='';
end

name= sprintf('Clicktrain, %ddB%.1fms/%sms/d=%d/%s', amplitude ,clickduration, icistring, trainduration,interleave_laserstr);
description=sprintf('Clicktrain, %d dB, %.1fms clickduration, %.1fms ramp, %d repeats, %d ms trainduration, %d ICIs (%sms), %ds duration per repeat, %s', amplitude, clickduration, ramp, nrepeats, trainduration, numicis, icistring, interleave_laserstr);
filename=sprintf('Clicktrain-%ddB-%.1fms-%sms-d%d%s.mat',amplitude,clickduration, icistring,trainduration,interleave_laserstr);

nn=0;
for rep=1:nrepeats
    for n=randperm(numicis)
        nn=nn+1;
        stimuli(nn).type='clicktrain';
        stimuli(nn).param.amplitude=amplitude;
        nclicks=ceil(trainduration/icis(n));
        stimuli(nn).param.nclicks=nclicks;
        stimuli(nn).param.ici=icis(n);
        stimuli(nn).param.clickduration=clickduration;
        stimuli(nn).param.next=next;
        stimuli(nn).param.ramp=ramp;
        stimuli(nn).param.duration=trainduration;
        stimuli(nn).param.laser=0;
        stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
        stimuli(nn).protocol_name=name;
        stimuli(nn).protocol_description=description;
        stimuli(nn).PlottingFunction='PlotClicktrain_PSTH';
        stimuli(nn).version='djmaus';
    end
end

if interleave_laser
    for rep=1:nrepeats
        for n=randperm(numicis)
            nn=nn+1;
            stimuli(nn).type='clicktrain';
            stimuli(nn).param.amplitude=amplitude;
            nclicks=ceil(trainduration/icis(n));
            stimuli(nn).param.nclicks=nclicks;
            stimuli(nn).param.ici=icis(n);
            stimuli(nn).param.clickduration=clickduration;
            stimuli(nn).param.next=next;
            stimuli(nn).param.ramp=ramp;
            stimuli(nn).param.duration=trainduration;
            stimuli(nn).param.laser=1;
            stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
            stimuli(nn).protocol_name=name;
            stimuli(nn).protocol_description=description;
            stimuli(nn).PlottingFunction='PlotClicktrain_PSTH';
            stimuli(nn).version='djmaus';
        end
    end

    %shuffle all stimuli
    rk=randperm(length(stimuli));
    for k=1:length(stimuli)
        newstim(k)=stimuli(rk(k));
    end
    stimuli=newstim;

end

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('ClicktrainProtocols')
cd ('ClicktrainProtocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);
