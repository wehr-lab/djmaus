function [filename,path]=MakeF1F2GapProtocol2(F1, F2, amplitude, ...
    include_whitenoise, duration, ramp, iti, nrepeats, gapdurs, include_silent_sound)
% same as MakeF1F2GapProtocol, but only presents the single F1-F2
% combination. This will allow efficient collection of many gap durs, but
% must be run in conjunction with a normal F1F2GapProtocol so that we get
% all the F1 and F2 combinations at 2 key gapdurs
%
%usage: MakeF1F2GapProtocol2(f1, f2, amplitude, ...
%    include_whitenoise, duration, ramp, iti, nrepeats, gapdurs, include_silent_sound)
% mw 3.6.2017
%
% creates a djmaus stimulus protocol file for a 2 tone gap stimulus
% with an F1 frequency tone followed by a gap and then a F2 frequency tone
% the combinations of F1-F1, F2-F2, F1-F2, and F2-F1 are used
% this is designed for a specific experiment in which F1 evokes ON response
% only, and F2 evokes OFF response only
%
% inputs:
% F1, F2: tone frequencies in Hz
% amplitude: we used one fixed amplitude in dB SPL
% include_whitenoise: 1 for yes and 0 for no
% duration: tone duration in ms
% ramp: on-off ramp duration in ms
% iti: inter trial interval in ms
% nrepeats: number of repetitions (different pseudorandom orders)
% gapdurs: the duration of gaps between F1 and F2
% include_silent_sound: 0 or 1 to include an extra stimulus that is silent,
%           useful for explicitly collecting spontaneous and/or laser-only
%           trials as a separate condition
%
%
% outputs:
% creates a suitably named stimulus protocol in
% home/lab/djmaus/stimuli/2Tone
%
%
%example calls:
%500 ms tones with 16 and 128 ms gaps, no WN or silent sound
% F1=4000; F2=10000; amplitude=70; include_whitenoise=0; duration=500; ramp=2; iti=1000; nrepeats=20;
% gapdurs=[2 4 8 16 32 64 128 256]; include_silent_sound=0;
% MakeF1F2GapProtocol2(F1, F2, amplitude, include_whitenoise, duration, ramp, iti, nrepeats, gapdurs, include_silent_sound)



if nargin==0; fprintf('\nno input');return;end
if nargin>10; error('\MakeF1F2GapProtocol: wrong number of arguments.'); end

numgapdurs=length(gapdurs);
StimPerRepeat=numgapdurs*(1+ include_whitenoise)+include_silent_sound; %


TotalNumStim=StimPerRepeat*nrepeats;
DurationPerRepeatSecs=StimPerRepeat*(2*duration+mean(gapdurs)+iti)/1000;%approx. duration per repeat
TotalDurationSecs=DurationPerRepeatSecs*nrepeats;


% filename and description


if include_silent_sound
    include_silent_soundstr='ISS';
else
    include_silent_soundstr='';
end
if include_whitenoise
    include_whitenoisestr='+WN';
else
    include_whitenoisestr='';
end
gapdursstring='';
for i=1:numgapdurs
    gapdursstring=[gapdursstring, sprintf('%g-', gapdurs(i))];
end
gapdursstring=gapdursstring(1:end-1); %remove trailing -



name= sprintf('F1F2Gap2 %s/%d & %d Hz/%d dB/tone dur %d ms/%s ms gaps/%s/%d ms iti/%d reps', ...
    include_whitenoisestr, F1,F2, amplitude,duration, gapdursstring, ...
    include_silent_soundstr, iti,nrepeats);
description= sprintf('F1F2Gap2 %s/%d & %d Hz/%d dB/tone dur %d ms/%s ms gaps/%s/%d ms iti/%d reps ,%d stim per rep, %d total stimuli, %ds per rep, %d s total dur', ...
    include_whitenoisestr, F1,F2, amplitude,duration, gapdursstring, ...
    include_silent_soundstr, iti,nrepeats, ...
    StimPerRepeat,TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs));

filename= sprintf('F1F2Gap2-%s-%d-%dHz%ddB-tonedur%dms-%s-msgaps-%s-%dmsiti-%dreps', ...
    include_whitenoisestr, F1,F2, amplitude,duration, gapdursstring, ...
    include_silent_soundstr, iti,nrepeats);

n=0;
reps=0:StimPerRepeat:(StimPerRepeat*nrepeats);
for rep=1:nrepeats
        for g=1:numgapdurs
        

        
   
        
         %F1-F2
        n=n+1;
        stimuli(n).type='2tone';
        stimuli(n).param.frequency=F1;
        stimuli(n).param.amplitude=amplitude;
        stimuli(n).param.probeamp=amplitude;
        stimuli(n).param.probefreq=F2;
        stimuli(n).param.duration=duration;
        stimuli(n).param.SOA=duration + gapdurs(g);
        stimuli(n).param.ramp=ramp;
        stimuli(n).param.next=iti;
        stimuli(n).param.laser=0;
        stimuli(n).stimulus_description=GetParamStr(stimuli(n));
        stimuli(n).protocol_name=name;
        stimuli(n).protocol_description=description;
        stimuli(n).PlottingFunction='Plot2Tone_PSTH';
        stimuli(n).version='djmaus';
        
        if include_whitenoise
            %WN-WN
            n=n+1;
            stimuli(n).type='2tone';
            stimuli(n).param.frequency=-1;
            stimuli(n).param.amplitude=amplitude;
            stimuli(n).param.probeamp=amplitude;
            stimuli(n).param.probefreq=-1;
            stimuli(n).param.duration=duration;
            stimuli(n).param.SOA=duration + gapdurs(g);
            stimuli(n).param.ramp=ramp;
            stimuli(n).param.next=iti;
            stimuli(n).param.laser=0;
            stimuli(n).stimulus_description=GetParamStr(stimuli(n));
            stimuli(n).protocol_name=name;
            stimuli(n).protocol_description=description;
            stimuli(n).PlottingFunction='Plot2Tone_PSTH';
            stimuli(n).version='djmaus';
        end
        
        % insert silent sounds
        if include_silent_sound
            n=n+1;
            stimuli(n).param.laser=0;
            stimuli(n).type='silentsound';
            stimuli(n).param.duration=duration;
            stimuli(n).param.ramp=0;
            stimuli(n).param.next=iti;
            stimuli(n).stimulus_description=GetParamStr(stimuli(n));
            stimuli(n).protocol_name=name;
            stimuli(n).stimulus_description=description;
            stimuli(n).PlottingFunction='Plot2Tone_PSTH';
            stimuli(n).version='djmaus';
        end
    end
end

%shuffle stimuli

stimorder=randperm(TotalNumStim); %order of single tone and 2Tone, random
for n=1:TotalNumStim
    shuffledstimuli(n)=stimuli(stimorder(n));
end
stimuli=shuffledstimuli;

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('2Tone')
cd ('2Tone')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);


