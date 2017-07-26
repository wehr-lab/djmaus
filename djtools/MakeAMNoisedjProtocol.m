function [filename,path]=MakeAMNoisejProtocol(modulation_rates, modulation_depth, amplitude, duration, next, ramp, nrepeats)
%usage:  MakeAMNoisejProtocol(modulation_rates, modulation_depth, amplitude, duration, isi, ramp, nrepeats)
%
% inputs:
%   modulation_rates      -   vector of modulation rates (Hz)
%   modulation_depth      -   modulation depth (0-100%)
%   duration              -   overall stimulus duration (ms)
%   amplitude             -   amplitude (dB)
%   isi                   -   inter stimulus interval (ms)
%   ramp                  -   inter stimulus interval (ms)
%   nrepeats              -   number of repeats 
% outputs:
%       - creates a suitably named stimulus protocol in djprefs.stimuli/AMProtocols
%
%
%example calls:
%   modulation_rates =[1 2 4 8 16 32 64 128 256];  modulation_depth=100; duration =2000; amplitude=80; next=1000; ramp=10; nrepeats=20;
%   MakeAMNoisedjProtocol(modulation_rates, modulation_depth, amplitude, duration, next, ramp, nrepeats)

nummodulation_rates=length(modulation_rates);
modratestring=sprintf('%d-', modulation_rates);modratestring=modratestring(1:end-1);
dur_per_repeat=nummodulation_rates*(duration+next)/1000;
totaldur=dur_per_repeat*nrepeats;

name= sprintf('AMNoise/%sHz/%d%%fdepth/%ddB/dur%dms/isi%dms/%dreps', modratestring, modulation_depth, amplitude, duration, next, nrepeats);
description=sprintf('AM Noise, %sHz, %d%% depth, %d dB, %d ms dur, %d ms isi, %.1f ms ramp, %d repeats,  %d sec duration per repeat, %d sec total duration', modratestring, modulation_depth, amplitude, duration, next, ramp, nrepeats, dur_per_repeat, totaldur)
filename= sprintf('AMNoise_%sHz_%d%%depth_%ddB_dur%dms_isi%dms_%dreps', modratestring, modulation_depth, amplitude, duration,next, nrepeats);


nn=0;
for rep=1:nrepeats
    for n=randperm(nummodulation_rates)
        nn=nn+1;
        stimuli(nn).type='AMNoise';
        stimuli(nn).param.amplitude=amplitude;
        stimuli(nn).param.modulation_frequency=modulation_rates(n);
        stimuli(nn).param.modulation_depth=modulation_depth/100;
        stimuli(nn).param.next=next;
        stimuli(nn).param.ramp=ramp;
        stimuli(nn).param.duration=duration;
        stimuli(nn).param.laser=0;
        stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
        stimuli(nn).protocol_name=name;
        stimuli(nn).protocol_description=description;
        stimuli(nn).PlottingFunction='PlotAMNoise_PSTH';
        stimuli(nn).version='djmaus';
    end
end

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('AM Protocols')
cd ('AM Protocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);
