function [filename,path]=MakeSquaredjProtocol(Amplitudes, Durations, Duties, isi, nrepeats)


% Usage: [filename,path]=MakeSquaredjProtocol(Amplitudes, Durations, Duties, isi, nrepeats;
%
%
% loosley based on MakeTonedjProtocol
%
% INPUTS:
% Amplitudes: vector of amplitudes in dB SPL 
%    (requires system to be calibrated) (can be single amplitude)
% Durations: vector of different square wave durations (in ms) (can be a single duration)
% Duties: vector of duty cycles (in percent positive) (can be a single value)
% isi: inter stimulus interval (onset-to-onset) in ms
% nrepeats: number of repetitions (different pseudorandom orders)
% OUTPUTS:
%       - creates a suitably named stimulus protocol in djprefs.stimuli\Tone Protocols
%       - returns name & path to protocol
% ------------------------------------------------------------------------
%
% example call: MakeSquaredjProtocol([20 40 60 80], [5 1], 50, 500, 10)
% gives 4 amplitudes, 2 durations, 1 duty and 10 repeats of these
%

if nargin==0; fprintf('\nno input');return;end

numamplitudes=length(Amplitudes);
numdurations=length(Durations);
numduties=length(Duties);

neworder=randperm( numamplitudes * numdurations * numduties);

amplitude=zeros(size(neworder));
dur=zeros(size(neworder));
duty=zeros(size(neworder));

StimPerRepeat=length(neworder);
TotalNumStim=StimPerRepeat*nrepeats;
DurationPerRepeatSecs=StimPerRepeat*(mean(Durations)+isi)/1000; %approx. duration per repeat
TotalDurationSecs=DurationPerRepeatSecs*nrepeats;

ampstring=sprintf('%d-', Amplitudes);ampstring=ampstring(1:end-1);
durstring=sprintf('%.2f-', Durations*1);durstring=durstring(1:end-1);    %%%%%%%%
dutystring=sprintf('%d-', Duties);dutystring=dutystring(1:end-1);

%put into stimuli structure
name= sprintf('Square-wave %da(%s),/%ddur(%s)/%dduty(%s)/%dmsisi/%d reps', ...
    numamplitudes,ampstring, numdurations, durstring, numduties,dutystring, isi, nrepeats);
description=sprintf('repeatingSquareWaves, %d amplitudes (%s dB SPL), %d durations (%sms),%d duties (%spercent), %d ms isi,%d stim per rep, %d repeats, %d total stimuli, %ds per rep, %d s total dur (%.1f min)',...
    numamplitudes,ampstring, numdurations, durstring, numduties, dutystring,...
    isi, StimPerRepeat, nrepeats,TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs), TotalDurationSecs/60);
filename=sprintf('repeatingSquareWaves_%da-%s_%ddur-%s_%dduty-%s_%dms-%dreps.mat',...
    numamplitudes,ampstring, numdurations, durstring, numduties, dutystring, isi, nrepeats);

%%% not well randomized here
nn=0;
for rep=1:nrepeats    
    neworder=randperm( numamplitudes);
    amplitudes = Amplitudes( neworder );
    for nAmp=1:length(amplitudes)
        neworder=randperm( numdurations);
        durs = Durations( neworder );
        for nDur=1:length(durs)
            neworder=randperm( numduties);
            duties = Duties( neworder );
            for nDuty=1:length(duties)
                nn=nn+1;
                stimuli(nn).type='square'; %use nn because stimuli(1) is name/description
                stimuli(nn).param.amplitude=amplitudes(nAmp);
                stimuli(nn).param.duration=durs(nDur);
                stimuli(nn).param.duty=duties(nDuty);
                stimuli(nn).param.next=isi;
                stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
                stimuli(nn).protocol_name=name;
                stimuli(nn).protocol_description=description;
                %%%% change this:
                stimuli(nn).PlottingFunction='PlotTC_PSTH';
                stimuli(nn).version='djmaus';
            end
        end
    end
end


global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('Tuning Curve protocols')
cd ('Tuning Curve protocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);

