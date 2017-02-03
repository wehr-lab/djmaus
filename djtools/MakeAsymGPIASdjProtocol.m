function MakeAsymGPIASdjProtocol(noiseamp, gapdurs, gapdelay, post_startle_duration, pulsedur, pulseamp, soa, soaflag, ...
    pulseramp,gapramps, isi, isi_var, interleave_laser, nrepeats)
% usage MakeAsymGPIASdjProtocol(noiseamp, gapdurs, gapdelay, post_startle_duration, 
%       pulsedur, pulseamp, soa, soaflag, ramp, iti, iti_var, interleave_laser, nrepeats)
%
%
%
% creates a djmaus stimulus protocol file for Asymmetrical GPIAS (gap-induced pre-pulse inhibition of acoustic startle
% response). can use multiple gap durations, gap is silent
% using variable ITI.
% What is an Asymmetric GPIAS?
% The gap ramps are specified individually. The offramp is the noise offset
% ramp (start of gap). The onramp is the noise onset (gap termination).
% Note that the conventional GPIAS stimuli always have gap ramp of zero.
% Note that the continuous background noise has ramps fixed at zero
% (which is required if they are to be continuous). 
% The ramps are centered on gap transitions. In other words, the background
% noise starts changing half-ramp before and stops changing half-ramp after
% each transition.
% This protocol is designed such that you provide a list of gap ramps (such
% as [0 5 10]), and the on and off ramps take each possible combination of
% these values. You could design it differently, all you need to do is
% specify what on & off ramps you want MakeAsymGPIAS to use.
%
%recent edits: 
%  -updated to djmaus version 9-2016
%  -added interleave_laser flag (0 or 1), so output can be already
%       laser-interleaved (no need for an extra ppalaser step) 9-2016
%  -changed to use whitenoise instead of band-passed noise. 9-2016
%  -added soaflag to specify whether soa is 'soa' or 'isi'
%  -changed gapdelay to specify time to gap offset instead of gap onset (so
%   that ppalaser comes on relative to gap offset in the 'isi' case) (note:
%   this is actually implemented in MakeGPIAS)
%  -introduced asymmetric ramps
%  -only one pulseamp allowed in this protocol builder
%   mw 10.11.2014
%
%NOTE: 
% inputs:
% gapramps: list of ramps for gap onset and termination, in ms
% noiseamp: amplitude of the continuous white noise, in dB SPL
% gapdurs: durations of the pre-pulse gap, in ms, in a vector, e.g. 50, or [0 50]
% gapdelay: delay from start of continuous noise to gap OFFSET, in ms
% post_startle_duration: duration of noise to play after the startle
%       stimulus has finished. We added this Oct 14, 2013 to allow extra time
%       for laser be on after the startle. 
% pulsedur: duration of the startle pulse in ms (can be 0 for no startle)
% pulseamp: amplitude of the startle pulse in dB SPL, only one value allowed here
% soa: Stimulus Onset Asynchrony in ms = time between gap onset and
%       startle pulse tone onset
% soaflag: can be either 'soa' (default), in which case soa value specifies the time
% between the onset of the gap and the onset of the startle, or else 'isi',
% in which case soa specifies the time between gap offset and startle
% onset. If anything other than 'isi' it will default to 'soa'.
% pulseramp: startle pulse on-off ramp duration in ms
% iti: inter trial interval (onset-to-onset) in ms
% iti_var: fractional variability of iti. Use 0 for fixed iti, or e.g. 0.1 to have iti vary by up to +-10%
% interleave_laser: 0 or 1 to duplicate all stimuli and interleave laser
%            and non-laser trials in random order
% nrepeats: number of repetitions (different pseudorandom orders)
%
% note: still using the variable isi for inter-trial interval, AKA iti
% outputs:
% creates a suitably named stimulus protocol in D:\lab\exper2.2\protocols\ASR Protocols
%
% Note that the ITI between GPIAS stimuli cannot be less than
% gapdelay+poststartle, and can only vary by increments of this duration.
%
%example calls:
% fixed iti of 10 seconds with interleaved laser:
%MakeGPIASdjProtocol(80, [0 2 4 6], 1000, 1000, 25, 100, 60, 'soa', 0, 10e3, 0, 1, 5)
%
% iti ranging from 10s to 20s (15 s on average)
%
%brief variable duration gaps, 60ms SOA
%MakeGPIASdjProtocol(80, [0 2 4 6], 1000, 1000, 25, 100, 60, 'soa', 0, 15e3, .33, 1, 15)
%
%brief gap, no startle, ability to deliver a long (1sec) laser pulse beyond
%startle offset time
%MakeGPIASdjProtocol(80, [10], 1000, 1000, 0, 100, 60, 'soa', 0, 15e3, .33, 1, 20)
%
%MakeGPIASdjProtocol(80, [0 1 2 4 8 16 32 64 128 256], 1000, 1000, 0, 100, 50, 'isi', 0, 1e3, 0, 0, 15)
%
%noiseamp=80; gapdurs=[16 128 ]; gapdelay=1000; poststartle=1000; gapramps=[0 4 8]
%pulsedur=0; pulseamp=0; soa=50; soaflag='isi'; pulseramp=0; isi=1000; isi_var=0; IL=0; nreps=50;
%MakeAsymGPIASdjProtocol(noiseamp, gapdurs, gapdelay, poststartle, pulsedur, pulseamp, soa, soaflag, pulseramp, gapramps, isi, isi_var, IL, nreps)

if ~strcmp(soaflag, 'isi')
    soaflag='soa';
    fprintf('\nusing soa of %d ms', soa)
else
    fprintf('\nusing isi of %d ms', soa)
end

if strcmp(soaflag, 'soa')
    if any(gapdurs>soa)
        fprintf('\n\n!!!!!!!\n\n')
        warning('at least one gap duration exceeds the soa, so that gap duration will be invalid (will be interrupted by startle during the gap)')
    end
end

if any(gapramps>min(gapdurs))
    warning('at least one gap ramp exceeds a gap duration... that won''t work')
end
if length(pulseamp)>1 error('only one pulseamp allowed in this protocol builder');end

%if post_startle_duration==0 error('please use a finite post_startle_duration');end

global pref
if isempty(pref) djPrefs;end

numgapdurs=length(gapdurs);
numgapramps=length(gapramps);
numgaprampcombos=numgapramps*factorial(numgapramps); %number of combinations of on and off ramps of each value

gapdursstring='';
for i=1:numgapdurs
    gapdursstring=[gapdursstring, sprintf('%g-', gapdurs(i))];
end
gapdursstring=gapdursstring(1:end-1); %remove trailing -

gaprampstring='';
for i=1:numgapramps
    gaprampstring=[gaprampstring, sprintf('%d-', gapramps(i))];
end
gaprampstring=gaprampstring(1:end-1); %remove trailing -

if interleave_laser==1
    [GapdurGrid,GapOnrampGrid, GapOfframpGrid, Lasers]=ndgrid( gapdurs , gapramps, gapramps, [0 1]);
    numlasers=2;
else
    [GapdurGrid,GapOnrampGrid, GapOfframpGrid, Lasers]=ndgrid( gapdurs ,gapramps, gapramps, 0);
    numlasers=1;
end


neworder=randperm( numgapdurs * numgapramps * numgapramps * numlasers);
rand_gapdurs=zeros(1, size(neworder, 2)*nrepeats);
rand_gapramps=zeros(1, size(neworder, 2)*nrepeats);
lasers=zeros(1, size(neworder, 2)*nrepeats);


for n=1:nrepeats
    neworder=randperm( numgapdurs * numgapramps * numgapramps * numlasers);
    rand_gapdurs( prod(size(GapdurGrid))*(n-1) + (1:prod(size(GapdurGrid))) ) = GapdurGrid( neworder );
    rand_gaponramps( prod(size(GapOnrampGrid))*(n-1) + (1:prod(size(GapOnrampGrid))) ) = GapOnrampGrid( neworder );
    rand_gapofframps( prod(size(GapOfframpGrid))*(n-1) + (1:prod(size(GapOfframpGrid))) ) = GapOfframpGrid( neworder );
    lasers( prod(size(Lasers))*(n-1) + (1:prod(size(Lasers))) ) = Lasers( neworder );
end

if interleave_laser
    interleave_laserstr='IL';
else
    interleave_laserstr='';
end

gpias_duration=gapdelay+max(rand_gapdurs)+soa+pulsedur+post_startle_duration; %actual duration

%note: for seamless playing of sounds, all buffers must be identical in
%length. So we are making short noise segments and using variable numbers
%of them

%next=-2000;%totally empirical value that allows psychportaudio rescheduling to work seamlessly
%was -1000, trying new values to get it working on Rig 2
next = -gapdelay/2;%testing mw 032410
%next = -.9*gapdelay;%testing mw 06.11.2014

this_isi_ms=round(isi+isi*isi_var*(2*rand(1)-1));
num_noises=round(this_isi_ms/gpias_duration);

GPIASPerRepeat=length(neworder);
StimPerRepeat=round(GPIASPerRepeat*isi/gpias_duration); %on average
DurationPerRepeatSecs=StimPerRepeat*(gpias_duration)/1000;%approx. duration per repeat

name= sprintf('AsymGPIAS-na%ddB-gd%sms-gr%sms-pd%dms-pa%ddb-pr%d-soa%dms(%s)-iti%d-itivar%d-%s-%dreps.mat',...
    noiseamp, gapdursstring, gaprampstring, round(pulsedur),pulseamp, round(pulseramp),soa,soaflag, isi,round(100*isi_var),interleave_laserstr, nrepeats);

description=''; %temp
filename=name;

n=0;

for noisenum=1:num_noises

    n=n+1;

    stimuli(n).type='whitenoise';
    stimuli(n).param.amplitude=noiseamp;
    stimuli(n).param.ramp=0;
    stimuli(n).param.loop_flg=0;
    stimuli(n).param.seamless=1;
    %     stimuli(n).param.duration=500;
    stimuli(n).param.duration=gpias_duration;
    stimuli(n).param.next=next; %totally empirical value that allows psychportaudio rescheduling to work seamlessly
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(nn).PlottingFunction='PlotAsymGPIAS_PSTH';
    stimuli(n).version='djmaus';
end

for kk=1:length(rand_gapdurs)

    n=n+1;
    stimuli(n).type='AsymGPIAS';
    stimuli(n).param.amplitude=noiseamp;
    stimuli(n).param.pulseramp=pulseramp;
    stimuli(n).param.soa=soa;
    stimuli(n).param.soaflag=soaflag;
    stimuli(n).param.loop_flg=0;
    stimuli(n).param.seamless=1;
    stimuli(n).param.duration=gpias_duration;
    stimuli(n).param.next=next; %totally empirical value that allows psychportaudio rescheduling to work seamlessly
    stimuli(n).param.gapdelay=gapdelay;
    stimuli(n).param.gapdur=rand_gapdurs(kk);
    stimuli(n).param.onramp=rand_gaponramps(kk);
    stimuli(n).param.offramp=rand_gapofframps(kk);
    stimuli(n).param.pulsedur=pulsedur;
    stimuli(n).param.pulseamp=pulseamp;
    stimuli(n).param.laser=lasers(kk);
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(nn).PlottingFunction='PlotAsymGPIAS_PSTH';
    stimuli(n).version='djmaus';
    
    %
    this_isi_ms=round(isi+isi*isi_var*(2*rand(1)-1));
    num_noises=round(this_isi_ms/gpias_duration);
    for noisenum=1:num_noises
        n=n+1;
        stimuli(n).type='whitenoise';
        stimuli(n).param.amplitude=noiseamp;
        stimuli(n).param.ramp=0;
        stimuli(n).param.loop_flg=0;
        stimuli(n).param.seamless=1;
        stimuli(n).param.duration=gpias_duration; %trying to set to same dur as gpias
        stimuli(n).param.next=next; %totally empirical value that allows psychportaudio rescheduling to work seamlessly
        stimuli(n).stimulus_description=GetParamStr(stimuli(n));
        stimuli(n).protocol_name=name;
        stimuli(n).protocol_description=description;
        stimuli(nn).PlottingFunction='PlotAsymGPIAS_PSTH';
        stimuli(n).version='djmaus';
    end

end

TotalNumStim=length(stimuli);
TotalDurationSecs=TotalNumStim*gpias_duration/1000;
description=sprintf('Asymmetric GPIAS protocol: noise amp:%ddB, gapdur: %sms, gapramps: %s ms, gapdelay: %dms, pulsedur%dms pulse amplitude:%ddb pulseramp:%dms SOA:%dms (%s) iti:%dms iti-var: %.1f %s %drepeats, %d stim per rep (avg), %d total stimuli, %ds per rep (avg), %d s total dur',...
    noiseamp, gapdursstring, gaprampstring, gapdelay, pulsedur, pulseamp, pulseramp, soa, soaflag, isi,round(100*isi_var),interleave_laserstr, nrepeats, StimPerRepeat, TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs));
for n=1:TotalNumStim
    stimuli(n).protocol_description=description;
end

cd(pref.stimuli) %where stimulus protocols are saved
warning off MATLAB:MKDIR:DirectoryExists
mkdir('GPIAS Protocols')
cd('GPIAS Protocols')
save(filename, 'stimuli')

fprintf('\nwrote file %s \n in directory %s', filename, pwd)

%  keyboard