function MakeGPIASdjProtocol(noiseamp, gapdurs, gapdelay, post_startle_duration, pulsedur, pulseamps, soa, soaflag, ...
    ramp, isi, isi_var, interleave_laser, nrepeats)
% usage MakeGPIASdjProtocol(noiseamp, gapdurs, gapdelay, post_startle_duration, 
%       pulsedur, pulseamps, soa, soaflag, ramp, iti, iti_var, interleave_laser, nrepeats)
%
%
%
% creates a djmaus stimulus protocol file for GPIAS (gap-induced pre-pulse inhibition of acoustic startle
% response). can use multiple gap durations, gap is silent
% using variable ITI.
%recent edits: 
%  -updated to djmaus version 9-2016
%  -added interleave_laser flag (0 or 1), so output can be already
%       laser-interleaved (no need for an extra ppalaser step) 9-2016
%  -changed to use whitenoise instead of band-passed noise. 9-2016
%  -added soaflag to specify whether soa is 'soa' or 'isi'
%  -changed gapdelay to specify time to gap offset instead of gap onset (so
%   that ppalaser comes on relative to gap offset in the 'isi' case) (note:
%   this is actually implemented in MakeGPIAS)
%   mw 06.09.2014
%
%NOTE: 
% inputs:
% noiseamp: amplitude of the continuous white noise, in dB SPL
% gapdurs: durations of the pre-pulse gap, in ms, in a vector, e.g. 50, or [0 50]
% gapdelay: delay from start of continuous noise to gap OFFSET, in ms
% post_startle_duration: duration of noise to play after the startle
%       stimulus has finished. We added this Oct 14, 2013 to allow extra time
%       for laser be on after the startle. 
% pulsedur: duration of the startle pulse in ms (can be 0 for no startle)
% pulseamps: amplitudes of the startle pulse in dB SPL, in a vector, e.g. 95, or [90 95 100]
% soa: Stimulus Onset Asynchrony in ms = time between gap onset and
%       startle pulse tone onset
% soaflag: can be either 'soa' (default), in which case soa value specifies the time
% between the onset of the gap and the onset of the startle, or else 'isi',
% in which case soa specifies the time between gap offset and startle
% onset. If anything other than 'isi' it will default to 'soa'.
% ramp: on-off ramp duration in ms
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
%noiseamp=80; gapdurs=[0 1 2 4 8 16 32 64 128 256]; gapdelay=1000; poststartle=1000;
%pulsedur=25; pulseamps=110; soa=50; soaflag='isi'; ramp=0; isi=15000; isi_var=.33; IL=1; nreps=10;
%MakeGPIASdjProtocol(noiseamp, gapdurs, gapdelay, poststartle, pulsedur, pulseamps, soa, soaflag, ramp, isi, isi_var, IL, nreps)

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

%if post_startle_duration==0 error('please use a finite post_startle_duration');end

global pref
if isempty(pref) djPrefs;end
if nargin~=13 error('\MakeGPIASdjProtocol: wrong number of arguments.'); end

numgapdurs=length(gapdurs);
numpulseamps=length(pulseamps);

gapdursstring='';
for i=1:numgapdurs
    gapdursstring=[gapdursstring, sprintf('%g-', gapdurs(i))];
end
gapdursstring=gapdursstring(1:end-1); %remove trailing -

pulseampsstring='';
for i=1:numpulseamps
    pulseampsstring=[pulseampsstring, sprintf('%d-', pulseamps(i))];
end
pulseampsstring=pulseampsstring(1:end-1); %remove trailing -

if interleave_laser==1
    [GapdurGrid,PulseampGrid, Lasers]=meshgrid( gapdurs , pulseamps, [0 1]);
    numlasers=2;
else
    [GapdurGrid,PulseampGrid, Lasers]=meshgrid( gapdurs , pulseamps, 0);
    numlasers=1;
end


neworder=randperm( numpulseamps * numgapdurs * numlasers);
rand_gapdurs=zeros(1, size(neworder, 2)*nrepeats);
rand_pulseamps=zeros(1, size(neworder, 2)*nrepeats);
lasers=zeros(1, size(neworder, 2)*nrepeats);


for n=1:nrepeats
    neworder=randperm( numpulseamps * numgapdurs * numlasers);
    rand_gapdurs( prod(size(GapdurGrid))*(n-1) + (1:prod(size(GapdurGrid))) ) = GapdurGrid( neworder );
    rand_pulseamps( prod(size(PulseampGrid))*(n-1) + (1:prod(size(PulseampGrid))) ) = PulseampGrid( neworder );
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

name= sprintf('GPIAS-na%ddB-gd%sms-pd%dms-pa%sdb-soa%dms(%s)-r%d-iti%d-itivar%d-%s-%dreps.mat',...
    noiseamp, gapdursstring, round(pulsedur), pulseampsstring, soa,soaflag, round(ramp), isi,round(100*isi_var),interleave_laserstr, nrepeats);

description=''; %temp
filename=name;

n=0;

for noisenum=1:num_noises

    n=n+1;

    stimuli(n).type='whitenoise';
    stimuli(n).param.amplitude=noiseamp;
    stimuli(n).param.ramp=ramp;
    stimuli(n).param.loop_flg=0;
    stimuli(n).param.seamless=1;
    %     stimuli(n).param.duration=500;
    stimuli(n).param.duration=gpias_duration;
    stimuli(n).param.next=next; %totally empirical value that allows psychportaudio rescheduling to work seamlessly
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(n).version='djmaus';
end

for kk=1:length(rand_gapdurs)

    n=n+1;
    stimuli(n).type='GPIAS';
    stimuli(n).param.amplitude=noiseamp;
    stimuli(n).param.ramp=ramp;
    stimuli(n).param.soa=soa;
    stimuli(n).param.soaflag=soaflag;
    stimuli(n).param.loop_flg=0;
    stimuli(n).param.seamless=1;
    stimuli(n).param.duration=gpias_duration;
    stimuli(n).param.next=next; %totally empirical value that allows psychportaudio rescheduling to work seamlessly
    stimuli(n).param.gapdelay=gapdelay;
    stimuli(n).param.gapdur=rand_gapdurs(kk);
    stimuli(n).param.pulsedur=pulsedur;
    stimuli(n).param.pulseamp=rand_pulseamps(kk);
    stimuli(n).param.laser=lasers(kk);
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(n).version='djmaus';
    
    %
    this_isi_ms=round(isi+isi*isi_var*(2*rand(1)-1));
    num_noises=round(this_isi_ms/gpias_duration);
    for noisenum=1:num_noises
        n=n+1;
        stimuli(n).type='whitenoise';
        stimuli(n).param.amplitude=noiseamp;
        stimuli(n).param.ramp=ramp;
        stimuli(n).param.loop_flg=0;
        stimuli(n).param.seamless=1;
        stimuli(n).param.duration=gpias_duration; %trying to set to same dur as gpias
        stimuli(n).param.next=next; %totally empirical value that allows psychportaudio rescheduling to work seamlessly
        stimuli(n).stimulus_description=GetParamStr(stimuli(n));
        stimuli(n).protocol_name=name;
        stimuli(n).protocol_description=description;
        stimuli(n).version='djmaus';
    end

end

TotalNumStim=length(stimuli);
TotalDurationSecs=TotalNumStim*gpias_duration/1000;
description=sprintf('GPIAS protocol: noise amp:%ddB, gapdur: %sms, gapdelay: %dms, pulsedur%dms pulse amplitude:%sdb SOA:%dms (%s) ramp:%dms iti:%dms iti-var: %.1f %s %drepeats, %d stim per rep (avg), %d total stimuli, %ds per rep (avg), %d s total dur',...
    noiseamp, gapdursstring, gapdelay, pulsedur, pulseampsstring, soa, soaflag, ramp, isi,round(100*isi_var),interleave_laserstr, nrepeats, StimPerRepeat, TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs));
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