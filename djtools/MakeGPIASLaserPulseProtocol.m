function MakeGPIASLaserPulseProtocol(Laserstarts, Laserpulsewidth, Lasernumpulses, Laserisi, noiseamp, gapdurs, gapdelay, post_startle_duration, pulsedur, pulseamp, soa, soaflag, ...
    ramp, isi, isi_var, nrepeats)


% usage  MakeGPIASLaserPulseProtocol(Laserstart, Laserpulsewidth, Lasernumpulses, Laserisi, noiseamp, gapdurs, gapdelay, post_startle_duration, pulsedur, pulseamp, soa, soaflag, ...
%     ramp, isi, isi_var, nrepeats)
%
%
% creates a djmaus stimulus protocol file for GPIAS (gap-induced pre-pulse inhibition of acoustic startle
% response). can use multiple gap durations, gap is silent, can use variable ITI.
%
%this version randomly interleaves trials with a laser pulse for each GPIAS trial,
% the parameters are set by the first arguments (Laserstarts, Laserpulsewidth, Lasernumpulses, Laserisi)
%Laserstarts can be an array, the other laser params should be single values
%
%
% example call
%  Laserstarts=[0 25]; Laserpulsewidth=5; Lasernumpulses=1; Laserisi=0; 
%  noiseamp=80; gapdurs=[0 32]; gapdelay=500; poststartle=500;
%  pulsedur=25; pulseamp=100; soa=50; soaflag='isi'; ramp=0; isi=12000; isi_var=0; nrepeats=20;
% MakeGPIASLaserPulseProtocol(Laserstarts, Laserpulsewidth, Lasernumpulses, Laserisi, noiseamp, gapdurs, gapdelay, poststartle, pulsedur, pulseamp, soa, soaflag, ramp, isi, isi_var, nrepeats)

if any(Laserpulsewidth==0)
    fprintf('\nremoving requested laserpulsewidth of 0 because this doesn''t make sense and is reserved for the laser-off condition');
    Laserpulsewidth=Laserpulsewidth(find(Laserpulsewidth));
end

%recent edits:
%  -updated to djmaus version 9-2016
%  -uses VarLaser params in stimuli to set laser params and ignore djmaus GUI, 9-2016
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
% outputs:
% creates a suitably named stimulus protocol in D:\lab\exper2.2\protocols\ASR Protocols
%
%example calls:
% fixed iti of 10 seconds:
%MakeVarGPIASdjProtocol(80, [0 2 4 6], 1000, 1000, 25, 100, 60, 'soa', 0, 10e3, 0, 1, 5)
%
% iti ranging from 10s to 20s (15 s on average)
%
%brief variable duration gaps, 60ms SOA
%MakeVarGPIASdjProtocol(80, [0 2 4 6], 1000, 1000, 25, 100, 60, 'soa', 0, 15e3, .33, 1, 15)
%
%brief gap, no startle, ability to deliver a long (1sec) laser pulse beyond
%startle offset time
%MakeVarGPIASdjProtocol(80, [10], 1000, 1000, 0, 100, 60, 'soa', 0, 15e3, .33, 1, 20)

%note: still using the variable isi for inter-trial interval, AKA iti

interleave_laser=1;
%fixing it to 1, otherwise you wouldn't eb using this function, you would
%use MakeGPIASdjProtocol instead

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
if nargin~=16 error('MakeGPIASLaserPulseProtocol: wrong number of arguments.'); end

numgapdurs=length(gapdurs);
numpulseamps=1;
numLaserstarts=length(Laserstarts);

gapdursstring='';
for i=1:numgapdurs
    gapdursstring=[gapdursstring, sprintf('%g-', gapdurs(i))];
end
gapdursstring=gapdursstring(1:end-1); %remove trailing -

Laserstartstring='';
for i=1:numLaserstarts
    Laserstartstring=[Laserstartstring, sprintf('%g-', Laserstarts(i))];
end
Laserstartstring=Laserstartstring(1:end-1); %remove trailing -

%original version did a straight interleave, i.e. half of trials had no
%laser, but instead we want only as many laser-off trials as one of the
%laserstarts
%so I'm removing the interleave and just adding in laser-off as a fake
%laserstart of -666

if interleave_laser==1
%     [GapdurGrid,VarLaserpulsewidthGrid, Lasers]=meshgrid( gapdurs , VarLaserpulsewidths, [0 1]);
    [GapdurGrid,LaserstartGrid]=meshgrid( gapdurs , [-666 Laserstarts]);
    numlasers=1;
else %cannot happen
%     [GapdurGrid,PulseampGrid, Lasers]=meshgrid( gapdurs , pulseamps, 0);
    numlasers=1;
end


neworder=randperm( (1+numLaserstarts) * numgapdurs * numlasers);
rand_gapdurs=zeros(1, size(neworder, 2)*nrepeats);
rand_Laserstarts=zeros(1, size(neworder, 2)*nrepeats);
lasers=zeros(1, size(neworder, 2)*nrepeats);


for n=1:nrepeats
    neworder=randperm( (1+numLaserstarts) * numgapdurs * numlasers);
    rand_gapdurs( prod(size(GapdurGrid))*(n-1) + (1:prod(size(GapdurGrid))) ) = GapdurGrid( neworder );
    rand_Laserstarts( prod(size(LaserstartGrid))*(n-1) + (1:prod(size(LaserstartGrid))) ) = LaserstartGrid( neworder );
   % lasers( prod(size(Lasers))*(n-1) + (1:prod(size(Lasers))) ) = Lasers( neworder );
end

if interleave_laser
    interleave_laserstr='IL-';
else
    interleave_laserstr='';
end

name= sprintf(['GPIAS-varlaserpulse-Lstarts%sms-Lpw%.1fms-%dp-isi%d-na%ddB', ...
    '-gd%sms-pd%dms-pa%ddb-soa%dms(%s)-r%d-iti%d-itivar%d-%s%dreps.mat'], ...
    Laserstartstring,Laserpulsewidth, Lasernumpulses, Laserisi, noiseamp, gapdursstring, round(pulsedur), pulseamp, soa,soaflag, round(ramp), isi,round(100*isi_var),interleave_laserstr, nrepeats);

description=sprintf(['GPIAS protocol with varlaserpulse, laserstarts: %s ms, ', ...
    'laserpulsewidth: %.1f ms, %d laser pulses, %.1f laser isi, noise amp:%ddB, gap duration: %sms, gapdelay: %dms, ',...
    'pulse duration: %dms pulse amplitude: %ddb SOA:%dms (%s) ramp: %dms iti: %dms iti-var: %.1f %s %drepeats'],...
    Laserstartstring,Laserpulsewidth, Lasernumpulses, Laserisi, noiseamp, gapdursstring, gapdelay, pulsedur, pulseamp, soa, soaflag, ramp, isi,round(100*isi_var),interleave_laserstr, nrepeats);
filename=name;


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
    stimuli(n).protocol_name=name;
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_description=description;
    stimuli(n).PlottingFunction='PlotGPIASflashtrain_PSTH';
    stimuli(n).version='djmaus';
    stimuli(n).param.VarLaser=0;
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
    stimuli(n).param.pulseamp=pulseamp;
    if rand_Laserstarts(kk)==-666 %using laserstart of -666 as the laser-off flag
        stimuli(n).param.laser=0;
    else
        stimuli(n).param.laser=1;
    end
    if stimuli(n).param.laser laserstr='Var laser ON'; else laserstr='laser OFF';end
    if stimuli(n).param.laser
        stimuli(n).param.VarLaser=1;
        %from hard-coded values at top
        stimuli(n).param.VarLaserstart=gapdelay+rand_Laserstarts(kk); %for example.  laser onset coincides w/ gap onset
        stimuli(n).param.VarLaserpulsewidth=Laserpulsewidth; 
        stimuli(n).param.VarLasernumpulses=Lasernumpulses;
        stimuli(n).param.VarLaserisi=Laserisi; %
    else
        stimuli(n).param.VarLaser=0;
    end
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(n).PlottingFunction='PlotGPIASLaserPulse_PSTH';
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
        stimuli(n).PlottingFunction='PlotGPIASflashtrain_PSTH';
        stimuli(n).param.VarLaser=0;
    end
    
end

TotalNumStim=length(stimuli);
TotalDurationSecs=TotalNumStim*gpias_duration/1000;
description=sprintf('%s total dur %d s (%d min)', description, round(TotalDurationSecs), round(TotalDurationSecs/60));
for n=1:TotalNumStim
    stimuli(n).protocol_description=description;
end


cd(pref.stimuli) %where stimulus protocols are saved
warning off MATLAB:MKDIR:DirectoryExists
mkdir('GPIAS Protocols')
cd('GPIAS Protocols')
save(filename, 'stimuli')

fprintf('\nwrote file %s \n in directory %s', filename, pwd)
