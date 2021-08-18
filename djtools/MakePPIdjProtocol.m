function [filename,path]=MakePPIdjProtocol(prepulsedurs, prepulseamps, pulsedur, pulseamp, soa, soaflag, ...
    ramp, iti, iti_var, interleave_laser, post_startle_duration, nrepeats)
% usage: MakePPIdjProtocol(prepulsedurs, prepulseamps, pulsedur, pulseamp, soa, soaflag, ...
%    ramp, iti, iti_var, interleave_laser, post_startle_duration, nrepeats)
%
%
%
% creates a djmaus stimulus protocol file for PPI (pre-pulse inhibition of acoustic startle
% response). Pre-pulse and startle pulse are both white noise.
% Can use multiple pre-pulse durations and amplitudes.
% You can only use a single startle pulse duration and amplitude.
%
% automatically includes a silent pre-pulse (-1000 dB) to provide a
% pure-startle condition for comparison. If you don't want this, you can
% set a flag below. If you request a -1000 dB or 0 ms prepulseamp it will not
% be duplicated.
%
%   mw 08.16.2021
%
%NOTE:
% inputs:
% prepulseamps: amplitudes of the pre-pulse, in dB, in a vector, e.g. 60, or [40 50 60]
% prepulsedurs: durations of the pre-pulse, in ms, in a vector, e.g. 50, or [0 50]
% pulsedur: duration of the startle pulse in ms (can be 0 for no startle)
% pulseamp: amplitude of the startle pulse in dB SPL
% soa: Stimulus Onset Asynchrony in ms = time between prepulse onset and
%       startle pulse onset
% soaflag: can be either 'soa' (default), in which case soa value specifies the time
% between the onset of the prepulse and the onset of the startle pulse, or else 'isi',
% in which case soa specifies the time between prepulse offset and startle
% onset. If anything other than 'isi' it will default to 'soa'.
% ramp: onset and offset ramp duration in ms, used for both prepulse and startle
% iti: inter trial interval (onset-to-onset) in ms
% iti_var: fractional variability of iti. Use 0 for fixed iti, or e.g. 0.1 to have iti vary by up to +-10%
% interleave_laser: 0 or 1 to duplicate all stimuli and interleave laser
%            and non-laser trials in random order
% post_startle_duration: duration of silence to play after the startle
%       stimulus has finished. We added this to allow extra time
%       for laser be on after the startle.
% nrepeats: number of repetitions (different pseudorandom orders)
%
% outputs:
% creates a suitably named stimulus protocol in djmaus/stimuli/PPI Protocols
% returns filename and path so you could programmatically create and use a protocol
%
%example calls:
% prepulsedurs=25; prepulseamps = 60; pulsedur=25; pulseamp=100; soa = 75;
% soaflag='soa'; ramp=0; iti=15000; iti_var=0; interleave_laser=0;
% nrepeats=10; post_startle_duration=0;
% MakePPIdjProtocol(prepulsedurs, prepulseamps, pulsedur, pulseamp, soa, soaflag, ramp, iti, iti_var, interleave_laser, post_startle_duration, nrepeats)

include_silent_prepulse=1; %set to 0 if you don't want to automatically include a silent prepulse condition.

if ismember(-1000, prepulseamps) | ismember(0, prepulsedurs)
    include_silent_prepulse=0;
    %no need to duplicate since there is already a silent prepulse
    %requested
end

if ~strcmp(soaflag, 'isi')
    soaflag='soa';
    fprintf('\nusing soa of %d ms', soa)
else
    fprintf('\nusing isi of %d ms', soa)
end

if strcmp(soaflag, 'soa')
    if any(prepulsedurs>soa)
        fprintf('\n\n!!!!!!!\n\n')
        warning('at least one prepulse duration exceeds the soa, so that prepulse duration will be invalid (will be interrupted by startle during the prepulse)')
    end
end


global pref
if isempty(pref) djPrefs;end
if nargin~=12 error('\MakePPIdjProtocol: wrong number of arguments.'); end

numprepulseamps=length(prepulseamps);
numprepulsedurs=length(prepulsedurs);

prepulsedursstring='';
for i=1:numprepulsedurs
    prepulsedursstring=[prepulsedursstring, sprintf('%g-', prepulsedurs(i))];
end
prepulsedursstring=prepulsedursstring(1:end-1); %remove trailing -

prepulseampsstring='';
for i=1:numprepulseamps
    prepulseampsstring=[prepulseampsstring, sprintf('%d-', prepulseamps(i))];
end
prepulseampsstring=prepulseampsstring(1:end-1); %remove trailing -


if interleave_laser==1
    [PrepulsedurGrid,PrepulseampGrid, Lasers]=meshgrid( prepulsedurs , prepulseamps, [0 1]);
    numlasers=2;
else
    [PrepulsedurGrid,PrepulseampGrid, Lasers]=meshgrid( prepulsedurs , prepulseamps, 0);
    numlasers=1;
end


neworder=randperm( numprepulseamps * numprepulsedurs  * numlasers);
rand_prepulsedurs=zeros(1, size(neworder, 2)*nrepeats);
rand_prepulseamps=zeros(1, size(neworder, 2)*nrepeats);
lasers=zeros(1, size(neworder, 2)*nrepeats);


for n=1:nrepeats
    neworder=randperm( numprepulseamps * numprepulsedurs * numlasers);
    rand_prepulsedurs( prod(size(PrepulsedurGrid))*(n-1) + (1:prod(size(PrepulsedurGrid))) ) = PrepulsedurGrid( neworder );
    rand_prepulseamps( prod(size(PrepulseampGrid))*(n-1) + (1:prod(size(PrepulseampGrid))) ) = PrepulseampGrid( neworder );
    lasers( prod(size(Lasers))*(n-1) + (1:prod(size(Lasers))) ) = Lasers( neworder );
end

if interleave_laser
    interleave_laserstr='IL';
else
    interleave_laserstr='';
end


name= sprintf('PPI-ppa%sdB-ppd%sms-sd%dms-sa%ddb-soa%dms(%s)-r%d-iti%d-itivar%d-%s-%dreps.mat',...
    prepulseampsstring, prepulsedursstring, pulsedur, pulseamp, soa,soaflag, round(ramp), iti,round(100*iti_var),interleave_laserstr, nrepeats);

description=''; %temp
filename=name;

n=0;


TotalDurationSecs=0;

if include_silent_prepulse 
% include a silent prepulse condition.
    for r=1:nrepeats
      n=n+1;
    stimuli(n).type='ASR';
    stimuli(n).param.ramp=ramp;
    stimuli(n).param.soa=soa;
    stimuli(n).param.soaflag=soaflag;
    stimuli(n).param.loop_flg=0;
    stimuli(n).param.seamless=0;
    PPI_duration= soa+pulsedur+post_startle_duration; %actual duration, note prepulse dur is zero here
    stimuli(n).param.duration=PPI_duration;
    this_iti=round(iti+iti*iti_var*(2*rand(1)-1));
    stimuli(n).param.next=this_iti;
    stimuli(n).param.pulsedur=pulsedur;
    stimuli(n).param.pulseamp=pulseamp;
    stimuli(n).param.prepulseamp=-1000;
    stimuli(n).param.prepulsedur=0;
    stimuli(n).param.laser=0;
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(n).PlottingFunction='PlotPPI_PSTH';
    stimuli(n).version='djmaus';
    
    TotalDurationSecs=TotalDurationSecs+(PPI_duration+this_iti)/1000;
    end
end

for kk=1:length(rand_prepulsedurs)
    
    switch soaflag
        case 'isi'
            PPI_duration=rand_prepulsedurs(kk) + soa+pulsedur+post_startle_duration; % duration
        case 'soa'
            PPI_duration= soa+pulsedur+post_startle_duration; %actual duration
    end
    n=n+1;
    stimuli(n).type='ASR';
    stimuli(n).param.ramp=ramp;
    stimuli(n).param.soa=soa;
    stimuli(n).param.soaflag=soaflag;
    stimuli(n).param.loop_flg=0;
    stimuli(n).param.seamless=0;
    stimuli(n).param.duration=PPI_duration;
    this_iti=round(iti+iti*iti_var*(2*rand(1)-1));
    stimuli(n).param.next=this_iti;
    stimuli(n).param.pulsedur=pulsedur;
    stimuli(n).param.pulseamp=pulseamp;
    stimuli(n).param.prepulseamp=rand_prepulseamps(kk);
    stimuli(n).param.prepulsedur=rand_prepulsedurs(kk);
    stimuli(n).param.laser=lasers(kk);
    stimuli(n).stimulus_description=GetParamStr(stimuli(n));
    stimuli(n).protocol_name=name;
    stimuli(n).protocol_description=description;
    stimuli(n).PlottingFunction='PlotPPI_PSTH';
    stimuli(n).version='djmaus';
    
    TotalDurationSecs=TotalDurationSecs+(PPI_duration+this_iti)/1000;
    
end

%add description with updated duration information
PPIPerRepeat=length(neworder);
DurationPerRepeatSecs=TotalDurationSecs/nrepeats;    %average duration per repeat
TotalNumStim=length(stimuli);
description=sprintf('PPI protocol, prepulsedurs: %s ms, prepulseamps: %s dB, startle pulse dur: %d ms, startle pulse amplitude: %d dB, SOA: %d ms (%s), ramp: %d ms, iti: %d ms, iti-var: %.1f, interleave laser: %d, post-startle dur %d ms, %d repeats, %d stim per rep, %d total stimuli, %ds per rep (avg), total dur %d s (%.1f min) ',...
    prepulsedursstring, prepulseampsstring, pulsedur, pulseamp, soa, soaflag, ramp, iti, iti_var, interleave_laser, post_startle_duration, nrepeats, PPIPerRepeat, TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs), TotalDurationSecs/60);
for n=1:TotalNumStim
    stimuli(n).protocol_description=description;
end

%shuffle so that the silent prepulse trials are interleaved
old_stim=stimuli;
shuf_order=randperm(TotalNumStim);
for n=1:TotalNumStim
    new_stimuli(n)=old_stim(shuf_order(n));
end
stimuli=new_stimuli;

cd(pref.stimuli) %where stimulus protocols are saved
warning off MATLAB:MKDIR:DirectoryExists
mkdir('PPI Protocols')
cd('PPI Protocols')
save(filename, 'stimuli')
path=pwd;
fprintf('\nwrote file %s \n in directory %s', filename, path)

