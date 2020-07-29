function ProcessTC_PupilMotion(varargin)
% processes continous traces from pupil and ball motion, similarly to LFP
% trace.
%make sure to analyze pupil first EyeVidDLC_Points ->
%post_DeepLab_pupil_correction to save pupil_long_axis.mat
%getBallMotion(datadir) to save moves_trace1.mat

djPrefs;
global pref

if nargin==0
    fprintf('\nno input\n')
    return
else
    datadir=varargin{1};
end
if nargin==1
    xlimits=[-300 900]; %x limits for axis
elseif nargin==2
    datadir=varargin{1};
    xlimits=varargin{2};
end

%cd to data dir and load Moves trace and Pupil trace
cd (datadir)
if ~exist('moves_trace1.mat')
    getBallMotion(datadir)
    load('moves_trace1.mat')
else
    load('moves_trace1.mat')
end
%normalize pupil
load('dirs.mat'); 
max_pupil=[]; min_pupil=[];
for d=1:length(dirs)
    try
        load([dirs{d} '\pupil_long_axis.mat'])
        l=smooth(laxis, 150);
        max_pupil=[max_pupil max(l)];
        min_pupil=[min_pupil min(l)];
    end
end
max_pupil=max(max_pupil);
min_pupil=min(min_pupil);
cont_files=dir('*.continuous');
[scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(cont_files(1).name); 
laxis(1:length(scaledtrace))=NaN;
if ~exist('pupil_long_axis.mat')
    try
        %EyeVidDLC_Points
    catch
        sprintf('\n No Pupil data was found. Make sure to analyze the video first.')
        cont_files=dir('*.continuous');
        [scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(cont_files(1).name); % load a sample cont file to know how long the recording was
        laxis(1:length(scaledtrace))=NaN;
    end
else
    load('pupil_long_axis.mat')
end

%fix delays in pupil trace
try
    pi_delay_calculator
catch
    pi_delay_on=NaN;
    pi_delay_off=NaN;
    save('pi_delay.mat', 'pi_delay_on', 'pi_delay_off')
end

load('pi_delay.mat')
on=floor(pi_delay_on/.030); %delay between OE and pi recording starts
off=floor(pi_delay_off/.030); %delay between OE and pi recording stops
if on>0
    try
        l=[zeros(on,1) laxis];
    catch
        try
            l=[zeros(on,1); laxis'];
        catch
            l=[zeros(on,1); laxis];
        end
    end
    laxis=l;
end
if off>0
laxis(end-off:end)=[]; %remove data points from the trace after OE turned off
end

laxis= fillmissing(laxis,'nearest');
%laxis=smooth(laxis,20);
laxis_normalized=laxis;
try
laxis_normalized=(laxis-min_pupil)/(max_pupil-min_pupil);
save('pupil_long_axis_normalized.mat', 'laxis_normalized')
end

[num, dem] = rat(length(scaledtrace)/length(laxis));
laxis1=resample(laxis_normalized, num, dem);
figure; plot(laxis_normalized);
TF=isoutlier(laxis1);
laxis1(TF)=NaN;
laxis1= fillmissing(laxis1,'nearest');
figure; plot(laxis1)

%fix sampling discrapances in moves trace
try
    load('moves_trace1.mat')
catch
    getBallMotion(pwd)
end

[num, dem] = rat(length(scaledtrace)/length(moves_trace1));
moves_trace=resample(moves_trace1, num, dem);
moves_trace=moves_trace-mode(moves_trace);

if exist('Events.mat')
    load('Events.mat')
else
    [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
end

j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'silentsound') ...
            |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise')| strcmp(Events(i).type, 'grating')
        j=j+1;
        alldurs(j)=Events(i).duration;
        if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, '2tone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).frequency;
        elseif strcmp(Events(i).type, 'whitenoise')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=-1000;
        elseif strcmp(Events(i).type, 'silentsound')
            allfreqs(j)=-2000;
            allamps(j)=nan; %flagging silent sound amp as nan
        elseif strcmp(Events(i).type, 'fmtone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).carrier_frequency;
        elseif strcmp(Events(i).type, 'grating')
            allfreqs(j)=Events(i).angle*1000;
            allamps(j)=Events(i).spatialfrequency;
        end
    end
end
allamps=allamps(~isnan(allamps)); %strip out nans from silent sound
freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
numfreqs=length(freqs);
numamps=length(amps);
numdurs=length(durs);

%check for laser in Events
%for now, setting AOPulseOn to 0 for all events
for i=1:length(Events)
    if isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
        LaserScheduled(i)=Events(i).laser; %whether the stim protocol scheduled a laser for this stim
        LaserOnOffButton(i)=Events(i).LaserOnOff; %whether the laser button was turned on
        LaserTrials(i)=LaserScheduled(i) & LaserOnOffButton(i);
    elseif isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %Not sure about this one. Assume no laser for now, but investigate.
        warning('ProcessTC_LFP: Cannot tell if laser button was turned on in djmaus GUI');
        LaserTrials(i)=0;
        Events(i).laser=0;%
        
    elseif ~isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %if neither of the right fields are there, assume no laser
        LaserTrials(i)=0;
        Events(i).laser=0;
    else
        error('wtf?')
    end
end
fprintf('\n%d laser pulses in this Events file', sum(LaserTrials))
try
    if sum(LaserOnOffButton)==0
        fprintf('\nLaser On/Off button remained off for entire file.')
    end
end
if sum(LaserTrials)>0
    IL=1;
else
    IL=0;
end

nreps=zeros(numfreqs, numamps, numdurs);
nrepsON=zeros(numfreqs, numamps, numdurs);
nrepsOFF=zeros(numfreqs, numamps, numdurs);
samprate=30e3;
j=0;
P1ON=[]; Moves1ON=[];
for i=1:length(Events)
    if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') |  strcmp(Events(i).type, 'silentsound') | ...
            strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating')
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate; %pos is in samples
        end
        laser=LaserTrials(i);
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        switch Events(i).type
            case {'tone', '2tone'}
                freq=Events(i).frequency;
                amp=Events(i).amplitude;
            case 'fmtone'
                freq=Events(i).carrier_frequency;%
                amp=Events(i).amplitude;
            case 'whitenoise'
                freq=-1000;
                amp=Events(i).amplitude;
            case 'silentsound'
                freq=-2000;
                amp=min(amps); %put silentsound in it's own column (freq=-2) in the lowest row
            case 'grating'
                amp=Events(i).spatialfrequency;
                freq=Events(i).angle*1000;
        end
        dur=Events(i).duration;
        findex= find(freqs==freq);
        aindex= find(amps==amp);
        dindex= find(durs==dur);
        nreps(findex, aindex, dindex)=nreps(findex, aindex, dindex)+1;
        P1(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=laxis1(region);
        Moves1(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=moves_trace(region);
        if laser
            nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
            P1ON(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=laxis1(region);
            Moves1ON(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=moves_trace(region);
        else
            nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
            P1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=laxis1(region);
            Moves1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=moves_trace(region);
        end
    end
end

outPM.P1ON=P1ON;
outPM.P1OFF=P1OFF;
outPM.Moves1ON=Moves1ON;
outPM.Moves1OFF=Moves1OFF;
outPM.P1=P1;
outPM.Moves1=Moves1;
outfilename=sprintf('outPupilMoves.mat')
save(outfilename, 'outPM')