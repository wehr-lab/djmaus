function ProcessTC_LFP(varargin)
%processes continuous tuning curve data from djmaus.
%sorts Open Ephys continuous data into a big response matrix.
%use for LFP, WC, or other continuous data (i.e. does not extract spikes)
%
%usage: ProcessTC_OE(datapath, [channel], [xlimits], [ylimits])
%channel is a number not a string
%saves output in an outfile
%
%notes: whitenoise plotted as freq=-1 kHz, silent sound as -2 kHz
%       data are in microvolts

djPrefs;
global pref

if nargin==0
    fprintf('\nno input\n')
    return
else
    datadir=varargin{1};
end
if nargin==1
    xlimits=[-20 100]; %x limits for axis
    ylimits=[-.1 .2];
    prompt=('Please enter channel number: ');
    channel=input(prompt) ;
elseif nargin==2
    channel=varargin{2};
    xlimits=[-20 100]; %x limits for axis
    ylimits=[-.1 .2];
elseif nargin==3
    channel=varargin{2};
    xlimits=varargin{3};
    ylimits=[-.1 .2];
elseif nargin==4
    channel=varargin{2};
    xlimits=varargin{3};
    ylimits=varargin{4};
else
    error('wrong number of arguments');
end
if ischar(channel) channel=str2num(channel);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

try
    load notebook.mat
catch
    warning('could not find notebook file')
end

%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);


%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

%get Events and soundcard trigger timestamps
[Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps

%check if this is an appropriate stimulus protocol
switch GetPlottingFunction(datadir)
    case {'PlotTC_LFP', 'PlotTC_PSTH'}
    otherwise
        error('This does not appear to be a tuning curve stimulus protcol');
end

node='';
nodes={};
j=0;
NodeIds=getNodes(pwd);
for i=1:length(NodeIds)
    filename=sprintf('%s_CH%d.continuous', NodeIds{i}, channel);
    
    if exist(filename,'file')
        node=NodeIds{i};
        j=j+1;
        nodes{j} =node;
    end
end
%note if there are multiple continuous recordings of the requested channel,
%this will only return the last node
%node=nodes{1};
if length(nodes)>1
    fprintf('\n%d nodes found for channel %d: ', length(nodes), channel)
    fprintf('\t%s', nodes{:})
end
fprintf('\n using node %s', node)


filename=sprintf('%s_CH%d.continuous', node, channel);
if exist(filename, 'file')~=2 %couldn't find it
    error(sprintf('could not find data file %s in datadir %s', filename, datadir))
end

[scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(filename);

%messages is a list of all network event, which includes the stimuli
%messages sent by djmaus, SCTchannel as well as the "ChangeDirectory" and
%"GetRecordingPath" messages sent by djmaus, as well as 2 initial system
%messages. I strip out the stimulus (sound) event and put them in "Events."
%Events is a list of sound event, which were sent by djmaus with the
%'TrialType' flag.

%sanity check, continued
% fid=fopen('temp.txt', 'a');
% for i=1:length(Events)
%     fprintf(fid, '%f, %s, freq %f\n', ...
%         Events(i).message_timestamp_sec, Events(i).type, Events(i).frequency);
% end
% fclose(fid);


%try to load laser and stimulus monitor
Lasertrace=0*scaledtrace;
stimtrace=0*scaledtrace;
try
    [Lasertrace, Lasertimestamps, Laserinfo] =load_open_ephys_data(getLaserfile('.'));
    %Lasertimestamps=Lasertimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
    Lasertimestamps=Lasertimestamps-Lasertimestamps(1);
    % Lasertrace=Lasertrace./max(abs(Lasertrace));
    fprintf('\nsuccessfully loaded laser trace')
end
try
    [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(getStimfile('.'));
    %Stimtimestamps=Stimtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
    Stimtimestamps=Stimtimestamps-Stimtimestamps(1);
    %  Stimtrace=Stimtrace./max(abs(Stimtrace));
    fprintf('\nsuccessfully loaded stim trace')
catch
    [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data('105_ADC3.continuous');
    fprintf('\ncouldnt find a stimtrace, loaded a different trace for plotting')
    
    if isempty(getLaserfile('.'))
        LaserRecorded=0;
    else
        LaserRecorded=1;
    end
    if isempty(getStimfile('.'))
        StimRecorded=0;
    else
        StimRecorded=1;
    end
    
    if LaserRecorded
        try
            [Lasertrace, Lasertimestamps, Laserinfo] =load_open_ephys_data(getLaserfile('.'));
            Lasertimestamps=Lasertimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
            Lasertrace=Lasertrace./max(abs(Lasertrace));
            fprintf('\nsuccessfully loaded laser trace\n')
        catch
            fprintf('\nfound laser file %s but could not load laser trace', getLaserfile('.'))
        end
    else
        fprintf('\nLaser trace not recorded')
    end
    if StimRecorded
        try
            [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(getStimfile('.'));
            Stimtimestamps=Stimtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
            Stimtrace=Stimtrace./max(abs(Stimtrace));
            fprintf('\nsuccessfully loaded stim trace')
        catch
            fprintf('\nfound stim file %s but could not load stim trace', getStimfile('.'))
        end
    else
        fprintf('\nSound stimulus trace not recorded')
    end
end


fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

%get freqs/amps
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
%if lasers were used, we'll un-interleave them and save ON and OFF data

M1=[];
M1ON=[];M1OFF=[];
M1ONLaser=[]; M1OFFLaser=[];
M1ONStim=[];M1OFFStim=[];

nreps=zeros(numfreqs, numamps, numdurs);
nrepsON=zeros(numfreqs, numamps, numdurs);
nrepsOFF=zeros(numfreqs, numamps, numdurs);

fprintf('\nextracting traces into BIG matrix...');

samprate=sampleRate;

%extract the traces into a big matrix M
j=0;
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
        
        if isempty(find(region<1)) %(disallow negative or zero start times)
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
            M1(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=scaledtrace(region);
            M1stim(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=Stimtrace(region);
            if laser
                nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
                M1ON(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=scaledtrace(region);
                M1ONLaser(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=Lasertrace(region);
                M1ONStim(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=Stimtrace(region);
            else
                nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
                M1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=scaledtrace(region);
                M1OFFLaser(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Lasertrace(region);
                M1OFFStim(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Stimtrace(region);
            end
        end
    end
end

%region=length(M1OFF);
traces_to_keep=[];
if ~isempty(traces_to_keep)
    fprintf('\n using only traces %d to %d, discarding others\n', traces_to_keep(1),traces_to_keep(end));
    mM1=mean(M1(:,:,:,traces_to_keep,:), 4);
    if ~isempty(M1ON)
        mM1ON=mean(M1ON(:,:,:,traces_to_keep,:), 4);
    else
        mM1ON = [];
        mM1stim = [];
        mM1ONStim = [];
        mM1ONLaser = [];
        mM1OFFStim = [];
        mM1OFFLaser = [];
    end
    mM1OFF=mean(M1OFF(:,:,:,traces_to_keep,:), 4);
else
    for aindex=1:numamps
        for findex=1:numfreqs
            for dindex=1:numdurs
                if nreps(findex, aindex, dindex)>0
                    mM1(findex, aindex, dindex,:)=mean(M1(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                    mM1stim(findex, aindex, dindex,:)=mean(M1stim(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1(findex, aindex, dindex,:)=zeros(size(region));
                    mM1stim(findex, aindex, dindex,:)=zeros(size(region));
                end
                if nrepsON(findex, aindex, dindex)>0
                    mM1ON(findex, aindex, dindex,:)=mean(M1ON(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                    mM1ONStim(findex, aindex, dindex,:)=mean(M1ONStim(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                    mM1ONLaser(findex, aindex, dindex,:)=mean(M1ONLaser(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1ON(findex, aindex, dindex,:)=zeros(size(region));
                    mM1ONStim(findex, aindex, dindex,:)=zeros(size(region));
                    mM1ONLaser(findex, aindex, dindex,:)=zeros(size(region));
                end
                if nrepsOFF(findex, aindex, dindex)>0
                    mM1OFF(findex, aindex, dindex,:)=mean(M1OFF(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    mM1OFFStim(findex, aindex, dindex,:)=mean(M1OFFStim(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    mM1OFFLaser(findex, aindex, dindex,:)=mean(M1OFFLaser(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1OFF(findex, aindex, dindex,:)=zeros(size(region));
                    mM1OFFStim(findex, aindex, dindex,:)=zeros(size(region));
                    mM1OFFLaser(findex, aindex, dindex,:)=zeros(size(region));
                end
            end
        end
    end
end


%find optimal axis ylimits
if ylimits<0
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            trace=mM1(findex, aindex, dindex,:);
            ylimits(1)=min(ylimits(1), min(trace));
            ylimits(2)=max(ylimits(2), max(trace));
        end
    end
end

fprintf('\nsaving out file...');

%assign outputs
out.scaledtrace=scaledtrace;
out.M1=M1;
out.M1ON=M1ON;
out.M1ONStim=M1ONStim;
out.M1ONLaser=M1ONLaser;
out.mM1ON=mM1ON;
out.mM1ONStim=mM1ONStim;
out.mM1ONLaser=mM1ONLaser;
out.M1OFF=M1OFF;
out.M1OFFStim=M1OFFStim;
out.M1OFFLaser=M1OFFLaser;
out.mM1OFF=mM1OFF;
out.mM1OFFStim=mM1OFFStim;
out.mM1OFFLaser=mM1OFFLaser;
out.M1stim=M1stim;
out.mM1stim=mM1stim;
out.mM1=mM1;
out.datadir=datadir;
out.freqs=freqs;
out.amps=amps;
out.durs=durs;
out.nreps=nreps;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.numfreqs=numfreqs;
out.numamps=numamps;
out.numdurs=numdurs;
out.traces_to_keep=traces_to_keep;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.IL=IL;
out.xlimits=xlimits;
out.ylimits=ylimits;
out.samprate=samprate;
% out.nstd=nstd;
out.channel=channel;
%should probably save header info and stuff like that



outfilename=sprintf('outLFP_ch%d.mat',channel);
save(outfilename, 'out')
outfilename=sprintf('outLFP_ch%d',channel);
save(outfilename, 'out', '-v7.3' )
fprintf('\n saved to %s', outfilename)

