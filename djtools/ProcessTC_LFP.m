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

djPrefs;
global pref

if nargin==0
    fprintf('\nno input\n')
    return
else
    datadir=varargin{1};
end
if nargin==1
    xlimits=[0 300]; %x limits for axis
    ylimits=[-.1 .2];
    prompt=('Please enter channel number: ');
    channel=input(prompt) ;
elseif nargin==2
    channel=varargin{2};
    xlimits=[0 300]; %x limits for axis
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

%fprintf('\nfound lostat file. using lostat %d %d',lostat );
%lostat=getlostat(expdate,session, filenum); %discard data after this position (in samples), [] to skip
%getgotat(expdate,session, filenum); %discard data before this position (in samples), [] to skip
%[datafile, Eventsfile, stimfile]=getfilenames(expdate, session, filenum);
%OEEventsfile=strrep(Eventsfile, 'AxopatchData1', 'OE');
%[D E S]=gogetdata(expdate,session,filenum);
fprintf('\nProcessing continuous data\nloading file: ')

%     fprintf('\ntrying to load %s...', datafile)
%     godatadir(expdate, session, filenum)
%     pathname='C:\Program Files\Open Ephys\ira_2014-02-14_13-43-56\'
%     cd(pathname)
%     filename=sprintf('102_CH/max(abs(Lasertrace))%s.continuous',channel)

%try to read OE filename from exper structure (only will work after
%02.14.14)

cd(pref.datapath)
cd(datadir)
filename=getContinuousFilename('.', channel);
if exist(filename, 'file')~=2 %couldn't find it
    filename=sprintf('114_CH%d.continuous', channel);
end
if exist(filename, 'file')~=2 %couldn't find it
    error(sprintf('could not find data file %s in datadir %s', filename, datadir))
end



[scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(filename);

%should record a copy of stim
stim=0*scaledtrace;

%might want to read settings.xml from this directory to grab info about
%filter settings, etc.

try
    load notebook.mat
catch
    warning('could not find notebook file')
end
%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

% ADC1: sound monitor
% ADC2: soundcard trigger monitor
% ADC3: laser monitor

%all_channels_data is the channel the Events are associated with
% the 0 channel Events are network Events
%all_channels_timestamps are in seconds
%all_channels_info is a struct with the following:
%Eventstype 5 are network events, according to https://open-ephys.atlassian.net/wiki/display/OEW/Network+Events
%   the content of those network event are stored in messages
%   Two network events are saved at the beginning of recording containing the system time and sampling information header.
%Eventstype 3 are TTL (digital input lines) one each for up & down
%eventId is 1 or 0 indicating up or down Events for the TTL signals

%FYI the nodeID corresponds to the processor that was the source of that
%Events or data. Which nodeIds correspond to which processor is listed in
%the settings.xml
%for example, in the test data file I'm working with, 100 is the Rhythm
%FPGA, and 102 is the bandpass filter.

%sanity check
fid=fopen('temp.txt', 'w');
for i=1:length(all_channels_timestamps)
    fprintf(fid, '%f, eventType %d, Id %d, channel %d\n', ...
        all_channels_timestamps(i), all_channels_info.eventType(i), ...
        all_channels_info.eventId(i), all_channels_data(i));
end
fclose(fid);


%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);


sound_index=0;
cont_files=dir('*.continuous');
[~, timestamps, ~] =load_open_ephys_data(cont_files(1).name);
StartAcquisitionSec=timestamps(1);
StartAcquisitionSamples=timestamps(1)*30e3;
for i=1:length(messages)
    str=messages{i};
    str2=strsplit(str);
    timestamp=str2num(str2{1});
    Events_type=str2{2};
    if strcmp(deblank(Events_type), 'StartAcquisition')
        %if present, a convenient way to find start acquisition time
        %for some reason not always present, though
%         StartAcquisitionSamples=timestamp;
%         StartAcquisitionSec=timestamp/sampleRate;
        check1=StartAcquisitionSamples;
    elseif strcmp(deblank(Events_type), 'Software')
%         StartAcquisitionSamples=timestamp;
%         StartAcquisitionSec=timestamp/sampleRate;
        check2=StartAcquisitionSamples;
    elseif strcmp(Events_type, 'TrialType')
        sound_index=sound_index+1;
        Events(sound_index).type=str2{3};
        for j=4:length(str2)
            str3=strsplit(str2{j}, ':');
            fieldname=str3{1};
            value=str2num(str3{2});
            Events(sound_index).(fieldname)= value;
        end
            Events(sound_index).message_timestamp_samples=timestamp - StartAcquisitionSamples;
            Events(sound_index).message_timestamp_sec=timestamp/sampleRate - StartAcquisitionSec;
        
        
        %get corresponding SCT TTL timestamp and assign to Event
        %first get all the soundcard triggers
        all_SCTs=[];
        for k=1:length(all_channels_timestamps)
            if all_channels_info.eventType(k)==3 & all_channels_info.eventId(k)==1 & all_channels_data(k)==2
                corrected_SCT=all_channels_timestamps(k)-StartAcquisitionSec;
                all_SCTs=[all_SCTs corrected_SCT];
            end
        end
        %find closest SCT by finding first before and first after, and
        %choosing whichever is closer
        [idx_after]=find(all_SCTs>Events(sound_index).message_timestamp_sec, 1); %find first SCT after the message timestamp
        [idx_before]=find(all_SCTs<=Events(sound_index).message_timestamp_sec, 1, 'last'); %find first SCT before the message timestamp
       
        if isempty(idx_after)
                        SCTtime_sec=all_SCTs(idx_before);
        elseif isempty(idx_before)
                        SCTtime_sec=all_SCTs(idx_after);
        elseif abs(diff([all_SCTs(idx_before), Events(sound_index).message_timestamp_sec])) <= abs(diff([all_SCTs(idx_after), Events(sound_index).message_timestamp_sec]))
            %before is closer
            SCTtime_sec=all_SCTs(idx_before);
        elseif abs(diff([all_SCTs(idx_before), Events(sound_index).message_timestamp_sec])) > abs(diff([all_SCTs(idx_after), Events(sound_index).message_timestamp_sec]))
                      %after is closer
            SCTtime_sec=all_SCTs(idx_after);
        else
        error('WTF how can this happen')
        end
        Events(sound_index).soundcard_trigger_timestamp_sec=SCTtime_sec;
        

    end
end

if ~exist('Events.mat')
save('Events.mat', 'Events')
elseif exist('Events.mat')
    load('Events.mat')
end
datatimestamps=datatimestamps-StartAcquisitionSec;

if exist('check1', 'var') & exist('check2', 'var')
    fprintf('start acquisition method agreement check: %d, %d', check1, check2);
end
fprintf('\nNumber of sound events (from network messages): %d', length(Events));
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d', length(all_SCTs));
if length(Events) ~=  length(all_SCTs)
    warning('ProcessTC_LFP: Number of sound events (from network messages) does not match Number of hardware triggers (soundcardtrig TTLs)')
    %[Events, all_SCTs, stimlog]=ResolveEventMismatch(Events, all_SCTs, stimlog  );
end

%messages is a list of all network event, which includes the stimuli
%messages sent by djmaus, SCTchannelas well as the "ChangeDirectory" and
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
    fprintf('\ncouldnt find a stimtrace, loaded a differnt trace for plotting')
end


monitor = 0;
if monitor
    figure
    set(gcf, 'pos', [-1853 555 1818 420]);

    %   I'm running the soundcard trigger (SCT) into ai1 as another sanity check.
    try
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(getSCTfile('.'));
        SCTtimestamps=SCTtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
    catch
        warning('cannot find soundcard trigger file')
    end
   

%     SCTtrace=SCTtrace./max(abs(SCTtrace));
    hold on
     plot(SCTtimestamps, SCTtrace./max(SCTtrace), 'k')
    plot(Lasertimestamps, Lasertrace./max(Lasertrace)+1, 'c')
    plot(Stimtimestamps, Stimtrace./max(Stimtrace)+2, 'm')
    plot(all_SCTs, .8+0*all_SCTs, 'g^')

    %plot "software trigs" i.e. network messages in red o's
    for i=1:length(Events)
        plot(Events(i).message_timestamp_sec, 0, 'ro');
        text(Events(i).message_timestamp_sec, .2, sprintf('network message #%d', i), 'color', 'r')
        plot(Events(i).soundcard_trigger_timestamp_sec, .5, 'g*');
        text(Events(i).soundcard_trigger_timestamp_sec, .7, sprintf('SCT #%d', i), 'color', 'g')
    end
    
    %all_channels_info.eventType(i) = 3 for digital line in (TTL), 5 for network Events
    %all_channels_info.eventId = 1 for rising edge, 0 for falling edge
    %all_channels_data(i) is the digital input line channel
    
    
    % plot TTL SCTs in green ^=on, v=off
    for i=1:length(all_channels_timestamps)
        if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==2
            plot(all_channels_timestamps(i)-StartAcquisitionSec, 1, 'g^')
            text(all_channels_timestamps(i)-StartAcquisitionSec, 2, 'TTL on/off')
        elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==2
            plot(all_channels_timestamps(i)-StartAcquisitionSec, 1, 'gv')
        end
    end
    
    scaledtrace=scaledtrace./max(abs(scaledtrace));
    plot(datatimestamps, scaledtrace, 'm')
    
    for i=1:length(Events)
        xlim([Events(i).message_timestamp_sec-.02 Events(i).message_timestamp_sec+.2])
        ylim([-4 2])
        pause(0.5)
    end
end %if monitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







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
        %region1=datatimestamps>start/samprate & datatimestamps<stop/samprate;
        %if length(scaledtrace(region1))~=length(region)
            
        %region=region1;
        %if isempty(find(region<1)) %(disallow negative or zero start times)
            %             if stop>lostat
            %                 fprintf('\ndiscarding trace (after lostat)')
            %             elseif start<gotat
            %                fprintf('\ndiscarding trace (before gotat)')
            %                %commented out by ira 09-05-2013
            %             else
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
%end

%region=length(M1OFF);
traces_to_keep=[];
if ~isempty(traces_to_keep)
    fprintf('\n using only traces %d, discarding others', traces_to_keep);
    mM1=mean(M1(:,:,:,traces_to_keep,:), 4);
    mM1ON=mean(M1ON(:,:,:,traces_to_keep,:), 4);
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
fprintf('\n saved to %s', outfilename)

