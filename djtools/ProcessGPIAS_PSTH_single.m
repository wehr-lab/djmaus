function ProcessGPIAS_PSTH_single(varargin)

%processes a single .t file of clustered spiking GPIAS data from djmaus
%
% usage: ProcessGPIAS_PSTH_single(datadir, t_filename, [xlimits],[ylimits])
% (xlimits, ylimits are optional)
% xlimits default to [-1.5*max(gapdurs) 2*soa]
% saves to outfile

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

try
    xlimits=varargin{3};
catch
    xlimits=[];
end
if isempty(xlimits)
    xlimits=[];
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end

filename=varargin{2};
[p,f,ext]=fileparts(filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

fprintf('\nchannel %d, cluster %d', channel, clust)

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

try
    load notebook.mat
    %check if nb and stimlog are actually there
    if ~exist('stimlog','var')
        stimlog=[];
        fprintf('\nfound notebook file but there was no stimlog in it!!!');
        fprintf('\n in principle it should be possible to reconstruct all stimuli\nfrom network events, but that would take a lot of coding, does this problem happen often enough to be worth the effort???');
        error('found notebook file but there was no stimlog in it!!!');
    end
    if ~exist('nb','var')
        nb=[];
        warning('found notebook file but there was no nb notebook information in it!!!');
    end
catch
    warning('could not find notebook file')
end

%read messages
messagesfilename='messages.events';
if exist(messagesfilename, 'file')~=2
    error('Could not find messages.events file. Is this the right data directory?')
end
[messages] = GetNetworkEvents(messagesfilename);


%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

%all_channels_data is the channel the Events are associated with
%the 0 channel Events are network Events
%all_channels_timestamps are in seconds
%all_channels_info is a struct with the following:
%Eventstype 5 are network events, according to https://open-ephys.atlassian.net/wiki/display/OEW/Network+Events
%the content of those network events are stored in messages
%Two network events are saved at the beginning of recording containing the system time and sampling information header.
%Eventstype 3 are TTL (digital input lines) one each for up & down
%eventId is 1 or 0 indicating up or down Events for the TTL signals
%FYI the nodeID corresponds to the processor that was the source of that
%Events or data. Which nodeIds correspond to which processor depends on
%your configuration, and is listed in the settings.xml. For example, 100
%might be the Rhythm FPGA, and 102 the bandpass filter (depending on your
%config).

% get all SCT timestamps
sound_index=0;
for i=1:length(messages)
    str=messages{i};
    str2=strsplit(str);
    timestamp=str2num(str2{1});
    Events_type=str2{2};
    if strcmp(deblank(Events_type), 'StartAcquisition')
        %if present, a convenient way to find start acquisition time
        %for some reason not always present, though
        StartAcquisitionSamples=timestamp;
        StartAcquisitionSec=timestamp/sampleRate;
        fprintf('\nStartAcquisitionSec=%g', StartAcquisitionSec)
        check1=StartAcquisitionSamples;
    elseif strcmp(deblank(Events_type), 'Software')
        StartAcquisitionSamples=timestamp;
        StartAcquisitionSec=timestamp/sampleRate;
        fprintf('\nStartAcquisitionSec=%g', StartAcquisitionSec)
        check2=StartAcquisitionSamples;
    elseif strcmp(Events_type, 'TrialType')
        sound_index=sound_index+1;
        Events(sound_index).type=str2{3};
        for j=4:length(str2)
            str3=strsplit(str2{j}, ':');
            fieldname=str3{1};
            value=str2num(str3{2});
            if isempty(value) %this happens if it's a real string (like soaflag), not a 'num'
                value=str3{2}; %just use the string. E.g. 'soa' or 'isi'
            end
            Events(sound_index).(fieldname)= value;
        end
        Events(sound_index).message_timestamp_samples=timestamp - StartAcquisitionSamples;
        Events(sound_index).message_timestamp_sec=timestamp/sampleRate - StartAcquisitionSec;
        all_SCTs=[];
        for k=1:length(all_channels_timestamps)
            if all_channels_info.eventType(k)==3 & all_channels_info.eventId(k)==1 & all_channels_data(k)==2
                corrected_SCT=all_channels_timestamps(k)-StartAcquisitionSec;
                all_SCTs=[all_SCTs corrected_SCT];
            end
        end
        %get corresponding SCT TTL timestamp and assign to Event
        %old way (from TC) won't work, since the network events get ahead of the SCTs.
        %another way (dumb and brittle) is to just use the corresponding
        %index, assuming they occur in the proper order with no drops or
        %extras
        try
            SCTtime_sec=all_SCTs(sound_index);
            Events(sound_index).soundcard_trigger_timestamp_sec=SCTtime_sec;
        catch
            %     one reason this won't work is if you stop recording during the
            %     play-out, i.e. stimuli were logged (to stimlog and network events)
            %     but never played before recording was stopped
            if sound_index>length(all_SCTs)
                Events(sound_index).soundcard_trigger_timestamp_sec=nan;
            end
        end
    end
end

fprintf('\nNumber of sound events (from network messages): %d', length(Events));
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d', length(all_SCTs));
try
    fprintf('\nNumber of logged stimuli in notebook: %d', length(stimlog));
catch
    fprintf('\nCould not find stimlog, no logged stimuli in notebook!!');
end
if length(Events) ~=  length(all_SCTs)
    warning('ProcessGPIAS_PSTH_single: Number of sound events (from network messages) does not match Number of hardware triggers (soundcardtrig TTLs)')
    [Events, all_SCTs, stimlog]=ResolveEventMismatch(Events, all_SCTs, stimlog  );
end

if exist('check1', 'var') & exist('check2', 'var')
    fprintf('start acquisition method agreement check: %d, %d', check1, check2);
end

%read MClust .t file
fprintf('\nreading MClust output file %s', filename)
spiketimes=read_MClust_output(filename)'/10000; %spiketimes now in seconds
%correct for OE start time, so that time starts at 0
spiketimes=spiketimes-StartAcquisitionSec;
totalnumspikes=length(spiketimes);
fprintf('\nsuccessfully loaded MClust spike data')
Nclusters=1;

monitor = 0; %should really clean this up / put in separate function or sub function or nested function
if monitor
    %   I'm running the soundcard trigger (SCT) into ai1 as another sanity check.
    SCTfname=getSCTfile(datadir);
    if isempty(SCTfname)
        warning('could not find soundcard trigger file')
    else
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
    end
    
    %     %here I'm loading a data channel to get good timestamps - the ADC timestamps are screwed up
    %     datafname=getContinuousFilename( pwd, 1 );
    %     [scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(datafname);
    
    SCTtimestamps=SCTtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
    
    
    %sanity check
    % fid=fopen('temp.txt', 'w');
    % for i=1:length(all_channels_timestamps)
    %     fprintf(fid, '%f, eventType %d, Id %d, channel %d\n', ...
    %         all_channels_timestamps(i), all_channels_info.eventType(i), ...
    %         all_channels_info.eventId(i), all_channels_data(i));
    % end
    % fclose(fid);
    
    
    
    %messages is a list of all network event, which includes the stimuli
    %messages sent by djmaus, as well as the "ChangeDirectory" and
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
    
    
    
    
    
    figure
    hold on
    %set(gcf, 'pos', [-1853 555 1818 420]);
    SCTtrace=SCTtrace./max(abs(SCTtrace));
    %     scaledtrace=scaledtrace./max(abs(scaledtrace));
    plot(SCTtimestamps, SCTtrace) %plotting with data timestamps to work around wierd bug in ADC timestamps
    % plot(datatimestamps, scaledtrace, 'm') %plotting with data timestamps to work around wierd bug in ADC timestamps
    
    hold on
    %plot "software trigs" i.e. network messages in red o's
    for i=1:length(Events)
        plot(Events(i).message_timestamp_sec, .25, 'ro');
        plot(Events(i).soundcard_trigger_timestamp_sec, 1, 'g*');
        text(Events(i).message_timestamp_sec, .5, sprintf('network message #%d', i))
        text(Events(i).message_timestamp_sec, .75, sprintf('%s', Events(i).type))
        text(Events(i).soundcard_trigger_timestamp_sec, 1.25, sprintf('SCT #%d', i))
    end
    
    %all_channels_info.eventType(i) = 3 for digital line in (TTL), 5 for network Events
    %all_channels_info.eventId = 1 for rising edge, 0 for falling edge
    %all_channels_data(i) is the digital input line channel
    
    
    % plot TTL SCTs in green ^=on, v=off
    for i=1:length(all_channels_timestamps)
        if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==2
            plot(all_channels_timestamps(i), 1, 'g^')
            text(all_channels_timestamps(i), 1, 'TTL on/off')
        elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==2
            plot(all_channels_timestamps(i), 1, 'gv')
        end
    end
    
    
    for i=1:length(Events)
        xlim([Events(i).message_timestamp_sec-2 Events(i).message_timestamp_sec+3])
        ylim([-5 2])
        pause(.01)
    end
end %if monitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\ncomputing tuning curve...');
samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS')
        j=j+1;
        allsoas(j)=Events(i).soa;
        allgapdurs(j)=Events(i).gapdur;
        allgapdelays(j)=Events(i).gapdelay;
        allpulseamps(j)=Events(i).pulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allnoiseamps(j)=Events(i).amplitude;
    elseif strcmp(Events(i).type, 'gapinnoise')
        allsoas(j)=Events(i).soa;
        allgapdurs(j)=Events(i).gapdur;
        allgapdelays(j)=Events(i).gapdelay;
        allpulseamps(j)=Events(i).pulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allnoiseamps(j)=Events(i).amplitude;
    end
    
end
gapdurs=unique(allgapdurs);
pulsedurs=unique(allpulsedurs);
soas=unique(allsoas);
gapdelays=unique(allgapdelays);
pulseamps=unique(allpulseamps);
pulsedurs=unique(allpulsedurs);
noiseamps=unique(allnoiseamps);
numgapdurs=length(gapdurs);
numpulseamps=length(pulseamps);
nrepsON=zeros( numgapdurs, numpulseamps);
nrepsOFF=zeros( numgapdurs, numpulseamps);

if length(noiseamps)~=1
    error('not able to handle multiple noiseamps')
end
if length(gapdelays)~=1
    error('not able to handle multiple gapdelays')
end
if length(pulsedurs)~=1
    error('not able to handle multiple pulsedurs')
end
if length(soas)~=1
    error('not able to handle multiple soas')
end
noiseamp=noiseamps;
soa=soas;
pulsedur=pulsedurs;
gapdelay=gapdelays;

%check for laser in Events
for i=1:length(Events)
    if isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
        LaserScheduled(i)=Events(i).laser; %whether the stim protocol scheduled a laser for this stim
        LaserOnOffButton(i)=Events(i).LaserOnOff; %whether the laser button was turned on
        LaserTrials(i)=LaserScheduled(i) & LaserOnOffButton(i);
        if isempty(stimlog(i).LaserStart)
            LaserStart(i)=nan;
            LaserWidth(i)=nan;
            LaserNumPulses(i)=nan;
            LaserISI(i)=nan;
        else
            LaserStart(i)=stimlog(i).LaserStart;
            LaserWidth(i)=stimlog(i).LaserWidth;
            LaserNumPulses(i)=stimlog(i).LaserNumPulses;
            LaserISI(i)=stimlog(i).LaserISI;
        end
        
    elseif isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %Not sure about this one. Assume no laser for now, but investigate.
        warning('ProcessGPIAS_PSTH_single: Cannot tell if laser button was turned on in djmaus GUI');
        LaserTrials(i)=0;
        Events(i).laser=0;
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
%try to load laser and stimulus monitor files
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

M1ON=[];M1OFF=[];
nrepsON=zeros(numgapdurs, numpulseamps);
nrepsOFF=zeros(numgapdurs, numpulseamps);

% %find optimal axis limits
if isempty(xlimits)
    xlimits(1)=-1.5*max(gapdurs);
    xlimits(2)=2*soa;
end
fprintf('\nprocessing with xlimits [%d-%d]', xlimits(1), xlimits(2))

%extract the traces into a big matrix M
j=0;
inRange=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') | strcmp(Events(i).type, 'gapinnoise')
        
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=LaserTrials(i);
        start=pos + gapdelay/1000 +xlimits(1)/1000; %start is in seconds
        stop=pos+ gapdelay/1000 + xlimits(2)/1000; %stop is in seconds
        region=round(start*samprate)+1:round(stop*samprate);
        if start>0 %(disallow negative or zero start times)
            gapdur=Events(i).gapdur;
            gdindex= find(gapdur==gapdurs);
            pulseamp=Events(i).pulseamp;
            paindex= find(pulseamp==pulseamps);
            st=spiketimes; %are in seconds
            st_inrange=st(st>start & st<stop); % spiketimes in region, in seconds relative to start of acquisition
            spikecount=length(st_inrange); % No. of spikes fired in response to this rep of this stim.
            inRange=inRange+ spikecount; %accumulate total spikecount in region
            spiketimes1=st_inrange*1000 - pos*1000 - gapdelay;%covert to ms after gap termination
            spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window
            if laser
                nrepsON(gdindex,paindex)=nrepsON(gdindex,paindex)+1;
                M1ON(gdindex,paindex, nrepsON(gdindex,paindex)).spiketimes=spiketimes1; % Spike times
                M1ONspikecounts(gdindex,paindex,nrepsON(gdindex,paindex))=spikecount; % No. of spikes
                M1spontON(gdindex,paindex, nrepsON(gdindex,paindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                if LaserRecorded
                    M1ONLaser(gdindex, paindex, nrepsON(gdindex,paindex),:)=Lasertrace(region);
                end
                if StimRecorded
                    M1ONStim(gdindex,paindex, nrepsON(gdindex,paindex),:)=Stimtrace(region);
                end
            else
                nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                
                M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex)).spiketimes=spiketimes1;
                M1OFFspikecounts(gdindex,paindex,nrepsOFF(gdindex,paindex))=spikecount;
                M1spontOFF(gdindex,paindex, nrepsOFF(gdindex,paindex))=spont_spikecount;
                if LaserRecorded
                    M1OFFLaser(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=Lasertrace(region);
                end
                if StimRecorded
                    M1OFFStim(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=Stimtrace(region);
                end
            end
        end
    end
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))
fprintf('\ntotal num spikes: %d', length(spiketimes))
fprintf('\nIn range: %d', inRange)

% Accumulate spiketimes across trials, for psth...
for gdindex=1:numgapdurs; % Hardcoded.
    for paindex=1:numpulseamps
        % on
        spiketimesON=[];
        for rep=1:nrepsON(gdindex,paindex)
            spiketimesON=[spiketimesON M1ON(gdindex,paindex, rep).spiketimes];
        end
        
        % All spiketimes for a given f/a/d combo, for psth:
        mM1ON(gdindex,paindex).spiketimes=spiketimesON;
        
        % off
        spiketimesOFF=[];
        for rep=1:nrepsOFF(gdindex,paindex)
            spiketimesOFF=[spiketimesOFF M1OFF(gdindex,paindex, rep).spiketimes];
        end
        mM1OFF(gdindex,paindex).spiketimes=spiketimesOFF;
    end
end


if ~IL %no laser pulses in this file
    mM1ONspikecount=[];
    sM1ONspikecount=[];
    semM1ONspikecount=[];
    M1ONspikecounts=[];
else
    mM1ONspikecount=mean(M1ONspikecounts,3); % Mean spike count
    sM1ONspikecount=std(M1ONspikecounts,[],3); % Std of the above
    semM1ONspikecount=sM1ONspikecount./sqrt(max(nrepsON(:))); % Sem of the above
    % Spont
    mM1spontON=mean(M1spontON,3);
    sM1spontON=std(M1spontON,[],3);
    for clust=1:Nclusters
        semM1spontON=sM1spontON./sqrt(max(nrepsON(:)));
    end
end
if isempty(M1OFF) %only laser pulses in this file
    mM1OFFspikecount=[];
    sM1OFFspikecount=[];
    semM1OFFspikecount=[];
    M1OFFspikecounts=[];
else
    mM1OFFspikecount=mean(M1OFFspikecounts,3); % Mean spike count
    sM1OFFspikecount=std(M1OFFspikecounts,[],3); % Std of the above
    semM1OFFspikecount(:,:)=sM1OFFspikecount./sqrt(max(nrepsOFF(:))); % Sem of the above
    % Spont
    mM1spontOFF=mean(M1spontOFF,3);
    sM1spontOFF=std(M1spontOFF,[],3);
    semM1spontOFF=sM1spontOFF./sqrt(max(nrepsOFF(:)));
end

%average laser and stimulus monitor M matrices across trials
if LaserRecorded
    for gdindex=1:numgapdurs
        for paindex =1:numpulseamps
            if nrepsON(gdindex,paindex)>0
                mM1ONLaser(gdindex, paindex, :)=mean(M1ONLaser(gdindex,paindex, 1:nrepsON(gdindex,paindex),:), 3);
            else %no reps for this stim, since rep=0
                mM1ONLaser(gdindex,paindex,:)=zeros(size(region));
            end
            if nrepsOFF(gdindex,paindex)>0
                mM1OFFLaser(gdindex,paindex,:)=mean(M1OFFLaser(gdindex,paindex, 1:nrepsOFF(gdindex,paindex),:), 3);
            else %no reps for this stim, since rep=0
                mM1OFFLaser(gdindex,paindex,:)=zeros(size(region));
            end
        end
    end
end
if StimRecorded
    for gdindex=1:numgapdurs
        for paindex =1:numpulseamps
            if nrepsON(gdindex,paindex)>0
                mM1ONStim(gdindex,paindex,:)=mean(M1ONStim(gdindex,paindex, 1:nrepsON(gdindex,paindex),:), 3);
            else %no reps for this stim, since rep=0
                mM1ONStim(gdindex,paindex,:)=zeros(size(region));
            end
            if nrepsOFF(gdindex,paindex)>0
                mM1OFFStim(gdindex,paindex,:)=mean(M1OFFStim(gdindex,paindex, 1:nrepsOFF(gdindex,paindex),:), 3);
            else %no reps for this stim, since rep=0
                mM1OFFStim(gdindex,paindex,:)=zeros(size(region));
            end
        end
    end
end

%save to outfiles
%one outfile for each cell

%after squeezing cluster, saves with the following dimensions:
% M1ON(numgapdurs, numpulseamps, nrepsON).spiketimes
% mM1ON(numgapdurs, numpulseamps).spiketimes
% mM1ONspikecount(numgapdurs, numpulseamps)

out.IL=IL;
out.Nclusters=Nclusters;
out.tetrode=channel;
out.channel=channel;
out.cluster=clust; %there are some redundant names here
out.cell=clust;
if IL
    out.M1ON=M1ON; %isn't this so much easier?
    out.mM1ON=mM1ON;
    out.mM1ONspikecount=mM1ONspikecount;
    out.sM1ONspikecount=sM1ONspikecount;
    out.semM1ONspikecount=semM1ONspikecount;
    out.mM1spontON=mM1spontON;
    out.sM1spontON=sM1spontON;
    out.semM1spontON=semM1spontON;
    out.M_LaserStart=M_LaserStart;
    out.M_LaserWidth=M_LaserWidth;
    out.M_LaserNumPulses=M_LaserNumPulses;
    out.M_LaserISI=M_LaserISI;
    out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);
    
else
    out.M1ON=[];
    out.mM1ONspikecount=[];
    out.sM1ONspikecount=[];
    out.semM1ONspikecount=[];
    out.mM1ON=[];
    out.mM1spontON=[];
    out.sM1spontON=[];
    out.semM1spontON=[];
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
    out.M_LaserStart=[];
    out.M_LaserWidth=[];
    out.M_LaserNumPulses=[];
    out.M_LaserISI=[];
end
if isempty(M1OFF)
    out.M1OFF=[];
    out.mM1OFF=[];
    out.mM1OFFspikecount=[];
    out.sM1OFFspikecount=[];
    out.semM1OFFspikecount=[];
    out.mM1spontOFF=[];
    out.sM1spontOFF=[];
    out.semM1spontOFF=[];
else
    out.M1OFF=M1OFF;
    out.mM1OFF=mM1OFF;
    out.mM1OFFspikecount=mM1OFFspikecount;
    out.sM1OFFspikecount=sM1OFFspikecount;
    out.semM1OFFspikecount=semM1OFFspikecount;
    out.mM1spontOFF=mM1spontOFF;
    out.sM1spontOFF=sM1spontOFF;
    out.semM1spontOFF=semM1spontOFF;
end
out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.xlimits=xlimits;
out.samprate=samprate;
out.datadir=datadir;
out.spiketimes=spiketimes;

out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel

if LaserRecorded
    if exist('M1ONLaser')
        out.M1ONLaser=M1ONLaser;
        out.mM1ONLaser=mM1ONLaser;
    else
        out.M1ONLaser=[];
        out.mM1ONLaser=[];
    end
    if exist('M1OFFLaser')
        out.M1OFFLaser=M1OFFLaser;
        out.mM1OFFLaser=mM1OFFLaser;
    else
        out.M1OFFLaser=[];
        out.mM1OFFLaser=[];
    end
else
    out.M1ONLaser=[];
    out.mM1ONLaser=[];
    out.M1OFFLaser=[];
    out.mM1OFFLaser=[];
end
if StimRecorded
    if exist('M1ONStim')
        out.M1ONStim=M1ONStim;
        out.mM1ONStim=mM1ONStim;
    else
        out.M1ONStim=[];
        out.mM1ONStim=[];
    end
    if exist('M1OFFStim')
        out.M1OFFStim=M1OFFStim;
        out.mM1OFFStim=mM1OFFStim;
    else
        out.M1OFFStim=[];
        out.mM1OFFStim=[];
    end
else
    out.M1ONStim=[];
    out.mM1ONStim=[];
    out.M1OFFStim=[];
    out.mM1OFFStim=[];
end

try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
end
outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
save (outfilename, 'out')
end




