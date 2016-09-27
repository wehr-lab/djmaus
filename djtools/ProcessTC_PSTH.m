function ProcessTC_PSTH(varargin)

%processes clustered spiking tuning curve data from djmaus
%
% usage: ProcessTC_PSTH(datadir, channel, [xlimits],[ylimits])
% (xlimits, ylimits are optional)
% xlimits default to [0 200]
% channel (tetrode) number should be an integer, if omitted you will be prompted for it
% automatically processes data for all cells clustered for that tetrode
% saves to outfile


if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

xlimits=[0 200];
try
    xlimits=varargin{4};
end
if isempty(xlimits)
    xlimits=[0 200];
end
try
    ylimits=varargin{5};
catch
    ylimits=[];
end
try
    channel=varargin{2};
catch
    prompt=('please enter channel number: ');
    channel=input(prompt);
end
if strcmp('char',class(channel))
    channel=str2num(channel);
end
try
    cell=varargin{3};
catch
    cell=[];
end
if strcmp('char',class(cell))
    cell=str2num(cell);
end


fprintf('\nprocessing with xlimits [%d-%d]', xlimits(1), xlimits(2))

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
        check1=StartAcquisitionSamples;
    elseif strcmp(deblank(Events_type), 'Software')
        StartAcquisitionSamples=timestamp;
        StartAcquisitionSec=timestamp/sampleRate;
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
        all_SCTs=[];
        for k=1:length(all_channels_timestamps)
            if all_channels_info.eventType(k)==3 & all_channels_info.eventId(k)==1 & all_channels_data(k)==2
                corrected_SCT=all_channels_timestamps(k)-StartAcquisitionSec;
                all_SCTs=[all_SCTs corrected_SCT];
            end
        end
        [idx]=find(all_SCTs>Events(sound_index).message_timestamp_sec, 1); %find first SCT after the message timestamp
        SCTtime_sec=all_SCTs(idx);
        %         SCTtime_sec=SCTtime_sec-StartAcquisitionSec; %correct for open-ephys not starting with time zero
        Events(sound_index).soundcard_trigger_timestamp_sec=SCTtime_sec;
        
    end
end

fprintf('\nNumber of sound events (from network messages): %d', length(Events));
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d', length(all_SCTs));
if length(Events) ~=  length(all_SCTs)
    error('ProcessTC_PSTH: Number of sound events (from network messages) does not match Number of hardware triggers (soundcardtrig TTLs)')
end

if exist('check1', 'var') & exist('check2', 'var')
    fprintf('start acquisition method agreement check: %d, %d', check1, check2);
if check1==check2    fprintf(' Good.'), end
end


basefn=sprintf('ch%d_simpleclust_*.t', channel);
d=dir(basefn);
numclusters=size(d, 1);
if numclusters==0 error('ProcessTC_PSTH: no cluster files found');
else fprintf('\n%d cluster files found - will process all of them', numclusters)
end
for clustnum=1:numclusters
    if clustnum<10
        fn=sprintf('ch%d_simpleclust_0%d.t', channel, clustnum);
    else
        fn=sprintf('ch%d_simpleclust_%d.t', channel, clustnum);
    end
    fprintf('\nreading MClust output file %s cluster %d', fn, clustnum)
    spiketimes(clustnum).spiketimes=read_MClust_output(fn)'/10000; %spiketimes now in seconds
    %correct for OE start time, so that time starts at 0
    spiketimes(clustnum).spiketimes=spiketimes(clustnum).spiketimes-StartAcquisitionSec;
    
    totalnumspikes(clustnum)=length(spiketimes(clustnum).spiketimes);
end
fprintf('\nsuccessfully loaded MClust spike data')
Nclusters=numclusters;


monitor = 0;
if monitor
    
    %   I'm running the soundcard trigger (SCT) into ai1 as another sanity check.
    SCTfname=getSCTfile(datadir);
    if isempty(SCTfname)
        warning('could not find soundcard trigger file')
    else
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
    end
    
    %here I'm loading a data channel to get good timestamps - the ADC timestamps are screwed up
    datafname=getContinuousFilename( pwd, 1 );
    [scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(datafname);
    
    SCTtimestamps=SCTtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
    datatimestamps=datatimestamps-StartAcquisitionSec;
    
    
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
    scaledtrace=scaledtrace./max(abs(scaledtrace));
    plot(datatimestamps, SCTtrace) %plotting with data timestamps to work around wierd bug in ADC timestamps
    plot(datatimestamps, scaledtrace, 'm') %plotting with data timestamps to work around wierd bug in ADC timestamps
    
    hold on
    %plot "software trigs" i.e. network messages in red o's
    for i=1:length(Events)
        plot(Events(i).message_timestamp_sec, .25, 'ro');
        plot(Events(i).soundcard_trigger_timestamp_sec, 1, 'g*');
        text(Events(i).message_timestamp_sec, .5, sprintf('network message #%d', i))
        text(Events(i).soundcard_trigger_timestamp_sec, .5, sprintf('SCT #%d', i))
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
        xlim([Events(i).message_timestamp_sec-.02 Events(i).message_timestamp_sec+.5])
        ylim([-5 2])
        pause(.1)
    end
end %if monitor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') ...
            |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise')| strcmp(Events(i).type, 'grating')
        j=j+1;
        alldurs(j)=Events(i).duration;
        if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, '2tone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).frequency;
        elseif strcmp(Events(i).type, 'whitenoise')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=-1;
        elseif strcmp(Events(i).type, 'fmtone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).carrier_frequency;
        elseif strcmp(Events(i).type, 'grating')
            allfreqs(j)=Events(i).angle*1000;
            allamps(j)=Events(i).spatialfrequency;
        end
    end
end
freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
numfreqs=length(freqs);
numamps=length(amps);
numdurs=length(durs);

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
        warning('ProcessTC_LFP: Cannot tell if laser button was turned on in djmaus GUI');
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

M1=[];M1ON=[];M1OFF=[];
nreps=zeros(numfreqs, numamps, numdurs);
nrepsON=zeros(numfreqs, numamps, numdurs);
nrepsOFF=zeros(numfreqs, numamps, numdurs);

%extract the traces into a big matrix M
j=0;
inRange=zeros(1, Nclusters);
for i=1:length(Events)
    if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') | ...
            strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating')
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        end
        laser=LaserTrials(i);
        start=(pos+xlimits(1)*1e-3);
        stop=(pos+xlimits(2)*1e-3);
        if start>0 %(disallow negative or zero start times)
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
                    freq=Events(i).carrier_frequency;
                    amp=Events(i).amplitude;
                case 'whitenoise'
                    freq=-1;
                    amp=Events(i).amplitude;
                case 'grating'
                    amp=Events(i).spatialfrequency;
                    freq=Events(i).angle*1000;
            end
            
            
            dur=Events(i).duration;
            findex= find(freqs==freq);
            aindex= find(amps==amp);
            dindex= find(durs==dur);
            nreps(findex, aindex, dindex)=nreps(findex, aindex, dindex)+1;
            for clust=1:Nclusters %could be multiple clusts (cells) per tetrode
                st=spiketimes(clust).spiketimes; %are in seconds
                spiketimes1=st(st>start & st<stop); % spiketimes in region
                spikecount=length(spiketimes1); % No. of spikes fired in response to this rep of this stim.
                inRange(clust)=inRange(clust)+ spikecount; %accumulate total spikecount in region
                spiketimes1=(spiketimes1-pos)*1000;%covert to ms after tone onset
                spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window
                
                if laser
                    if clust==1
                        nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
                    end
                    M1ON(clust, findex,aindex,dindex, nrepsON(findex, aindex, dindex)).spiketimes=spiketimes1; % Spike times
                    M1ONspikecounts(clust, findex,aindex,dindex,nrepsON(findex, aindex, dindex))=spikecount; % No. of spikes
                    M1spontON(clust, findex,aindex,dindex, nrepsON(findex, aindex, dindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                    M_LaserStart(clust, findex,aindex,dindex, nrepsON(findex, aindex, dindex))=LaserStart(i);
                    M_LaserWidth(clust, findex,aindex,dindex, nrepsON(findex, aindex, dindex))= LaserWidth(i);
                    M_LaserNumPulses(clust, findex,aindex,dindex, nrepsON(findex, aindex, dindex))= LaserNumPulses(i);
                    M_LaserISI(clust, findex,aindex,dindex, nrepsON(findex, aindex, dindex))= LaserISI(i);
                else
                    if clust==1
                        nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
                    end
                    M1OFF(clust, findex,aindex,dindex, nrepsOFF(findex, aindex, dindex)).spiketimes=spiketimes1;
                    M1OFFspikecounts(clust, findex,aindex,dindex,nrepsOFF(findex, aindex, dindex))=spikecount;
                    M1spontOFF(clust, findex,aindex,dindex, nrepsOFF(findex, aindex, dindex))=spont_spikecount;
                end
            end
        end
    end
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))
for clust=1:Nclusters %could be multiple clusts (cells) per tetrode
    fprintf('\ncell %d:', clust)
    fprintf('\ntotal num spikes: %d', length(spiketimes(clust).spiketimes))
    fprintf('\nIn range: %d', inRange(clust))
end


% Accumulate spiketimes across trials, for psth...
for dindex=1:length(durs); % Hardcoded.
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            for clust=1:Nclusters
                
                % on
                spiketimesON=[];
                spikecountsON=[];
                for rep=1:nrepsON(findex, aindex, dindex)
                    spiketimesON=[spiketimesON M1ON(clust, findex, aindex, dindex, rep).spiketimes];
                end
                
                % All spiketimes for a given f/a/d combo, for psth:
                mM1ON(clust, findex, aindex, dindex).spiketimes=spiketimesON;
                
                % off
                spiketimesOFF=[];
                spikecountsOFF=[];
                for rep=1:nrepsOFF(findex, aindex, dindex)
                    spiketimesOFF=[spiketimesOFF M1OFF(clust, findex, aindex, dindex, rep).spiketimes];
                end
                mM1OFF(clust, findex, aindex, dindex).spiketimes=spiketimesOFF;
            end
        end
    end
end

if ~IL %no laser pulses in this file
    mM1ONspikecount=[];
    sM1ONspikecount=[];
    semM1ONspikecount=[];
    M1ONspikecounts=[];
else
    mM1ONspikecount=mean(M1ONspikecounts,5); % Mean spike count
    sM1ONspikecount=std(M1ONspikecounts,[],5); % Std of the above
    for clust=1:Nclusters
        semM1ONspikecount(clust, :,:)=squeeze(sM1ONspikecount(clust, :,:))./sqrt(max(max(nrepsON(:,:,1)))); % Sem of the above
    end
    % Spont
    mM1spontON=mean(M1spontON,5);
    sM1spontON=std(M1spontON,[],5);
    for clust=1:Nclusters
        semM1spontON(clust, :,:)=squeeze(sM1spontON(clust, :,:))./sqrt(max(max((nrepsON(:,:,1)))));
    end
end
if isempty(mM1OFF) %no laser pulses in this file
    mM1OFFspikecount=[];
    sM1OFFspikecount=[];
    semM1OFFspikecount=[];
    M1OFFspikecounts=[];
else
    mM1OFFspikecount=mean(M1OFFspikecounts,5); % Mean spike count
    sM1OFFspikecount=std(M1OFFspikecounts,[],5); % Std of the above
    for clust=1:Nclusters
        semM1OFFspikecount(clust, :,:)=squeeze(sM1OFFspikecount(clust, :,:))./sqrt(max(max(nrepsOFF(:,:,1)))); % Sem of the above
    end
    % Spont
    mM1spontOFF=mean(M1spontOFF,5);
    sM1spontOFF=std(M1spontOFF,[],5);
    for clust=1:Nclusters
        semM1spontOFF(clust, :,:)=squeeze(sM1spontOFF(clust, :,:))./sqrt(max(max(nrepsOFF(:,:,1))));
    end
end

%save to outfiles
%one outfile for each cell
%all previously existing outfiles for this tetrode will be deleted, to
%avoid duplicated or mismatched outfiles when reclustering

%after squeezing cluster, saves with the following dimensions:
% M1ON(findex,aindex,dindex, nrepsON).spiketimes
% mM1ON(findex,aindex,dindex).spiketimes
% mM1ONspikecount(findex,aindex,dindex)

for clust=1:Nclusters
    out.IL=IL;
    out.Nclusters=Nclusters;
    out.tetrode=channel;
    out.channel=channel;
    out.cluster=clust; %there are some redundant names here
    out.cell=clust;
    if IL
        sz=size(M1ON);
        out.M1ON=reshape(M1ON(clust,:,:,:,:), sz(2:end));
        out.mM1ONspikecount=(mM1ONspikecount(clust,:,:,:)); % Mean spikecount for each laser/f/a combo.
        out.sM1ONspikecount=(sM1ONspikecount(clust,:,:,:));
        out.semM1ONspikecount=(semM1ONspikecount(clust,:,:,:));
        %I am not sure if I have this right:
        if numfreqs>numamps
            out.mM1ON=squeeze(mM1ON(clust,:,:,:))';
        else
            out.mM1ON=squeeze(mM1ON(clust,:,:,:));
        end
        out.mM1spontON=(mM1spontON(clust,:,:,:)); % Spont spikes.
        out.sM1spontON=(sM1spontON(clust,:,:,:));
        out.semM1spontON=(semM1spontON(clust,:,:,:));
        out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
        out.LaserWidth=unique(LaserWidth);
        out.LaserNumPulses=unique(LaserNumPulses);
        out.LaserISI=unique(LaserISI);
        
        sz=size(M_LaserStart);
        out.M_LaserStart=reshape(M_LaserStart(clust,:,:,:,:), sz(2:end));
        out.M_LaserWidth=reshape(M_LaserWidth(clust,:,:,:,:), sz(2:end));
        out.M_LaserNumPulses=reshape(M_LaserNumPulses(clust,:,:,:,:), sz(2:end));
        out.M_LaserISI=reshape(M_LaserISI(clust,:,:,:,:), sz(2:end));
        
        
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
    sz=size(M1OFF);
    out.M1OFF=reshape(M1OFF(clust,:,:,:,:), sz(2:end)); % All spiketimes, trial-by-trial.
    %I am not sure if I have this right:
    if numfreqs>numamps
        out.mM1OFF=squeeze(mM1OFF(clust,:,:,:))'; % Accumulated spike times for *all* presentations of each laser/f/a combo.
    else
        out.mM1OFF=squeeze(mM1OFF(clust,:,:,:)); % Accumulated spike times for *all* presentations of each laser/f/a combo.
    end
    out.mM1OFFspikecount=(mM1OFFspikecount(clust,:,:,:));
    out.sM1OFFspikecount=(sM1OFFspikecount(clust,:,:,:));
    out.semM1OFFspikecount=(semM1OFFspikecount(clust,:,:,:));
    out.mM1spontOFF=(mM1spontOFF(clust,:,:,:));
    out.sM1spontOFF=(sM1spontOFF(clust,:,:,:));
    out.semM1spontOFF=(semM1spontOFF(clust,:,:,:));
    out.amps=amps;
    out.freqs=freqs;
    out.durs=durs;
    out.numamps=numamps;
    out.numfreqs=numfreqs;
    out.numdurs=numdurs;
    out.nreps=nreps;
    out.nrepsON=nrepsON;
    out.nrepsOFF=nrepsOFF;
    out.xlimits=xlimits;
    out.samprate=samprate;
    out.datadir=datadir;
    out.spiketimes=spiketimes;
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




