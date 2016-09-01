function ProcessLFP_TC(varargin)
%sorts processed Open Ephys continuous data into a big response matrix
%for LFP, WC, or other continuous data (i.e. does not extract spikes)
%usage: ProcessTC_OE(datafile, [channel], [xlimits], [ylimits])
%channel is a number not a string
%saves output in an outfile
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
fprintf('\nload file 1: ')

%     fprintf('\ntrying to load %s...', datafile)
%     godatadir(expdate, session, filenum)
%     pathname='C:\Program Files\Open Ephys\ira_2014-02-14_13-43-56\'
%     cd(pathname)
%     filename=sprintf('102_CH%s.continuous',channel)

%try to read OE filename from exper structure (only will work after
%02.14.14)

    cd(pref.datapath)
    cd(datadir)
filename=sprintf('102_CH%d.continuous', channel);
if exist(filename, 'file')~=2 %couldn't find it
    filename=sprintf('100_CH%d.continuous', channel);
end
if exist(filename, 'file')~=2 %couldn't find it
    error(sprintf('could not find data file %s in datadir %s', filename, datadir))
end

 
 
  [scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(filename);

  %should record a copy of stim
  stim=0*scaledtrace;

%   I'm running the soundcard trigger (SCT) into ai1 as another sanity check.
[SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data('100_ADC1.continuous');

    %might want to read settings.xml from this directory to grab info about
    %filter settings, etc.
load notebook.mat

%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

%all_channels_data (I think) is the channel the Events is associated with
% the 0 channel Events are network Events?
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
% fid=fopen('temp.txt', 'w');
% for i=1:length(all_channels_timestamps)
%     fprintf(fid, '%f, eventType %d, Id %d, channel %d\n', ...
%         all_channels_timestamps(i), all_channels_info.eventType(i), ...
%         all_channels_info.eventId(i), all_channels_data(i));
% end
% fclose(fid);
    
    
%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);


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
            Events(sound_index).message_timestamp_samples=timestamp - StartAcquisitionSamples;
            Events(sound_index).message_timestamp_sec=timestamp/sampleRate - StartAcquisitionSec;
        end
        
        %get corresponding SCT TTL timestamp and assign to Event
        all_SCTs=[];
        for k=1:length(all_channels_timestamps)
            if all_channels_info.eventType(k)==3 & all_channels_info.eventId(k)==1 & all_channels_data(k)==2
                all_SCTs=[all_SCTs all_channels_timestamps(k)];
            end
        end
        [idx]=find(all_SCTs>Events(sound_index).message_timestamp_sec, 1); %find first SCT after the message timestamp
        SCTtime_sec=all_SCTs(idx);
        SCTtime_sec=SCTtime_sec-StartAcquisitionSec; %correct for open-ephys not starting with time zero
        Events(sound_index).soundcard_trigger_timestamp_sec=SCTtime_sec;
        
    end
end

if exist('check1', 'var') & exist('check2', 'var')
    fprintf('start acquisition method agreement check: %d, %d', check1, check2);
end
fprintf('\nreminder to cross-check network mesasges with stimlog')
            
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





monitor = 0;
if monitor
    figure
    set(gcf, 'pos', [-1853 555 1818 420]);
    plot(SCTtimestamps, SCTtrace)
    hold on
    %plot "software trigs" i.e. network messages in red o's
    for i=1:length(Events)
        plot(Events(i).message_timestamp_sec, 0, 'ro');
        plot(Events(i).soundcard_trigger_timestamp_sec, 1e4, 'g*');
        text(Events(i).message_timestamp_sec, 1e4, sprintf('network message #%d', i))
        text(Events(i).soundcard_trigger_timestamp_sec, .5e4, sprintf('SCT #%d', i))
    end
    
    %all_channels_info.eventType(i) = 3 for digital line in (TTL), 5 for network Events
    %all_channels_info.eventId = 1 for rising edge, 0 for falling edge
    %all_channels_data(i) is the digital input line channel
    
    
    % plot TTL SCTs in green ^=on, v=off
    for i=1:length(all_channels_timestamps)
        if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==2
            plot(all_channels_timestamps(i), 1e3, 'g^')
            text(all_channels_timestamps(i), 1e4, 'TTL on/off')
        elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==2
            plot(all_channels_timestamps(i), 1e3, 'gv')
        end
    end
    
    plot(datatimestamps, scaledtrace, 'm')
    
    for i=1:length(Events)
        xlim([Events(i).message_timestamp_sec-.02 Events(i).message_timestamp_sec+.08])
        ylim([-4 2]*1e4)
        pause(.5)
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
%for now, setting AOPulseOn to 0 for all events
for i=1:length(Events)
    if 0%isfield(Events(i), 'AOPulseOn')
        aopulseon(i)=Events(i).AOPulseOn;
    else
        aopulseon(i)=0;
        Events(i).AOPulseOn=0;
    end
end
fprintf('\n%d laser pulses in this Events file', sum(aopulseon))
if sum(aopulseon)>0
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
for i=1:length(Events)
    if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') | ...
            strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating')
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate; %pos is in samples
        end
        aopulseon=Events(i).AOPulseOn;
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        if isempty(find(region<1)) %(disallow negative or zero start times)
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
                M1(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=scaledtrace(region);
                M1stim(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=stim(region);
                if aopulseon
                    nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
                    M1ON(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=scaledtrace(region);
                else
                    nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
                    M1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=scaledtrace(region);
                end
%             end
        end
    end
end

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
                else %no reps for this stim, since rep=0
                    mM1ON(findex, aindex, dindex,:)=zeros(size(region));
                end
                 if nrepsOFF(findex, aindex, dindex)>0
                    mM1OFF(findex, aindex, dindex,:)=mean(M1OFF(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1OFF(findex, aindex, dindex,:)=zeros(size(region));
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

%%
% added by ira 09-17-13
% this way there is only one outfile that can be used for PlotTC and
% PlotTC_psth
% high_pass_cutoff=300; %Hz
%     fprintf('\nhigh-pass filtering at %d Hz', high_pass_cutoff);
%     [b,a]=butter(1, high_pass_cutoff/(samprate/2), 'high');
%     filteredtrace=filtfilt(b,a,scaledtrace);
%             
%             nstd=thresh/std(filteredtrace);
%             fprintf('\nusing absolute spike detection threshold of %.2f mV (%.2f sd)', thresh, nstd);
%         thresh=nstd*std(filteredtrace);
%         if thresh>1
%             fprintf('\nusing spike detection threshold of %.2f mV (%.2f sd)', thresh, nstd);
%         elseif thresh<=1
%             fprintf('\nusing spike detection threshold of %.2f mV (%.2f sd)', thresh, nstd);
%         end
%         refract=15;
%     fprintf('\nusing refractory period of %.1f ms (%d samples)', 1000*refract/samprate,refract );
%     spikes=find(abs(filteredtrace)>thresh);
%     dspikes=spikes(1+find(diff(spikes)>refract));
%     try dspikes=[spikes(1) dspikes'];
%     catch
%         fprintf('\n\ndspikes is empty; either the cell never spiked or the nstd is set too high\n. Ingore if plotting LFPs');
%         
%     end
    

%%
%assign outputs
out.scaledtrace=scaledtrace;
out.M1=M1;
out.M1ON=M1ON;
out.M1OFF=M1OFF;
out.mM1ON=mM1ON;
out.mM1OFF=mM1OFF;
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
out.xlimits=xlimits;
out.ylimits=ylimits;
out.samprate=samprate;
% out.nstd=nstd;
out.channel=channel;
%should probably save header info and stuff like that


outfilename=sprintf('outLFP_ch%d',channel);
save(outfilename, 'out')
fprintf('\n saved to %s', outfilename)

