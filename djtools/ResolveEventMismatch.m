function [ Events, all_SCTs, stimlog ] = ResolveEventMismatch(Events, all_SCTs, stimlog  )
%try to resolve the situation when the number of sound events (from network
%messages) does not match Number of hardware triggers (soundcardtrig TTLs

nstimlog=length(stimlog);
nEvents=length(Events);
nSCTs=length(all_SCTs);
fprintf('\nNumber of stimuli logged to notebook (stimlog): %d', nstimlog);
fprintf('\nNumber of sound events (from network messages): %d', nEvents);
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d', nSCTs);

if nstimlog~=nEvents
    error('ResolveEventMismatch: number of logged stim does not match number of Events')
end

%first double-check that the stimlog and Events match
if nstimlog==nEvents
    for i=1:nstimlog
        if ~strcmp(stimlog(i).type, Events(i).type)
            error('ResolveEventMismatch: stimlog does not match Events')
        end
        if ~stimlog(i).param.duration==Events(i).duration
            error('ResolveEventMismatch: stimlog does not match Events')
        end
        if strcmp(stimlog(i).type,'tone')
            if ~stimlog(i).param.frequency==Events(i).frequency
                error('ResolveEventMismatch: stimlog does not match Events')
            end
            if ~stimlog(i).param.amplitude==Events(i).amplitude
                error('ResolveEventMismatch: stimlog does not match Events')
            end
        end
    end
end
fprintf('\nGood. Logged stimuli match network Events.')
    
if nstimlog==nEvents & nSCTs==1+ nEvents
    %there's one extra soundcard trigger, find it and remove it
    
    %if the extra one comes before any events, don't need to do anything
    if all_SCTs(1) < Events(1).message_timestamp_sec
        %we're fine, because extra SCTs before any message is automatically ignored
        fprintf('\nMismatched resolved. Extra (initial) soundcard trigger discarded') 
    elseif  all_SCTs(end) > Events(end).message_timestamp_sec
        warning(sprintf('\nMismatched possibly resolved by discarding extra (last) soundcard trigger. \nFurther investigation is highly recommended!!!!!!!')) 
    else
        error('ResolveEventMismatch: this case is not handled yet')
    end
else
    error('ResolveEventMismatch: this case is not handled yet')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if you want to investigate the situation, here are some things to try:
if 0
    SCTfname=getSCTfile(pwd);
    if isempty(SCTfname)
        warning('could not find soundcard trigger file')
    else
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
    end
    
    %here I'm loading a data channel to get good timestamps - the ADC timestamps are screwed up
    %datafname=getContinuousFilename( pwd, 1 );
    %[scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(datafname);
    
    SCTtimestamps=SCTtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
    %datatimestamps=datatimestamps-StartAcquisitionSec;
    
     figure
    hold on
    %set(gcf, 'pos', [-1853 555 1818 420]);
    SCTtrace=SCTtrace./max(abs(SCTtrace));
    %scaledtrace=scaledtrace./max(abs(scaledtrace));
    plot(SCTtimestamps, SCTtrace) %plotting with data timestamps to work around wierd bug in ADC timestamps
    %plot(datatimestamps, scaledtrace, 'm') %plotting with data timestamps to work around wierd bug in ADC timestamps
    
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
%     for i=1:length(all_channels_timestamps)
%         if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==2
%             plot(all_channels_timestamps(i), 1, 'g^')
%             text(all_channels_timestamps(i), 1, 'TTL on/off')
%         elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==2
%             plot(all_channels_timestamps(i), 1, 'gv')
%         end
%     end
    
    
%    for i=1:length(Events)
numEvents=length(Events);
for i=numEvents-5:numEvents
        xlim([Events(i).message_timestamp_sec-.02 Events(i).message_timestamp_sec+2])
        ylim([-5 2])
        drawnow
        pause(.1)
    end
end

