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
    warning('ResolveEventMismatch: number of logged stim does not match number of Events')
    %try to match them up
    %situation: AsymGPIAS stimlog missed the first stimulus
    switch stimlog(1).type
        case 'AsymGPIAS'
            for i=1:10
                events(i,:)=[Events(i).gapdur Events(i).onramp Events(i).offramp];
                stim(i,:)=[stimlog(i).param.gapdur stimlog(i).param.onramp stimlog(i).param.offramp];
            end
            C=corr(stim', events');
            offset=find(C(1,:)==1)-1; %num stimuli that stimlog is missing
            if offset>=1
                for i=2:length(Events)
                    temp(i)=stimlog(i-1);
                end
                stimlog=temp;
            end
            %implicitly, stimlog(1) has the fields but they are all []
            warning(sprintf('\nLooks like the first %d stimuli were not logged to notebook. Mismatch possibly resolved by shifting stimlog by %d. \nFurther investigation is highly recommended!!!!!!!'), offset, offset)
            return
    end
    if nEvents~=nstimlog & nEvents==nSCTs
        nMismatch = nEvents-nstimlog;
        if Events(1).message_timestamp_sec>Events(nMismatch+1).message_timestamp_sec
            %this is a specific case where (we think) Kat might have pressed record in
            %OE before clicking record (again) in djmaus. Must investigate further.
            %the first nMismatch events are spurious (or at least don't match stimlog) and have incorrect timestamps
            Events=Events(nMismatch+1:end);
            all_SCTs=all_SCTs(nMismatch+1:end);
            warning(sprintf('\nLooks like the first nMismatch stimuli were not logged to notebook. Mismatch possibly resolved by discarding first 2 events. \nFurther investigation is highly recommended!!!!!!!'))
            return
        end
        %a more general approach would be to look where Events.timestamps
        %get restarted. We think the first n are timed relative to clicking
        %OE record button, then n+1 onwards are relative to the
        %StartAcquisition signal that happens shortly afterwards when the user pressed record
        %in djmaus. Finding this reset point should be doable.
    end
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
        
        %old comment:
        %we're fine, because extra SCTs before any message is automatically ignored
        
        %update (2-8-2017)
        %I don't think we can ignore it. removing it is required. not sure
        %what is different about this case (/Volumes/d/lab/djmaus/Data/kip/2017-01-09_13-07-45_mouse-3B13A)
        
        %strip out first SCT which is extraneous
        all_SCTs=all_SCTs(2:end);
        for i=1:length(all_SCTs)
            Events(i).soundcard_trigger_timestamp_sec=all_SCTs(i);
        end
        fprintf('\nMismatched resolved. Extra (initial) soundcard trigger discarded')
    elseif  all_SCTs(end) - Events(end).message_timestamp_sec > 0.1
        warning(sprintf('\nMismatch possibly resolved by discarding extra (last) soundcard trigger. \nFurther investigation is highly recommended!!!!!!!'))
    else
        check_for_timestamp_discontinuities
        error('ResolveEventMismatch: this case is not handled yet')
%%%%%%%%%Special exception added by nick
%         for i=1:nEvents
%         temp(i) = Events(i).soundcard_trigger_timestamp_sec;
%         end
% 
%         temp = temp(1,1:167);
%         temp(1,168:180) = all_SCTs(1,169:181);
% 
%         for i = 1:length(Events)
%         Events(i).soundcard_trigger_timestamp_sec = temp(1,i);
%         end
%         all_SCTs = temp;
    end
elseif strcmp(stimlog(1).protocol_name(1:5), 'GPIAS')
    %     one problem with GPIAS stimuli is if you stop recording during the
    %     play-out, i.e. stimuli were logged (to stimlog and network events)
    %     but never played before recording was stopped
    if length(Events)>length(all_SCTs)
        warning(sprintf('\nIt looks like Recording was stopped before all GPIAS stimuli finished playing. Mismatch possibly resolved by discarding extra Events. \nFurther investigation is highly recommended!!!!!!!'))
        Events=Events(1:length(all_SCTs));
    end
elseif  nstimlog==nEvents & nSCTs== nEvents -1
    if strcmp(stimlog(1).protocol_name, 'ArchPVRev1')
        %initially with this stimulus type we had some problems with java heap memory and file not
        %finishing or something like that?
        warning(sprintf('\nMismatch possibly resolved by discarding extra (last) Event. \nFurther investigation is highly recommended!!!!!!!'))
        Events=Events(1:nSCTs);
        stimlog=stimlog(1:nSCTs);
    end
elseif nstimlog==nEvents & nSCTs < nEvents
    %maybe user stopped open-ephys recording manually but djmaus was still running
    %this happened spuriously at least once - when Aldis came in to check,
    %open-ephys had stopped itself but djmaus was still running. Not sure
    %how this could happen, but investigation confirms everything is fine up until the last soundcard trigger
    
    %adding special case
    [p,n,e]=fileparts(pwd);
    switch n
        case '2016-11-10_09-14-26_mouse-7093'
            Events=Events(1:nSCTs);
            stimlog=stimlog(1:nSCTs);
        case '2017-08-16_11-01-22_mouse-7765'
            Events=Events(1:nSCTs);
            stimlog=stimlog(1:nSCTs);
        otherwise
            warning('There''s a problem.')
            fprintf('\nnumber of soundcardtriggers is less than the number of Events & stimuli in stimlog')
            fprintf('\n maybe user stopped open-ephys recording manually but djmaus was still running')
            fprintf('\n investigate, and if this appears to be so, add a special case')
            check_for_timestamp_discontinuities
            error('ResolveEventMismatch: this case is not handled yet')
    end
elseif nstimlog==nEvents & nSCTs > nEvents
    for i=1:nEvents
        temp(i) = Events(i).soundcard_trigger_timestamp_sec;
    end
    [~,ind] = intersect(all_SCTs,temp);
    all_SCTs = all_SCTs(ind);
    nSCTs = length(ind);
    
else
    check_for_timestamp_discontinuities
    error('ResolveEventMismatch: this case is not handled yet')
end


function check_for_timestamp_discontinuities
        warning('checking for timestamp discontinuities...')
    SCTfname=getSCTfile(pwd);
    if isempty(SCTfname)
        warning('could not find soundcard trigger file')
    else
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
            sampleRate=SCTinfo.header.sampleRate; %in Hz
    end
        % are there discontinuities in the timestamps?
      d=  (sampleRate*unique(diff(SCTtimestamps)));
    fprintf('\nhere are the unique values of dt (in units of 1/sampleRate) ')
    fprintf('\nthese should all be 1.0000 to within machine precision ')
    fprintf('\n%.4f', d)
    fprintf('\n')
    fprintf('\nif there are some values >> 1.0000, this means there are discontinuities ')
    fprintf('\nin the timestamps. This is fatal (non-recoverable) and indicates drop-outs by ')
    fprintf('\nopen-ephys. On rigs where this happens a lot, we have (mostly?) solved it by ')
    fprintf('\nswitching to a 2-computer set-up. \n\n')
    if any(d>1.1)
        error('ResolveEventMismatch: there apppear to be discontinuities in the timestamps. This is fatal (non-recoverable) and indicates drop-outs by open-ephys. See above message');
    else
        fprintf('\nThere don''t appear to be any timestamp discontinuities \n\n')

    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%if you want to investigate the situation, here are some things to try:
if 0
    
    
    
    
    %look for StartAcquisitionSec = xxx in command window and execute to set StartAcquisitionSec
    SCTfname=getSCTfile(pwd);
    if isempty(SCTfname)
        warning('could not find soundcard trigger file')
    else
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
    end
    
    %all this is just to get StartAcquisitionSec
    messagesfilename='messages.events';
    [messages] = GetNetworkEvents(messagesfilename);
    Eventsfilename='all_channels.events';
    [all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
    sampleRate=all_channels_info.header.sampleRate; %in Hz
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
            break
        elseif strcmp(deblank(Events_type), 'Software')
            StartAcquisitionSamples=timestamp;
            StartAcquisitionSec=timestamp/sampleRate;
            fprintf('\nStartAcquisitionSec=%g', StartAcquisitionSec)
            break
        end
    end
        
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
        text(Events(i).message_timestamp_sec, .35, sprintf('network message #%d', i),'color', 'r')
        text(Events(i).soundcard_trigger_timestamp_sec, .95, sprintf('SCT #%d', i), 'color','g')
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

