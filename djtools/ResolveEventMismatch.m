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
            %we're fine, because extra SCTs before any message is automatically ignored
            fprintf('\nMismatched resolved. Extra (initial) soundcard trigger discarded')
        elseif  all_SCTs(end) > Events(end).message_timestamp_sec
            warning(sprintf('\nMismatch possibly resolved by discarding extra (last) soundcard trigger. \nFurther investigation is highly recommended!!!!!!!'))
        else
            error('ResolveEventMismatch: this case is not handled yet')
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
            otherwise
                fprintf('\n maybe user stopped open-ephys recording manually but djmaus was still running')
                fprintf('\n investigate, and if this appears to be so, add a special case')
                error('ResolveEventMismatch: this case is not handled yet')
        end
        
    else
        error('ResolveEventMismatch: this case is not handled yet')
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
    
