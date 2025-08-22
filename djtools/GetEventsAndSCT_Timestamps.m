function [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog)
global pref
if isempty(pref); djPrefs;end
% Here are some general notes on the format of Events and messages
%
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
er=0; % manually fix SCT if there are too many
% get all SCT timestamps
sound_index=0;

cont_files=dir('*.continuous');
[~, timestamps, ~] =load_open_ephys_data(cont_files(1).name); %grab any continuous file in the recording
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
            if all_channels_info.eventType(k)==3 & all_channels_info.eventId(k)==1 & all_channels_data(k)==pref.SCT_digital_line_in
                corrected_SCT=all_channels_timestamps(k)-StartAcquisitionSec;
                all_SCTs=[all_SCTs corrected_SCT];
            end
        end
%         figure; plot(diff(all_SCTs), 'ko')
%         if er
%             save('all_SCTs.mat','all_SCTs')
%         end
        %get corresponding SCT TTL timestamp and assign to Event
        %old way (from TC) won't work, since the network events get ahead of the SCTs.
        %another way (dumb and brittle) is to just use the corresponding
        %index, assuming they occur in the proper order with no drops or
        %extras
        try
            if er
            load('all_SCTs.mat')
            end
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
%  figure; plot(diff(all_SCTs), 'ko')
if length(Events) ~=  length(all_SCTs)
     warning('GetEventsAndSCT_Timestamps: Number of sound events (from network messages) does not match Number of hardware triggers (soundcardtrig TTLs)')
   % THERE_IS_A_PROBLEM
   fprintf('\n%d Events but %d SCTs, calling ResolveEventMismatch...',length(Events),length(all_SCTs)) 
warning('\ncommenting out ResolveEventMismatch this is a HACK and needs to be investigated!!!!') 

%     [Events, all_SCTs, stimlog]=ResolveEventMismatch(Events, all_SCTs, stimlog);
end

if exist('check1', 'var') & exist('check2', 'var')
    fprintf('start acquisition method agreement check: %d, %d', check1, check2);
end

fprintf('\nsuccessfully processed network messages and soundcard triggers into Events.')
fprintf('\nNumber of sound events (from network messages): %d', length(Events));
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d', length(all_SCTs));
