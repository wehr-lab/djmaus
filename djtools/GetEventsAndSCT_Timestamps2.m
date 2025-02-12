function [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps2(BonsaiPath)
global pref
if isempty(pref); djPrefs;end

%updated to work with new OE file formats and file hierarchy with version 0.6 and open-ephys-matlab-tools
% -mike 9.2023

%I stored OE info in Sky, but sometimes we want to run without a camera, so
%I'm changin it to pull from OEinfo instead
%-mike 2.12.25'
[~,BonsaiFolder,~]=fileparts(BonsaiPath);
OEinfofilename=sprintf('OEinfo-%s', BonsaiFolder);
load(OEinfofilename)

try
    if startsWith(OEversion, 'old')
        error('GetEventsAndSCT_Timestamps2: this data is from old version of open ephys, use GetEventsAndSCT_Timestamps instead')
    elseif startsWith(OEversion, '0.6')
        % OK, do nothing
        fprintf('\nOE %s', OEversion)
    else
        OEversion=split(OEversion, '.');
        if str2num(OEversion{1})==0 & str2num(OEversion{2})>=6
            % OK, do nothing
            fprintf('\nOE %s', OEversion)
        else
            error ('GetEventsAndSCT_Timestamps2: is this data from OpenEphys 0.6 or greater? I cant tell')
        end
    end
catch
    fprintf('GetEventsAndSCT_Timestamps2: is this data from OpenEphys 0.6 or greater? I cant tell');
end

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

%messages=Sky.messages;
%TTL=Sky.TTL;

% this works with OE 06.4, but you have to locate the continuous files:
% cont_files=dir('*.continuous');
% [~, timestamps, ~] =load_open_ephys_data(cont_files(1).name); %grab any continuous file in the recording
% StartAcquisitionSamples=timestamps(1)*30e3;

StartAcquisitionSamples=messages{1, 'sample_number'};
StartAcquisitionSec=messages{1, 'timestamps'};
if strcmp(messages{1, 'text'}, 'GetRecordingPath')
    %good, this is what we expect in messages(1,:). Do nothing.
elseif startsWith(messages{1, 'text'}, 'StartRecord')
    %also good, this is another valid possibility we expect in messages(1,:). Do nothing.
elseif startsWith(messages{1, 'text'}, 'StartAcquisition')
    %I think this is also OK, not sure why it didn't start occuring until 8.10.23 but should be fine. Do nothing.
else
    error('GetEventsAndSCT_Timestamps2: unexpected first messsage text. Not sure how to proceed.')
end

%I had no idea how to get actual values for these in OE 0.6.4, but figured
%it out.
%using 0 is wrong. But you can use messages(1,:) as confirmed by the
%following check:
%for    '/Volumes/Projects/AldisLateralSeptumPhys/2023-07-19_11-00-18_mouse-1934/2023-07-19_11-00-20mouse-1934/Record Node 109'
%the first timestamp is 1.925708800000000e+03
%and the first sample is 57771264
% and messages(1,:) is the same
% timestamps    sample_number           text       
%     __________    _____________    __________________
% 
%     1925.7088       57771264       "GetRecordingPath"

all_SCTs=[];
for k=1:height(TTL)
    if TTL.state(k) ... %is true, means rising
            & TTL.line(k)==pref.SCT_digital_line_in %should be line 2
        corrected_SCT=TTL.timestamp(k)-StartAcquisitionSec;
        all_SCTs=[all_SCTs corrected_SCT];
    end
end

%Here's where we can recover from a situation where the TTLs were not
%recorded (e.g. the BNC only went to an ADC channel)
if ~height(TTL) %this means TTL is empty (although it doesn't show up as empty(TTL) bc it's a dataframe)
    str=sprintf('\n\n%s\n\n%s\n\n%s\n\n', '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!', ...
        'TTLs are empty. Maybe the soundcard triggers weren''t plugged into the Digital In channel? Attempting to recover soundcard triggers from messages and/or analog soundcardtrig', ...
        '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    warning(str)
end

for i=1:height(messages)
    
    timestamp=messages.timestamps(i);
    sample_number=messages.sample_number(i);
    str2=split(messages.text(i));
    
    Events_type=str2{1};
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
        
        
        Events(sound_index).type=str2{2};
        for j=3:length(str2)
            str3=strsplit(str2{j}, ':');
            fieldname=str3{1};
            value=str2num(str3{2});
            if isempty(value) %this happens if it's a real string (like soaflag), not a 'num'
                value=str3{2}; %just use the string. E.g. 'soa' or 'isi'
            end
            Events(sound_index).(fieldname)= value;
        end
        Events(sound_index).message_timestamp_samples=sample_number - StartAcquisitionSamples;
        Events(sound_index).message_timestamp_sec=timestamp - StartAcquisitionSec;
        if ~isempty(all_SCTs) %the normal case, we have digital TTL soundcard triggers
        SCTtime_sec=all_SCTs(sound_index);
        else
            %the digital TTL soundcard triggers are missing, so we try to
            %recover using the message timestamp (very suboptimal!)
            %-mike 04.25.2024
            SCTtime_sec=Events(sound_index).message_timestamp_sec;
            Events(sound_index).soundcardtrig_is_missing='using message timestamp as workaround';
        end
        Events(sound_index).soundcard_trigger_timestamp_sec=SCTtime_sec;



%         figure; plot(diff(all_SCTs), 'ko')
%         if er
%             save('all_SCTs.mat','all_SCTs')
%         end
        %get corresponding SCT TTL timestamp and assign to Event
        %old way (from TC) won't work, since the network events get ahead of the SCTs.
        %another way (dumb and brittle) is to just use the corresponding
        %index, assuming they occur in the proper order with no drops or
        %extras
%         try
%             if er
%             load('all_SCTs.mat')
%             end
%             SCTtime_sec=all_SCTs(sound_index);
%             Events(sound_index).soundcard_trigger_timestamp_sec=SCTtime_sec;
%         catch
%             %     one reason this won't work is if you stop recording during the
%             %     play-out, i.e. stimuli were logged (to stimlog and network events)
%             %     but never played before recording was stopped
%             if sound_index>length(all_SCTs)
%                 Events(sound_index).soundcard_trigger_timestamp_sec=nan;
%             end
%         end
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
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d\n', length(all_SCTs));
