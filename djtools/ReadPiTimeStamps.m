function [pupil_data]=ReadPiTimeStamps(file)
% Reads json .txt file with timepstamps (received during video recording),
%converts them to numbers to sync with recorded neural data
% + 2 triggers: first trigger is start of recording, last trigger is the end of recording
%ira 08.29.17

%% get video timestramps
TS=loadjson(file);

for i=1:length(TS)
    ts0=char(TS{(i)}); %pull timestamps out and convert them to char
    ts(i).ts=datevec(ts0, 'yyyy-mm-dd HH:MM:SS.FFF'); %convert ts to numetic array
    if i==1
        ts(1).ts_zero=0; %subtract starting time and start at zero
    else
        ts(i).ts_zero=etime(ts(i).ts,ts(1).ts);
        
    end
    %sprintf('%.4f',ts(i).ts_zero)
end
duration=etime(ts(end).ts,ts(1).ts);


sprintf('\nVideo recorded for %.4f seconds\n', duration)


%% sync with neural data
sampleRate=30e3;
% read messages, to get timestramps from OE recording first

if exist('Events.mat')
    load('Events.mat')
else
    Eventsfilename='all_channels.events';
    [all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
    messagesfilename='messages.events';
    [messages] = GetNetworkEvents(messagesfilename);
    sound_index=0;
    figure(10)
    hold on
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
            plot(sound_index,SCTtime_sec,'ro')
            
        end
    end
    
    save('Events.mat', 'Events') %save Events so we dont have to do this again
end

if length(ts)-2~=length(Events)
    warning('\nNumber of SCTs recorded on Pi and in OE do not match. Pi=%d, OE=%d',length(ts)-2, length(Events))
end

er=Events(1).soundcard_trigger_timestamp_sec-ts(2).ts_zero;  %difference between when OE started recording and when Pi started recording
for i=1:length(ts)
    ts(i).ts_fixed=ts(i).ts_zero+er;
    ts(i).ts_fixed_samprate=ts(i).ts_fixed*30e3;
    ts(i).ts_zero_samprate=ts(i).ts_zero*30e3;
    if i~=1 && i~=length(ts)
    dif(i)=Events(i-1).soundcard_trigger_timestamp_sec-ts(i).ts_fixed;
    end
end
fprintf('\nMean Difference between OE and Pi SCTs is %.5fsec, error= %.5fsec \n', mean(dif), er)
pupil_file=sprintf('pupil_%s', file);

try
    pupil_size=loadjson(pupil_file);
catch
    error('Could not load pupil size data, %s', pupil_file)
end
x=(0:length(pupil_size)-1)';
f=fit(x, pupil_size', 'pchip');

smooth_data=smooth(pupil_size);
figure;
hold on
plot(x, pupil_size, 'k.')
plot(x, smooth_data, 'r-')
% save('ts.mat','ts')
% save('smooth_data.mat', 'smooth_data')

video_fps=length(smooth_data)/ts(end).ts_zero; %frames per second used for analysis (different from recorded fps

for i=1:length(ts)
    plot(ts(i).ts_zero*video_fps, (mean(smooth_data)), 'go')
end
legend('pupil size', 'smoothed pupil size', 'SCTs')
resampled_data=resample(smooth_data,x, 30000/video_fps);
fprintf( 'Assuming OE recorded at 30k')

figure;
plot((0:length(resampled_data)-1),resampled_data,'b-')
hold on
for i=1:length(ts)
    plot(ts(i).ts_zero_samprate, (mean(resampled_data)), 'go')
end
legend('resampled pupil size', 'SCTs')
pupil_data.ts=ts;
pupil_data.smooth_data=smooth_data;
pupil_data.resampled_data=resampled_data;
pupil_data.pupil_size=pupil_size;
pupil_data.er=er;
pupil_data.dur=duration;

save('pupil_data.mat', 'pupil_data')




