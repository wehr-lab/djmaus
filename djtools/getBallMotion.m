function [moves_trace1, pulses_ms, motion_on_ms, motion_off_ms]=getBallMotion(varargin)

% Read events from raspberry Pi and decode running speed.

% How it works:
%       loads 'all_channels.events' in specified folder; uses pwd if folder is not specified.
%       uses moving median to smooth data.

% Note: pulse with is hardcoded to be 10 ms right now. This value is
% defined in ball_motion.py -> GPIO.PWM rate. If you are not sure what the pulse width is, check
% rapsberry pi, otherwise speed will be wrong. Currently the max speed that
% it can record is 32 cm/sec which has been sufficient for Rig3 set-up, but
% will have to be reevaluated if the set-up is changed. Speed is estimated
% based on 12 inch ball diameter.

% Duty cycles of PWM is 500 Hz. 50% - stationary, >50% - moving
% forward, <50% moving backwards. Pi changes duty cycle every 100ms (delay
% period) 09.18.2017 ira.

smoothN = 20; % smothing over N, where timetrace is N*PulseWidth (10 ms)
pulseWidth = 10; %ms

% cd to directory
if nargin==0
    directory = pwd;
elseif nargin==1
    directory=varargin{1};
end
cd(directory)

% read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);

% sanity check
fid = fopen('temp.txt', 'w');
for i = 1:length(all_channels_timestamps)
    fprintf(fid, '%f, eventType %d, Id %d, channel %d\n', ...
        all_channels_timestamps(i), all_channels_info.eventType(i), ...
        all_channels_info.eventId(i), all_channels_data(i));
end
fclose(fid);

%% get timestamps of on and off events
on_indx = 0;  off_indx = 0; indx = 0;
pChannel = 3;

% figure; hold on; % check events
for i=1:length(all_channels_timestamps)
    
    if all_channels_info.eventType(i) == pChannel && all_channels_info.eventId(i) == 1 && all_channels_data(i)== pChannel % on event
        on_indx = on_indx + 1;
        indx = indx + 1;
        % plot(all_channels_timestamps(i), 1, 'g^')
        motion_on(on_indx) = all_channels_timestamps(i);
        pulses(indx) = all_channels_timestamps(i);
        
    elseif all_channels_info.eventType(i) == pChannel && all_channels_info.eventId(i) == 0 && all_channels_data(i)== pChannel % off event
        off_indx = off_indx + 1;
        indx = indx + 1;
        % plot(all_channels_timestamps(i), 1, 'gv')
        motion_off(off_indx) = all_channels_timestamps(i);
        pulses(indx) = all_channels_timestamps(i);
    end % if event is on the pChannel
    
end % for events

start_record=all_channels_timestamps(1); % first timestamp of recording
motion_on = motion_on - start_record; % subtract to start at 0
motion_off = motion_off - start_record;
pulses = pulses - start_record;

%this needs investigating
if motion_off(2) < motion_on(1)
    motion_off(1)=[];  % drop first off timestamp
    pulses(1)=[];
    fprintf('\n Something went wrong, timestamps are out of order \n')
end
if motion_on(2) < motion_off(1)
    motion_on(1)=[]; % drop first on timestamp
    pulses(1)=[];
    fprintf('\n Something went wrong, timestamps are out of order \n')
end
num_pulses = max(length(motion_on), length(motion_off)); % each pulse has an on and an off event

% convert timestamps to ms
pulses_ms = pulses * 1000;
motion_on_ms = motion_on * 1000;
motion_off_ms = motion_off * 1000;

%% calculate widths
for i = 1 : length(pulses_ms)
    if sum(pulses_ms(i) == motion_off_ms)
        ind(i)=0;
    else
        ind(i)=1;
    end
end
pulse_matrix = [pulses' ind']; % creat a matrix of all pulses and their id

if motion_off(1) < motion_on(1) % start with an on event
    for i = 1 : num_pulses - 2
        M(i) = (motion_on_ms(i) - motion_off_ms(i+1)) / (motion_on_ms(i+1) - motion_on_ms(i));
    end
else
    for i = 1 : num_pulses - 1
        M(i) = (motion_on_ms(i) - motion_off_ms(i)) / (motion_on_ms(i+1) - motion_on_ms(i));
    end
end

%% preprocess the trace, smooth spikes, resample if events were dropped, etimate speed (cm/sec)
M = M - median(M); % center at 0
moves_trace = movmedian(M, smoothN);

% this is not the best way to estimate speed, but it works. Speed is only
% an approximation from regression from known recordings calculated by
% hand.
load('ball_motion_test.mat'); % load test points from ball (add more for accuraccy if needed)
[r,m,b] = regression(X,Y);
moves_trace = b + m * moves_trace; % convert it to cm / sec
ts_error = length(motion_on_ms) - length(moves_trace);

if ts_error ~= 0
    fprintf('\n some on and off events were missed. This resulted in time error of %d ms \n', ts_error*pulseWidth)
end

% not all pulses are 10 ms, the sampling of ball motion is not even
% resample the trace at 10 ms to make sure it aligned with neural data
% (sometimes pi GPIO pins had delays and jitter resulting in speed spikes,
% this smoothies out the spikes)
try
    [scaledtrace, ~, ~] =load_open_ephys_data(sprintf('114_CH%d.continuous', 11));
catch
    [scaledtrace, ~, ~] =load_open_ephys_data(sprintf('114_CH%d.continuous', 11));
end
    
recording_length=length(scaledtrace)/30e3;
trace_length=round(length(scaledtrace)/30e3*100);

t = round(motion_on(1:length(moves_trace))*1000); % irregular timestamps

%interpolate new running_trace now regularly at pulseWidth(10 ms)
[moves_trace1, t2] = resample(moves_trace,t,pulseWidth/100);
%t2 - regular timestamps

% plot the result
figure;
plot(moves1_trace)
ylabel('cm/s')

%%
save('moves_trace1.mat', 'moves_trace1');
save('moves.mat', 'pulses_ms', 'motion_on_ms',  'motion_off_ms',  'pulse_matrix');


