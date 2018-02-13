function [pulses_ms,puls_ms, motion_on_ms, motion_off_ms]=getBallMotion(varargin)

%leads'all_channels.events'in current folder (or specified folder)
%finds ball movement pwm signal and calculates the velocity
%uses 'rlowess' to smooth data - A robust version of 'lowess' that assigns lower weight to outliers in the regression.
%                 The method assigns zero weight to data outside six mean absolute deviations.
% 'lowess'- Local regression using weighted linear least squares and a 1st degree polynomial model

% 09.18.2017 duty cycles of PWM is 500 Hz. 50%- stationary, >50%- moving
% forward, <50% moving backwards. Pi changes duty cycle every 100ms (delay
% period)

%djPrefs;
%global pref
if nargin==0
    fprintf('\nno input\n')
else
    directory=varargin{1};
end
if nargin==1
    directory=varargin{1};
    cd(directory)
end
cd(directory)
%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

ch=3; %

%sanity check
fid=fopen('temp.txt', 'w');
for i=1:length(all_channels_timestamps)
    fprintf(fid, '%f, eventType %d, Id %d, channel %d\n', ...
        all_channels_timestamps(i), all_channels_info.eventType(i), ...
        all_channels_info.eventId(i), all_channels_data(i));
end
fclose(fid);

pulse_length=0.010; %s, assuming freq= 100 Hz
on=0;off=0;
j=0;
%eventType=3, eventId=1 On and 0 Off, channel =3 (0 start)
for i=1:length(all_channels_timestamps)
    if all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==1 & all_channels_data(i)==3
        on=on+1;
        j=j+1;
        %plot(all_channels_timestamps(i), 1, 'g^')
        motion_on(on)=all_channels_timestamps(i);
        pulses(j)=all_channels_timestamps(i);
    elseif all_channels_info.eventType(i)==3 & all_channels_info.eventId(i)==0 & all_channels_data(i)==3
        off=off+1;
        j=j+1;
        %plot(all_channels_timestamps(i), 1, 'gv')
        motion_off(off)=all_channels_timestamps(i);
        pulses(j)=all_channels_timestamps(i);
    end
end
start_record=all_channels_timestamps(1);%getFirstTimestamp;
motion_on=motion_on-start_record; %start at 0
motion_off=motion_off-start_record;
pulses=pulses-start_record;
    
if motion_off(2)<motion_on(1)
    motion_off(1)=[]; %this needs investigating
    pulses(1)=[];
end
if motion_on(2)<motion_off(1)
    motion_on(1)=[]; %this needs investigating
    pulses(1)=[];
end
num_pulses=max(length(motion_on), length(motion_off));
j=0;

pulses_ms=pulses*1000;
motion_on_ms=motion_on*1000;
motion_off_ms=motion_off*1000;
for i=1:length(pulses_ms)
    if sum(pulses_ms(i)==motion_off_ms)
        ind(i)=0;
    else
        ind(i)=1;
    end
end
puls_ms=[pulses_ms' ind']; %creat a matrix of all pulses and their id



% if motion_off(1)<motion_on(1)
%     moves=pul*-1; %flip sign so that positive # would mean forward
% else
%     moves=pul;
% end
% moves=moves-mean(moves); %center at 0
% moves=(moves*6.4)*1000; %max speed of the ball is about 32 cm/s
%figure; plot(1:length(velocity), velocity)
% tic
% moves=smooth(moves, 'rlowess');
% toc
% figure; plot(1:length(moves), moves)
% yticks1=[-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5]*6.4;
% yticks(yticks1)
% for i=1:length(yticks1)
% yticklabel1{i}=sprintf('%.1f', yticks1(i));
% end
% yticklabels(yticklabel1);
% ylabel('cm/s')
figure; plot(diff(pulses_ms));
save('moves.mat','pulses_ms','puls_ms', 'motion_on_ms', 'motion_off_ms')








