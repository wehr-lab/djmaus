function sweep=MakeFMSweep(varargin)

% Note: now we use absolute SPL instead of attenuation

% NOTE: this was the previous header
% function sweep=MakeFMSweep(start_frequency,stop_frequency,speed,duration,loglin,attenuation,samplerate,ramp)

global exper pref

% creates an FM sweep starting from start_frequency (Hz) to stop_frequency
% (Hz). Start>stop for downward sweeps and vice versa. Either speed (oct/s) or duration (ms) must be supplied, if both
% then duration is used.
% Input:
%  start_frequency            -   start frequency (Hz)
%  stop_frequency             -   stop frequency (Hz)
%  speed                      -   sweep speed (oct/s)
%  duration                   -   duration of the sweep (ms)
%  loglin                     -   1=logarithmic sweep
%                                 0=linear sweep
% % % % % % % % % % % % % % % % %  attenuation                -   attenuation (dB)
%  amplitude                  -   sound pressure level of the sound (dB)
%  samplerate                 -   sampling rate (Hz)
%  ramp                       -   length of ascending/descending ramp (ms)
%  
%  Output:
%  sweep                      -   sweep of the requested parameters; empty if unsuccessful
%  
 
sweep=[];
if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
start_frequency=params.start_frequency;
stop_frequency=params.stop_frequency;
speed=params.speed;
duration=params.duration;
loglin=params.loglin;
% attenuation=params.attenuation;
amplitude=params.amplitude;
ramp=params.ramp;

if isempty(speed) & isempty(duration)
    return;
end

if stop_frequency<start_frequency       % downward sweep
    down=1;
    x=stop_frequency;
    start_frequency=stop_frequency;
    stop_frequency=x;
else
    down=0;
end
     
if isempty(duration)
    % compute the duration given the speed
    octaves=log2(stop_frequency)-log2(start_frequency);
    duration=abs(octaves/speed);        % duration is now in s
    duration=round(duration*1000)/1000; % rounded to ms
else    
    % we have the duration
    duration=duration/1000;
end

% amplitude=10*10.^(-attenuation/20);
amplitude=10*(10.^((amplitude-pref.maxSPL)/20));

t=0:1/samplerate:duration;
if loglin
    style='logarithmic';
else
    style='linear';
end

sweep=amplitude*chirp(t,start_frequency,duration,stop_frequency,style);

if down
    sweep=fliplr(sweep);
end

[edge,ledge]=MakeEdge(ramp,samplerate);
sweep(1:ledge)=sweep(1:ledge).*fliplr(edge);
sweep((end-ledge+1):end)=sweep((end-ledge+1):end).*edge;

