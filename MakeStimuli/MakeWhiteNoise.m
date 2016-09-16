function noise=MakeWhiteNoise(varargin)


% 
% Creates a gaussian white-noise sample with the given parameters
%
%Usage: noise=MakeWhiteNoise(param, samplerate)
%
% Inputs: 
% param should have the following fields:
%  duration           -   duration of the stimulus (ms)
%  amplitude          -   sound pressure level of the sound (dB)
%  ramp               -   length of an ascending/descending edge in ms
%
% samplerate         -   required sampling rate (Hz)
%  
%  Output:
%  noise              -   the requested sample; empty if unsuccessful
%  
global pref

noise=[];
if nargin<2
    return;
end  
  
params=varargin{1};
samplerate=varargin{2};
duration=params.duration;
% attenuation=params.attenuation;
if isfield(params,'amplitude')
    amplitude=params.amplitude;
else
    amplitude=params.attenuation;
end
ramp=params.ramp;
if duration<2*ramp warning('duration too short for requested ramp');end
    
duration=duration/1000;                     % switch to seconds

noise=randn(1,round(duration*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
[edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
noise(1:ledge)=noise(1:ledge).*fliplr(edge);
noise((end-ledge+1):end)=noise((end-ledge+1):end).*edge;

noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V

%amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); 
amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %mw 080107
noise=amplitude.*noise;
