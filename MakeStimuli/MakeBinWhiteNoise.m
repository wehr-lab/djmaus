function noise=MakeBinWhiteNoise(varargin)

global pref

% Note: now we use absolute SPL instead of attenuation!!!

% NOTE: this was the original header
% function noise=MakeWhiteNoise(duration,attenuation,samplerate,ramp)

% creates a gaussian white-noise sample with the given parameters
% Input:
%  duration           -   duration of the stimulus (ms)
% % % % % % % % % % % %  attenuation        -   attenuation (dB)
%  amplitude          -   sound pressure level of the sound (dB)
%  samplerate         -   required sampling rate (Hz)
%  ramp               -   length of an ascending/descending edge MS???
%  
%  Output:
%  noise              -   the requested sample; empty if unsuccessful
%  

noise=[];
if nargin<2
    return;
end  
  
params=varargin{1};
samplerate=varargin{2};
duration=params.duration;
% attenuation=params.attenuation;
Ramplitude=params.Ramplitude;
Lamplitude=params.Lamplitude;
ramp=params.ramp;

Ramplitude=1*(10.^((Ramplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
Lamplitude=1*(10.^((Lamplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V

duration=duration/1000;                     % switch to seconds
Rnoise=randn(1,round(duration*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
Lnoise=randn(1,round(duration*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
[edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
Rnoise(1:ledge)=Rnoise(1:ledge).*fliplr(edge);
Rnoise((end-ledge+1):end)=Rnoise((end-ledge+1):end).*edge;

Lnoise(1:ledge)=Lnoise(1:ledge).*fliplr(edge);
Lnoise((end-ledge+1):end)=Lnoise((end-ledge+1):end).*edge;

Rnoise=Rnoise./(max(abs(Rnoise)));             % normalize, so we could fit to +/-10V
Lnoise=Lnoise./(max(abs(Lnoise)));             % normalize, so we could fit to +/-10V

%amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); 
%amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %mw 080107
Rnoise=Ramplitude.*Rnoise;
Lnoise=Lamplitude.*Lnoise;
noise(1,:)=Rnoise;
noise(2,:)=Lnoise;
