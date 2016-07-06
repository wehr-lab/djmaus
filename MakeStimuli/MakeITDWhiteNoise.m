function noise=MakeITDWhiteNoise(varargin)

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
amplitude=params.amplitude;
ramp=params.ramp;

amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V

duration=duration/1000;                     % switch to seconds
Rnoise=randn(1,round(duration*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;

%temp:
% Rnoise=0*Rnoise;
% Rnoise(50:60)=1;

[edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
Rnoise(1:ledge)=Rnoise(1:ledge).*fliplr(edge);
Rnoise((end-ledge+1):end)=Rnoise((end-ledge+1):end).*edge;
Rnoise=Rnoise./(max(abs(Rnoise)));             % normalize, so we could fit to +/-10V

%amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); 
%amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %mw 080107
Rnoise=amplitude.*Rnoise;

itd=params.itd; %in microsec
%positive ITD is defined as right ear leading
itd_samp=round(itd*samplerate/1000000); %itd in samples

Rnoise=[Rnoise zeros(1, abs(itd_samp))];

if itd_samp>0
    Lnoise=circshift(Rnoise, [1 itd_samp]);  
elseif itd_samp==0
    Lnoise=Rnoise;
elseif itd_samp<0
    Lnoise=Rnoise;
    Rnoise=circshift(Lnoise, [1 -itd_samp]);       
end

noise(1,:)=Rnoise;
noise(2,:)=Lnoise;
