function tone=Make2Tone(varargin)

% Note: now we use absolute SPL instead of attenuation!!!

% Note: this was the original header
%function tone=MakeTone(frequency,attenuation,duration,samplerate,ramp)
%last update: mw 10.30.2015 to allow use of whitenoise when freq==-1

global pref

% Creates a two tone sequence consisting of two pure tones of a 
% given frequency, attenuation, duration, at a
% given sample rate, with an ascending/descending ramp of a given length,
% at a given SOA (stimulus onset asynchrony)
% All lengths (duration, ramp) are in ms
% frequency is in Hz, attenuation in dB
% Input:
%  frequency          -   frequency of the tone (Hz)
% % % % % % % % % % % %  attenuation        -   attenuation of the tone (dB) 
% % % % % % % % % % % %                         !!!!!attenuation is relative to the max. sound pressure level 
% % % % % % % % % % % %                         as specified in Prefs.m (pref.maxSPL)
%  amplitude          -   sound pressure level of the tone (dB)            
%  duration           -   duration of the tone (ms)
%  samplerate         -   required sampling rate
%  ramp               -   length of an rising/falling edge (ascending/descending ramp) (ms)
%  probefreq          -   frequency of the probe tone (Hz)
%  probeamp           -   sound pressure level of the probe tone (dB) 
%  SOA                -   stimulus onset asynchrony (time in ms between onset of tone and probe tone)
%  
%  Output:
%  tone               -   the specified tone (empty if unsuccessful)
%  
 
tone=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
frequency=params.frequency;
% attenuation=params.attenuation;
if isfield(params,'amplitude')
    amplitude=params.amplitude;
else
    amplitude=params.attenuation;
end
duration=params.duration;
ramp=params.ramp;
probefreq=params.probefreq;
probeamp=params.probeamp;
SOA=params.SOA;



%     amplitude=10*(10.^((-attenuation)/20));
%    amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-10<x<10), i.e. pref.maxSPL=+_10V

%masker 
amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
duration_s=duration/1000;                     % adjust the duration to seconds
t=0:1/samplerate:duration_s;                  % length of the sampled trial
if frequency==-1 %white noise
    noise=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
    noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V
    masker=amplitude.*noise;
else %tone
    masker=amplitude*sin(frequency*2*pi*t);       % the new tone itself
end
if ramp>0
    [edge,ledge]=MakeEdge(ramp,samplerate);     % prepare the edges
    masker(1:ledge)=masker(1:ledge).*fliplr(edge);
    masker((end-ledge+1):end)=masker((end-ledge+1):end).*edge;
end

%SOA interval
SOA_silence=zeros(1, samplerate*(SOA-duration)/1000);

%probe 
probeamp=1*(10.^((probeamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
t=0:1/samplerate:duration_s;                  % length of the sampled trial

if frequency==-1 %white noise
        noise=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
        noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V
        probe=probeamp.*noise;       % corresponds to t=0:1/samplerate:duration;
else %tone
    probe=probeamp*sin(probefreq*2*pi*t);       % the new tone itself
end
if ramp>0
    [edge,ledge]=MakeEdge(ramp,samplerate);     % prepare the edges
    probe(1:ledge)=probe(1:ledge).*fliplr(edge);
    probe((end-ledge+1):end)=probe((end-ledge+1):end).*edge;
end

tone=[masker SOA_silence probe];

