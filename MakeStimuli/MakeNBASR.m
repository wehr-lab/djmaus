function tone=MakeNBASR(varargin)

%makes an ASR (acoustic startle response) stimulus with a 
%prepulse and white noise startle pulse
%modified from makeASR 

% Note: now we use absolute SPL instead of attenuation!!!


global pref

% Creates a two tone sequence consisting of two noise bursts, the first is 
% narrow band, the second white noise, attenuation, duration, at a
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
ramp=params.ramp;
prepulsedur=params.prepulsedur;
prepulseamp=params.prepulseamp;
prepulsefreq=params.prepulsefreq;
prepulsebandwidth=params.prepulsebandwidth;
pulsedur=params.pulsedur;
pulseamp=params.pulseamp;
SOA=params.soa;

%params for narrow band pre pulse;
noise_params.type='noise';
noise_params.amplitude=prepulseamp;
noise_params.filter_operation='bandpass';
noise_params.center_frequency=prepulsefreq;
noise_params.lower_frequency=prepulsefreq-prepulsebandwidth/2;
noise_params.upper_frequency=prepulsefreq+prepulsebandwidth/2;
noise_params.ramp=ramp;
noise_params.duration=prepulsedur;
    
%narrow band prepulse
prepulse=MakeNoise(noise_params, samplerate);
% 
% amplitude=1*(10.^((prepulseamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
% duration_s=prepulsedur/1000;                     % adjust the duration to seconds
% t=0:1/samplerate:duration_s;                  % length of the sampled trial
% noise=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
% [edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
% noise(1:ledge)=noise(1:ledge).*fliplr(edge);
% noise((end-ledge+1):end)=noise((end-ledge+1):end).*edge;
% noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V
% prepulse=amplitude.*noise;
%%%%%%

%SOA interval
SOA_silence=zeros(1, samplerate*(SOA-prepulsedur)/1000);

%pulse
amplitude=1*(10.^((pulseamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
duration_s=pulsedur/1000;                     % adjust the duration to secondst=0:1/samplerate:duration_s;                  % length of the sampled trial

t=0:1/samplerate:duration_s;                  % length of the sampled trial
noise=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
[edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
noise(1:ledge)=noise(1:ledge).*fliplr(edge);
noise((end-ledge+1):end)=noise((end-ledge+1):end).*edge;
noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V
pulse=amplitude.*noise;


tone=[prepulse SOA_silence pulse];

