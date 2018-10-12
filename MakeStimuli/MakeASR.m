function tone=MakeASR(varargin)

%makes an ASR (acoustic startle response) stimulus with a narrow band
%prepulse and white noise startle pulse
%modified from makeASR
% Note: now we use absolute SPL instead of attenuation!!!

% Note: this was the original header
%function tone=MakeTone(frequency,attenuation,duration,samplerate,ramp)

global pref

% Creates a two tone sequence consisting of two noise bursts,
% attenuation, duration, at a
% given sample rate, with an ascending/descending ramp of a given length,
% at a given SOA (stimulus onset asynchrony)
% All lengths (duration, ramp) are in ms
% frequency is in Hz, attenuation in dB
% Input:
%  frequency          -   frequency of the tone (Hz)
% % % % % % % % % % % %  attenuation        -   attenuation of the tone (dB)
% % % % % % % % % % % %                         !!!!!attenuation is relative to the max. sound pressure level
% % % % % % % % % % % %                         as specified in Prefs.m (pref.maxSPL)
% prepulsefreq  - center frequency for narrow band pre pulse
% prepulsebandwidth - bandwidth for narrow band pre pulse (octaves)
%  amplitude          -   sound pressure level of the tone (dB)
%  duration           -   duration of the tone (ms)
%  samplerate         -   required sampling rate
%  ramp               -   length of an rising/falling edge (ascending/descending ramp) (ms)
%  probefreq          -   frequency of the probe tone (Hz)
%  probeamp           -   sound pressure level of the probe tone (dB)
%  SOA                -   stimulus onset asynchrony (time in ms between onset of tone and probe tone)
%  soaflag            -   % soaflag: can be either 'soa' (default), in which case soa value specifies the time
%                           between the onset of the gap and the onset of the startle, or else 'isi',
%                           in which case soa specifies the time between gap offset and startle
%                           onset. If anything other than 'isi' it will default to 'soa'.
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
pulsedur=params.pulsedur;
pulseamp=params.pulseamp;
SOA=params.soa;
if isfield(params, 'soaflag')
    soaflag=params.soaflag;
else
    soaflag='soa';
end
prepulsedur=params.prepulsedur;
prepulseamp=params.prepulseamp;

if prepulsedur>0 %you can use 0 ms dur to indivate no prepulse - mw 09.24.15
    
    %narrowband noise prepulse now handled by MakeNBASR
    % % narrowband noise prepulse
    % noiseparam.amplitude=params.prepulseamp;
    % noiseparam.filter_operation='bandpass';
    % noiseparam.center_frequency=params.prepulsefreq;
    % noiseparam.lower_frequency=params.prepulsefreq/2^(params.prepulsebandwidth/2);
    % noiseparam.upper_frequency=params.prepulsefreq*2^(params.prepulsebandwidth/2);
    % noiseparam.ramp=ramp;
    % noiseparam.duration=prepulsedur;
    % prepulse=MakeNoise(noiseparam, samplerate);
    %
    % amplitude=1*(10.^((prepulseamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
    % noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V
    % prepulse=amplitude.*noise;
    
    
    %white noise prepulse
    
    amplitude=1*(10.^((prepulseamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
    duration_s=prepulsedur/1000;                     % adjust the duration to seconds
    t=0:1/samplerate:duration_s;                  % length of the sampled trial
    noise=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
    [edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
    noise(1:ledge)=noise(1:ledge).*fliplr(edge);
    noise((end-ledge+1):end)=noise((end-ledge+1):end).*edge;
    noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V
    prepulse=amplitude.*noise;
    
else
    prepulse=[];
end
%%%%%%

%SOA interval

switch soaflag
    case 'isi'
        SOA_silence=zeros(1, samplerate*(SOA)/1000);
    case 'soa'
        SOA_silence=zeros(1, samplerate*(SOA-prepulsedur)/1000);
end

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

