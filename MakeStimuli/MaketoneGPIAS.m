function tone=MakeGPIAS(varargin)

% Note: now we use absolute SPL instead of attenuation!!!


global pref

% Creates a GPIAS consisting of a continuous background noise, a gap, and a startle pulse. 
%
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
%noiseamp=params.noiseamp;
%noisefreq=params.noisefreq;
%noisebandwidth=params.noisebandwidth;
gapdur=params.gapdur;
gapdelay=params.gapdelay;
pulsedur=params.pulsedur;
pulseamp=params.pulseamp;
ramp=params.ramp;
SOA=params.soa;
if isfield(params, 'soaflag')
    soaflag=params.soaflag;
else
    soaflag='soa';
end
gapdur_samples=gapdur*samplerate/1000;
gapdelay_samples=gapdelay*samplerate/1000;
SOA_samples=SOA*samplerate/1000;

noise_params=params;
%noise_params.duration=params.duration+5000; %make it a second longer to cover up isi
noise_params.duration=params.duration; %don't make it a second longer to cover up isi
continuous_noise=MakeWhiteNoise(noise_params, samplerate);
%continuous noise
% amplitude=1*(10.^((noiseamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
% duration_s=continuous_noise_dur/1000;                     % adjust the duration to seconds
% t=0:1/samplerate:duration_s;                  % length of the sampled trial
% noise=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
% [edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
% noise(1:ledge)=noise(1:ledge).*fliplr(edge);
% noise((end-ledge+1):end)=noise((end-ledge+1):end).*edge;
% noise=noise./(max(abs(noise)));             % normalize, so we could fit to +/-10V

%insert Gap 
% the following line is for gapdelay specifying time to gap onset
%continuous_noise(gapdelay_samples:gapdelay_samples+gapdur_samples)=zeros(size(gapdelay_samples:gapdelay_samples+gapdur_samples));
% the following line is for gapdelay specifying time to gap offset
continuous_noise(gapdelay_samples-gapdur_samples:gapdelay_samples)=zeros(size(gapdelay_samples-gapdur_samples:gapdelay_samples));

%pulse
if pulsedur>0 %allow for no pulse if duration ==0
    amplitude=1*(10.^((pulseamp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
    duration_s=pulsedur/1000;                     % adjust the duration to secondst=0:1/samplerate:duration_s;                  % length of the sampled trial
    t=0:1/samplerate:duration_s;                  % length of the sampled trial
    noisepulse=randn(1,round(duration_s*samplerate)+1);       % corresponds to t=0:1/samplerate:duration;
    [edge,ledge]=MakeEdge(ramp,samplerate);     % and add the edge
    noisepulse(1:ledge)=noisepulse(1:ledge).*fliplr(edge);
    noisepulse((end-ledge+1):end)=noisepulse((end-ledge+1):end).*edge;
    noisepulse=noisepulse./(max(abs(noisepulse)));             % normalize, so we could fit to +/-10V
    pulse=amplitude.*noisepulse;
    
    % insert pulse after SOA interval
    switch soaflag
        case 'isi'
            pulse_start=gapdelay_samples+SOA_samples;
            pulse_stop=pulse_start+length(pulse);
            continuous_noise(pulse_start+1:pulse_stop)=pulse;
        case 'soa'
            pulse_start=gapdelay_samples-gapdur_samples+SOA_samples;
            pulse_stop=pulse_start+length(pulse);
            continuous_noise(pulse_start+1:pulse_stop)=pulse;
    end
end

tone=continuous_noise;
% figure(2)
% plot(tone)


