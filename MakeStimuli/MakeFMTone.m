function tone=MakeFMTone(varargin)

global pref

% NOTE: this was the original header
% function tone=MakeFMTone(Fc, carrier_phase, Fm, modulation_phase, modulation_index, amplitude, duration, samplerate, ramp)

% [This function generates frequency-modulated tones]
%
% function y=MakeFMTone(Fc, carrier_phase, Fm, mod_phase, mf, A, stimlength, SR)
% Fc            -   carrier frequency (in Hz)
% carrier_phase -   carrier phase
% Fm            -   modulation frequency (in Hz)
% modulation_phase  -   modulation phase
% modulation_index  -   modulation index (frequency deviation in Hz)
% amplitude     -   amplitude (dB)
% duration      -   duration (in msec)
% samplerate    -   Sampling rate (in Hz)
% ramp          - rising/falling edge (in ms)

% orig: sFM.m Wang Lab, Johns Hopkins University (Edited on January 13, 2004 by Tom Lu)

tone=[];

if nargin<2
    return;
end

params          =varargin{1};
samplerate      =varargin{2};
Fc              =params.carrier_frequency;
carrier_phase   =params.carrier_phase;
Fm              =params.modulation_frequency;
modulation_phase=params.modulation_phase;
modulation_index=params.modulation_index;
amplitude       =params.amplitude;
duration        =params.duration;
ramp            =params.ramp;

%    amplitude=10*(10.^((amplitude-pref.maxSPL)/20));
    amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %mw 01-19-2012 !

    npts=samplerate*(duration/1000);
    x=(1:npts)/samplerate;
    tone=sin(2*pi*Fc*x - modulation_index/Fm*(cos(2*pi*Fm*x + modulation_phase - pi/2)) + carrier_phase);
    tone=amplitude*tone;

    [edge,ledge]=MakeEdge(ramp,samplerate);     % prepare the edges
    tone(1:ledge)=tone(1:ledge).*fliplr(edge);
    tone((end-ledge+1):end)=tone((end-ledge+1):end).*edge;


% figure
% subplot(211)
% plot((1:length(tone))/samplerate*1000, tone)
% xlabel('Time (ms)')
% ylabel('amplitude')
% 
% subplot(212)
% plot((1:length(fft(tone)))/length(fft(tone))*(samplerate/1000), abs(fft(tone)))
% specgram(tone,[],samplerate/1000)
% ylabel('Frequency (kHz)');
% 
