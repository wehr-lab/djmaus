function amnoise = MakeAMNoise(varargin) 

global pref

% NOTE: this was the original header
% function y = MakeAMNoise(modulation_frequency, modulation_depth, duration, amplitude, samplerate, ramp) 

% [This function generates amplitude-modulated noises]
%
% function y = MakeAMNoise(mod_freq, am_depth, duration, SR) 
% modulation_frequency      -   modulation frequency (Hz)
% modulation_depth          -   modulation depth (0-1)
% duration                  -   duration (ms)
% amplitude                 -   amplitude (dB)
% samplerate                -   sampling rate (Hz)
% ramp                      -   rising/falling edge (in ms)

% orig: Wang Lab, Johns Hopkins University (Edited on Janyary 13, 2004 by Tom Lu)

amnoise=[];

if nargin<2
    return;
end

params              =varargin{1};
samplerate          =varargin{2};
modulation_frequency=params.modulation_frequency;
modulation_depth    =params.modulation_depth;
duration            =params.duration;
amplitude           =params.amplitude;
ramp                =params.ramp;

    yc = randn(1,duration/1000*samplerate); % generate noise signal  
    ymod =(1+modulation_depth*cos(2*pi*modulation_frequency*(1:length(yc))/samplerate + pi));
    y = ymod .* yc;  %apply modulation to signal
    zmax = 2.576;  % 99-prctile (two tailed)
    y(find(y > zmax*max(ymod))) = zmax*max(ymod);
    y(find(y < -zmax*max(ymod))) = -zmax*max(ymod);
    y(end-1:end)=[];
    y = y/(zmax*max(abs(ymod)));

% This was the original Wang's ramp    
%         stim=y;
%         RFT = ramp; %rise fall time
%         onset = RFT * samplerate/1000;
%         dur = length(stim);
%         out = stim;
% 
%         for i=1:onset
%           out(i) = stim(i) *i/onset;
%           out(dur-i+1) = stim(dur-i+1) *i/onset;
%         end
%         y=out;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% our standard ramp
[edge,ledge]=MakeEdge(ramp,samplerate);     % prepare the edges
y(1:ledge)=y(1:ledge).*fliplr(edge);
y((end-ledge+1):end)=y((end-ledge+1):end).*edge;

%amplitude=10*(10.^((amplitude-pref.maxSPL)/20));
amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %mw 01-19-2012

y=amplitude.*y;    
amnoise=y;    

% figure
% subplot(211)
% plot((1:length(y))/samplerate*1000, y)
% xlabel('Time (ms)')
% ylabel('amplitude')
% 
% subplot(212)
% temp=fft(y);
% plot((1:length(temp))/length(temp)*(samplerate/1000), abs(temp))
% specgram(y,[],samplerate/1000)
% ylabel('Frequency (kHz)')
