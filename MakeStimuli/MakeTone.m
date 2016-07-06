function tone=MakeTone(varargin)
% usage: tone=MakeTone(params, samplerate)
%inputs: 
% param (should have the following fields:)
%   frequency          -   frequency of the tone (Hz)
%   amplitude          -   sound pressure level of the tone (dB)            
%   duration           -   duration of the tone (ms)
%   ramp               -   length of an rising/falling edge (ascending/descending ramp) (ms)
% samplerate         -    sampling rate in hz

% Note: now we use absolute SPL instead of attenuation!!!

% Note: this was the original header
%function tone=MakeTone(frequency,attenuation,duration,samplerate,ramp)

global pref

% Creates a pure tone of a given frequency, attenuation, duration, at a
% given sample rate, with an ascending/descending ramp of a given length.
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



%     amplitude=10*(10.^((-attenuation)/20));
%    amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-10<x<10), i.e. pref.maxSPL=+_10V
    amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
    duration=duration/1000;                     % adjust the duration to seconds
    t=0:1/samplerate:duration;                  % length of the sampled trial
    tone=amplitude*sin(frequency*2*pi*t);       % the new tone itself
    if ramp>0
        [edge,ledge]=MakeEdge(ramp,samplerate);     % prepare the edges
        tone(1:ledge)=tone(1:ledge).*fliplr(edge);
        tone((end-ledge+1):end)=tone((end-ledge+1):end).*edge;
    end

% figure
% subplot(211)
% plot((1:length(tone))/samplerate*1000, tone)
% xlabel('Time (ms)')
% ylabel('amplitude')
