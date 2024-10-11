function sqWave=MakeSquare(varargin)
% usage: sqWave=MakesqWave(params, samplerate)
%inputs:
% params (should have the following fields:)
%   duty            -   percent of cycle positive
%   amplitude          -   sound pressure level of the tone (dB)
%   duration           -   duration of the square wave (ms)
% samplerate         -    sampling rate in hz

% Note: we use absolute SPL instead of attenuation!!!
% patterned after: tone=MakeTone(params, samplerate)
%%%% NOT SURE IF sqWave HAS TO START AND END WITH 0 ???

global pref

% Creates a squarewave of a given duty cycle, attenuation, duration, at a
% given samplerate
% duty is percent positive
% duration is in ms
% amplitude in dB
%
%  Output:
%  sqWave               -   the specified squarewave (empty if unsuccessful)
%

sqWave=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
duty=params.duty;
amplitude=params.amplitude;
duration=params.duration;


amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
duration=duration/1000;                     % adjust the duration to seconds
t = 0:1/samplerate:duration-1/samplerate;   % length of the sampled trial
sqWave=amplitude*square(2*pi/duration*t,duty);       % the new squareWave itself

% figure
% plot(0:1/30:params.duration-1/30,sqWave)
% xlabel('Time (ms)')
% ylabel('amplitude')

% pad to minimum length of 1 ms (to accomodate soundcardtrigger) mw.9.6.24
if duration<1
    silence=zeros(1, samplerate);
    silence(1:length(sqWave))=sqWave;
    sqWave=silence;
end
