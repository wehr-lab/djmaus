function samples=LoadSoundfile(varargin)
% usage: tone=LoadSoundfile(params, samplerate)
%inputs: 
% param (should have the following fields:)
%   frequency          -   frequency of the tone (Hz)
%   amplitude          -   sound pressure level of the tone (dB)            
%   duration           -   duration of the tone (ms)
%   ramp               -   length of an rising/falling edge (ascending/descending ramp) (ms)
% samplerate         -    sampling rate in hz

% Note: now we use absolute SPL instead of attenuation!!!


global pref

 
 
samples=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
    amplitude=params.amplitude;
    amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V

duration=params.duration;
ramp=params.ramp;

sourcefile=params.sourcefile;
sourcepath=params.sourcepath;
cd(pref.stimuli)
cd('Soundfile Protocols')
cd(sourcepath)
load(sourcefile)
samples=sample.sample;
samples=amplitude*samples;

%NO: mike 4.23.2024 % amplitude already applied during soundfile creation
