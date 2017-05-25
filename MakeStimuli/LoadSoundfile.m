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

duration=params.duration;
ramp=params.ramp;

sourcefile=params.sourcefile;
sourcepath=params.sourcepath;
cd(pref.stimuli)
cd('Soundfile Protocols')
cd(sourcepath)
load(sourcefile)
samples=sample.sample;

% amplitude already applied during soundfile creation
