function silentsound=MakeSilentSound(varargin)


% Creates a "silent sound", which is a token of white-noise at -1000 dB SPL
% this is useful for explicitly collecting spontaneous activity or
% laser-only as its own stimulus condition 
% 
% Usage: silentsound=MakeSilentSound(param, samplerate)
% Inputs: 
% param should have the following fields:
%   duration           -   duration of the stimulus (ms)
%   amplitude          -   sound pressure level of the sound (dB SPL)
%   ramp               -   length of an ascending/descending edge MS???
% 
% samplerate         -   required sampling rate (Hz)
%   
%  Output:
%  noise              -   the requested sample; empty if unsuccessful
%
global pref

params=varargin{1};
samplerate=varargin{2};
params.amplitude=-1000;

silentsound=[];
if nargin<2
    return;
end
silentsound=MakeWhiteNoise(params, samplerate);


