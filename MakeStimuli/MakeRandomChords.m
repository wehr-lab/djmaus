function [chords, params]=MakeRandomChords(varargin)

% Note: now we use absolute SPL instead of attenuation!!!

global pref

% Creates a random chords stimulus (ala deCharms)
% Input:
%   frequency       -   [min max] frequencies
%   toneduration    -   duration of individual tones
%   amplitude       -   amplitude of individual tone components
%   octavesteps     -   how many frequencies per octave
%   density         -   number of frequency components in each chord
%   ramp            -   ramp of each individual tone
%   duration        -   total required duration
% Output:
%   chords          -   requested stimulus (empty if unsuccessful)
%   params          -   input params plus some additional:
%                           duration (may be different than requested)
%                           spectrum (matrix containing frequency
%                                     components of individual chords)

chords=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};

frequency=      params.frequency;
toneduration=   params.toneduration;
octavesteps=    params.octavesteps;
density=        params.density;
amplitude=      params.amplitude;
duration=       params.duration;
ramp=           params.ramp;

maxfreq=max(frequency);
minfreq=min(frequency);
octaves=log2(maxfreq/minfreq);
steps=0:(1/octavesteps):octaves;
freqs=minfreq*(2.^steps);
nfreqs=length(freqs);

if density>nfreqs   % we wouldn't have enough frequency components
    return;
end

toneparam.duration=toneduration;
% toneparam.ramp=ramp;
toneparam.ramp=0;       % ramp will be added later to the final sound
toneparam.amplitude=amplitude;
toneparam.frequency=freqs(1);

tone=MakeTone(toneparam,samplerate);
sampledtoneduration=length(tone);
tones=zeros(nfreqs,sampledtoneduration);
tones(1,:)=tone;
for k=2:nfreqs
    toneparam.frequency=freqs(k);
    tones(k,:)=MakeTone(toneparam,samplerate);
end

amplitude=10*(10.^((amplitude-pref.maxSPL)/20));
duration=(duration-mod(duration,toneduration));    % duration in seconds

nbins=duration/toneduration;

chords=zeros(1,sampledtoneduration*nbins);

chordparams=zeros(nbins,density);

ramp=params.ramp;   % requested ramp
[edge,ledge]=MakeEdge(ramp,samplerate);
for c=1:nbins
    chordfreq=randperm(nfreqs);
    chordfreq=chordfreq(1:density);
    onechord=sum(tones(chordfreq,:));
    % now add the ramp
     if ramp>0
        onechord(1:ledge)=onechord(1:ledge).*fliplr(edge);
        onechord((end-ledge+1):end)=onechord((end-ledge+1):end).*edge;
    end    
end

    chords(((c-1)*sampledtoneduration+1):c*sampledtoneduration)=onechord;
    chordparams(c,:)=freqs(chordfreq);
end

params.spectrum=chordparams;