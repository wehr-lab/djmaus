function [edge, ledge]=MakeEdge(ramp, samplerate);

% generates rising/falling edge for the stimuli.
% Input:
% ramp          -   length of the edge (ms)
% samplerate    -   sampling rate to use for the new edge(Hz)
% Output:
% edge          -   the brand new edge itself
% ledge         -   length of the new edge
%
% Returns empty edge and ledge=0 if unsuccessful
%
% Usage:
% The following adds edges to 'sample'
%   [edge,ledge]=MakeEdge(5,samplerate)
%   sample(1:ledge)=sample(1:ledge).*fliplr(edge);
%   sample((end-ledge+1):end)=sample((end-ledge+1):end).*edge;
%

edge=[];
ledge=0;
if nargin<2
    return;
end

omega=(1e3/ramp)*(acos(sqrt(0.1))-acos(sqrt(0.9)));
t=0:1/samplerate:pi/2/omega;
t=t(1:(end-1));
edge=(cos(omega*t)).^2;
ledge=length(edge);
