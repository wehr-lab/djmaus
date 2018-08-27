function [t, fr]=GaussSmooth(spiketimes, sigma, xlimits)
%usage: [t, fr]=GaussSmooth(spiketimes, sigma, xlimits)
%smooths a spiketrain by gaussian convolution
%inputs:
%  spiketimes: a list of spiketimes, in ms, relative to stimulus onset
%       (e.g. [75.2406 -15.4977 176.7524 49.1524 ... ])
%  sigma: standard deviation of gaussian kernel, in ms
%  xlimits: xlimits for time vector, in ms
%outputs:
%  t: time vector for smoothed firing rate, in ms
%  fr: smoothed firing rate, in Hz
%for example, you could plot(t, fr); xlabel('time, ms'); ylabel('firing rate, Hz')

%gaussian kernel 
winsize=1000; %in 10x ms
gk=gausswin(winsize, .1*winsize/sigma);
gk=gk./sum(gk); %->unit area

% create 10x upsampled spiketrain with hist
[spiketrain, t]=hist(spiketimes, xlimits(1):.1:xlimits(2)); %

%convolve
c=conv(spiketrain, gk, 'same'); %c is in spikes/10xms
fr=10000*c; %convert back to Hz

