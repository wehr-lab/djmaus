function JustPlotIt( datapath, varargin )
% tries to figure out the appropriate way to plot the data in a data
% directory, and does it. 
% For djmaus/open-ephys data.
%
% usage: JustPlotIt(datapath)
if nargin==0
    datapath=pwd;
end

stimparams=GetStimParams(datapath);
stimparams.stimtypes={'tone', 'whitenoise'}

if strcmp(stimparams.stimtypes,  'AsymGPIAS')
        PlotAsymGPIAS_PSTH(datapath, varargin{:})
elseif any(strcmp(stimparams.stimtypes,  'GPIAS'))
        PlotAsymGPIAS_PSTH(datapath, varargin{:})
elseif any(strcmp(stimparams.stimtypes,  'tone'))
            PlotTC_PSTH(datapath, varargin{:})
else
        fprintf('\nsorry, could not figure out the appropriate plotting function for the data in %s', datapath)
end
