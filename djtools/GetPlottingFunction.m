function PlottingFunction=GetPlottingFunction(varargin)

%tries to figure out which plotting function is appropriate for some djmaus data
%Assume that we want the PSTH version (should be straightforward to plot
%LFPs instead if desired)
%
% usage: PlottingFunction=GetPlottingFunction([datadir])
% defaults to current directory

try
    datadir=varargin{1};
catch
    datadir=pwd; %default to pwd
end

PlottingFunction=[];
match=0;
cd(datadir)
stimparams=GetStimParams(datadir); %uses stimlog

%firsst check if stimulus protocol is new enough to contain the
%PlottingFunction field
if isfield(stimparams(1), 'PlottingFunction')
    PlottingFunction=stimparams(1).PlottingFunction;
    return
end

if any(strcmp(stimparams.stimtypes, 'clicktrain'))
    PlottingFunction='PlotClicktrain_PSTH';
    match=match+1;
end

if strfind(stimparams(1).protocol_description, 'Flashtrain')
    PlottingFunction='PlotFlashtrain_PSTH';
    match=match+1;
end

if any(strcmp(stimparams.stimtypes, 'AsymGPIAS'))
    PlottingFunction='PlotAsymGPIAS_PSTH';
    match=match+1;
end

if any(strcmp(stimparams.stimtypes, 'GPIAS'))
    PlottingFunction='PlotGPIAS_PSTH';
    match=match+1;
end

%tuning curve (tones/WN) by process of elimination
if ~match
    if all(strcmp(stimparams.stimtypes, 'whitenoise')) ...
            | ...
            all(strcmp(stimparams.stimtypes, 'tone')) ...
            | ...
            all(strcmp(stimparams.stimtypes, {'tone', 'whitenoise'})) ... %works because it's sorted
            | ...
            all(strcmp(stimparams.stimtypes, sort({'tone', 'whitenoise', 'silentsound'})))  %works because it's sorted
        PlottingFunction='PlotTC_PSTH';
    end
end