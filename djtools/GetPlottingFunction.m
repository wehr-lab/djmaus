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

if strfind(stimparams(1).protocol_description, 'GPIAS protocol with flashtrain')
    PlottingFunction='PlotGPIASflashtrain_PSTH';
    match=match+1;
end

%tuning curve (tones/WN) by process of elimination
if ~match
    switch length(stimparams.stimtypes)
        case 1
            if all(strcmp(stimparams.stimtypes, 'whitenoise')) | ...
                    all(strcmp(stimparams.stimtypes, 'tone')) | ...
                    all(strcmp(stimparams.stimtypes, 'silentsound')) %I think PlotTC_PSTH would be appropriate for silentsound only
                PlottingFunction='PlotTC_PSTH';
            end
        case 2
            if all(strcmp(stimparams.stimtypes, {'tone', 'whitenoise'})) | ...
                    all(strcmp(stimparams.stimtypes, {'silentsound', 'whitenoise'})) | ...
                    all(strcmp(stimparams.stimtypes, {'silentsound', 'tone'})) 
                    PlottingFunction='PlotTC_PSTH';
            end
        case 3
            if  all(strcmp(stimparams.stimtypes, sort({'silentsound', 'tone', 'whitenoise'}))) ...  %works because it's sorted
                    PlottingFunction='PlotTC_PSTH';
            end
    end
end