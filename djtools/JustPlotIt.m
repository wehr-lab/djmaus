function JustPlotIt(varargin )
% tries to figure out the appropriate way to plot the data in a data
% directory, and does it. 
% For djmaus/open-ephys data.
%
% usage: JustPlotIt([datapath], [other parameters])
%  defaults to current directory
%  if you include a list of other parameters (like xlimits, tetrode number,
%  etc) they will get passed on to the appropriate function. They have to
%  be in the right order for that function.

if nargin==0
    datapath=pwd;
end


try
    PlottingFunction=GetPlottingFunction
    feval(PlottingFunction, varargin{:})
catch
    fprintf('\nsorry, could not figure out the appropriate plotting function for the data in %s', datapath)
end
