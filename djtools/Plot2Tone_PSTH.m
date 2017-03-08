function Plot2Tone_PSTH(varargin)

%plots clustered spiking 2 tone data from djmaus
%
% usage: Plot2Tone_PSTH([datadir], [tetrode], [clust], [xlimits],[ylimits], [binwidth])
% Everything is optional.
% if omitted, the following defaults are used:
%     datadir defaults to the current directory
%     tetrode defaults to all tetrodes in data directory
%     cluster defaults to all clusters
%     xlimits defaults to [-.5*dur 1.5*dur]
%     ylimits defaults to an autoscaled value
%
% channel (tetrode) number should be an integer
% clust can be an integer or an array of integers
%
%Processes data if outfile is not found;

rasters=1;

try
    datadir=varargin{1};
catch
    datadir=pwd; %default to pwd
end
try
    channel=varargin{2};
    if ischar(channel)
        channel=str2num(channel);
    end
catch
    channel=[];     %default to all tetrodes
end
try
    clust=varargin{3};
catch
    clust=[]; %to plot all clusts
end
if strcmp('char',class(clust))
    clust=str2num(clust);
end
try
    xlimits=varargin{4};
catch
    xlimits=[];
end
try
    ylimits=varargin{5};
catch
    ylimits=[];
end
try
    binwidth=varargin{6};
catch
    binwidth=5;
end

cd(datadir)
if isempty(channel)     %default to all tetrodes
    d=dir('*.t');
    for i=1:length(d)
        fn=d(i).name;
        Plot2Tone_PSTH_single(datadir, fn, xlimits, ylimits)
    end
else %user specified a channel
    if isempty(clust) % default to all clusters
        d=dir(sprintf('ch%d*.t', channel));
        for i=1:length(d)
            fn=d(i).name;
            Plot2Tone_PSTH_single(datadir, fn, xlimits, ylimits)
        end
    else %user specified a channel and a cluster
        if clust<10
            fn=sprintf('ch%d_simpleclust_0%d.t', channel, clust);
            Plot2Tone_PSTH_single(datadir, fn, xlimits, ylimits)
        else
            fn=sprintf('ch%d_simpleclust_%d.t', channel, clust);
            Plot2Tone_PSTH_single(datadir, fn, xlimits, ylimits)
            
        end
    end
end






