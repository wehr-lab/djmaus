function PlotAsymGPIAS_PSTH(varargin)

%plots clustered spiking AsymGPIAS data from djmaus
%
% usage: PlotAsymGPIAS_PSTH([datapath], [tetrode], [clust], [xlimits],[ylimits], [binwidth])
% (all inputs are optional)
%     datadir defaults to the current directory
%     tetrode defaults to all tetrodes in data directory
%     cluster defaults to all clusters
%     xlimits defaults to [-1.5*max(gapdurs) 2*soa];
%     ylimits defaults to an autoscaled value
%
% tetrode number should be an integer
% clust can be an integer or an array of integers
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;


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
        PlotAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
    end
else %user specified a channel
    if isempty(cluster) % default to all clusters
        d=dir(sprintf('ch%d*.t', channel));
        for i=1:length(d)
            fn=d(i).name;
            PlotAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
        end
    else %user specified a channel and a cluster
        if cluster<10
            fn=sprintf('ch%d_simpleclust_%0d.t', channel, cluster);
            PlotAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
                    else
            fn=sprintf('ch%d_simpleclust_%d.t', channel, cluster);
            PlotAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
            
        end
    end
end






