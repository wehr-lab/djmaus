function PlotGPIAS_PSTH(varargin)

%plots clustered spiking GPIAS data from djmaus
%
% usage: PlotGPIAS_PSTH([datapath], [tetrode], [clust], [xlimits],[ylimits], [binwidth], ['force-reprocess'])
% (all inputs are optional)
%     datadir defaults to the current directory
%     tetrode defaults to all tetrodes in data directory
%     cluster defaults to all clusters
%     xlimits defaults to [-1.5*max(gapdurs) 2*soa];
%     ylimits defaults to an autoscaled value
%     'force-reprocess' can be in any position
%
% tetrode number should be an integer
% clust can be an integer or an array of integers
%
%Processes data if outfile is not found;

force_reprocess=0;
for i=1:nargin
    try
        switch varargin{i}
            case {'force-reprocess', ...
                    'force_reprocess', ...
                    'reprocess', ...
                    'forcereprocess'}
                force_reprocess=1;
                varargin=varargin(setdiff(1:nargin, i));
        end
    end
end

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
if ischar(clust)
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
        if force_reprocess
            PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, 'force-reprocess')
        else
            PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth)
        end
    end
    if isempty(d)
        fprintf('\nNo clustered data found (no .t files in this directory)')
    end
    
else %user specified a channel
    if isempty(clust) % default to all clusters
        d=dir(sprintf('ch%d*.t', channel));
        for i=1:length(d)
            fn=d(i).name;
            if force_reprocess
                PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, 'force-reprocess')
            else
                PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth)
            end
        end
    else %user specified a channel and a cluster
        if clust<10
            fn=sprintf('ch%d_simpleclust_0%d.t', channel, clust);
            if force_reprocess
                PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, 'force-reprocess')
            else
                PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth)
            end
        else
            for j=1:length(clust)
                fn=sprintf('ch%d_simpleclust_%d.t', channel, clust(j));
                if force_reprocess
                    PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, 'force-reprocess')
                else
                    PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth)
                end
            end
            
        end
    end
end

