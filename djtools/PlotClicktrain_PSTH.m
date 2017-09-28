function PlotClicktrain_PSTH(varargin)

%plots clustered spiking Clicktrain data from djmaus
%
% usage: PlotClicktrain_PSTH([datapath], [tetrode], [clust], [xlimits],[ylimits], [binwidth])
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
try
    % Do we want to show every plot as we make it?
    show_plots = varargin{7};
catch
    % Default to yes
    show_plots = 1;
end


cd(datadir)
if isempty(channel)     %default to all tetrodes
    d=dir('*.t');
    for i=1:length(d)
        fn=d(i).name;
        PlotClicktrain_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, show_plots)
    end
else %user specified a channel
    if isempty(clust) % default to all clusters
        d=dir(sprintf('ch%d*.t', channel));
        for i=1:length(d)
            fn=d(i).name;
            PlotClicktrain_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, show_plots)
        end
    else %user specified a channel and a cluster
        d=0; %to avoid throwing "empty" error below
        if clust<10
            fn=sprintf('ch%d_simpleclust_0%d.t', channel, clust);
            PlotClicktrain_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, show_plots)
        else
            fn=sprintf('ch%d_simpleclust_%d.t', channel, clust);
            PlotClicktrain_PSTH_single(datadir, fn, xlimits, ylimits, binwidth, show_plots)
        end
    end
end

if isempty(d)
    fprintf('\nNo clustered data found (no .t files in this directory)')
end