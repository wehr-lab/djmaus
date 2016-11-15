function ProcessAsymGPIAS_PSTH(varargin)

%Processes clustered spiking AsymGPIAS data from djmaus
%
% usage: ProcessAsymGPIAS_PSTH([datadir], [tetrode], [cluster], [xlimits], [ylimits])
% Everything is optional.
% if omitted, the following defaults are used:
%     datadir defaults to the current directory
%     tetrode defaults to all tetrodes in data directory
%     cluster defaults to all clusters
%     xlimits defaults to [0 200]
%     ylimits defaults to an autoscaled value
%
% channel (tetrode) number should be an integer
% saves to outfile


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
    %    d=dir('*.t');
    %    for i=1:length(d)
    %        fn=d(i).name;
    %        ch=strsplit(fn , '_');
    %        c=strsplit(ch{1} , 'ch');
    %        chs(i)=str2num(c{2});
    %    end
    %    channel=unique(chs);
    %
end

try
    cluster=varargin{3};
catch
    cluster=[]; %if isempty(cluster), process all clusters
end

try
    xlimits=varargin{4};
catch
    xlimits=[];
end
if isempty(xlimits)
    xlimits=[];
end
try
    ylimits=varargin{5};
catch
    ylimits=[];
end


if isempty(channel)     %default to all tetrodes
    d=dir('*.t');
    for i=1:length(d)
        fn=d(i).name;
        ProcessAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
    end
else %user specified a channel
    if isempty(cluster) % default to all clusters
        d=dir(sprintf('ch%d*.t', channel));
        for i=1:length(d)
            fn=d(i).name;
            ProcessAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
        end
    else %user specified a channel and a cluster
        if cluster<10
            fn=sprintf('ch%d_simpleclust_%0d.t', channel, cluster);
            ProcessAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
            
        else
            fn=sprintf('ch%d_simpleclust_%d.t', channel, cluster);
            ProcessAsymGPIAS_PSTH_single(datadir, fn, xlimits, ylimits)
            
        end
    end
end


