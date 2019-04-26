function KS_PlotGPIAS_PSTH(varargin)

%plots clustered spiking GPIAS data from djmaus
%
% usage: PlotGPIAS_PSTH([datapath], [tetrode], [clust], [xlimits],[ylimits], [binwidth])
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

cd(datadir)
if exist('dirs.mat','file')
    load('dirs.mat')
    try
        masterdir=dirs{1};
    catch
        masterdir=pwd;
    end
    cd(masterdir);
    sp = loadKSdir(pwd);
    
    %These two lines adapted from: https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/phy-users/Ydu1etOXwF0/-vEM9Rx_BgAJ
    tempChannelAmps = squeeze(max(sp.temps,[],2))-squeeze(min(sp.temps,[],2)); % amplitude of each template on each channel, size nTemplates x nChannels
    [~,maxChannel] = max(tempChannelAmps,[],2); % index of the largest amplitude channel for each template, size nTemplates x 1
    cd(datadir);
    
    [ChanClustCellnum] = getChanClustCellnum(sp,maxChannel);
    
    if isempty(channel)     %default to all tetrodes
        for i = 1:size(ChanClustCellnum,1)
            KS_PlotGPIAS_PSTH_single(datadir, ChanClustCellnum(i,:), xlimits, ylimits, binwidth)
        end
    else    %user specified a channel
        if isempty(clust)   %default to all clusters
            d = find(ChanClustCellnum(:,1) == channel);
            if isempty(d)
                fprintf('\nNo clustered data found for this channel')
            else
                for i=1:length(d)
                    KS_PlotGPIAS_PSTH_single(datadir, ChanClustCellnum(d(i),:), xlimits, ylimits, binwidth)
                end
            end
        else %user specified a channel and a cluster
            d = find(ChanClustCellnum(:,2) == clust,1);
            KS_PlotGPIAS_PSTH_single(datadir, ChanClustCellnum(d,:), xlimits, ylimits, binwidth)
        end
        
    end

else
    fprintf('No dirs, so Im doing nothing')
end

end

function [ChanClustCellnum] = getChanClustCellnum(sp,maxChannel)
    ChanClustCellnum = []; %
    %sp.cids are the kilosort cell IDs which you didn't label as noise. Found in 'cluster_groups.csv'
    for i=1:length(sp.cids) %each iteration here therefore corresponds to the 'cellnum' used in readKiloSortOutput.
        template = []; Utemplate = []; numUtemplate = [];
        template = sp.spikeTemplates(sp.clu==sp.cids(i)); %identifies the template for every spike in this cluster/cell_ID
        Utemplate = unique(template); %identify all the unique templates contained in this cluster/cell_ID

        if length(Utemplate) == 1 %If the cluster contains spikes from only one template
            template = ((template(1,1))+1); %then use this template (+1 because they are indexed from 0:nTemplates-1 and maxChannel is 1:nTemplates)
            chan = maxChannel(template,1); %retun the channel that was identified to have the highest amplitude for this template
        else
            for k = 1:length(Utemplate) %If the cluster contains spikes from more than one template
                numUtemplate(k,1) = sum(template(:) == Utemplate(k,1)); %find how many spikes in this cluster/cell_ID there are, for each template 
                [val, idx] = max(numUtemplate); %find which template is most prevalent
                template = (Utemplate(idx,1)+1); %use this template (+1 because they are indexed from 0 and maxChannel is indexed 1:nTemplates)
                chan = maxChannel(template,1); %retun the channel that was identified to have the highest amplitude for this template
            end
        end
        chan = chan - 1;
        
        ChanClustCellnum(i,:) = [chan,sp.cids(i),i];
        %fn is normally the name of the .t file
        %In the case of Kilosort data however, readKiloSortOutput is used in the Process_single
        %function instead, and instead of using any .t files, it requires the
        %'cellnum' which it uses to iterate through the non-noise clusters in 'cluster_groups.csv'
        %So for now at least, I have fn = [clust, channel, cellnum]
                    %fn = [sp.cids(i),chan, i]; 
        %So for kilosort data, fn input to PlotClicktrain_PSTH_single will no longer be a string.
        %I made some small changes to accomodate for this:
        %   lines 44-55 in PlotClicktrain_PSTH_single
        %   lines 30-43 in ProcessClicktrain_PSTH_single
                    %KS_PlotGPIAS_PSTH_single(datadir, fn, xlimits, ylimits, binwidth)
    end
    %end of Nick addition 8/31/18. Some slight rearrangements below.
end