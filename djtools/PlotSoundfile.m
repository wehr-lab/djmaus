function PlotSoundfile(varargin)

%plots clustered spiking tuning curve data from djmaus
%
% usage: PlotSoundfile([datadir], [tetrode], [clust], [xlimits],[ylimits], [binwidth])
% Everything is optional.
% if omitted, the following defaults are used:
%     datadir defaults to the current directory
%     tetrode defaults to all tetrodes in data directory
%     cluster defaults to all clusters
%     xlimits defaults to [-100 200]
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
    %first check for outfiles and use those if possible
    d=dir('outPSTH*.mat');
    if ~isempty(d)
        for i=1:length(d)
            PlotSoundfile_single(datadir, d(i).name, xlimits, ylimits, binwidth)
        end
    elseif exist('dirs.mat','file') %Nick addition 8/31/18 - defaults to kilosort output if present. Otherwise defaults to .t files.

        load('dirs.mat')
        masterdir=dirs{1};
        try
        cd(masterdir);
        catch
            remain = masterdir;
            while (remain ~= "")
                [token,remain] = strtok(remain, '\');
            end
                cd ..
                cd(token)
        end
        sp = loadKSdir(pwd);
        %The rest of this addition is to find the channel your cell is on.
        %It looks through the templates associated with the cell, finds the
        %most prevalent template within the cell, and returns the channel with the highest
        %amplitude for that template.
        %These two lines adapted from: https://groups.google.com/forum/?utm_medium=email&utm_source=footer#!msg/phy-users/Ydu1etOXwF0/-vEM9Rx_BgAJ
        tempChannelAmps = squeeze(max(sp.temps,[],2))-squeeze(min(sp.temps,[],2)); % amplitude of each template on each channel, size nTemplates x nChannels
        [~,maxChannel] = max(tempChannelAmps,[],2); % index of the largest amplitude channel for each template, size nTemplates x 1
        
        cd(datadir);
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
            %fn is normally the name of the .t file
            %In the case of Kilosort data however, readKiloSortOutput is used in the Process_single
            %function instead, and instead of using any .t files, it requires the
            %'cellnum' which it uses to iterate through the non-noise clusters in 'cluster_groups.csv'
            %So for now at least, I have fn = [clust, channel, cellnum]
            fn= [sp.cids(i),chan, i]; 
            %So for kilosort data, fn input to PlotClicktrain_PSTH_single will no longer be a string.
            %I made some small changes to accomodate for this:
            %   lines 44-55 in PlotClicktrain_PSTH_single
            %   lines 30-43 in ProcessClicktrain_PSTH_single
            PlotSoundfile_single(datadir, fn, xlimits, ylimits, binwidth)
        end
        %end of Nick addition 8/31/18. Some slight rearrangements below.
    else
        d=dir('*.t');
        if isempty(d)
            fprintf('\nNo clustered data found (no .t files in this directory)')
            PlotSoundfile_single(datadir, fn, xlimits, ylimits, binwidth)
        end
        for i=1:length(d)
            fn=d(i).name;
            PlotSoundfile_single(datadir, fn, xlimits, ylimits, binwidth)
        end
    end
else %user specified a channel
    if isempty(clust) % default to all clusters
        d=dir(sprintf('ch%d*.t', channel));
        for i=1:length(d)
            fn=d(i).name;
            PlotSoundfile_single(datadir, fn, xlimits, ylimits)
        end
    else %user specified a channel and a cluster
        if clust<10
            fn=sprintf('ch%d_simpleclust_0%d.t', channel, clust);
            PlotSoundfile_single(datadir, fn, xlimits, ylimits)
        else
            fn=sprintf('ch%d_simpleclust_%d.t', channel, clust);
            PlotSoundfile_single(datadir, fn, xlimits, ylimits)
            
        end
    end
end






