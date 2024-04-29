function PlotTC_PSTH(varargin)

%plots clustered spiking tuning curve data from djmaus
%
% usage: PlotTC_PSTH([datadir], [tetrode], [clust], [xlimits],[ylimits], [binwidth])
% usage: PlotTC_PSTH([datadir], [cellnum], [], [xlimits],[ylimits], [binwidth])
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

%first try for new OE format -mike 04.28.2024
if exist([datadir, filesep, 'dirs.mat'],'file') | exist([datadir, filesep, 'Bdirs.mat'],'file')
    try
        load([datadir, filesep, 'dirs.mat'])
    catch
        load([datadir, filesep, 'Bdirs.mat'])
    end
    masterdir=Bdirs{1};
    cd(masterdir);
    d=dir('OEinfo*.mat');
    if ~isempty(d) %if OEinfo exists, we're using new OE format -mike 04.28.2024
        load(d.name)
        cd (dirs{1})
        outfilenames=dir('outPSTH_*.mat');
        if isempty(outfilenames)
            sorted_units_filename=dir('SortedUnits*.mat');
            load(sorted_units_filename.name)
            ProcessTC_PSTH_single2([], SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, [],[])
            outfilenames=dir('outPSTH_*.mat');
        end
        for i=1:length(outfilenames) %each iteration here therefore corresponds to the 'cellnum' used in readKiloSortOutput.
            PlotTC_PSTH_single2(pwd, outfilenames(i).name, xlimits,ylimits, binwidth)
        end
    end
else %try old OE and kilosort format

    if isempty(channel)     %default to all tetrodes
        %Nick addition 8/31/18 - defaults to kilosort output if present. Otherwise defaults to .t files.
        if exist([datadir, filesep, 'dirs.mat'],'file')
            load([datadir, filesep, 'dirs.mat'])
            masterdir=dirs{1};
            cd(masterdir);
            try
                cd('Spikes') %mw 04.26.2024 should rewrite to load SortedUnits instead
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
                PlotTC_PSTH_single(datadir, fn, xlimits, ylimits)
            end
            %end of Nick addition 8/31/18. Some slight rearrangements below.
        else
            d=dir('*.t');
            if isempty(d)
                fprintf('\nNo clustered data found (no .t files in this directory)')
                PlotTC_PSTH_single(datadir, fn, xlimits, ylimits)
            end
            for i=1:length(d)
                fn=d(i).name;
                PlotTC_PSTH_single(datadir, fn, xlimits, ylimits)
            end
        end
    else %user specified a channel
        if isempty(clust) % default to all clusters
            d=dir(sprintf('ch%d*.t', channel));
            if isempty(d) || channel == -1
                load('dirs.mat')
                sp = loadKSdir(dirs{1});
                cell_indx=find(sp.cgs==2);
                for c = 1:length(cell_indx)
                    t_filename = sprintf('ch-1_simpleclust_%d.t', cell_indx(c));
                    PlotTC_PSTH_single(datadir, t_filename, xlimits, ylimits)
                end
            else
                for i=1:length(d)
                    fn=d(i).name;
                    PlotTC_PSTH_single(datadir, fn, xlimits, ylimits)
                end
            end
        else %user specified a channel and a cluster
            if clust<10
                fn=sprintf('ch%d_simpleclust_0%d.t', channel, clust);
                PlotTC_PSTH_single(datadir, fn, xlimits, ylimits)
            else
                fn=sprintf('ch%d_simpleclust_%d.t', channel, clust);
                PlotTC_PSTH_single(datadir, fn, xlimits, ylimits)

            end
        end
    end

end




