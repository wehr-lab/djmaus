function PlotAsymGPIAS_PSTH(varargin)

%plots clustered spiking AsymGPIAS data from djmaus
%
% usage: PlotAsymGPIAS_PSTH(datapath, [channel], [clust], [xlimits],[ylimits], [binwidth])
% (channel, clust, xlimits, ylimits are optional)
% channel number should be an integer
% clust can be an integer or an array of integers, or defaults to all clusts
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

try
    channel=varargin{2};
catch
    prompt=('please enter channel number: ');
    channel=input(prompt);
end
if strcmp('char',class(channel))
    channel=str2num(channel);
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


t_filename=sprintf('ch%s_simpleclust_0%s.t', channel, clust);


djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessAsymGPIAS_PSTH(datadir,  channel, xlimits, ylimits);
end

if isempty(clust)
    basefn=sprintf('outPSTH_ch%dc*.mat',channel);
    d=dir(basefn);
    numclusters=size(d, 1);
    if numclusters==0
        ProcessAsymGPIAS_PSTH(datadir,  channel, xlimits, ylimits);
        basefn=sprintf('outPSTH_ch%dc*.mat',channel);
        d=dir(basefn);
        numclusters=size(d, 1);
        if numclusters==0
            error('ProcessAsymGPIAS_PSTH: no cluster files found');
        end
    else fprintf('\nno cluster specified\n%d outfiles found', numclusters)
        if numclusters>1 fprintf(' -  will plot all of them');end
    end
    for clustnum=1:numclusters
        
        outfilename{clustnum}=d(clustnum).name;
    end
else %a clust was specified
    outfilename{1}=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
end

for clustindex=1:length(outfilename) %main cluster loop
    fprintf('\nclustindex=%d', clustindex)
    if exist(outfilename{clustindex},'file')
        load(outfilename{clustindex})
    else
        ProcessAsymGPIAS_PSTH(datadir,  channel, xlimits, ylimits);
        load(outfilename{clustindex});
    end
    
    IL=out.IL; %whether there are any interleaved laser trials
    numpulseamps=out.numpulseamps;
    numgapdurs=out.numgapdurs;
    pulseamps=out.pulseamps;
    gapdurs=out.gapdurs;
    gapdelay=out.gapdelay;
    samprate=out.samprate; %in Hz
    mM1ON=out.mM1ON;
    mM1OFF=out.mM1OFF;
    M1ON=out.M1ON;
    M1OFF=out.M1OFF;
    nrepsON=out.nrepsON;
    nrepsOFF=out.nrepsOFF;
    soa=out.soa;
    onramps=out.onramps;
    offramps=out.offramps;
    numonramps=out.numonramps;
    numofframps=out.numofframps;
    
    
    fs=10; %fontsize
    
    % %find optimal axis limits
    if isempty(xlimits)
        xlimits(1)=-1.5*max(gapdurs);
        xlimits(2)=2*soa;
    end
    fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))
    
    if isempty(ylimits)
        ymax=0;
        for gdindex=1:numgapdurs
            for onrampindex=1:numonramps
                for offrampindex=1:numofframps
                    if isempty(M1OFF)
                        
                        st=mM1ON(gdindex,onrampindex, offrampindex).spiketimes;
                        nr=nrepsON(gdindex,onrampindex, offrampindex);
                    else
                        st=mM1OFF(gdindex,onrampindex, offrampindex).spiketimes;
                        nr=nrepsOFF(gdindex,onrampindex, offrampindex);
                    end
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(st, X);
                    N=N./nr; %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));
                end
            end
            ylimits=[-.3 ymax];
        end
    end
    
    if ~isempty(M1OFF)
        
        %plot the mean tuning curve OFF
        for gdindex=1:numgapdurs
            figure
            p=0;
            subplot1(numonramps, numofframps, 'Max', [.95 .9], 'Gap', [.01 .01])
            for onrampindex=1:numonramps
                for offrampindex=1:numofframps
                    
                    p=p+1;
                    subplot1(p)
                    hold on
                    spiketimes1=mM1OFF(gdindex,onrampindex, offrampindex).spiketimes; %spiketimes are in ms relative to gap termination
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nrepsOFF(gdindex,onrampindex, offrampindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=0;
                    yl=ylimits;
                    inc=(yl(2))/max(nrepsOFF(:));
                    if rasters==1
                        for n=1:nrepsOFF(gdindex,onrampindex, offrampindex)
                            spiketimes2=M1OFF(gdindex,onrampindex, offrampindex, n).spiketimes;
                            offset=offset+inc;
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    bar(x, N,1,'facecolor','k','edgecolor','k');
                    
                    if gapdurs(gdindex)>0
                        line([0 0],[ylim],'color','m')
                        line(-[(gapdurs(gdindex)) (gapdurs(gdindex))],[ylim],'color','m')
                    end
                    line(xlimits, [0 0], 'color', 'k')
                    ylimits2(2)=ylimits(2)*3;
                    ylimits2(1)=-2;
                    ylim(ylimits2)
                    
                    xlim(xlimits)
                    set(gca, 'fontsize', fs)
                    %set(gca, 'xticklabel', '')
                    %set(gca, 'yticklabel', '')
                    %title(sprintf('onramp: %d, offramp: %d', onramps(onrampindex), offramps(offrampindex)))

                end
            end
            
            subplot1(1)
            h=title(sprintf('%s: \ntetrode%d cell %d, nreps: %d-%d, OFF',datadir,channel,out.cluster,min(nrepsOFF(:)),max(nrepsOFF(:))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            
            %label amps and freqs
            p=0;
            for onrampindex=1:numonramps
                for offrampindex=1:numofframps
                    
                    p=p+1;
                    subplot1(p)
                    vpos=ylimits(2);
                    text(xlimits(1), vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r')
                    vpos=ylimits2(2)*.9;
                    text(0, vpos, sprintf('%d', offramps(offrampindex)), 'color', 'r')
                    gapdur=gapdurs(gdindex);
                    text(-gapdur, vpos, sprintf('%d', onramps(onrampindex)), 'color', 'r')
                    set(gca, 'yticklabel', '');
                end
            end
            %turn on ytick for bottom-most plot
            set(gca, 'yticklabelmode', 'auto');
            
        end
    end
    
    
    if IL
        %plot the mean tuning curve ON
        for gdindex=1:numgapdurs
            figure
            p=0;
            subplot1(numonramps, numofframps, 'Max', [.95 .9], 'Gap', [.01 .01])
            for onrampindex=1:numonramps
                for offrampindex=1:numofframps
                    p=p+1;
                    subplot1(p)
                    hold on
                    spiketimes1=mM1ON(gdindex,onrampindex, offrampindex).spiketimes; %spiketimes are in ms relative to gap termination
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nrepsON(gdindex,onrampindex, offrampindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=0;
                    yl=ylimits;
                    inc=(yl(2))/max(nrepsON(:));
                    if rasters==1
                        for n=1:nrepsON(gdindex,onrampindex, offrampindex)
                            spiketimes2=M1ON(gdindex,onrampindex, offrampindex, n).spiketimes;
                            offset=offset+inc;
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    bar(x, N,1,'facecolor','g','edgecolor','k');
                    
                    if gapdurs(gdindex)>0
                        line([0 0],[ylim],'color','m')
                        line(-[(gapdurs(gdindex)) (gapdurs(gdindex))],[ylim],'color','m')
                    end
                    line(xlimits, [0 0], 'color', 'k')
                    ylimits2(2)=ylimits(2)*3;
                    ylimits2(1)=-2;
                    ylim(ylimits2)
                    
                    xlim(xlimits)
                    set(gca, 'fontsize', fs)
                    %set(gca, 'xticklabel', '')
                    %set(gca, 'yticklabel', '')
                    
                  %  title(sprintf('onramp: %d, offramp: %d', onramps(onrampindex), offramps(offrampindex)))
                end
            end
            
            subplot1(1)
            h=title(sprintf('%s: \ntetrode%d cell%d, nreps: %d-%d, ON',datadir,channel,out.cluster,min(nrepsON(:)),max(nrepsON(:))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            
            %label amps and freqs
            p=0;
            for onrampindex=1:numonramps
                for offrampindex=1:numofframps
                    p=p+1;
                    subplot1(p)
                    vpos=ylimits(2);
                    text(xlimits(1), vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r')
                    set(gca, 'yticklabel', '');
                    vpos=ylimits2(2)*.9;
                    text(0, vpos, sprintf('%d', offramps(offrampindex)), 'color', 'r')
                    gapdur=gapdurs(gdindex);
                    text(-gapdur, vpos, sprintf('%d', onramps(onrampindex)), 'color', 'r')
                end
            end
        end
        %turn on ytick for bottom-most plot
        set(gca, 'yticklabelmode', 'auto');
        
    end %plot ON
    
    
    
end %main cluster loop
