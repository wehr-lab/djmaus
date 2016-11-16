function PlotAsymGPIAS_PSTH_single(varargin)

%plots a single file of clustered spiking AsymGPIAS data from djmaus
%
% usage: PlotAsymGPIAS_PSTH(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% (xlimits, ylimits, binwidth are optional)
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

t_filename=varargin{2};

try
    xlimits=varargin{3};
catch
    xlimits=[];
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end
try
    binwidth=varargin{5};
catch
    binwidth=5;
end

if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessAsymGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits);
end

[p,f,ext]=fileparts(t_filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\n%s', t_filename)
fprintf('\n%s', outfilename)

cd(datadir)

if exist(outfilename,'file')
    load(outfilename)
else
    ProcessAsymGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
end
    
%if xlimits don't match, force preprocess
if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
    ProcessAsymGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
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
