function PlotTC_PSTH(varargin)

%plots clustered spiking tuning curve data from djmaus
%
% usage: PlotTC_PSTH(datapath, [channel], [clust], [xlimits],[ylimits], [binwidth])
% (channel, clust, xlimits, ylimits are optional)
% xlimits default to [0 200]
% channel number should be an integer
% clust can be an integer or an array of integers, or defaults to all clusts
%
%Processes data if outfile is not found;

rasters=1;

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
    xlimits=[0 200];
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

fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

if isempty(clust)
    basefn=sprintf('outPSTH_ch%dc*.mat',channel);
    d=dir(basefn);
    numclusters=size(d, 1);
    if numclusters==0
        ProcessTC_PSTH(datadir,  channel, xlimits, ylimits);
        basefn=sprintf('outPSTH_ch%dc*.mat',channel);
        d=dir(basefn);
        numclusters=size(d, 1);
        if numclusters==0
            error('ProcessTC_PSTH: no cluster files found');
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
        ProcessTC_PSTH(datadir,  channel, xlimits, ylimits);
        load(outfilename{clustindex});
    end
    
    IL=out.IL; %whether there are any interleaved laser trials
    freqs=out.freqs;
    amps=out.amps;
    durs=out.durs;
    nreps=out.nreps;
    numfreqs=out.numfreqs;
    numamps=out.numamps;
    numdurs=out.numdurs;
    samprate=out.samprate; %in Hz
    mM1ON=out.mM1ON;
    mM1OFF=out.mM1OFF;
    M1ON=out.M1ON;
    M1OFF=out.M1OFF;
    nrepsON=out.nrepsON;
    nrepsOFF=out.nrepsOFF;
    
    fs=10; %fontsize
    
    % %find optimal axis limits
    if isempty(ylimits)
        ymax=0;
        for aindex=[numamps:-1:1]
            for findex=1:numfreqs
                for dindex=1:numdurs
                    st=mM1OFF(findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(st, X);
                    N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));
                end
            end
        end
        ylimits=[-.3 ymax];
    end
    
    
    
    %plot the mean tuning curve OFF
    for dindex=1:numdurs
        figure
        p=0;
        subplot1(numamps,numfreqs, 'Max', [.95 .9])
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                hold on
                spiketimes1=mM1OFF(findex, aindex, dindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(spiketimes1, X);
                N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                offset=0;
                yl=ylimits;
                inc=(yl(2))/max(max(max(nreps)));
                if rasters==1
                    for n=1:nrepsOFF(findex, aindex, dindex)
                        spiketimes2=M1OFF(findex, aindex, dindex, n).spiketimes;
                        offset=offset+inc;
                        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    end
                end
                bar(x, N,1);
                line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
                line(xlimits, [0 0], 'color', 'k')
                ylimits2(2)=ylimits(2)*3;
                ylimits2(1)=-2;
                ylim(ylimits2)
                
                xlim(xlimits)
                set(gca, 'fontsize', fs)
                %set(gca, 'xticklabel', '')
                %set(gca, 'yticklabel', '')
                
            end
        end
        subplot1(1)
        h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, OFF',datadir,channel,out.cluster,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
        
        %label amps and freqs
        p=0;
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                if findex==1
                    text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                end
            end
        end
    end
    
    if IL
        %plot the mean tuning curve ON
        for dindex=1:numdurs
            figure
            p=0;
        subplot1(numamps,numfreqs, 'Max', [.95 .9])
        for aindex=numamps:-1:1
                for findex=1:numfreqs
                    p=p+1;
                    subplot1(p)
                    hold on
                    spiketimes1=mM1ON(findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=0;
                    yl=ylimits;
                    inc=(yl(2))/max(max(max(nreps)));
                    if rasters==1
                        for n=1:nrepsON(findex, aindex, dindex)
                            spiketimes2=M1ON(findex, aindex, dindex, n).spiketimes;
                            offset=offset+inc;
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    bar(x, N,1);
                    line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
                    line(xlimits, [0 0], 'color', 'k')
                    ylimits2(2)=ylimits(2)*3;
                    ylimits2(1)=-2;
                    ylim(ylimits2(:))
                    
                    xlim(xlimits)
                    set(gca, 'fontsize', fs)
                    set(gca, 'xticklabel', '')
                    set(gca, 'yticklabel', '')
                    
                end
            end
            subplot1(1)
            h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, ON',datadir,channel,out.cluster,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            

            %label amps and freqs
            p=0;
            for aindex=numamps:-1:1
                for findex=1:numfreqs
                    p=p+1;
                    subplot1(p)
                    if findex==1
                        text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                    end
                    if aindex==1
                        if mod(findex,2) %odd freq
                            vpos=ylimits(1)-mean(ylimits);
                        else
                            vpos=ylimits(1)-mean(ylimits);
                        end
                        text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                    end
                end
            end
        end
    end
    
end %main cluster loop
