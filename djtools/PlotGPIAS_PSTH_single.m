function PlotGPIAS_PSTH_single(varargin)

%plots a single file of clustered spiking GPIAS data from djmaus
%
% usage: PlotGPIAS_PSTH(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% (xlimits, ylimits, binwidth are optional)
%
%Processes data if outfile is not found;

plotOFFON=1;
rasters=1;
force_reprocess=1;

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
    ProcessGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
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
    fprintf('\nloaded outfile')
else
    fprintf('\ncould not find outfile, calling ProcessGPIAS_PSTH_single...')
    ProcessGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
    load(outfilename);
end

%if xlimits don't match, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%d %d] but xlimits in outfile are [%d %d], calling ProcessAsymGPIAS_PSTH_single...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        ProcessGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
        load(outfilename);
    end
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
if isfield(out, 'LaserRecorded')
    LaserRecorded=out.LaserRecorded;
    M1ONLaser=out.M1ONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
    mM1ONLaser=out.mM1ONLaser;
    M1OFFLaser=out.M1OFFLaser;
    mM1OFFLaser=out.mM1OFFLaser;
else
    LaserRecorded=0;
end
if isfield(out, 'StimRecorded')
    StimRecorded=out.StimRecorded;
    M1ONStim=out.M1ONStim;
    mM1ONStim=out.mM1ONStim;
    M1OFFStim=out.M1OFFStim;
    mM1OFFStim=out.mM1OFFStim;
else
    StimRecorded=0;
end

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
        for paindex=1:numpulseamps
            if isempty(M1OFF)
                
                st=mM1ON(gdindex,paindex).spiketimes;
                nr=nrepsON(gdindex,paindex);
            else
                st=mM1OFF(gdindex,paindex).spiketimes;
                nr=nrepsOFF(gdindex,paindex);
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

if ~isempty(M1OFF)
    
    %plot the mean tuning curve OFF
    figure('position',[600 200 600 750])
    
    p=0;
    subplot1(numgapdurs, numpulseamps, 'Max', [.95 .96])
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            
            p=p+1;
            subplot1(p)
            hold on
            spiketimes1=mM1OFF(gdindex,paindex).spiketimes; %spiketimes are in ms relative to gap termination
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(spiketimes1, X);
            N=N./nrepsOFF(gdindex,paindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            offset=0;
            yl=ylimits;
            inc=(yl(2))/max(nrepsOFF(:));
            if rasters==1
                for n=1:nrepsOFF(gdindex,paindex)
                    spiketimes2=M1OFF(gdindex,paindex, n).spiketimes;
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
            ylimits2(2)=ylimits(2)*2.2;
            ylimits2(1)=-2;
            ylim(ylimits2)
            
            xlim(xlimits)
            set(gca, 'fontsize', fs)
            %set(gca, 'xticklabel', '')
            %set(gca, 'yticklabel', '')
            if StimRecorded
                Stimtrace=squeeze(mM1OFFStim(gdindex,paindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace);
                height=.1*diff(ylimits2); %reasonable height
                %height=.5*diff(ylimits2); %magnified height for scrutinizing stimulus
                Stimtrace=height*Stimtrace;
                t=1:length(Stimtrace);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                offset=-range(Stimtrace);
                ylimits2(1)=ylimits2(1)+2*offset;
                %  offset=ylimits(1)+.1*diff(ylimits);
                plot(t, Stimtrace+offset, 'm')
                ylim(ylimits2)
            else
                %do nothing
            end
            if LaserRecorded
                height=.05*diff(ylimits2);
                offset=-height;
                for rep=1:nrepsOFF(gdindex,paindex)
                    Lasertrace=squeeze(M1OFFLaser(gdindex,paindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=height*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end
                ylimits2(1)=ylimits2(1)-2*range(Lasertrace);
                ylim(ylimits2)
            end
            
            text(xlimits(1)-10, ylimits(2), sprintf('%d', gapdurs(gdindex)), 'color', 'r')
            if gdindex<numgapdurs
                set(gca, 'yticklabel', '');
            end
        end
    end
    
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell %d, nreps: %d-%d, OFF',datadir,channel,out.cluster,min(nrepsOFF(:)),max(nrepsOFF(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %print to pdf
    print('-dpdf', sprintf('tet%d-cell%d',channel,out.cluster))

    %label amps and freqs
%     p=0;
%     for paindex=1:numpulseamps
%         p=p+1;
%         subplot1(p)
%         vpos=ylimits(2);
%         text(xlimits(1), vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r')
%         set(gca, 'yticklabel', '');
%     end

    %turn on ytick for bottom-most plot
    %set(gca, 'yticklabelmode', 'auto');
    
end           %plot the mean tuning curve OFF


if IL
    %plot the mean tuning curve ON
    figure('position',[1200 200 600 750])
    p=0;
    subplot1(numgapdurs, numpulseamps, 'Max', [.95 .9])
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            hold on
            spiketimes1=mM1ON(gdindex,paindex).spiketimes; %spiketimes are in ms relative to gap termination
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(spiketimes1, X);
            N=N./nrepsON(gdindex,paindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            offset=0;
            yl=ylimits;
            inc=(yl(2))/max(nrepsON(:));
            if rasters==1
                for n=1:nrepsON(gdindex,paindex)
                    spiketimes2=M1ON(gdindex,paindex, n).spiketimes;
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
                        
            if StimRecorded
                Stimtrace=squeeze(mM1ONStim(gdindex,paindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.05*diff(ylimits)*Stimtrace;
                t=1:length(Stimtrace);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                offset=ylimits(1)+.1*diff(ylimits);
                plot(t, Stimtrace+offset, 'm')
            else
                %do nothing
            end
            if LaserRecorded
                for rep=1:nrepsON(gdindex,paindex)
                    Lasertrace=squeeze(M1ONLaser(gdindex,paindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=.05*diff(ylimits)*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end
            end
            
        end
    end
    
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell%d, nreps: %d-%d, ON',datadir,channel,out.cluster,min(nrepsON(:)),max(nrepsON(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %label amps and freqs
    p=0;
        for gdindex=1:numgapdurs
        p=p+1;
        subplot1(p)
        vpos=ylimits(2);
        text(xlimits(1), vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r')
            if gdindex<numgapdurs
                set(gca, 'yticklabel', '');
            end
    end
    %turn on ytick for bottom-most plot
    set(gca, 'yticklabelmode', 'auto');
    

end %            %plot the mean tuning curve ON

if IL % plot the mean tuning curve OFF/ON
    if plotOFFON==1
            
    figure('position',[1200 200 600 750])
    p=0;
    subplot1(numgapdurs, numpulseamps, 'Max', [.95 .9])
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            hold on
            spiketimesON=mM1ON(gdindex,paindex).spiketimes; %spiketimes are in ms relative to gap termination
            spiketimesOFF=mM1OFF(gdindex,paindex).spiketimes; %spiketimes are in ms relative to gap termination
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [NON, xON]=hist(spiketimesON, X);
            [NOFF, xOFF]=hist(spiketimesOFF, X);
            NON=NON./nrepsON(gdindex,paindex); %normalize to spike rate (averaged across trials)
            NOFF=NOFF./nrepsOFF(gdindex,paindex); 
            NON=1000*NON./binwidth; %normalize to spike rate in Hz
            NOFF=1000*NOFF./binwidth; %normalize to spike rate in Hz
            offset=0;
            yl=ylimits;
            max_rep=max(max(nrepsON(:)), max(nrepsOFF(:)));
            inc=(yl(2))/max_rep;
            if rasters==1
                for n=1:nrepsOFF(gdindex,paindex)
                    spiketimes2=M1OFF(gdindex,paindex, n).spiketimes;
                    offset=offset+inc;
                    h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                end
                for n=1:nrepsON(gdindex,paindex)
                    spiketimes2=M1ON(gdindex,paindex, n).spiketimes;
                    offset=offset+inc;
                    h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.g');
                end
                        

            end
            bar(xON, NON,1,'facecolor','g','edgecolor','g');
            bar(xOFF, NOFF,1,'facecolor','none','edgecolor','k');
            
            if gapdurs(gdindex)>0
                line([0 0],[ylim],'color','m')
                line(-[(gapdurs(gdindex)) (gapdurs(gdindex))],[ylim],'color','m')
            end
            line(xlimits, [0 0], 'color', 'k')
            ylimits2(2)=ylimits(2)*2+offset;
            ylimits2(1)=-2;
            ylim(ylimits2)
            
            xlim(xlimits)
            set(gca, 'fontsize', fs)
            %set(gca, 'xticklabel', '')
            %set(gca, 'yticklabel', '')
                        
            if StimRecorded
                Stimtrace=squeeze(mM1ONStim(gdindex,paindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.05*diff(ylimits)*Stimtrace;
                t=1:length(Stimtrace);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                offset=ylimits(1)+.1*diff(ylimits);
                plot(t, Stimtrace+offset, 'm')
            else
                %do nothing
            end
            if LaserRecorded
                for rep=1:nrepsON(gdindex,paindex)
                    Lasertrace=squeeze(M1ONLaser(gdindex,paindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=.05*diff(ylimits)*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end
            end
            
        end
    end
    
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell%d, nreps: %d-%d, OFF/ON',datadir,channel,out.cluster,min(nrepsON(:)),max(nrepsON(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %label amps and freqs
    p=0;
        for gdindex=1:numgapdurs
        p=p+1;
        subplot1(p)
        vpos=ylimits(2);
        text(xlimits(1), vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r')
            if gdindex<numgapdurs
                set(gca, 'yticklabel', '');
            end
    end
    %turn on ytick for bottom-most plot
    set(gca, 'yticklabelmode', 'auto');
    end
    
end %plot ON OFF





