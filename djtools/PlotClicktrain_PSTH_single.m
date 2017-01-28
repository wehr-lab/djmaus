function PlotClicktrain_PSTH_single(varargin)

%plots a single file of clustered spiking Clicktrain data from djmaus
%
% usage: PlotClicktrain_PSTH_single(datapath, t_filename, [xlimits],[ylimits], [binwidth])
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
    ProcessClicktrain_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
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
    fprintf('\ncould not find outfile, calling ProcessClicktrain_PSTH_single...')
    ProcessClicktrain_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
    load(outfilename);
end

%if xlimits don't match, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%d %d] but xlimits in outfile are [%d %d], calling ProcessClicktrain_PSTH_single...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        ProcessClicktrain_PSTH_single(datadir,  t_filename, xlimits, ylimits, binwidth);
        load(outfilename);
    end
end

IL=out.IL; %whether there are any interleaved laser trials
numicis=out.numicis;
icis=out.icis;
numdurs=out.numdurs;
durs=out.durs;
nclicks=out.nclicks;
numnclicks=out.numnclicks;
samprate=out.samprate; %in Hz
mMtON=out.mMtON;
mMtOFF=out.mMtOFF;
MtON=out.MtON;
MtOFF=out.MtOFF;
nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;
if isfield(out, 'LaserRecorded')
    LaserRecorded=out.LaserRecorded;
    MtONLaser=out.MtONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
    mMtONLaser=out.mMtONLaser;
    MtOFFLaser=out.MtOFFLaser;
    mMtOFFLaser=out.mMtOFFLaser;
else
    LaserRecorded=0;
end
if isfield(out, 'StimRecorded')
    StimRecorded=out.StimRecorded;
    MtONStim=out.MtONStim;
    mMtONStim=out.mMtONStim;
    MtOFFStim=out.MtOFFStim;
    mMtOFFStim=out.mMtOFFStim;
else
    StimRecorded=0;
end

fs=10; %fontsize

% %find optimal axis limits
if isempty(xlimits)
    xlimits(1)=-.5*max(durs);
    xlimits(2)=1.5*max(durs);
end
fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

if isempty(ylimits)
    ymax=0;
    for iciindex=1:numicis
        if isempty(MtOFF)
            st=mMtON(iciindex).spiketimes;
            nr=nrepsON(iciindex);
        else
            st=mMtOFF(iciindex).spiketimes;
            nr=nrepsOFF(iciindex);
        end
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(st, X);
        N=N./nr; %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        ymax= max(ymax,max(N));
    end
    ylimits=[-.3 ymax];
end

if ~isempty(MtOFF)
    
    %plot the mean tuning curve OFF
    figure
    p=0;
    subplot1(numicis,1, 'Max', [.95 .9])
    for iciindex=1:numicis
        p=p+1;
        subplot1(p)
        hold on
        spiketimes1=mMtOFF(iciindex).spiketimes; %spiketimes are in ms relative to gap termination
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nrepsOFF(iciindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        offset=0;
        yl=ylimits;
        inc=(yl(2))/max(nrepsOFF(:));
        if rasters==1
            for n=1:nrepsOFF(iciindex)
                spiketimes2=MtOFF(iciindex, n).spiketimes;
                offset=offset+inc;
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
            end
        end
        bar(x, N,1,'facecolor','k','edgecolor','k');
        line(xlimits, [0 0], 'color', 'k')
        ylimits2(2)=ylimits(2)*2.2;
        ylimits2(1)=-2;
        ylim(ylimits2)
                xlim(xlimits)
        set(gca, 'fontsize', fs)

        if StimRecorded %plot actual recorded stimulus trace in magenta
            Stimtrace=squeeze(mMtOFFStim(iciindex, :));
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
        end
        
        %draw lines where clicks should be
        L=[];
        for k=1:nclicks(iciindex)
            clickonset=(k-1)*(icis(iciindex));
            height=.08*diff(ylimits2); %reasonable height
            offset=ylimits2(1);
            c=       [.5 .5 1] ;
            l=line(clickonset*[1 1], [offset offset+height], 'color', c, 'linewidth', 2);
            L=[L l];
        end
        %you could comment out this line to save time -- it takes ~1 sec
        uistack(L, 'bottom')
        
        if LaserRecorded & IL %plot laser square pulses
            height=.05*diff(ylimits2);
            offset=-height;
            for rep=1:nrepsOFF(iciindex)
                Lasertrace=squeeze(MtOFFLaser(iciindex,rep, :));
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=height*Lasertrace;
                plot( t, Lasertrace+offset, 'c')
            end
            ylimits2(1)=ylimits2(1)-2*range(Lasertrace);
            ylim(ylimits2)
        end
    end
    ymin=inf;
    ymax=0;
    for iciindex=1:numicis
        subplot1(iciindex)
        yl=ylim;
        if yl(1)<ymin ymin=yl(1);end
        if yl(2)>ymax ymax=yl(2);end
    end
    for iciindex=1:numicis
        subplot1(iciindex)
        ylim([ymin ymax])
    end
    
    
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell %d, nreps: %d-%d, OFF',datadir,channel,out.cluster,min(nrepsOFF(:)),max(nrepsOFF(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %label amps and freqs
    p=0;
    for iciindex=1:numicis
        p=p+1;
        subplot1(p)
        vpos=.8*ylimits(2);
        text(xlimits(1)+100, vpos, sprintf('%d', icis(iciindex)), 'color', 'r', 'fontsiz', 14)
        set(gca, 'yticklabel', '');
    end
    
    %turn on ytick for bottom-most plot
    set(gca, 'yticklabelmode', 'auto');
    
end           %plot the mean tuning curve OFF


if IL
    %plot the mean tuning curve ON
    figure
    p=0;
    subplot1(numicis, 'Max', [.95 .9])
    for iciindex=1:numicis
        p=p+1;
        subplot1(p)
        hold on
        spiketimes1=mMtON(iciindex).spiketimes; %spiketimes are in ms relative to gap termination
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nrepsON(iciindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        offset=0;
        yl=ylimits;
        inc=(yl(2))/max(nrepsON(:));
        if rasters==1
            for n=1:nrepsON(iciindex)
                spiketimes2=MtON(iciindex, n).spiketimes;
                offset=offset+inc;
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
            end
        end
        bar(x, N,1,'facecolor','g','edgecolor','k');
        line(xlimits, [0 0], 'color', 'k')
        ylimits2(2)=ylimits(2)*3;
        ylimits2(1)=-2;
        ylim(ylimits2)
        
        xlim(xlimits)
        set(gca, 'fontsize', fs)
        
        if StimRecorded
            Stimtrace=squeeze(mMtOFFStim(iciindex, :));
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=.05*diff(ylimits)*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            offset=ylimits(1)+.1*diff(ylimits);
            plot(t, Stimtrace+offset, 'm')
        else
            %do nothing
            keyboard
        end
        if LaserRecorded
            for rep=1:nrepsOFF(iciindex)
                Lasertrace=squeeze(MtOFFLaser(iciindex,rep, :));
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=.05*diff(ylimits)*Lasertrace;
                plot( t, Lasertrace+offset, 'c')
            end
        end
    end
    
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell%d, nreps: %d-%d, ON',datadir,channel,out.cluster,min(nrepsON(:)),max(nrepsON(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %label amps and freqs
    p=0;
    for iciindex=1:numicis
        p=p+1;
        subplot1(p)
        vpos=ylimits(2);
        text(xlimits(1), vpos, sprintf('%d', icis(iciindex)), 'color', 'r')
    end
    %turn on ytick for bottom-most plot
    set(gca, 'yticklabelmode', 'auto');
    
end %            %plot the mean tuning curve ON




