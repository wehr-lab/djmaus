function PlotClicktrain_PSTH2(varargin)

%plots clustered spiking Clicktrain data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%
% usage: PlotClicktrain_PSTH2([datapath], [cells], [xlimits],[ylimits], [binwidth])
% all inputs are optional
% defaults to datapath=pwd, xlimits and ylimits autoscaled, binwidth=5ms
% cells defaults to all cells in data directory (or specify from [1:numcells])
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;
windowpos=[200 100 1381 680];
if ismac  windowpos=[ -1415 479 1381 680];end %mike's pref
maxwindows=20; %raise a "continue?" box after this many windows to avoid crashing
fs=12; %fontsize
printtofile=1; %print figures to postscript file
closewindows=1; %close windows as soon as you print them

if nargin==0
    datadir=pwd;
else
    datadir=varargin{1};
end
if isempty(datadir) datadir=pwd;end

try
    cells=varargin{2};
catch
    cells=[];
end
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
    if isempty(binwidth), binwidth=5;end

catch
    binwidth=5;
end

outfilename='outPSTH.mat';
cd(datadir)

if force_reprocess
    fprintf('\nForce re-process\n')
    warning('Force re-process')
    if exist('./notebook.mat')==2
        %we're in Ephys folder, so go up one
        cd ..
    end
    ProcessSession
    if ~exist('SortedUnitsFile') fprintf('\nNo Sorted Units File, probably because there is no kilosort data. \nPlotTC_PSTH2 will fail.'); end
    load(SortedUnitsFile)
    ProcessClicktrain_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits)
    load(outfilename);
end

if exist(outfilename,'file')
    load(outfilename)
    fprintf('\nloaded outfile.')
else
    fprintf('\ncould not find outfile, calling ProcessSession...')
    if exist('./notebook.mat')==2
        %we're in Ephys folder, so go up one
        cd ..
    end
    ProcessSession
    if ~exist('SortedUnitsFile') fprintf('\nNo Sorted Units File, probably because there is no kilosort data. PlotTC_PSTH2 will fail.'); end
    load(SortedUnitsFile)
    ProcessClicktrain_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits)
    load(outfilename);
end

if isempty(ylimits)
    autoscale_ylimits=1;
else
    autoscale_ylimits=0;
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%.1f %.1f] but xlimits in outfile are [%.1f %.1f], re-processing...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        ProcessClicktrain_PSTH2(out.SortedUnits, out.BonsaiPath, out.EphysPath, out.EphysPath_KS, xlimits, ylimits)
        load(outfilename);
    end
else %xlimits is empty
    xlimits=out.xlimits;
end
fprintf('\nusing xlimits [%d %d]', xlimits)
fprintf('\nusing binwidth %d', binwidth)
fprintf('\n%d cells',length(out.SortedUnits))
if isempty(cells)
    cells=1:length(out.SortedUnits);
    fprintf('\nplotting all cells')
else
    fprintf('\nplotting cells ')
    fprintf('%d ', cells)
end

if printtofile
    pdffilename=sprintf('%s-figs.pdf', out.BonsaiFolder);
    delete(pdffilename)
    close all
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
mMtONspikecount=out.mMtONspikecount;
sMtONspikecount=out.sMtONspikecount;
mMtOFFspikecount=out.mMtOFFspikecount;
sMtOFFspikecount=out.sMtOFFspikecount;
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


% %find optimal axis limits
if isempty(xlimits)
    xlimits(1)=-.5*max(durs);
    xlimits(2)=1.5*max(durs);
end
fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
try
    M_LaserStart=out.M_LaserStart;
    M_LaserWidth=out.M_LaserWidth;
    M_LaserNumPulses=out.M_LaserNumPulses;
    M_LaserISI=out.M_LaserISI;
end
try
    LaserStart=out.LaserStart;
    LaserWidth=out.LaserWidth;
    LaserNumPulses=out.LaserNumPulses;
    LaserISI=out.LaserISI;
end
MtONStim=out.MtONStim;
MtONLaser=out.MtONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
mMtONStim=out.mMtONStim;
% kluge below Kip 4/21
try
    mM1ONLaser=out.mM1ONLaser;
    M1OFFLaser=out.M1OFFLaser;
    mM1OFFLaser=out.mM1OFFLaser;
catch
    LaserRecorded=0;
end
% kluge below Kip 4/21
try
    mM1OFFStim=out.mM1OFFStim;
    M1OFFStim=out.M1OFFStim;
catch
    StimRecorded = 0;
end

%hardcoding for probe P128-2
distance=20;
chans_per_shank=64;

for cellnum=cells
    tic
    fprintf('\ncell %d/%d', cellnum, cells(end))
    fprintf('\n%d spikes',  length(out.SortedUnits(cellnum).spiketimes))

    if ismac %you could check this on windows, too, in case of excessive figure windows
        f=findobj('type', 'figure');
        if length(f)>maxwindows
            windowcheck=questdlg(sprintf('you''ve exceeded %d figure windows, continue?', maxwindows), 'lots of windows', 'continue', 'abort', 'continue');
            switch windowcheck
                case 'abort'
                    fprintf('\nuser aborted due to too many figure windows')
                    return
            end
        end
    end

    chan=out.SortedUnits(cellnum).channel;
    if chan<=chans_per_shank
        shank=1;
        depth=distance*chan;
    else
        shank=2;
        depth=distance*chan-chans_per_shank;
    end
    fprintf('\ncell %d, chan %d, shank %d, raw depth %d', cellnum, chan, shank, depth)

    % %find optimal axis limits
    if autoscale_ylimits
        ymax=0;
        for iciindex=1:numicis
            st=mMtOFF(cellnum,iciindex).spiketimes;
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(st, X);
            N=N./nrepsOFF(cellnum, iciindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            ymax= max(ymax,max(N));

            if IL
                st=mMtON(cellnum,iciindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(st, X);
                N=N./nrepsON(cellnum, iciindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));
            end

        end
        if isempty(ymax)
            ymax=0;
            fprintf('\nmaybe this cell never fired a spike ?');
        end
        ylimits=[-.3 ymax];
    end
    fprintf('\nusing ylimits [%.1f %.1f]', ylimits);

    %plot the mean tuning curve OFF


    figure('position',windowpos, 'PaperOrientation', 'landscape')
    p=0;
    subplot1(numicis,1, 'Max', [.95 .9])
    for iciindex=1:numicis
        p=p+1;
        subplot1(p)
        hold on
        spiketimes1=mMtOFF(cellnum,iciindex).spiketimes; %spiketimes are in ms relative to gap termination
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nrepsOFF(cellnum,iciindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        offset=0;
        yl=ylimits;
        inc=(yl(2))/max(nrepsOFF(:));
        if rasters==1
            for n=1:nrepsOFF(cellnum,iciindex)
                spiketimes2=MtOFF(cellnum,iciindex, n).spiketimes;
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
        if nclicks(iciindex) < 1500 %past a certain number of clicks there's no point
            for k=1:nclicks(iciindex)
                clickonset=(k-1)*(icis(iciindex));
                height=.05*diff(ylimits2); %reasonable height
                offset=ylimits2(1);
                c=       [.5 .5 1] ;
                l=line(clickonset*[1 1], [offset offset+height], 'color', 'm', 'linewidth', 2);
                L=[L l];
            end
        end


        if LaserRecorded & IL %plot laser square pulses
            height=.05*diff(ylimits2);
            offset=-height;
            for rep=1:nrepsOFF(cellnum,iciindex)
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
    h=title(sprintf('%s: \ncell %d, chan %d, shank %d, raw depth %d um, %d spikes, nreps: %d-%d, OFF ',datadir,cellnum, chan, shank, depth, length(out.SortedUnits(cellnum).spiketimes),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
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

    % %lengthen figure
    % pos=get(gcf, 'pos');
    % pos(4)=900;
    % pos(2)=78;
    % set(gcf, 'pos',pos)


    if IL
        %plot the mean tuning curve ON
        if show_plots == 0
            figure('Visible','off')
        else
            figure
        end
        p=0;
        subplot1(numicis, 1, 'Max', [.95 .9])
        for iciindex=1:numicis
            p=p+1;
            subplot1(p)
            hold on
            spiketimes1=mMtON(cellnum,iciindex).spiketimes; %spiketimes are in ms relative to gap termination
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(spiketimes1, X);
            N=N./nrepsON(cellnum,iciindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            offset=0;
            yl=ylimits;
            inc=(yl(2))/max(nrepsON(:));
            if rasters==1
                for n=1:nrepsON(cellnum,iciindex)
                    spiketimes2=MtON(cellnum,iciindex, n).spiketimes;
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
                Stimtrace=squeeze(mMtONStim(iciindex, :));
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
          

            if LaserRecorded
                for rep=1:nrepsOFF(cellnum,iciindex)
                    Lasertrace=squeeze(MtONLaser(iciindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=.05*diff(ylimits)*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end
            end
        end

        subplot1(1)
        h=title(sprintf('%s: \nchannel%d cell%d, nreps: %d-%d, ON',datadir,channel,out.cluster,min(nrepsON(:)),max(nrepsON(:))));
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

        %lengthen figure
        pos=get(gcf, 'pos');
        pos(4)=900;
        pos(2)=78;
        set(gcf, 'pos',pos)

    end %            %plot the mean tuning curve ON

    % plot cycle histograms
        figure('position',windowpos, 'PaperOrientation', 'landscape')
    yl=0;
    subplot1(numicis, 2)
    p=0;
    for iciindex=[1:numicis]
        phaseON=out.PhaseON(cellnum,iciindex).phase;
        phaseOFF=out.PhaseOFF(cellnum,iciindex).phase;
        [NON xON]=hist(phaseON, [0:pi/10:2*pi]); hold on;
        [NOFF xOFF]=hist(phaseOFF, [0:pi/10:2*pi]);
        p=p+1; subplot1(p)
        bOFF=bar(xOFF, NOFF,1);
        xlim([0 2*pi])
        yl=max(yl, ylim);
        ylabel(sprintf('ici %dms',icis(iciindex)));

        p=p+1; subplot1(p)
        bON=bar(xON, NON,1);
        set(bON, 'facecolor', 'c','edgecolor', 'c');
        set(bOFF, 'facecolor', 'k','edgecolor', 'k');
        xlim([0 2*pi])
        yl=max(yl, ylim);
        %ylabel(sprintf('ici %dms',icis(iciindex)));
    end

    for n=1:p
        subplot1(p)
        ylim(yl)
    end
    xlabel('phase')
    subplot1(1)
    h=title(sprintf('cycle spike histogram, %s: \ncell %d, chan %d, shank %d, raw depth %d um, %d spikes, nreps: %d-%d, OFF ',datadir,cellnum, chan, shank, depth, length(out.SortedUnits(cellnum).spiketimes),min(nrepsOFF(:)),max(nrepsOFF(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    subplot1([1,2])
    h=title('ON');
set(h, 'color', 'c')

% pos=get(gcf, 'pos');
    % pos(4)=900;
    % pos(2)=78; %moved figure down
    % set(gcf, 'pos',pos)

    % plot vector strength, rayleigh statistic, etc
    figure('position',windowpos, 'PaperOrientation', 'landscape')
    subplot1(3, 1)
    subplot1(1)
    plot(1:numicis, out.VsOFF(cellnum, :), 'k-o', 1:numicis, out.VsON(cellnum, :), 'c-o')
    ylabel('vector strength')
    set(gca, 'xticklabel', icis)
    xlabel('ici, ms')
    h=title(sprintf('phase-locking statistics, %s: \ncell %d, chan %d, shank %d, raw depth %d um, %d spikes, nreps: %d-%d, OFF ',datadir,cellnum, chan, shank, depth, length(out.SortedUnits(cellnum).spiketimes),min(nrepsOFF(:)),max(nrepsOFF(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')


    subplot1(2)
    plot(1:numicis, out.RZOFF(cellnum, :), 'k-o', 1:numicis, out.RZON(cellnum, :), 'c-o')
    ylabel('Rayleigh statistic')
    set(gca, 'xticklabel', icis)
    xlabel('ici, ms')
    line(xlim, 13*[1 1], 'linestyle', '--')

    subplot1(3)
    semilogy(1:numicis, out.p_OFF(cellnum, :), 'k-o', 1:numicis, out.p_ON(cellnum, :), 'c-o')
    set(gca, 'yscale', 'log')
    ylabel('p-value of vector strength')
    set(gca, 'xticklabel', icis)
    xlabel('ici, ms')
    ylim([0 1])
    line(xlim, .001*[1 1], 'linestyle', '--')

    % pos=get(gcf, 'pos');
    % pos(2)=420;
    % pos(4)=595;
    % set(gcf, 'pos',pos)


    figure('position',windowpos, 'PaperOrientation', 'landscape')
    hold on
    if ~isempty(mMtONspikecount)
        e=errorbar(1:numicis, mMtONspikecount(cellnum, :), sMtONspikecount(cellnum, :), 'c-o');
    end
    if ~isempty(mMtOFFspikecount)
        e=errorbar(1:numicis, mMtOFFspikecount(cellnum, :), sMtOFFspikecount(cellnum, :), 'k-o');
    end
    set(gca, 'xtick', 1:numicis, 'xticklabel', icis)
    xlabel('ici, ms')
    ylabel('spike count')
    h=title(sprintf('spike count, %s: \ncell %d, chan %d, shank %d, raw depth %d um, %d spikes, nreps: %d-%d, OFF ',datadir,cellnum, chan, shank, depth, length(out.SortedUnits(cellnum).spiketimes),min(nrepsOFF(:)),max(nrepsOFF(:))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

    if printtofile
        %print figures to postscript file
        pdffilename=sprintf('%s-figs.pdf', out.BonsaiFolder);
        f=findobj('type', 'figure');
        for idx=1:length(f)
            exportgraphics(f(idx),pdffilename,'Append',true)
        end
        if closewindows %doing this in a separate step to ensure matlab prints all figs first
            for idx=1:length(f)
                close(f(idx))
            end
        end
    end
    fprintf('\telapsed time %.0fs', toc)

end %for cellnum








