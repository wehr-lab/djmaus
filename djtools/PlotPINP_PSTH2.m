function PlotPINP_PSTH2(varargin)

%plots clustered spiking PINP data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%
% usage: PlotPINP_PSTH2([datapath], [cells], [xlimits],[ylimits], [binwidth])
%
% Everything is optional.
% if omitted, the following defaults are used:
%     datadir defaults to the current directory
%     cells defaults to all cells in data directory (or specify from [1:numcells])
%     xlimits defaults to [0 200]
%     ylimits defaults to an autoscaled value
%
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;
windowpos=[200 100 1381 680];
maxwindows=20; %raise a "continue?" box after this many windows to avoid crashing
if ismac  windowpos=[ -1415 479 1381 680];end %mike's pref
fs=12;
printtofile=1; %print figures to postscript file
closewindows=1; %close windows as soon as you print them

if nargin==0
    datadir=pwd;
else
    datadir=varargin{1};
    if isempty(datadir) datadir=pwd;end
end

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
    load(SortedUnitsFile)
    ProcessPINP_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits)
    load(outfilename);
end

if exist(outfilename,'file')
    load(outfilename)
    fprintf('\nloaded outfile.')
    fprintf('\n%d cells',length(out.SortedUnits))
    if isempty(cells)
        cells=1:length(out.SortedUnits);
        fprintf('\nplotting all cells')
    else
        fprintf('\nplotting cells ')
        fprintf('%d ', cells)
    end

else
    fprintf('\ncould not find outfile, calling ProcessSession...')
    if exist('./notebook.mat')==2
        %we're in Ephys folder, so go up one
        cd ..
    end
    ProcessSession
    load(SortedUnitsFile)
    ProcessPINP_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits)
    load(outfilename);
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%.1f %.1f] but xlimits in outfile are [%.1f %.1f], re-processing...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        ProcessPINP_PSTH2(out.SortedUnits, out.BonsaiPath, out.EphysPath, out.EphysPath_KS, xlimits, ylimits)
        load(outfilename);
    end
end
fprintf('\nusing xlimits [%d %d]', out.xlimits)

if isempty(ylimits)
    autoscale_ylimits=1;
else
    autoscale_ylimits=0;
end

if printtofile
    delete figs.ps
    delete figs.pdf
    close all
end

MPulse=out.MPulse;
mMPulse=out.mMPulse;
mMSilentSoundOFF=out.mMSilentSoundOFF;
MSilentSoundOFF=out.MSilentSoundOFF;
mMTrain=out.mMTrain;
MTrain=out.MTrain;
mMPulseLasertrace=out.mMPulseLasertrace;
MPulseLasertrace=out.MPulseLasertrace;
mMSilentSoundOFFLasertrace=out.mMSilentSoundOFFLasertrace;
MSilentSoundOFFLasertrace=out.MSilentSoundOFFLasertrace;
mMTrainLasertrace=out.mMTrainLasertrace;
MTrainLasertrace=out.MTrainLasertrace;
mMPulseStimtrace=out.mMPulseStimtrace;
MPulseStimtrace=out.MPulseStimtrace;
mMSilentSoundOFFStimtrace=out.mMSilentSoundOFFStimtrace;
MSilentSoundOFFStimtrace=out.MSilentSoundOFFStimtrace;
mMTrainStimtrace=out.mMTrainStimtrace;
MTrainStimtrace=out.MTrainStimtrace;

silentsounddurs=out.silentsounddurs;
trainnumpulses=out.trainnumpulses;
pulsewidths=out.pulsewidths;
trainisis=out.trainisis;
trainpulsewidths=out.trainpulsewidths;
numsilentsounddurs=out.numsilentsounddurs;
numtrainnumpulses=out.numtrainnumpulses;
numpulsewidths=out.numpulsewidths;
numtrainisis=out.numtrainisis;
numtrainpulsewidths=out.numtrainpulsewidths;

nrepsOFF=out.nrepsOFF;
nrepsPulse=out.nrepsPulse;
nrepsTrain=out.nrepsTrain;

numnexts=out.numnexts;
nexts=out.nexts;
next=min(nexts);
laserstarts=out.laserstarts;
numlaserstarts=out.numlaserstarts;


LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
samprate=out.samprate; %in Hz
if isempty(xlimits) xlimits=out.xlimits;end

for cellnum=cells
        fprintf('\ncell %d/%d', cellnum, length(cells))

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

    % %find optimal axis limits
    if autoscale_ylimits
        ymax=0;
        for dindex=1:numsilentsounddurs
            st=mMSilentSoundOFF(cellnum, dindex).spiketimes;
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(st, X);
            N=N./nrepsOFF(cellnum, dindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            ymax= max(ymax,max(N));
        end
        for pwindex=1:numpulsewidths
            st=mMPulse(cellnum, pwindex).spiketimes;
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(st, X);
            N=N./nrepsPulse(cellnum, dindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            ymax= max(ymax,max(N));
        end
        for tnpindex=1:numtrainnumpulses
            for tpwindex=1:numtrainpulsewidths
                for tiindex=1:numtrainisis
                    st=mMTrain(cellnum, tnpindex,tpwindex,tiindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(st, X);
                    N=N./nrepsOFF(cellnum, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));
                end
            end
        end


        ylimits=[-.3 ymax*2.5];
    end
    if cellnum==1 fprintf('\nusing ylimits [%.1f %.1f]', ylimits); end


    %plot the mean tuning curve silent sound no laser
    figure('position',windowpos)

    subplot1(numsilentsounddurs,1, 'Max', [.95 .9])
    p=0;
    for dindex=1:numsilentsounddurs
        p=p+1;
        subplot1(p)

        hold on
        spiketimes1=mMSilentSoundOFF(cellnum, dindex).spiketimes;
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nrepsOFF(cellnum, dindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        bar(x,N,1);
        %offset=0;
        %yl=ylim;
        yl=ylimits;
        offset=-yl(2)/2;
        inc=(yl(2)/2.5)/max(nrepsOFF);
        if rasters==1
            for n=1:nrepsOFF
                spiketimes2=MSilentSoundOFF(cellnum, dindex, n).spiketimes;
                offset=offset+inc;
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
            end
        end

        if StimRecorded
            yl(1)=yl(1)-.05*diff(yl);
            offsetS=yl(1);
            for rep=1:nrepsOFF
                Stimtrace=squeeze(MSilentSoundOFFStimtrace(dindex, rep, :))';
                Stimtrace= Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.25*diff(yl)*Stimtrace;
                t=1:length(Stimtrace);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                plot(t, Stimtrace+offsetS, 'm')
            end
        else
            line([0 0+silentsounddurs(dindex)], ylimits(1)+[0 0], 'color', [.5 .5 .5], 'linewidth', 5)
        end
        if LaserRecorded
            yl(1)=yl(1)-.05*diff(yl);
            offsetL=yl(1);
            for rep=1:nrepsOFF
                Lasertrace=squeeze(MSilentSoundOFFLasertrace(dindex, rep, :))';
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=.05*diff(yl)*Lasertrace;
                plot( t, Lasertrace+offsetL, 'c')
            end
        end
        if ~isempty(yl) ylim(yl); end

        %     xlim(xlimits) %(makes less sense for silent sound)
        xlim([-next/2 silentsounddurs(dindex)+next/2]);
        set(gca, 'fontsize', fs)
        h=title(sprintf('%s, cell %d\nSilent Sound no laser, nreps: %d-%d',out.BonsaiFolder, cellnum,min(nrepsOFF(:)),max(nrepsOFF(:))));
        set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    end
    

    %plot the psth for silent sound with single laser pulse
    figure('position',windowpos)
    subplot1(numpulsewidths,1, 'Max', [.95 .9])
    p=0;
    if numlaserstarts>1 warning('\nusing only minimum laserstart for t-test');end
    laserstart=min(laserstarts);
    if numsilentsounddurs>1 warning('\nusing only minimum silentsound duration for t-test');end
    ssdindex=1;

    for pwindex=1:numpulsewidths
        p=p+1;
        subplot1(p)

        hold on
        spiketimes1=mMPulse(cellnum, pwindex).spiketimes;
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nrepsPulse(cellnum, pwindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        bar(x,N,1);
        %offset=0;
        %yl=ylim;
        yl=ylimits;
        offset=-yl(2)/2;
        inc=(yl(2)/2.5)/max(nrepsPulse);
        if rasters==1
            for n=1:nrepsPulse(cellnum, pwindex)
                spiketimes2=MPulse(cellnum, pwindex, n).spiketimes;
                offset=offset+inc;
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
            end
        end

        if StimRecorded
            yl(1)=yl(1)-.05*diff(yl);
            offsetS=yl(1);
            for rep=1:nrepsPulse(cellnum, pwindex)
                Stimtrace=squeeze(MPulseStimtrace(pwindex, rep, :))';
                Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.25*diff(yl)*Stimtrace;
                t=1:length(Stimtrace);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                plot(t, Stimtrace+offsetS, 'm')
            end
        else
            % do nothing
        end
        if LaserRecorded
            yl(1)=yl(1)-.05*diff(yl);
            offsetL=yl(1);
            for rep=1:nrepsPulse(cellnum, pwindex)
                Lasertrace=squeeze(MPulseLasertrace(pwindex, rep, :))';
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=.05*diff(yl)*Lasertrace;
                plot( t, Lasertrace+offsetL, 'c')
            end
        else
            line([laserstart laserstart+pulsewidths(pwindex)], ylimits(1)+[0 0], 'color', 'c', 'linewidth', 5)
        end
        if ~isempty(yl) ylim(yl); end

        %xlim(xlimits)
        xlim([-200 pulsewidths(pwindex)+200]);
        set(gca, 'fontsize', fs)
        h=title(sprintf('%s, cell %d\nSilent Sound with single laser pulse, nreps: %d-%d',out.BonsaiFolder,cellnum,min(nrepsOFF(:)),max(nrepsOFF(:))));
        set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')


        %How about a t-test for effect of Laser in first 75 ms
        for rep=1:nrepsPulse(cellnum, pwindex)
            stop=laserstart+75;
            st=MPulse(cellnum, pwindex,rep).spiketimes;
            spiketimes=st(st>laserstart & st<stop); % spiketimes in region
            ON(rep)=length(spiketimes);
        end
        for rep=1:nrepsOFF(cellnum, ssdindex)
            stop=laserstart+75;
            st=MSilentSoundOFF(rep).spiketimes;
            spiketimes=st(st>laserstart & st<stop); % spiketimes in region
            OFF(rep)=length(spiketimes);
        end
        [h,p_pulse]=ttest2(ON, OFF, 'tail', 'right');
        yl=ylim;
        xl=xlim;
        text(xl(1)+25, .95*yl(2), sprintf('h=%d, p=%.4f effect of laser (t-test) 0-75ms', h,p_pulse), 'fontsize', 14)
        fprintf('\ncell %d: h=%d, p=%.4f effect of laser (1-tailed t-test) 0-75ms',cellnum, h,p_pulse)
    end

    %plot the psth for silent sound with laser train
    for tpwindex=1:numtrainpulsewidths

        figure('position',windowpos)
        subplot1(numtrainnumpulses,numtrainisis, 'Max', [.95 .9])
        p=0;
        for tnpindex=1:numtrainnumpulses
            for tiindex=1:numtrainisis
                p=p+1;
                subplot1(p)

                hold on
                spiketimes1=mMTrain(cellnum, tnpindex,tpwindex,tiindex).spiketimes;

                %use this code to plot histograms
                %             X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                %             [N, x]=hist(spiketimes1, X);
                %             N=N./nrepsTrain(tnpindex,tpwindex,tiindex); %normalize to spike rate (averaged across trials)
                %             N=1000*N./binwidth; %normalize to spike rate in Hz
                %             bar(x,N,1);

                %use this code to plot smoothed firing rate
                my_xlimits=xlimits;
                [t, fr]=GaussSmooth(spiketimes1, 5, my_xlimits);
                fr=fr./nrepsTrain(cellnum, tnpindex,tpwindex,tiindex);
                plot(t, fr, 'k')
                xlabel('time, ms'); ylabel('firing rate, Hz')



                %offset=0;
                %yl=ylim;
                yl=ylimits;
                offset=-yl(2)/2;
                inc=(yl(2)/2.5)/max(nrepsTrain(cellnum, tnpindex,tpwindex,tiindex));
                if rasters==1
                    for rep=1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex)
                        spiketimes2=MTrain(cellnum, tnpindex,tpwindex,tiindex, rep).spiketimes;
                        offset=offset+inc;
                        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    end
                end

                %compute reliability
                npulses=trainnumpulses(tnpindex);
                trainisi=trainisis(tiindex);
                P=[];
                for pnum=1:npulses
                    offset=-yl(2)/2;
                    start=0+(pnum-1)*trainisi;
                    stop=start+50;
                    for rep=1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex)
                        spiketimes2=MTrain(cellnum, tnpindex,tpwindex,tiindex, rep).spiketimes;
                        st=spiketimes2(find(spiketimes2>start & spiketimes2<stop));
                        offset=offset+inc;
                        h=plot(st, yl(2)+ones(size(st))+offset, '.r');
                        p_rep(rep)=~isempty(st);
                    end
                    P(pnum)=mean(p_rep);
                end
                str= [sprintf('reliability p=') sprintf('%.2f, ', P)];str=str(1:end-2);
                text(xl(1)+25, .95*yl(2),str, 'fontsize', 12)
                fprintf('\ncell %d: %s',cellnum, str)

                %compute latencies
                L=[];
                for pnum=1:npulses
                    %fr is 10 samples/ms
                    start=1+round(10*(0+(pnum-1)*trainisi)-10*my_xlimits(1));
                    stop=start+100*10;
                    if stop>length(fr)
                        warning('latencies clipped beacuse axis limits of data are shorter than pulse train')
                        L(pnum)=nan;
                    else
                        frwin=fr(start:stop);
                        [pk, pki]=find(frwin==max(frwin));
                    %pki is index of peak, in 10x ms relative to laser onset
                    pki=pki(1);
                    L(pnum)=pki/10;
                    end
                end

                if StimRecorded
                    yl(1)=yl(1)-.05*diff(yl);
                    offsetS=yl(1);
                    for rep=1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex)
                        Stimtrace=squeeze(MTrainStimtrace(tnpindex,tpwindex,tiindex, rep, :))';
                        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                        Stimtrace=.25*diff(yl)*Stimtrace;
                        t=1:length(Stimtrace);
                        t=1000*t/out.samprate; %convert to ms
                        t=t+out.xlimits(1); %correct for xlim in original processing call
                        plot(t, Stimtrace+offsetS, 'm')
                    end
                else
                    %do nothing
                    %line([0 0+silentsounddurs(dindex)], ylimits(1)+[0 0], 'color', [.5 .5 .5], 'linewidth', 5)
                end
                if LaserRecorded
                    yl(1)=yl(1)-.05*diff(yl);
                    offsetL=yl(1);
                    for rep=1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex)
                        Lasertrace=squeeze(MTrainLasertrace(tnpindex,tpwindex,tiindex, rep, :))';
                        Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                        Lasertrace=.05*diff(yl)*Lasertrace;
                        plot( t, Lasertrace+offsetL, 'c')
                    end
                else
                    trainisi=trainisis(tiindex);
                    trainpulsewidth=trainpulsewidths(tpwindex);
                    for pnum=1:trainnumpulses(tnpindex)
                        %                     line([laserstart laserstart+trainpulsewidth], ylimits(1)+[0 0], 'color', 'c', 'linewidth', 5)
                        line([laserstart+trainisi*(pnum-1) laserstart+trainisi*(pnum-1)+trainpulsewidth], ylimits(1)+[0 0], 'color', 'c', 'linewidth', 5)
                    end
                end
                if ~isempty(yl) ylim(yl); end

                xlim(xlimits)
                %xlim([-next/2 silentsounddurs(dindex)+next/2]);
                set(gca, 'fontsize', fs)
                h=title(sprintf('%s, cell %d\nSilent Sound with laser train, nreps: %d-%d',out.BonsaiFolder,cellnum,min(nrepsOFF(:)),max(nrepsOFF(:))));
                set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            end
        end
        for rep=1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex)
            stop=0+1000;
            st=MTrain(cellnum, :,:,:,rep).spiketimes;
            spiketimes=st(st>laserstart & st<stop); % spiketimes in region
            ON(rep)=length(spiketimes);
        end
        for rep=1:nrepsOFF(cellnum, ssdindex)
            stop=laserstart+1000;
            st=MSilentSoundOFF(cellnum, rep).spiketimes;
            spiketimes=st(st>0 & st<stop); % spiketimes in region
            OFF(rep)=length(spiketimes);
        end
        [h,p_train]=ttest2(ON, OFF, 'tail', 'right');

        %text(xl(1)+25, .95*yl(2), sprintf('h=%d, p=%.4f effect of laser (t-test) 0-1000ms', h,p_train), 'fontsize', 14)
        fprintf('\ncell %d: h=%d, p=%.4f effect of laser (1-tailed t-test) 0-1000ms',cellnum, h,p_train)

    end

    %plot reliability
    figure('position',windowpos)
    %subplot1(2,2, 'XTickL', 'All', 'YTickL', 'All')
    subplot(221)
    h=plot(1:npulses, P, 'o-');
    ylim([0 1])
    xlim([0 npulses+1])
    xlabel('pulse number')
    ylabel('reliability')
    set(h, 'markersize', 10, 'markerfacecolor', 'k')
    % str= sprintf('%.2f, ', P);str=str(1:end-2);
    % text(1, .5,'reliability: p=', 'fontsize', 10)
    % text(1, .4, str, 'fontsize', 10)
    title( sprintf('cell%d mean reliability: p=%.2f +- %.2f', cellnum, mean(P), std(P)))

    %plot latency
    subplot(222)
    h=plot(1:npulses, L, 'o-');
    ylim([0 max(L)])
    xlim([0 npulses+1])
    xlabel('pulse number')
    ylabel('latency')
    set(h, 'markersize', 10, 'markerfacecolor', 'k')
    % str= sprintf('%.2f, ', L);str=str(1:end-2);
    % text(1, 2,'latency: L=', 'fontsize', 10)
    % text(1, 1, str, 'fontsize', 10)
    title( sprintf('mean latency: L=%.2f +- %.1f', mean(L), std(L)))

    %plot latency vs. reliability
    subplot(223)
    h=plot(L,P, 'o');
    xlim([0 max(L)+1])
    ylim([0 max(P)+1])
    ylabel('reliability')
    xlabel('latency')
    set(h, 'markersize', 10, 'markerfacecolor', 'k')
        h=title(sprintf('%s, cell %d',out.BonsaiFolder, cellnum));
        set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

        if printtofile
            %print figures to postscript file
            f=findobj('type', 'figure');
            for idx=1:length(f)
                
                figure(f(idx))
                % orient landscape
                % % print figs -dpsc2 -append -bestfit
                exportgraphics(f(idx),'figs.pdf','Append',true)

                if closewindows
                    close
                end
            end
        end

end %for cellnum



