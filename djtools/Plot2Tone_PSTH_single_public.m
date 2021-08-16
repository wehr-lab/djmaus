function Plot2Tone_PSTH_single_public(datadir, clust)
%example: Plot2Tone_PSTH_single_public(pwd, 21) to plot cell 21 (outfile name) from the
%outfile located in working directory

channel=-1; %siliprobes
outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\n%s', outfilename)

cd(datadir)
ylimits=[];
binwidth=5;
rasters=1;
if exist(outfilename,'file')
    load(outfilename)
    fprintf('\nloaded outfile.')
end
SOAs=out.SOAs; % SOAs and LaserISIs are the same
%if xlimits are requested but don't match those in outfile, force preprocess
xlimits=out.xlimits;

IL=out.IL; %whether there are any interleaved laser trials
freqs=out.freqs; % only one -1 WN
probefreqs=out.probefreqs; % probe freq -1 WN or 0 silent sound
durs=out.durs; % 25 ms one dur for WN and silent sound

nreps=out.nreps;
numfreqs=out.numfreqs;
numprobefreqs=out.numprobefreqs;
numamps=out.numamps;
numdurs=out.numdurs;
numSOAs=out.numSOAs;
samprate=out.samprate; %in Hz
M1=out.M1;
mM1=out.mM1;
M1spikecounts=out.M1spikecounts;
mM1spikecount=out.mM1spikecount;
sM1spikecount=out.sM1spikecount;
semM1spikecount=out.semM1spikecount;

numLaserNumPulses=out.numLaserNumPulses;
numLaserISIs=out.numLaserISIs;
LaserNumPulses=out.LaserNumPulses;
LaserISIs=out.LaserISIs;

if isempty(xlimits) xlimits=out.xlimits;end
if numdurs>1 error('cannot handle multiple durations'), end
dindex=1;
LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
M1Laser=out.M1Laser;
mM1Laser=out.mM1Laser;
M1Stim=out.M1Stim;
mM1Stim=out.mM1Stim;

fs=10; %fontsize
% %find optimal axis limits
if isempty(ylimits)
    ymax=0;
    for findex=1:numfreqs
        for pfindex=1:numprobefreqs
            for SOAindex=1:numSOAs
                st=mM1(findex, pfindex, SOAindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(st, X);
                N=N./nreps(findex, pfindex, SOAindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));
                
            end
        end
    end
    
    ylimits=[-.3 ymax];
end

numy=numfreqs*numprobefreqs;

%plot the mean tuning curve OFF
figure('position',[650 100 600 600])
p=0;

%since SOA and laserISIs are the same, those are the two paramaters varying
%in this stimulus, i will plot them as one


% Plot OFF WN and ON silent sound laser pulses on the same plot, 2 tone
% only
subplot1(numSOAs-1,1, 'Max', [.95 .9])
for findex=[numfreqs:-1:1]
    pfindex=1;
    p=0;
    for SOAindex=2:numSOAs
        p=p+1;
        subplot1(p)
        hold on
        spiketimes=mM1(findex, pfindex, SOAindex).spiketimes; %pfindex should be 2, or -1 (WN)
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes, X);
        N=N./nreps(findex, pfindex, SOAindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        if findex==numfreqs
            offset=0;
        else
            offset=offset1;
        end
        yl=ylimits;
        inc=(yl(2))/max(max(max(nreps)));
        if rasters==1
            for n=1:nreps(findex, pfindex, SOAindex)
                spiketimes2=M1(findex, pfindex, SOAindex, n).spiketimes;
                offset=offset+inc;
                if findex==1
                    h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k'); %findex 1 = WN
                else
                    h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.g'); %findex 2= laser
                end
            end
        end
        if findex==1
            bar(x,N,1,'facecolor','none','edgecolor','k');
        else
            bar(x,N,1,'facecolor','g','edgecolor','g');
        end
        vpos=ylimits(2)-2*inc;
        
        offsetS=ylimits(1)-.1*diff(ylimits);
        if StimRecorded
            Stimtrace=squeeze(mM1Stim(findex, pfindex, SOAindex, :));
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=1.25*diff(ylimits)*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            stimh=plot(t, Stimtrace+offsetS, 'm');
            uistack(stimh, 'bottom')
            ylimits2(1)=2*offsetS;
        else
            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            ylimits2(1)=-2;
        end
        if LaserRecorded
            for rep=1:nreps(findex, pfindex, SOAindex)
                Lasertrace=squeeze(M1Laser(findex, pfindex, SOAindex,rep, :));
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=.05*diff(ylimits)*Lasertrace;
                plot( t, Lasertrace+offsetS, 'c')
            end
        end
        line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
        line(xlimits, [0 0], 'color', 'k')
        ylimits2(2)=ylimits(2)+offset;
        ylim(ylimits2)
        
        xlim(xlimits)
        set(gca, 'fontsize', fs)
        set(gcf, 'Position',  [403 74 868 911])
        if p==1
            h=title(sprintf('%s: \n cell%d dur: %dms, nreps: %d, OFF KS id: %d, uQ: %.3f, depth: %d',datadir,out.cluster,out.durs,...
                max(max(max(nreps))), out.KiloSort_ID,out.uQ, out.depth));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            
        end
        if p==1
            ylabel('FR Hz')
        end
        xl = xlim; yl = ylim;
        ylabel(sprintf('%d ms',SOAs(SOAindex)));
        %         h=text(xlimits(1), range(yl),sprintf('%d ms',SOAs(SOAindex)), 'color', 'r');
        %         set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','position',[xlimits(1) range(yl)*.8 0])
    end
    offset1=offset;
end

%plot single pulse or single WN only
figure;
for findex=[numfreqs:-1:1] %plot laser first
    pfindex=2; %single stimulus
    SOAindex=1; % 0 SOA
    
    hold on
    spiketimes=mM1(findex, pfindex, SOAindex).spiketimes; %pfindex should be 2, or -1 (WN)
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes, X);
    N=N./nreps(findex, pfindex, SOAindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    if findex==numfreqs
        offset=0;
    end
    yl=ylimits;
    inc=(yl(2))/max(max(max(nreps)));
    if rasters==1
        for n=1:nreps(findex, pfindex, SOAindex)
            spiketimes2=M1(findex, pfindex, SOAindex, n).spiketimes;
            offset=offset+inc;
            if findex==1
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k'); %findex 1 = WN
            else
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.g'); %findex 2= laser
            end
        end
    end
    if findex==1
        bar(x,N,1,'facecolor','none','edgecolor','k');
    else
        bar(x,N,1,'facecolor','g','edgecolor','g');
    end
    %vpos=mean(ylimits);
    vpos=ylimits(2)-2*inc;
    
    offsetS=ylimits(1)-.1*diff(ylimits);
    if StimRecorded
        Stimtrace=squeeze(mM1Stim(findex, pfindex, SOAindex, :));
        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
        Stimtrace=1.25*diff(ylimits)*Stimtrace;
        t=1:length(Stimtrace);
        t=1000*t/out.samprate; %convert to ms
        t=t+out.xlimits(1); %correct for xlim in original processing call
        stimh=plot(t, Stimtrace+offsetS, 'm');
        uistack(stimh, 'bottom')
        ylimits2(1)=2*offsetS;
    else
        line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
        ylimits2(1)=-2;
    end
    if LaserRecorded
        for rep=1:nreps(findex, pfindex, SOAindex)
            Lasertrace=squeeze(M1Laser(findex, pfindex, SOAindex,rep, :));
            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
            Lasertrace=.05*diff(ylimits)*Lasertrace;
            plot( t, Lasertrace+offsetS, 'c')
        end
    end
    line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
    line(xlimits, [0 0], 'color', 'k')
    ylimits2(2)=ylimits(2)+offset;
    ylim(ylimits2)
    
    xlim(xlimits)
    set(gca, 'fontsize', fs)
    h=title('WN burst and laser pulse alone');
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    ylabel('FR Hz')
    xl = xlim; yl = ylim;
%     h=text(0, range(yl),sprintf('%d ms',SOAs(SOAindex)), 'color', 'r');
%     set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','position',[0 range(yl)*.8 0])
    
end
