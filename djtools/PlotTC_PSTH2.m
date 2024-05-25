function PlotTC_PSTH2(varargin)

%plots a  clustered spiking tuning curve data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%
% usage: PlotTC_PSTH2([datapath], [cells], [xlimits],[ylimits], [binwidth])
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
    load(SortedUnitsFile)
    ProcessTC_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits)
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
    load(SortedUnitsFile)
    ProcessTC_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits)
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
        ProcessTC_PSTH2(out.SortedUnits, out.BonsaiPath, out.EphysPath, out.EphysPath_KS, xlimits, ylimits)
        load(outfilename);
    end
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
    delete figs.ps
    delete figs.pdf
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
if isempty(xlimits) xlimits=out.xlimits;end

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
M1ONStim=out.M1ONStim;
M1ONLaser=out.M1ONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
mM1ONStim=out.mM1ONStim;
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

%silent sound
SilentSoundON=out.SilentSoundON;
SilentSoundOFF=out.SilentSoundOFF;
SilentSoundONspikecount=out.SilentSoundONspikecount;
SilentSoundOFFspikecount=out.SilentSoundOFFspikecount;
mSilentSoundON=out.mSilentSoundON;
mSilentSoundOFF=out.mSilentSoundOFF;
SilentSoundONStim=out.SilentSoundONStim;
SilentSoundOFFStim=out.SilentSoundOFFStim;
SilentSoundONLaser=out.SilentSoundONLaser;
SilentSoundOFFLaser=out.SilentSoundOFFLaser;
nreps_ssON=out.nreps_ssON;
nreps_ssOFF=out.nreps_ssOFF;


for cellnum=cells
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

    % %find optimal axis limits
    if autoscale_ylimits
        ymax=0;
        for aindex=[numamps:-1:1]
            for findex=1:numfreqs
                for dindex=1:numdurs
                    st=mM1OFF(cellnum,findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(st, X);
                    N=N./nreps(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));

                    if IL
                        st=mM1ON(cellnum,findex, aindex, dindex).spiketimes;
                        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                        [N, x]=hist(st, X);
                        N=N./nreps(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                        N=1000*N./binwidth; %normalize to spike rate in Hz
                        ymax= max(ymax,max(N));
                    end
                end
            end
        end
        if isempty(ymax) 
            ymax=0;
            fprintf('\nmaybe this cell never fired a spike ?');
        end
        ylimits=[-.3 ymax];
    end
    fprintf('\nusing ylimits [%.1f %.1f]', ylimits);
    


    if numamps>=numdurs
        %plot freqs x amps, with each dur in a separate window
        numy=numamps;
        ystep=-1;
        ystart=numy;
        yend=1;
        numw=numdurs;
    elseif numdurs>numamps
        %plot freqs x durs, with each amp in a separate window
        numw=numamps;
        numy=numdurs;
        ystart=1;
        ystep=1;
        yend=numy;
    end

    if numdurs==1 & numamps==1
        numy=numfreqs;
        numx=1;
    else
        numx=numfreqs;
    end

    %plot the mean tuning curve OFF
    for windex=1:1                                                              % Normally windex=1:numw, hardcoded here for efficiency
        figure('position',windowpos, 'PaperOrientation', 'landscape')
        p=0;
        subplot1(numy,numx, 'Max', [.95 .9])
        for yindex=ystart:ystep:yend
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                hold on
                if numamps>=numdurs, aindex=yindex; dindex=windex;
                else aindex=windex; dindex=yindex;
                end
                spiketimes1=mM1OFF(cellnum, findex, aindex, dindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(spiketimes1, X);
                N=N./nrepsOFF(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                offset=ylimits(2);
                yl=ylimits;
                inc=(yl(2))/max(max(max(nreps)));
                if rasters==1
                    for n=1:nrepsOFF(cellnum, findex, aindex, dindex)
                        spiketimes2=M1OFF(cellnum, findex, aindex, dindex, n).spiketimes;
                        offset=offset+inc;
                        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    end
                end
                b=bar(x,N,1);
                set(b, 'facecolor', 'k');

                offsetS=ylimits(1)+.05*diff(ylimits);
                if StimRecorded
                    Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
                    Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                    Stimtrace=.25*diff(ylimits)*Stimtrace;
                    t=1:length(Stimtrace);
                    t=1000*t/out.samprate; %convert to ms
                    t=t+out.xlimits(1); %correct for xlim in original processing call
                    stimh=plot(t, Stimtrace+offsetS, 'm');
                    uistack(stimh, 'bottom')
                    ylimits2(1)=0;%2*offsetS;
                else
                    line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                    ylimits2(1)=-2;
                end
                if LaserRecorded
                    for rep=1:nrepsOFF(cellnum, findex, aindex, dindex)
                        Lasertrace=squeeze(M1OFFLaser(findex, aindex, dindex,rep, :));
                        Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                        Lasertrace=.05*diff(ylimits)*Lasertrace;
                        t=1:length(Lasertrace);
                        t=1000*t/out.samprate; %convert to ms
                        t=t+out.xlimits(1); %correct for xlim in original processing call
                        plot( t, Lasertrace+offsetS, 'c')
                    end
                end
                ylimits2(2)=ylimits(2)+offset;
                %ylimits2(2)=2*ylimits(2);
                try
                    ylimits(2) = max([ylimits(2) 1]);
                    ylim(ylimits2)
                end
                xlim(xlimits)
                %xlim([-50 100])
                set(gca, 'fontsize', fs)
                %set(gca, 'xticklabel', '')
                %set(gca, 'yticklabel', '')
                if p==1
                    h=title(sprintf('%s: \ncell%d %d spikes, %dms, nreps: %d-%d, OFF ',datadir,cellnum, length(out.SortedUnits(cellnum).spiketimes), durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
                    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
                end
                xl = xlim; yl = ylim;
                %             h=text(0, range(yl),sprintf('%d ms',durs(dindex)));
                %             set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','position',[0 range(yl)*.8 0])



                %label amps and freqs
                subplot1(p)
                if findex==1
                    if numamps>=numdurs
                        %text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(yindex)))
                        ylabel(sprintf('%ddB',amps(yindex)))
                    else
                        text(xlimits(1)+diff(xlimits)/20, mean(ylimits), int2str(durs(yindex)))
                    end
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=mean(ylimits);
                        %                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=mean(ylimits);
                    end

                    if freqs(findex)==-1
                        xlabel('WN')
                    else
                        xlabel(sprintf('%.1f', freqs(findex)/1000))
                    end

                end

                if ~isempty(mSilentSoundOFF(cellnum).spiketimes) && yindex==1 && findex==1
                    %p=p+1;
                    subplot1(p)
                    spiketimes1=mSilentSoundOFF(cellnum).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nreps_ssOFF(cellnum); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=ylimits(2);
                    yl=ylimits;
                    inc=(yl(2))/max(max(max(nreps)));
                    if rasters==1
                        for n=1:nreps_ssOFF
                            spiketimes2=SilentSoundOFF(cellnum, n).spiketimes;
                            offset=offset+inc;
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    b=bar(x,N,1);
                    set(b, 'facecolor', 'k');
                    vpos=mean(ylimits);
                    text(xlimits(1), vpos, 'SS')
                    ylimits2(2)=ylimits(2)+offset;
                    %ylimits2(2)=2*ylimits(2);
                    ylimits(2) = max([ylimits(2) 1]);
                    ylim(ylimits2)
                    xlim(xlimits)
                    %xlim([-50 100])
                    set(gca, 'fontsize', fs)

                    offsetS=ylimits(1)+.05*diff(ylimits);
                    if StimRecorded
                        Stimtrace=mean(SilentSoundOFFStim,1);
                        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                        Stimtrace=.25*diff(ylimits)*Stimtrace;
                        t=1:length(Stimtrace);
                        t=1000*t/out.samprate; %convert to ms
                        t=t+out.xlimits(1); %correct for xlim in original processing call
                        offset=ylimits(1)+.1*diff(ylimits);
                        stimh=plot(t, Stimtrace+offsetS, 'm');
                        uistack(stimh, 'bottom')
                        ylimits2(1)=0;%2*offsetS;
                    else
                        line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                        ylimits2(1)=-2;
                    end
                    if LaserRecorded
                        for rep=1:nreps_ssOFF
                            Lasertrace=squeeze(SilentSoundOFFLaser(rep, :));
                            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                            Lasertrace=.05*diff(ylimits)*Lasertrace;
                            plot( t, Lasertrace+offsetS, 'c')
                        end
                    end
                end
            end
        end
    end

    if IL
        %plot the mean tuning curve ON
        for windex=1:numw
        figure('position',windowpos, 'PaperOrientation', 'landscape')
                
            p=0;
            subplot1(numy,numx, 'Max', [.95 .9])
            for yindex=ystart:ystep:yend
                for findex=1:numfreqs
                    p=p+1;
                    subplot1(p)
                    hold on
                    if numamps>=numdurs, aindex=yindex; dindex=windex;
                    else aindex=windex; dindex=yindex;
                    end
                    spiketimes1=mM1ON(cellnum, findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nrepsON(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=ylimits(2);
                    yl=ylimits;
                    inc=(yl(2))/max(nrepsON(:));
                    if rasters==1
                        for n=1:nrepsON(cellnum, findex, aindex, dindex)
                            spiketimes2=M1ON(cellnum, findex, aindex, dindex, n).spiketimes;
                            offset=offset+inc;
                            %this should plot a cyan line for every trial among the rasters
                            %it should accomodate trial-by-trial changes to
                            %Laser params
                            MLaserStart=M_LaserStart(findex,aindex,dindex, n);
                            MLaserWidth=M_LaserWidth(findex,aindex,dindex, n);
                            MLaserNumPulses=M_LaserNumPulses(findex,aindex,dindex, n);
                            MLaserISI=M_LaserISI(findex,aindex,dindex, n);
                            for np=1:MLaserNumPulses
                                plot([MLaserStart+(np-1)*(MLaserWidth+MLaserISI) MLaserStart+(np-1)*(MLaserWidth+MLaserISI)+MLaserWidth], [1 1]+yl(2)+offset, 'c')
                            end
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    b=bar(x, N,1);
                    set(b, 'facecolor', 'k');
                    line(xlimits, [0 0], 'color', 'k')
                    ylimits2(2)=ylimits(2)+offset;
                    if StimRecorded
                        Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
                        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                        Stimtrace=.25*diff(ylimits)*Stimtrace;
                        t=1:length(Stimtrace);
                        t=1000*t/out.samprate; %convert to ms
                        t=t+out.xlimits(1); %correct for xlim in original processing call
                        offset=ylimits(1)+.1*diff(ylimits);
                        stimh=plot(t, Stimtrace+offsetS, 'm');
                        uistack(stimh, 'bottom')
                        ylimits2(1)=0;%2*offsetS;
                    else
                        line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                        ylimits2(1)=-2;
                    end


                    %ylimits2(2)=2*ylimits(2);
                    try
                        ylimits(2) = max([ylimits(2) 1]);
                        ylim(ylimits2)
                    end

                    if LaserRecorded
                        for rep=1:nrepsON(findex, aindex, dindex)
                            Lasertrace=squeeze(M1ONLaser(findex, aindex, dindex,rep, :));
                            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                            Lasertrace=.05*diff(ylimits)*Lasertrace;
                            plot( t, Lasertrace+offset, 'c')
                        end
                    else

                        X=[xlimits(1), LaserStart, LaserStart, LaserStart+LaserWidth, LaserStart+LaserWidth, xlimits(2)];
                        height=diff(ylimits)/5;
                        Y=[0 0 height height 0 0];
                        line(X,Y, 'color', 'c', 'linewidth', 2)
                    end

                    xlim(xlimits)
                    %xlim([-50 100])
                    set(gca, 'fontsize', fs)
                    % set(gca, 'xticklabel', '')
                    % set(gca, 'yticklabel', '')

                    subplot1(1)
                    h=title(sprintf('%s\ncell%d %dms, nreps: %d-%d, ON',datadir,cellnum,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
                    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

                    %label amps and freqs
                    subplot1(p)
                    if findex==1
                        if numamps>=numdurs
                            %text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(yindex)))
                            ylabel(sprintf('%ddB',amps(yindex)))
                        else
                            text(xlimits(1)+diff(xlimits)/20, mean(ylimits), int2str(durs(yindex)))
                        end
                    end
                    if aindex==1
                        if mod(findex,2) %odd freq
                            vpos=mean(ylimits);
                            %                        vpos=ylimits(1)-mean(ylimits);
                        else
                            vpos=mean(ylimits);
                        end
                        %text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                        if freqs(findex)==-1
                            xlabel('WN')
                        else
                            xlabel(sprintf('%.1f', freqs(findex)/1000))
                        end
                    end

                    %plot response to silent sound laser pulse
                    if ~isempty(mSilentSoundON(cellnum).spiketimes) && yindex==1 && findex==1
                        subplot1(p)
                        hold on
                        spiketimes1=mSilentSoundON(cellnum).spiketimes;
                        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                        [N, x]=hist(spiketimes1, X);
                        N=N./nreps_ssON(cellnum); %normalize to spike rate (averaged across trials)
                        N=1000*N./binwidth; %normalize to spike rate in Hz
                        offset=max(N);
                        yl=ylimits;
                        inc=(offset)/nreps_ssON(cellnum);
                        if rasters==1
                            for n=1:nreps_ssON(cellnum)
                                spiketimes2=SilentSoundON(cellnum, n).spiketimes;
                                offset=offset+inc;
                                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                            end
                        end
                        b=bar(x,N,1);
                        set(b, 'facecolor', 'k');
                        vpos=mean(ylimits);
                        text(xlimits(1), vpos, 'SS')
                        ylimits2(2)=ylimits(2)+2+ offset;
                        %ylimits2(2)=2*ylimits(2);
                        try
                            %ylimits(2) = max([ylimits(2) 1]);
                            ylim(ylimits2)
                        end
                        xlim(xlimits)
                        set(gca, 'fontsize', fs)
                        offsetS=ylimits(1)+.05*diff(ylimits);
                        if StimRecorded
                            Stimtrace=squeeze(mean(SilentSoundONStim,1));
                            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                            Stimtrace=.25*diff(ylimits)*Stimtrace;
                            t=1:length(Stimtrace);
                            t=1000*t/out.samprate; %convert to ms
                            t=t+out.xlimits(1); %correct for xlim in original processing call
                            offset=ylimits(1)+.1*diff(ylimits);
                            stimh=plot(t, Stimtrace+offsetS, 'm');
                            uistack(stimh, 'bottom')
                            ylimits2(1)=0;%2*offsetS;
                        else
                            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                            ylimits2(1)=-2;
                        end
                        if LaserRecorded
                            for rep=1:nreps_ssON(cellnum)
                                Lasertrace=squeeze(SilentSoundONLaser(rep,:));
                                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                                Lasertrace=.05*diff(ylimits)*Lasertrace;
                                plot( t, Lasertrace+offsetS, 'c')
                            end
                        end
                    end
                end
            end
        end
    end

      if printtofile
            %print figures to postscript file
            f=findobj('type', 'figure');
            for idx=1:length(f)
                %figure(f(idx))
                % orient landscape
                % % print figs -dpsc2 -append -bestfit
                exportgraphics(f(idx),'figs.pdf','Append',true)

                if closewindows
                    close
                end
            end
      end

end