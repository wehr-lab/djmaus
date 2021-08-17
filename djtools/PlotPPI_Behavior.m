function PlotPPI_Behavior(datadir)

%plots PPI behavioral data from djmaus using the tilt platform (pressure
%transducer)
%
%usage: PlotPPI_Behavior(datadir)
%Processes data if outfile is not found;
% adds on GTR and burst data from all cells with outfiles in pwd

% modified from PlotGPIAS_Behavior_kip mw 8-16-21

if nargin==0 datadir=pwd;end

flag.plot = 0;  % plot all startles for each prepulsedur?

PreStartleWindowms=[-100 0]; % in ms relative to onset of startle pulse
PostStartleWindowms=[0 100]; % in ms relative to onset of startle-pulse
ISIWindowms=[0 60]; % in ms relative to onset of pre-pulse    %added by APW 3_31_14


[~, FNtemp, ~] = fileparts(datadir);
FNtemp=FNtemp(1:19); %trim mouseID from dirname

force_reprocess=0;
if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessPPI_BehaviorTilt(datadir,4);
end

cd(datadir)
% d=dir('outPPI_Behavior_combined.mat');
d=dir('outPPI_Behavior.mat');
if isempty(d)
    ProcessPPI_BehaviorTilt(datadir,4);
end

% d=dir('outPPI_Behavior_combined.mat');
d=dir('outPPI_Behavior*.mat');
fprintf('\nfound %d outfiles.', length(d))

for outindex=1:length(d)
    % outfilename=sprintf('outPPI_Behavior.mat');
    outfilename=d(outindex).name;
    fprintf('\nloading %s', outfilename)
    load(outfilename)
    
    
    
    IL=out.IL; %whether there are any interleaved laser trials
    numpulseamps=out.numpulseamps;
    numprepulsedurs=out.numprepulsedurs;
    pulseamps=out.pulseamps;
    prepulsedurs=out.prepulsedurs;
    samprate=out.samprate; %in Hz
    mM1ON=out.mM1ON;
    mM1OFF=out.mM1OFF;
    M1ON=out.M1ON;
    M1OFF=out.M1OFF;
    nrepsON=out.nrepsON;
    nrepsOFF=out.nrepsOFF;
    soa=out.soa;
    isi=out.isi;
    soaflag=out.soaflag;
    PeakON=out.PeakON;
    PeakOFF=out.PeakOFF;
    mPeakON=out.mPeakON;
    mPeakOFF=out.mPeakOFF;
    semPeakON=out.semPeakON;
    semPeakOFF=out.semPeakOFF;
    percentPPI_ON=out.percentPPI_ON;
    percentPPI_OFF=out.percentPPI_OFF;
    pON = out.pON;
    pOFF = out.pOFF;
    samprate=out.samprate;
    M1ONstim=out.M1ONstim;
    M1OFFstim=out.M1OFFstim;
    mM1ONstim=out.mM1ONstim;
    mM1OFFstim=out.mM1OFFstim;
    xlimits=out.xlimits;

    if isfield(out, 'mouseID') & strcmp(out.outfilename, 'outPPI_Behavior.mat')
        mouseID=out.mouseID;
    elseif isfield(out, 'mouse2ID') & strcmp(out.outfilename, 'outPPI_BehaviorMouse2.mat')  
        mouseID=out.mouse2ID;
    else mouseID='???';
    end
    % %find optimal axis limits
    if ~isempty(mM1OFF)
        ylimits(1)=min(mM1OFF(:));
        ylimits(2)=max(mM1OFF(:));
    elseif ~isempty(mM1ON)
        ylimits(1)=min(mM1ON(:));
        ylimits(2)=max(mM1ON(:));
    end
    
    fs=10;
    
    %plot the tuning curves
    
    % Plot the actual trace with mean trace overlayed
    % Separated into 2 figures with laser OFF/ON
    if true
        if ~isempty (mM1OFF)
            
            %plot the mean Laser OFF tuning curve
            for paindex=1:numpulseamps
                figure;
                p=0;
                subplot1(numprepulsedurs,1)
                for gdindex=1:numprepulsedurs
                    p=p+1;
                    subplot1(p)
                    hold on
                    
                    offset=10*std(M1OFF(:));
                    
                    % add the stimulus in magenta
                    stimtrace=squeeze(mM1OFFstim(gdindex, paindex,:));
                    stimtrace=stimtrace-median(stimtrace(1:1000));
                    t=1:length(stimtrace);
                    t=1e3*t/samprate;
                    t=t+xlimits(1);
                    plot(t, stimtrace-offset, 'm')
                    
                    % plot each trial in blue
                    for i=1:nrepsOFF(gdindex, paindex)
                        trace1=squeeze(M1OFF(gdindex, paindex,i,:));
                        trace1=trace1-median(trace1(1:1000));
                        plot(t, trace1+i*offset, 'b');
                    end
                    
                    % plot the mean trace in red
                    trace1=squeeze(mM1OFF(gdindex, paindex,:));
                    trace1=trace1-median(trace1(1:1000));
                    
                    plot(t, trace1, 'r')
                    
                    %ylim([ylimits(1)-3*offset ylimits(2)])
                    xlim(xlimits)
                    ylabel(sprintf('%d ms',prepulsedurs(gdindex)));
                    text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nrepsOFF(gdindex, paindex)))
                    %axis off
                    
                end
                
                subplot1(1)
                h=title(sprintf('%s mouse %s:\nnreps: %d-%d, OFF',FNtemp,mouseID, min(nrepsOFF(:)),max(nrepsOFF(:))));
                set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','interpreter','none')
                
                subplot1(numprepulsedurs)
                xlabel('Time (ms)');
                
                pos=get(gcf, 'pos');
                pos(4)=pos(4)+160; %make taller 
                pos(2)=pos(4)-160; %make taller 
                set(gcf, 'pos', pos)

            end
            
            %plot the mean Laser ON tuning curve
            if ~isempty (mM1ON)
                figure;
                p=0;
                subplot1(numprepulsedurs,1)
                for gdindex=1:numprepulsedurs
                    p=p+1;
                    subplot1(p)
                    hold on
                    
                    offset=10*std(M1ON(:));
                    
                    % add the stimulus in magenta
                    stimtrace=squeeze(mM1ONstim(gdindex, paindex,:));
                    stimtrace=stimtrace-median(stimtrace(1:1000));
                    t=1:length(stimtrace);
                    t=1e3*t/samprate;
                    t=t+xlimits(1);
                    plot(t, stimtrace-offset, 'm')
                    
                    % plot each trial in blue
                    for i=1:nrepsON(gdindex, paindex)
                        trace1=squeeze(M1ON(gdindex, paindex,i,:));
                        trace1=trace1-median(trace1(1:1000));
                        plot(t, trace1+i*offset, 'b');
                    end
                    
                    
                    % plot the mean trace in red
                    trace1=squeeze(mM1ON(gdindex, paindex,:));
                    trace1=trace1-median(trace1(1:1000));
                    plot(t, trace1, 'r')
                    
                    % add the laser in cyan
                    try
                        line([laserstart laserstart+laserwidth],-2.5*offset*[1 1],'color','c','linewidth',5);   % AW commented out to enable no laser prepulse detection plots
                    end
                    
                    xlim(xlimits)
                    ylabel(sprintf('%d ms',prepulsedurs(gdindex)));
                    text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nrepsON(gdindex, paindex)))
                    %axis off
                    
                end
                subplot1(1)
                h=title(sprintf('%s:\nnreps: %d-%d, ON',FNtemp,min(nrepsON(:)),max(nrepsON(:))));
                set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
                subplot1(numprepulsedurs)
                xlabel('Time (ms)');
                
                pos=get(gcf, 'pos');
                pos(1)=pos(1)+pos(3); %shift ON to right
                set(gcf, 'pos', pos)
            end
        end
        
    end
    
    %plot the percent PPI
    %note: the percentPPI assumes only a single pulseamp, so if we ever do
    %more than one pulse amp we need to update that.
    Hfig_percentPPI = figure;
    hold on
    gd=1:numprepulsedurs;
    if ~isempty(percentPPI_OFF)
        plot(gd, percentPPI_OFF, 'k-o')
    end
    if ~isempty(percentPPI_ON)
        plot(gd, percentPPI_ON, 'c-o')
    end
    set(gca, 'xtick', 1:numprepulsedurs)
    set(gca, 'xticklabel', prepulsedurs)
    h = title ([FNtemp, ' mouse ', mouseID,' percent PPI']);
    set(h,'interpreter','none')
    xlabel('prepulse duration')
    ylabel('percent PPI')
    if isfield(out, 'generated_by')
        if strcmp(out.generated_by, 'Outfile_Combiner')
            text(.1, .9, 'combined outfile', 'units', 'normal')
        end
    end
    
    %plot the mean peak rectified startle
    for paindex=1:numpulseamps
        Hfig_meanStartle(paindex) = figure;
        hold on
        gd=1:numprepulsedurs;
        if ~isempty(mPeakOFF)
            errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), 'k-o')
        end
        if ~isempty(mPeakON)
            errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), 'c-o')
        end
        set(gca, 'xtick', 1:numprepulsedurs)
        set(gca, 'xticklabel', prepulsedurs)
        h = title ([FNtemp, ' mouse ', mouseID, ' startle']);
        set(h,'interpreter','none')
        xlabel('prepulse duration')
        ylabel('startle response +- sem')
        if isfield(out, 'generated_by')
            if strcmp(out.generated_by, 'Outfile_Combiner')
                text(.1, .9, 'combined outfile', 'units', 'normal')
            end
        end
    end
    
    %plot all trials peak rectified startle
    if flag.plot
        for paindex=1:numpulseamps
            Hfig_allStartle(paindex) = figure;
            hold on
            gd=1:numprepulsedurs;
            if ~isempty(mPeakOFF)
                errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), 'k-o')
                plot(gd, squeeze(PeakOFF(:,paindex,:)), 'ko')
            end
            if ~isempty(mPeakON)
                errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), 'c-o')
                plot(gd, squeeze(PeakON(:,paindex,:)), 'co')
            end
            
            for iprepulse = 2:numprepulsedurs
                str = {['PPI ' num2str(round(percentPPI_OFF(iprepulse),1)) ' %']; ['p ~ ' num2str(round(pOFF(iprepulse),3))]};
                text(iprepulse-.25,mPeakOFF(iprepulse)*.95,str)
            end
            
            set(gca, 'xtick', 1:numprepulsedurs)
            set(gca, 'xticklabel', prepulsedurs)
            h = title ([FNtemp, ' mouse ', mouseID, ' startle']);
            set(h,'interpreter','none')
            xlabel('prepulse duration')
            ylabel('startle responses all trials')
        end
    end
    
end
return
%not sure why it's looking for cellular data here

% now process GTR and burst responses for all cells in pwd
tempDir = dir([pwd '\out*']);

if length(tempDir)>1
    % Hfig_percentPPI2 = figure; hold on
    % xlabel('percent PPI')
    % ylabel('percent burst or GTR')
    X1 = [];
    Y1 = [];
    Y2 = [];
    
    for icell = 2:length(tempDir)
        analysisWindow = 50;
        load(tempDir(icell).name)
        
        nrepsOFF=min(out.nrepsOFF);
        M1OFF = out.M1OFF;
        prepulsedurs = out.prepulsedurs;
        [minprepulse, ind_minprepulse] =min(prepulsedurs);
        numprepulsedurs = length(prepulsedurs);
        burstDelay = out.soa;
        
        % make arrays (last value is SPONT)
        PGIspikesPerRep = nan(numprepulsedurs+1,nrepsOFF);
        %     H_PGI = nan(1,numprepulsedurs);
        %     P_PGI = nan(1,numprepulsedurs);
        deltaGTR = nan(1,numprepulsedurs);
        BURSTspikesPerRep = nan(numprepulsedurs+1,nrepsOFF);
        %     H_burst = nan(1,numprepulsedurs);
        %     P_burst = nan(1,numprepulsedurs);
        deltaBurst = nan(1,numprepulsedurs);
        
        if minprepulse           % if there is no "pure startle" condition, then use data before minprepulse
            for irep = 1:nrepsOFF
                temp = M1OFF(ind_minprepulse,1,irep).spiketimes;
                PGIspikesPerRep(numprepulsedurs+1,irep) = length(find(temp>-(minprepulse+analysisWindow) & temp<=-minprepulse));
                BURSTspikesPerRep(numprepulsedurs+1,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
            end
        else        % minprepulse==0
            for irep = 1:nrepsOFF
                temp = M1OFF(ind_minprepulse,1,irep).spiketimes;
                PGIspikesPerRep(numprepulsedurs+1,irep) = length(find(temp>0 & temp<=analysisWindow));
                BURSTspikesPerRep(numprepulsedurs+1,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
            end
        end
        
        for idur = 1:numprepulsedurs
            for irep = 1:nrepsOFF
                temp = M1OFF(idur,1,irep).spiketimes;
                PGIspikesPerRep(idur,irep) = length(find(temp>0 & temp<=analysisWindow));
                BURSTspikesPerRep(idur,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
            end
            % hypothesis and P-values
            %[H_PGI(idur),P_PGI(idur)] = ttest(PGIspikesPerRep(idur,:), PGIspikesPerRep(numprepulsedurs+1,:),'tail','both');
            %[H_burst(idur),P_burst(idur)] = ttest(BURSTspikesPerRep(idur,:), BURSTspikesPerRep(numprepulsedurs+1,:),'tail','both');
        end
        
        meanGTR = mean(PGIspikesPerRep,2);
        meanBurst = mean(BURSTspikesPerRep,2);
        figure(Hfig_allStartle)
        for idur = 1:numprepulsedurs
            deltaGTR(idur) = mean(PGIspikesPerRep(idur,:)) / max1(meanGTR);
            deltaBurst(idur) = mean(BURSTspikesPerRep(idur,:)) / max1(meanBurst);
            plot(idur*ones(1,nrepsOFF),PGIspikesPerRep(idur,:)/ max1(meanGTR),'ro')
        end
        %ylim([0 10]);
        
        figure(Hfig_percentPPI)
        plot(2:numprepulsedurs, deltaGTR(2:numprepulsedurs)*100,'ro-')
        plot(2:numprepulsedurs, deltaBurst(2:numprepulsedurs)*100,'bo-')
        
        figure(Hfig_meanStartle)
        plot(2:numprepulsedurs, mean(PGIspikesPerRep(2:numprepulsedurs,:),2),'ro-')
        plot(2:numprepulsedurs, mean(BURSTspikesPerRep(2:numprepulsedurs,:),2),'bo-')
        
        figure; hold on
        plot(percentPPI_OFF, deltaGTR*100,'ro')
        plot(percentPPI_OFF, deltaBurst*100,'bo')
        xlabel('percent PPI')
        ylabel('percent burst or GTR')
        [B1,~,~,~,STATS1] = regress(deltaGTR' *100,[ones(size(percentPPI_OFF')) percentPPI_OFF']);
        Hline(1) = plot(percentPPI_OFF,B1(1)+B1(2).*percentPPI_OFF,'r');
        [B2,~,~,~,STATS2] = regress(deltaBurst' *100,[ones(size(percentPPI_OFF')) percentPPI_OFF']);
        Hline(2) = plot(percentPPI_OFF,B2(1)+B2(2).*percentPPI_OFF,'b');
        legend(Hline,['GTR   %var=' num2str(round(STATS1(1)*100)) '  p=' num2str(round(STATS1(3),2))], ...
            ['burst   %var=' num2str(round(STATS2(1)*100))  '  p=' num2str(round(STATS2(3),2))]);
        title(tempDir(icell).name)
        
        X1 = [X1 percentPPI_OFF];
        Y1 = [Y1 deltaGTR*100];
        Y2 = [Y2 deltaBurst*100];
        
    end
    
    figure; hold on
    plot(X1,Y1,'ro')
    plot(X1,Y2,'bo')
    xlabel('percent PPI')
    ylabel('percent burst or GTR')
    [B1,BINT1,R1,RINT1,STATS1] = regress(Y1',[ones(size(X1')) X1']);
    Hline(1) = plot(X1,B1(1)+B1(2).*X1,'r');
    [B2,BINT2,R2,RINT2,STATS2] = regress(Y2',[ones(size(X1')) X1']);
    Hline(2) = plot(X1,B2(1)+B2(2).*X1,'b');
    legend(Hline,['GTR   %var=' num2str(round(STATS1(1)*100)) '  p=' num2str(round(STATS1(3),2))], ...
        ['burst   %var=' num2str(round(STATS2(1)*100))  '  p=' num2str(round(STATS2(3),2))]);
end
