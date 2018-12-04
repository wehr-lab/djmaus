function PlotGPIAS_Behavior_kip(datadir)

%plots GPIAS behavioral data from djmaus and open-ephys accelerometers
%
%usage: PlotGPIAS_Behavior(datadir)
%Processes data if outfile is not found;
% adds on GTR and burst data from all cells with outfiles in pwd
% _kip version adds on GTR and burst data from all cells with outfiles in pwd
% should be cleaned up

if nargin==0 datadir=pwd;end

flag.plot = 0;  % plot all startles for each gapdur?

PreStartleWindowms=[-100 0]; % in ms relative to onset of startle pulse
PostStartleWindowms=[0 100]; % in ms relative to onset of startle-pulse
ISIWindowms=[0 60]; % in ms relative to onset of pre-pulse    %added by APW 3_31_14


[~, FNtemp, ~] = fileparts(datadir);
force_reprocess=1;
if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessGPIAS_BehaviorTilt(datadir,4);
end

cd(datadir)
outfilename=sprintf('outGPIAS_Behavior.mat');
d=dir(outfilename);
if ~isempty(d)
    load(outfilename)
else
    ProcessGPIAS_Behavior(datadir)
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
isi=out.isi;
soaflag=out.soaflag;
PeakON=out.PeakON;
PeakOFF=out.PeakOFF;
mPeakON=out.mPeakON;
mPeakOFF=out.mPeakOFF;
semPeakON=out.semPeakON;
semPeakOFF=out.semPeakOFF;
percentGPIAS_ON=out.percentGPIAS_ON;
percentGPIAS_OFF=out.percentGPIAS_OFF;
pON = out.pON;
pOFF = out.pOFF;
samprate=out.samprate;
M1ONstim=out.M1ONstim;
M1OFFstim=out.M1OFFstim;
mM1ONstim=out.mM1ONstim;
mM1OFFstim=out.mM1OFFstim;
xlimits=out.xlimits;

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
            subplot1(numgapdurs,1)
            for gdindex=1:numgapdurs
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
                ylabel(sprintf('%d ms',gapdurs(gdindex)));
                text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nrepsOFF(gdindex, paindex)))
                %axis off
                
            end
            subplot1(1)
            h=title(sprintf('%s:\nnreps: %d-%d, OFF',FNtemp,min(nrepsOFF(:)),max(nrepsOFF(:))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','interpreter','none')
            
            subplot1(numgapdurs)
            xlabel('Time (ms)');
        end
        
        %plot the mean Laser ON tuning curve
        if ~isempty (mM1ON)
            figure;
            p=0;
            subplot1(numgapdurs,1)
            for gdindex=1:numgapdurs
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
                    line([laserstart-gapdelay laserstart-gapdelay+laserwidth],-2.5*offset*[1 1],'color','c','linewidth',5);   % AW commented out to enable no laser gap detection plots
                end
                
                xlim(xlimits)
                ylabel(sprintf('%d ms',gapdurs(gdindex)));
                text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nrepsON(gdindex, paindex)))
                %axis off
                
            end
            subplot1(1)
            h=title(sprintf('%s:\nnreps: %d-%d, ON',FNtemp,min(nrepsON(:)),max(nrepsON(:))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            subplot1(numgapdurs)
            xlabel('Time (ms)');
            
            pos=get(gcf, 'pos');
            pos(1)=pos(1)+pos(3); %shift ON to right
            set(gcf, 'pos', pos)
        end
    end
    
end

%plot the percent GPIAS
%note: the percentGPIAS assumes only a single pulseamp, so if we ever do
%more than one pulse amp we need to update that.
Hfig_percentGPIAS = figure;
hold on
gd=1:numgapdurs;
if ~isempty(percentGPIAS_OFF)
    plot(gd, percentGPIAS_OFF, 'k-o')
end
if ~isempty(percentGPIAS_ON)
    plot(gd, percentGPIAS_ON, 'c-o')
end
set(gca, 'xtick', 1:numgapdurs)
set(gca, 'xticklabel', gapdurs)
h = title ([FNtemp ' percent GPIAS']);
set(h,'interpreter','none')
xlabel('gap duration')
ylabel('percent GPIAS')


%plot the mean peak rectified startle
for paindex=1:numpulseamps
    Hfig_meanStartle(paindex) = figure;
    hold on
    gd=1:numgapdurs;
    if ~isempty(mPeakOFF)
        errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), 'k-o')
    end
    if ~isempty(mPeakON)
        errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), 'c-o')
    end
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    h = title ([FNtemp ' startle']);
    set(h,'interpreter','none')
    xlabel('gap duration')
    ylabel('startle response +- sem')
end

%plot all trials peak rectified startle
if flag.plot
for paindex=1:numpulseamps
    Hfig_allStartle(paindex) = figure;
    hold on
    gd=1:numgapdurs;
    if ~isempty(mPeakOFF)
        errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), 'k-o')
        plot(gd, squeeze(PeakOFF(:,paindex,:)), 'ko')
    end
    if ~isempty(mPeakON)
        errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), 'c-o')
        plot(gd, squeeze(PeakON(:,paindex,:)), 'co')
    end
    
    for igap = 2:numgapdurs
        str = {['GPIAS ' num2str(round(percentGPIAS_OFF(igap),1)) ' %']; ['p ~ ' num2str(round(pOFF(igap),3))]};
        text(igap-.25,mPeakOFF(igap)*.95,str)
    end
    
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    h = title ([FNtemp ' startle']);
    set(h,'interpreter','none')
    xlabel('gap duration')
    ylabel('startle responses all trials')
end
end

% now process GTR and burst responses for all cells in pwd
tempDir = dir([pwd '\out*']);

if length(tempDir)>1
% Hfig_percentGPIAS2 = figure; hold on
% xlabel('percent GPIAS')
% ylabel('percent burst or GTR')
X1 = [];
Y1 = [];
Y2 = [];

for icell = 2:length(tempDir)
    analysisWindow = 50;
    load(tempDir(icell).name)
    
    nrepsOFF=min1(out.nrepsOFF);
    M1OFF = out.M1OFF;
    gapdurs = out.gapdurs;
    [minGap, ind_minGap] =min(gapdurs);
    numgapdurs = length(gapdurs);
    burstDelay = out.soa;
    
    % make arrays (last value is SPONT)
    PGIspikesPerRep = nan(numgapdurs+1,nrepsOFF);
%     H_PGI = nan(1,numgapdurs);
%     P_PGI = nan(1,numgapdurs);
    deltaGTR = nan(1,numgapdurs);
    BURSTspikesPerRep = nan(numgapdurs+1,nrepsOFF);
%     H_burst = nan(1,numgapdurs);
%     P_burst = nan(1,numgapdurs);
    deltaBurst = nan(1,numgapdurs);
    
    if minGap           % if there is no 0 gap, then use data before minGap
        for irep = 1:nrepsOFF
            temp = M1OFF(ind_minGap,1,irep).spiketimes;
            PGIspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap+analysisWindow) & temp<=-minGap));
            BURSTspikesPerRep(numgapdurs+1,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
        end
    else        % minGap==0
        for irep = 1:nrepsOFF
            temp = M1OFF(ind_minGap,1,irep).spiketimes;
            PGIspikesPerRep(numgapdurs+1,irep) = length(find(temp>0 & temp<=analysisWindow));
            BURSTspikesPerRep(numgapdurs+1,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
        end
    end
    
    for idur = 1:numgapdurs
        for irep = 1:nrepsOFF
            temp = M1OFF(idur,1,irep).spiketimes;
            PGIspikesPerRep(idur,irep) = length(find(temp>0 & temp<=analysisWindow));
            BURSTspikesPerRep(idur,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
        end
         % hypothesis and P-values
        %[H_PGI(idur),P_PGI(idur)] = ttest(PGIspikesPerRep(idur,:), PGIspikesPerRep(numgapdurs+1,:),'tail','both');
        %[H_burst(idur),P_burst(idur)] = ttest(BURSTspikesPerRep(idur,:), BURSTspikesPerRep(numgapdurs+1,:),'tail','both');
   end
   
   meanGTR = mean(PGIspikesPerRep,2);
   meanBurst = mean(BURSTspikesPerRep,2);
   figure(Hfig_allStartle)
   for idur = 1:numgapdurs
       deltaGTR(idur) = mean(PGIspikesPerRep(idur,:)) / max1(meanGTR);
       deltaBurst(idur) = mean(BURSTspikesPerRep(idur,:)) / max1(meanBurst);
       plot(idur*ones(1,nrepsOFF),PGIspikesPerRep(idur,:)/ max1(meanGTR),'ro')
   end
    %ylim([0 10]);
    
    figure(Hfig_percentGPIAS)
    plot(2:numgapdurs, deltaGTR(2:numgapdurs)*100,'ro-')
    plot(2:numgapdurs, deltaBurst(2:numgapdurs)*100,'bo-')
    
    figure(Hfig_meanStartle)
    plot(2:numgapdurs, mean(PGIspikesPerRep(2:numgapdurs,:),2),'ro-')
    plot(2:numgapdurs, mean(BURSTspikesPerRep(2:numgapdurs,:),2),'bo-')
    
    figure; hold on
    plot(percentGPIAS_OFF, deltaGTR*100,'ro')
    plot(percentGPIAS_OFF, deltaBurst*100,'bo')
    xlabel('percent GPIAS')
    ylabel('percent burst or GTR')
    [B1,~,~,~,STATS1] = regress(deltaGTR' *100,[ones(size(percentGPIAS_OFF')) percentGPIAS_OFF']);
    Hline(1) = plot(percentGPIAS_OFF,B1(1)+B1(2).*percentGPIAS_OFF,'r');
    [B2,~,~,~,STATS2] = regress(deltaBurst' *100,[ones(size(percentGPIAS_OFF')) percentGPIAS_OFF']);
    Hline(2) = plot(percentGPIAS_OFF,B2(1)+B2(2).*percentGPIAS_OFF,'b');
    legend(Hline,['GTR   %var=' num2str(round(STATS1(1)*100)) '  p=' num2str(round(STATS1(3),2))], ...
        ['burst   %var=' num2str(round(STATS2(1)*100))  '  p=' num2str(round(STATS2(3),2))]);
    title(tempDir(icell).name)
    
    X1 = [X1 percentGPIAS_OFF];
    Y1 = [Y1 deltaGTR*100];
    Y2 = [Y2 deltaBurst*100];
    
end

figure; hold on
plot(X1,Y1,'ro')
plot(X1,Y2,'bo')
xlabel('percent GPIAS')
ylabel('percent burst or GTR')
[B1,BINT1,R1,RINT1,STATS1] = regress(Y1',[ones(size(X1')) X1']);
Hline(1) = plot(X1,B1(1)+B1(2).*X1,'r');
[B2,BINT2,R2,RINT2,STATS2] = regress(Y2',[ones(size(X1')) X1']);
Hline(2) = plot(X1,B2(1)+B2(2).*X1,'b');
legend(Hline,['GTR   %var=' num2str(round(STATS1(1)*100)) '  p=' num2str(round(STATS1(3),2))], ...
    ['burst   %var=' num2str(round(STATS2(1)*100))  '  p=' num2str(round(STATS2(3),2))]);
end