function PlotGPIAS_Behavior_HeadFixed(varargin)
dbstop if error
%plots GPIAS behavioral data from djmaus and pressure sensor device using
%one of ADC inputs
%usage: PlotGPIAS_Behavior(datadir1, datadir2)
% plots one plot if diven one dir, plots pre and post plots if given two
% dirs


if nargin==0 datadir{1}=pwd;end
if nargin==1
    datadir{1}=varargin{1};
elseif nargin==2
    datadir{1}=varargin{1};
    datadir{2}=varargin{2};
end

num=length(datadir);

colorOFF=['b','g'];
colorON=['c','r'];
PreStartleWindowms=[-100 0]; % in ms relative to onset of startle pulse
PostStartleWindowms=[0 100]; % in ms relative to onset of startle-pulse
% ISIWindowms=[0 60]; % in ms relative to onset of pre-pulse    %added by APW 3_31_14

for N=1:num
    force_reprocess=0;
    if force_reprocess
        fprintf('\nForce re-process\n')
        ProcessGPIAS_Behavior_HeadFixed(datadir{N})
    end
    
    cd(datadir{N})
    outfilename=sprintf('outGPIAS_Behavior.mat');
    d=dir(outfilename);
    if ~isempty(d)
        load(outfilename)
    else
        ProcessGPIAS_Behavior_HeadFixed(datadir{N})
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
    samprate=out.samprate;
    M1ONstim=out.M1ONstim;
    M1OFFstim=out.M1OFFstim;
    mM1ONstim=out.mM1ONstim;
    mM1OFFstim=out.mM1OFFstim;
    xlimits=out.xlimits;
    
    % %find optimal axis limits
    if ~isempty(mM1OFF)
        ylimits(1,:,N)=min(mM1OFF(:));
        ylimits(:,2,N)=max(mM1OFF(:));
    elseif ~isempty(mM1ON)
        ylimits(1,:,N)=min(mM1ON(:));
        ylimits(:,2,N)=max(mM1ON(:));
    end
    
    fs=10;
    
    %plot the tuning curves
    
    % Plot the actual trace with mean trace overlayed
    % Separated into 2 figures with laser OFF/ON
    if true
        if ~isempty (mM1OFF)
            
            %plot the mean Laser OFF tuning curve
            for paindex=1:numpulseamps
                figure(N);
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
                    
                    plot(t, trace1*10, 'r')
                    
                    %ylim([ylimits(1)-3*offset ylimits(2)])
                    xlim(xlimits)
                    ylabel(sprintf('%d ms',gapdurs(gdindex)));
                    text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nrepsOFF(gdindex, paindex)))
                    %axis off
                    
                end
                subplot1(1)
                h=title(sprintf('%s:\nnreps: %d-%d, OFF',datadir{N},min(nrepsOFF(:)),max(nrepsOFF(:))));
                set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
                
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
                h=title(sprintf('%s:\nnreps: %d-%d, ON',datadir{N},min(nrepsON(:)),max(nrepsON(:))));
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
    figure(10);
    hold on
    gd=1:numgapdurs;
    if ~isempty(percentGPIAS_OFF)
        plot(gd, percentGPIAS_OFF, sprintf('%s-o', colorOFF(N)))
    end
    if ~isempty(percentGPIAS_ON)
        plot(gd, percentGPIAS_ON, sprintf('%s-o', colorON(N)))
    end
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    title ( 'percent GPIAS')
    xlabel('gap duration')
    ylabel('percent GPIAS')
    if N==2
        legend('pre','post')
    end
    
    
    %plot the mean peak rectified startle
    for paindex=1:numpulseamps
        figure(11);
        hold on
        gd=1:numgapdurs;
        if ~isempty(mPeakOFF)
            errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), sprintf('%s-o', colorOFF(N)))
        end
        if ~isempty(mPeakON)
            errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), sprintf('%s-o', colorON(N)))
        end
        set(gca, 'xtick', 1:numgapdurs)
        set(gca, 'xticklabel', gapdurs)
        title ('LaserON/OFF startle')
        xlabel('gap duration')
        ylabel('startle response +- sem')
        if N==2
            legend('pre','post')
        end
    end
    
    
    %plot all trials peak rectified startle
    for paindex=1:numpulseamps
        figure(12);
        hold on
        gd=1:numgapdurs;
        if ~isempty(mPeakOFF)
            errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), sprintf('%s-o', colorOFF(N)))
            plot(gd, squeeze(PeakOFF(:,paindex,:)), sprintf('%so', colorOFF(N)))
        end
        if ~isempty(mPeakON)
            errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), sprintf('%s-o', colorON(N)))
            plot(gd, squeeze(PeakON(:,paindex,:)), sprintf('%so', colorON(N)))
        end
        set(gca, 'xtick', 1:numgapdurs)
        set(gca, 'xticklabel', gapdurs)
        title ('LaserON/OFF startle')
        xlabel('gap duration')
        ylabel('startle responses all trials')
        if N==2
            legend('pre','post')
        end
    end
    
    
    
    
end


return
% code from here on is a mess, a mixture of re-computing
% integrated startle (already done in ProcessGPIAS_Behavior) and writing
% stuff to text files.

% Plot the integration of the abs(trace)
% Seperated into 2 figures with laser OFF/ON
if true
    
    LaserOFFstartle=nan(numgapdurs, numpulseamps,max(nrepsOFF));
    LaserONstartle=nan(numgapdurs, numpulseamps,max(nrepsON));
    
    %plot the mean Laser OFF tuning curve
    if ~isempty (mM1OFF)
        
        figure;
        p=0;
        subplot1(numgapdurs,1)
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            hold on
            
            %convert Pre/PostStartleWindow to samples:
            switch soaflag
                case 'soa' %soa actually refers to soa
                    PreStartleWindow=1+(PreStartleWindowms-xlimits(1))*samprate/1000;
                    PostStartleWindow=(PostStartleWindowms-xlimits(1)+soa)*samprate/1000;
                case 'isi' %soa actually refers to isi
                    PreStartleWindow=1+(PreStartleWindowms-xlimits(1)-gapdurs(gdindex))*samprate/1000;
                    PostStartleWindow=(PostStartleWindowms-xlimits(1)+soa)*samprate/1000;
            end
            fprintf('\ngap %d: PreStartleWindow: %d %d', gapdurs(gdindex), PreStartleWindow)
            fprintf('\ngap %d: PostStartleWindow: %d %d', gapdurs(gdindex), PostStartleWindow)
            
            
            % plot each trial in blue
            prestartle=[];
            poststartle=[];
            ISIamplitude=[];
            %note: PreStartleWindow, PostStartleWindow and ISIWindow set at top
            sumtrace=0;
            for i=1:nrepsOFF(gdindex, paindex)
                trace1=squeeze(M1OFF(gdindex, paindex,i,:));
                trace1=trace1-median(trace1(1:1000));
                trace1=abs(trace1);
                sumtrace=sumtrace+trace1;
                sumprestartle=sum(trace1(PreStartleWindow(1):PreStartleWindow(2)));
                prestartle=[prestartle sumprestartle];
                sumstartle=sum(trace1(PostStartleWindow(1):PostStartleWindow(2)));
                poststartle=[poststartle sumstartle];
                sumISIamplitude=sum(trace1(ISIWindow(1):ISIWindow(2)));%added by APW 3_31_14
                ISIamplitude=[ISIamplitude sumISIamplitude];%added by APW 3_31_14
                clear t sumprestartle sumstartle sumISIamplitude %added by APW 3_31_14
            end
            if nrepsOFF(gdindex, paindex)>0
                % add the gap stimulus in magenta
                switch soaflag{:}
                    case 'soa' %soa actually refers to soa
                        line([0 gapdurs(gdindex)],[0 0],'color','m','linewidth',10);
                    case 'isi' %soa actually refers to isi
                        line([-gapdurs(gdindex) 0],[0 0],'color','m','linewidth',10);
                end
                % add the startle stimulus in magenta
                line([soa soa+pulsedur],[0 0],'color','g','linewidth',10);
                
                sumtrace=sumtrace./nrepsOFF(gdindex, paindex);
                trace2=[0 sumtrace' 0];
                
                t=1:length(trace1);
                t=t/10;
                t=t+xlimits(1);
                t=[t(1) t t(end)];
                patch(t,trace2,'b','edgecolor','b')
                
                
                if length(prestartle)==length(poststartle)
                    [H,P,CI,STATS] = ttest(prestartle,poststartle,[],'left');
                else
                    [H,P,CI,STATS] = ttest2(prestartle,poststartle,[],'left');
                end
                LaserOFFstartle(gdindex, paindex,1:length(poststartle))=poststartle;
                LaserOFFprestartle(gdindex, paindex,1:length(prestartle))=prestartle;
                LaserOFFISIamplitude(gdindex, paindex,1:length(ISIamplitude))=ISIamplitude; %added by APW 3_31_14
                
                % for gdindex=1
                %    gap0LaserOFFstartle=LaserOFFstartle
                % end
                
                % plot the mean trace in red
                trace1=squeeze(mM1OFF(gdindex, paindex,:));
                trace1=trace1-median(trace1(1:1000));
                t=1:length(trace1);
                t=t/10;
                t=t+xlimits(1);
                plot(t, abs(trace1), 'r');
                
                % add the actual integration windows in grey
                %                    plot(t([PreStartleWindow(1) PreStartleWindow(2)]), -.1*diff(ylimits)+0*t([PreStartleWindow(1) PreStartleWindow(2)]), 'color',[.8 .8 .8], 'linewidth', 4);
                %                    plot(t([PostStartleWindow(1) PostStartleWindow(2)]), -.1*diff(ylimits)+0*t([PostStartleWindow(1) PostStartleWindow(2)]), 'color',[.8 .8 .8], 'linewidth', 4);
                %         ch=get(gca, 'children');
                %         set(gca, 'children', ch([3:length(ch)-1 1 2 length(ch)]));
                
                %ylim(ylimits)
                xlim(xlimits)
                ylabel(sprintf('%d ms',gapdurs(gdindex)));
                yl=ylim;
                text(xlimits(1)+10,yl(2)/2,sprintf('n=%d\np = %.3f',nrepsOFF(gdindex, paindex),P))
                text(xlimits(1)+10,yl(1)/2,sprintf('Pre-startle mean = %.1f +/- %.1f\n(%d:%d samples)',mean(prestartle),std(prestartle), PreStartleWindow(1),PreStartleWindow(2)))
                text(50,yl(1)/2,sprintf('Post-startle mean = %.1f +/- %.1f\n(%d:%d samples)',mean(poststartle),std(poststartle), PostStartleWindow(1),PostStartleWindow(2)))
                
                clear prestartle poststartle
            else %nreps==0, missing data for this gapdur
                LaserOFFstartle(gdindex, paindex,1:max(nrepsOFF(:)))=nan;
                LaserOFFprestartle(gdindex, paindex,1:max(nrepsOFF(:)))=nan;
                LaserOFFISIamplitude(gdindex, paindex,1:max(nrepsOFF(:)))=nan;%added by APW 3_31_14
                
            end
        end
        subplot1(1)
        h=title(sprintf('Laser OFF -- Integration, %s-%s-%s', expdate, session, filenum));
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
            
            % plot each trial in blue
            prestartle=[];
            poststartle=[];
            sumtrace=0;
            for i=1:nrepsON(gdindex, paindex)
                trace1=squeeze(M1ON(gdindex, paindex,i,:));
                trace1=trace1-median(trace1(1:1000));
                trace1=abs(trace1);
                sumtrace=sumtrace+trace1;
                sumprestartle=sum(trace1(PreStartleWindow(1):PreStartleWindow(2)));
                prestartle=[prestartle sumprestartle];
                sumstartle=sum(trace1(PostStartleWindow(1):PostStartleWindow(2)));
                poststartle=[poststartle sumstartle];
                clear t sumprestartle sumstartle
            end
            if nrepsON(gdindex, paindex)>0
                try
                    line([laserstart-gapdelay laserstart-gapdelay+laserwidth],-2.5*offset*[1 1],'color','c','linewidth',5);   % AW commented out to enable no laser gap detection plots
                end
                
                % add the gap stimulus in magenta
                switch soaflag{:}
                    case 'soa' %soa actually refers to soa
                        line([0 gapdurs(gdindex)],[0 0],'color','m','linewidth',10);
                    case 'isi' %soa actually refers to isi
                        line([-gapdurs(gdindex) 0],[0 0],'color','m','linewidth',10);
                end
                
                
                % add the startle stimulus in red
                if length(pulsedurs)~=1
                    warning('This function only knows how to plot one pulsedur') %#ok
                elseif length(pulsedurs)==1
                    line([soa soa+pulsedurs],[0 0],'color','g','linewidth',10);
                end
                
                sumtrace=sumtrace./nrepsON(gdindex, paindex);
                trace2=[0 sumtrace' 0];
                
                t=1:length(trace1);
                t=t/10;
                t=t+xlimits(1);
                t=[t(1) t t(end)];
                patch(t,trace2,'b','edgecolor','b')
                if length(prestartle)==length(poststartle)
                    [H,P,CI,STATS] = ttest(prestartle,poststartle,[],'left');
                else
                    [H,P,CI,STATS] = ttest2(prestartle,poststartle,[],'left');
                end
                LaserONstartle(gdindex, paindex,1:length(poststartle))=poststartle;
                LaserONprestartle(gdindex, paindex,1:length(prestartle))=prestartle;
                LaserONISIamplitude(gdindex, paindex,1:length(ISIamplitude))=ISIamplitude; %added by APW 3_31_14
                
                % plot the mean trace in red
                trace1=squeeze(mM1ON(gdindex, paindex,:));
                trace1=trace1-median(trace1(1:1000));
                t=1:length(trace1);
                t=t/10;
                t=t+xlimits(1);
                plot(t, abs(trace1), 'r');
                
                % add the actual integration windows in grey
                plot(t([PreStartleWindow(1) PreStartleWindow(2)]), -.1*diff(ylimits)+0*t([PreStartleWindow(1) PreStartleWindow(2)]), 'color',[.8 .8 .8], 'linewidth', 4);
                plot(t([PostStartleWindow(1) PostStartleWindow(2)]), -.1*diff(ylimits)+0*t([PostStartleWindow(1) PostStartleWindow(2)]), 'color',[.8 .8 .8], 'linewidth', 4);
                
                %ylim(ylimits)
                xlim(xlimits)
                ylabel(sprintf('%d dB',gapdurs(gdindex)));
                
                yl=ylim;
                text(xlimits(1)+10,yl(2)/2,sprintf('n=%d\np = %.3f',nrepsON(gdindex, paindex),P))
                text(xlimits(1)+10,yl(1)/2,sprintf('Pre-startle mean = %.1f +/- %.1f\n(%d:%d samples)',mean(prestartle),std(prestartle), PreStartleWindow(1),PreStartleWindow(2)))
                text(50,yl(1)/2,sprintf('Post-startle mean = %.1f +/- %.1f\n(%d:%d samples)',mean(poststartle),std(poststartle), PostStartleWindow(1),PostStartleWindow(2)))
                
            else %nreps==0, missing data for this gapdur
                LaserONstartle(gdindex, paindex,1:max(nrepsON(:)))=nan;
                LaserONprestartle(gdindex, paindex,1:max(nrepsON(:)))=nan;
                LaserONISIamplitude(gdindex, paindex,1:max(nrepsON(:)))=nan; %added by APW 3_31_14
            end
            
        end
        subplot1(1)
        h=title(sprintf('Laser ON -- Integration, %s-%s-%s', expdate, session, filenum));
        subplot1(numgapdurs)
        xlabel('Time (ms)');
        
        pos=get(gcf, 'pos');
        pos(1)=pos(1)+pos(3); %shift ON to right
        set(gcf, 'pos', pos)
        
        
        for gdindex=1:numgapdurs
            
            if sum(~isnan(LaserONstartle(gdindex, paindex,:)))==sum(~isnan(LaserOFFstartle(gdindex, paindex,:)))
                [H,P,CI,STATS] = ttest(LaserONstartle(gdindex, paindex,:),LaserOFFstartle(gdindex, paindex,:));
                if H==0
                    fprintf('\nAt %d ms, the laser did not affect the startle (ttest: p = %.3f)',gapdurs(gdindex),P);
                elseif H==1
                    fprintf('\nAt %d ms, the laser SIGNIFICANTLY changed the startle (ttest: p = %.3f)',gapdurs(gdindex),P);
                end
                
            else
                [H,P,CI,STATS] = ttest2(LaserONstartle(gdindex, paindex,:),LaserOFFstartle(gdindex, paindex,:));
                if H==0
                    fprintf('\nAt %d ms, the laser did not affect the startle (ttest2: p = %.3f)',gapdurs(gdindex),P);
                elseif H==1
                    fprintf('\nAt %d ms, the laser SIGNIFICANTLY changed the startle (ttest2: p = %.3f)',gapdurs(gdindex),P);
                end
            end
            
        end
        
        for gdindex=2:numgapdurs
            [H,P,CI,STATS] = ttest(LaserOFFstartle(1,paindex,:),LaserOFFstartle(gdindex, paindex,:));
            if H==0
                fprintf('\nGap %d: no different from Gap 0 (laser off) (ttest: p = %.3f)',gapdurs(gdindex),P);
            elseif H==1
                fprintf('\nGap %d: significantly different from Gap 0 (laser off) (ttest: p = %.3f)',gapdurs(gdindex),P);
            end
        end
        
    end
    
    
    figure;hold on
    %    errorbar(nanmean(squeeze(LaserOFFstartle(:,paindex,:))'),nanstd(squeeze(LaserOFFstartle(:,paindex,:))')/sqrt(size(LaserOFFstartle,3)), 'k-o')
    errorbar(nanmedian(squeeze(LaserOFFstartle(:,paindex,:))'),nanstd(squeeze(LaserOFFstartle(:,paindex,:))')/sqrt(size(LaserOFFstartle,3)), 'k-o')
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    %    errorbar(nanmean(squeeze(LaserONstartle(:,paindex,:))'),nanstd(squeeze(LaserONstartle(:,paindex,:))')/sqrt(size(LaserONstartle,3)), 'c-o')
    errorbar(nanmedian(squeeze(LaserONstartle(:,paindex,:))'),nanstd(squeeze(LaserONstartle(:,paindex,:))')/sqrt(size(LaserONstartle,3)), 'c-o')
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    title ('LaserON/OFF startle')
    xlabel('gap duration')
    ylabel('startle response +- sem')
    legend('Laser OFF', 'Laser ON')
    fprintf('\n\n')
    
    mean(squeeze(LaserOFFstartle(:,paindex,:))')
    mean(squeeze(LaserONstartle(:,paindex,:))')
    
    median(squeeze(LaserOFFstartle(:,paindex,:))')
    median(squeeze(LaserONstartle(:,paindex,:))')
end

txtfilename=sprintf('%s-%s-%sout.txt', expdate, session, filenum);
fid=fopen(txtfilename, 'wt');
for paindex=1:numpulseamps
    fprintf(fid, '\npulse amp: %d', pulseamps(paindex));
    
    fprintf(fid, '\ngap durs:');
    for gdindex=1:numgapdurs
        fprintf(fid, '\t%d', gapdurs(gdindex));
    end
    fprintf(fid, '\n');
    for gdindex=1:numgapdurs
        fprintf(fid, 'pre\tISI\tpost\t', gapdurs(gdindex));
    end
    
    fprintf(fid, '\nLaserON\n');
    
    for rep= 1:size(LaserONstartle, 3)
        for gdindex=1:numgapdurs
            fprintf(fid, '%f\t', LaserONprestartle(gdindex, paindex, rep));
            fprintf(fid, '%f\t', LaserONISIamplitude(gdindex, paindex, rep));
            fprintf(fid, '%f\t', LaserONstartle(gdindex, paindex, rep));
        end
        fprintf(fid, '\n');
    end
    
    fprintf(fid, '\nLaserOFF\n');
    
    for rep= 1:size(LaserOFFstartle, 3)
        for gdindex=1:numgapdurs
            fprintf(fid, '%f\t', LaserOFFprestartle(gdindex, paindex, rep));
            fprintf(fid, '%f\t', LaserOFFISIamplitude(gdindex, paindex, rep));
            fprintf(fid, '%f\t', LaserOFFstartle(gdindex, paindex, rep));
        end
        fprintf(fid, '\n');
    end
    fprintf(fid, '\n');
    fprintf(fid, '\n');
    
end
fclose(fid);






