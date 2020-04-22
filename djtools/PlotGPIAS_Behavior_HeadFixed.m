function PlotGPIAS_Behavior_HeadFixed(varargin)
dbstop if error
%plots GPIAS behavioral data from djmaus and pressure sensor device using
%one of ADC inputs
%usage: PlotGPIAS_Behavior(datadir1, datadir2)
% plots one plot if diven one dir, plots pre and post plots if given two
% dirs (pre dur, post dur, ...)


if nargin==0 datadir{1}=pwd;end
if nargin==1
    datadir{1}=varargin{1};
elseif nargin==2
    datadir{1}=varargin{1};
    datadir{2}=varargin{2};
elseif nargin==3
    datadir{1}=varargin{1};
    datadir{2}=varargin{2};
    datadir{3}=varargin{3};
elseif nargin==4
    datadir{1}=varargin{1};
    datadir{2}=varargin{2};
    datadir{3}=varargin{3};
    datadir{4}=varargin{4};
elseif nargin==5
    datadir{1}=varargin{1};
    datadir{2}=varargin{2};
    datadir{3}=varargin{3};
    datadir{4}=varargin{4};
    datadir{5}=varargin{5};
end

num=length(datadir);
colorOFF=['k','b','r','c','y','m'];

colorON=['c','r','y'];
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
    load('notebook.mat')
    lines{N}=nb.notes;
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
    startle_window=[0 75];
    %plot the tuning curves
    
    % Plot the actual trace with mean trace overlayed
    % Separated into 2 figures with laser OFF/ON
    if true
        if ~isempty (mM1OFF)
            start=(startle_window(1)-xlimits(1))*samprate/1000;
            stop=start+diff(startle_window)*samprate/1000;
            %plot the mean Laser OFF tuning curve
            for paindex=1:numpulseamps
                %figure(N);
                p=0;
                %subplot1(numgapdurs,1)
                for gdindex=1:numgapdurs
                    p=p+1;
                    %subplot1(p)
                    figure; %plots traces in separate windows for better visibility
                    
                    hold on
                    
                    offset=10*std(M1OFF(:));
                    
                    % add the stimulus in magenta
                    stimtrace=squeeze(mM1OFFstim(gdindex, paindex,:));
                    stimtrace=stimtrace-median(stimtrace(1:1000));
                    t=1:length(stimtrace);
                    t=1e3*t/samprate;
                    t=t+xlimits(1);
                    plot(t, .2*stimtrace-offset*2, 'm')
                    
                    % plot each trial in blue
                    for i=1:nrepsOFF(gdindex, paindex)
                        trace1=squeeze(M1OFF(gdindex, paindex,i,:));
                        trace1=trace1-median(trace1(1:500));
                        plot(t, trace1);
                    end
                    
                    % plot the mean trace in red
                    trace1=squeeze(mM1OFF(gdindex, paindex,:));
                    trace1=trace1-median(trace1(1:1000));
                    
                    plot(t, trace1-offset, 'r')
                    
                    %ylim([ylimits(1)-3*offset ylimits(2)])
                    xlim(xlimits)
                    ylabel(sprintf('%d ms',gapdurs(gdindex)));
                    text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nrepsOFF(gdindex, paindex)))
                    %axis off
                    h=title(sprintf('%s:\nnreps: %d-%d, OFF',datadir{N},min(nrepsOFF(:)),max(nrepsOFF(:))));
                set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
                
%                 subplot1(numgapdurs)
                xlabel('Time (ms)');
                end
            end
            %close all
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
    figure(200);
    hold on
    gd=1:numgapdurs;
    if num>1 % calculate %off based on initial startle to 0 
        for p=1:numgapdurs
            m1=abs(mPeakOFF(1, paindex)); %first data should be pre session, post sessions will be calculated from that
            m2=abs(mPeakOFF(p, paindex));
            percentGPIAS_OFF(p)=((m1-m2)/m1)*100;
            percentGPIAS_OFF(p)=percentGPIAS_OFF(p)-percentGPIAS_OFF(1);
        end
    end
    if ~isempty(percentGPIAS_OFF)
        plot(gd, percentGPIAS_OFF, sprintf('%s-o', colorOFF(N)))
    end
    if ~isempty(percentGPIAS_ON)
        plot(gd, percentGPIAS_ON, sprintf('%s-o', colorON(N)))
    end
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    title ( 'percent GPIAS attenuation')
    xlabel('gap duration')
    ylabel('percent GPIAS (normalized to pre 0 dur)')
    if N==num
        legend(lines)
    end
    
    
    %plot the mean peak rectified startle
    for paindex=1:numpulseamps
        figure(201);
        hold on
        gd=1:numgapdurs;
        if ~isempty(mPeakOFF)
            errorbar(gd, abs(mPeakOFF(:,paindex)),semPeakOFF(:,paindex), sprintf('%s-o', colorOFF(N)))
            plot(gd, abs(squeeze(PeakOFF(:,paindex,:))), sprintf('%so', colorOFF(N)))
        end
        if ~isempty(mPeakON)
            errorbar(gd, abs(mPeakON(:,paindex)),semPeakON(:,paindex), sprintf('%s-o', colorON(N)))
        end
        set(gca, 'xtick', 1:numgapdurs)
        set(gca, 'xticklabel', gapdurs)
        title ('Laser OFF startle')
        xlabel('gap duration')
        ylabel('startle response +- sem')
%         if N==num
%             legend(lines{1},lines{2}, lines{3})
%         end
    end
    
    
    %plot all trials peak rectified startle
%     for paindex=1:numpulseamps
%         figure(202);
%         hold on
%         gd=1:numgapdurs;
%         if ~isempty(mPeakOFF)
%             errorbar(gd, mPeakOFF(:,paindex),semPeakOFF(:,paindex), sprintf('%s-o', colorOFF(N)))
%             plot(gd, squeeze(PeakOFF(:,paindex,:)), sprintf('%so', colorOFF(N)))
%         end
%         if ~isempty(mPeakON)
%             errorbar(gd, mPeakON(:,paindex),semPeakON(:,paindex), sprintf('%s-o', colorON(N)))
%             plot(gd, squeeze(PeakON(:,paindex,:)), sprintf('%so', colorON(N)))
%         end
%         set(gca, 'xtick', 1:numgapdurs)
%         set(gca, 'xticklabel', gapdurs)
%         title ('LaserON/OFF startle')
%         xlabel('gap duration')
%         ylabel('startle responses all trials')
%         if N==num
%             legend(lines{1},lines{2}, lines{3})
%         end
    %end
    
    %plot mean traces for all gap dur for pre and post
    figure(300);

    hold on
    for paindex=1:numpulseamps
        p=0;
        offset=0;
        for gdindex=1:numgapdurs
            p=p+1;
            hold on
            offset1=10*std(M1OFF(:));
            
            % add the stimulus in magenta
%             if N==num
%             stimtrace=squeeze(mM1OFFstim(gdindex, paindex,:));
%             stimtrace=stimtrace-median(stimtrace(1:1000));
%             t=1:length(stimtrace);
%             t=1e3*t/samprate;
%             t=t+xlimits(1);
%             stimtrace=stimtrace*.5; %make it smaller
%             plot(t, stimtrace-offset1+offset, 'm')
%             end
            ticks1(p)=offset;
            
            % calculate and plot SEM of a trace
            for i=1:nrepsOFF(gdindex, paindex)
                trace1=squeeze(M1OFF(gdindex, paindex,i,:));
                tracesOFF(i,:)=trace1-median(trace1(1:1000));
            end
            semTraceOFF(gdindex, paindex,:)=std(tracesOFF)/sqrt(nrepsOFF(gdindex, paindex));
            % plot the mean trace in red
            trace1=squeeze(mM1OFF(gdindex, paindex,:));
            trace1=trace1-median(trace1(1:1000));
            color1=sprintf('-%s',colorOFF(N)); %plot pre, poset 0 and post 1 in different colors
            b=shadedErrorBar(t, trace1+offset, squeeze(semTraceOFF(gdindex, paindex,:)), 'lineProps',color1 );
            
            hold on
           % plot(t, ones(length(trace1))*offset, 'Color', [.5 .5 .5]);
            offset=offset-5;
            
            xlim([-100 200])
            if N==3
            
            %ylim([ylimits1(1)-offset ylimits1(2)])
            end
            
            %ylabels(p,:)=sprintf('%d ms',gapdurs(gdindex));
            %axis off
        end
        
    end
    %legend(b,lines{N});
    if N==num 
        load('notebook.mat')
        xlabel('Time (ms)');
        h=title(sprintf('Mean traces for all gap durations mouse %s',nb.mouseID));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            yticks(fliplr(ticks1))
            yticklabels(fliplr(gapdurs))
            ylabel('gap durations (ms)')
            %legend(lines)
            plot(zeros(60), -54:5,  'Color', [.5 .5 .5]);
            ylim([-30 5])
            h=gcf;
            set(h,'Position', [680 42 814 954])
    end
end


return
