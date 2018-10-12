function PlotLFPs_depth(varargin)
% plots LFPs from tone or WN stimulus according to depth
if nargin==0
    fprintf('\nno input');
    return;
end

djPrefs;
global pref
datadir=pwd;
datadir=varargin{1};
cd(datadir)
depth=varargin{2};
xlimits=[-100 300];
samprate=30e3;
low_pass_cutoff=300;
ylimits=[0 0];

depth=sqrt(depth^2+depth^2); %find diagonal line


tetrodes=dir('*.spikes');
if length(tetrodes)==8
    if exist('outLFP_ch32.mat') %to avoid processing twice
        load('outLFP_ch32.mat')
        if xlimits~=out.xlimits
            for i=1:32
                ProcessTC_LFP(pwd, i, xlimits, [])
            end
        else
            fprintf('\nfound an outfile, will skip processing')
        end
    else
        for i=1:32
            ProcessTC_LFP(datadir, i, xlimits, [])
        end
    end
    order=[8, 24, 2, 29, 7, 26, 15, 21, 11, 23, 12, 28, 6, 18, 13, 22, 5, 27, 4, 31, 10, 20, 9, 25, 14, 30, 3, 19, 16, 32, 1, 17];
    figure;
    %subplot1(32,1);
    offset=0;
    for i=1:32
        j=order(i);
        d{i}=sprintf('%.0f',depth-25*i-50);
        filename=sprintf('outLFP_ch%d.mat', j);
        load(filename);
        mM1OFF=out.mM1OFF;
        freqs=out.freqs;
        WN=find(freqs==-1000);
        trace=squeeze(mM1OFF(WN,1,1,:)); %use only WN stimulus for now
        amps=out.amps;
        durs=out.durs;
        %filter the trace in case it's raw data
        [b,a]=butter(1, low_pass_cutoff/(samprate/2), 'low');
        trace1=filtfilt(b,a,trace);
        trace1=trace1 -mean(trace1(1:100));
        t=1:length(trace1);
        t=1000*t/out.samprate; %convert to ms
        t=t+xlimits(1); %correct for xlim in original processing call
        %     line([0 0+durs], ylimits(1)+[0 0]+offset, 'color', 'm', 'linewidth', 5)
        hold on; plot(t, trace1-offset, 'k');
        
        %find ylimits
        if min([trace1])<ylimits(1); ylimits(1)=min([trace1]);end
        if max([trace1])>ylimits(2); ylimits(2)=max([trace1]);end
        yticks1(i)=mean(trace1(1:100)-offset);
        
        
        offset=offset+250;
        ylimits1(i,:,:)=ylimits; %use max and min of all traces later to adjust axis
        xlim(xlimits)
        
        if i==1
            subplot1(1)
            h=title(sprintf('%s, depth=%.2f',datadir, depth));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')
        end
        axis on
        
    end
elseif length(tetrodes)==4
    order= [13 29 4 20 15 31 3 19 16 32 1 17 2 18 14 30];
    if exist('outLFP_ch16.mat') %to avoid processing twice
        load('outLFP_ch16.mat')
        if xlimits~=out.xlimits
            for i=1:16
                ProcessTC_LFP(datadir, i, xlimits, [])
            end
        else
            fprintf('\nfound an outfile, will skip processing')
        end
    else
        for i=1:16
            ProcessTC_LFP(datadir, i, xlimits, [])
        end
    end
    for i=1:16
        j=order(i);
        d{i}=sprintf('%.0f',depth-50*i-25);
        filename=sprintf('outLFP_ch%d.mat', j);
        load(filename);
        mM1OFF=out.mM1OFF;
        freqs=out.freqs;
        WN=find(freqs==-1000);
        trace=squeeze(mM1OFF(WN,1,1,:)); %use only WN stimulus for now
        amps=out.amps;
        durs=out.durs;
        %filter the trace in case it's raw data
        [b,a]=butter(1, low_pass_cutoff/(samprate/2), 'low');
        trace1=filtfilt(b,a,trace);
        trace1=trace1 -mean(trace1(1:100));
        t=1:length(trace1);
        t=1000*t/out.samprate; %convert to ms
        t=t+xlimits(1); %correct for xlim in original processing call
        %     line([0 0+durs], ylimits(1)+[0 0]+offset, 'color', 'm', 'linewidth', 5)
        hold on; plot(t, trace1-offset, 'k');
        
        %find ylimits
        if min([trace1])<ylimits(1); ylimits(1)=min([trace1]);end
        if max([trace1])>ylimits(2); ylimits(2)=max([trace1]);end
        yticks1(i)=mean(trace1(1:100)-offset);
        
        
        offset=offset+250;
        ylimits1(i,:,:)=ylimits; %use max and min of all traces later to adjust axis
        xlim(xlimits)
        
        if i==1
            subplot1(1)
            h=title(sprintf('%s, depth=%.2f',datadir, depth));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')
        end
        axis on
    end
end

yticks1=fliplr(yticks1); %reverse to start at lowest and increase
yticks([yticks1])
yticklabels(d)
%plot(zeros(int16(abs(yticks1(1)))),yticks1(1):-1,'r') %plot a line of sound onset
ylim([yticks1(1)-500, 250])
set(gcf, 'Position', [680 22 546 956])
