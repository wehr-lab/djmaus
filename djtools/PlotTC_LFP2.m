function PlotTC_LFP2(varargin)

% plots continuous tuning curve data from djmaus
% using new OpenEphys and kilosort file formats and hierarchy
%
% usage: PlotTC_LFP([datapath], [xlimits],[ylimits])
% all inputs are optional
% defaults to datapath=pwd, xlimits [0 200], ylimits autoscaled
% Plots all continuous channels
% Processes data if outfile is not found;

% set filtering flag and cutoffs (if isempty(lo_pass_cutoff), then hipass created)
flag.filt = 0;
hi_pass_cutoff=3000;
lo_pass_cutoff=300;
printtofile=1; %print figures to postscript file
closewindows=1; %close windows as soon as you print them
write_depth_textfile=1; %create or edit depth.txt file
force_reprocess=0;

if nargin==0
    datadir=pwd;
else
    datadir=varargin{1};
end
if isempty(datadir) datadir=pwd;end

try
    xlimits=varargin{2};
catch
    xlimits=[];
end
try
    ylimits=varargin{3};
catch
    ylimits=[];
end


if flag.filt
    error('filtering not implemented for multichannel probe data')
    if isempty(lo_pass_cutoff)
        [a,b]=butter(1, hi_pass_cutoff/(30e3/2), 'high');
        fprintf('\hi-pass filtering [%d]\n', hi_pass_cutoff)
    else
        [b,a]=butter(2, [lo_pass_cutoff hi_pass_cutoff]/(30e3/2));
        fprintf('\nband-pass filtering [%d-%d]\n', lo_pass_cutoff, hi_pass_cutoff)
    end
end

outfilename='outLFP.mat';
cd(datadir)

try
    dOEinfo=dir('OEinfo*.mat'); %if we're in bonsai dir
    load(fullfile(dOEinfo.folder, dOEinfo.name))
catch %if we're in ephys dir
    dOEinfo=dir('../OEinfo*.mat');
    load(fullfile(dOEinfo.folder, dOEinfo.name))
end

if printtofile
    pdffilename=sprintf('%s-LFP-figs.pdf', BonsaiFolder);
    delete(pdffilename)
end

if force_reprocess
    fprintf('\nForce Re-process')
    fprintf('\ncalling ProcessTC_LFP2')
    ProcessTC_LFP2(datadir, xlimits, ylimits);
end

cd(BonsaiPath)
cd(EphysPath)

if printtofile
    pdffilename=sprintf('%s-LFP-figs.pdf', BonsaiFolder);
    delete(pdffilename)
end

d=dir(outfilename);
if ~isempty(d)
    tic
    fprintf('\nloading outfile %s ...', outfilename);
    load(outfilename)
    fprintf('\tdone');
    toc
else
    fprintf('\ndid not find outfile, calling ProcessTC_LFP2')
    ProcessTC_LFP2(datadir, xlimits, ylimits);
    load(outfilename);
end

%hard coded for P128-2 distance=20 um, 64 ch on each shank (1260 um total), shanks are 500 um apart
%CSD(2 shanks, 64 ch, duration)
%should rewrite this to flexibly use various probes
pitch=20;
chans_per_shank=64;


M1stim=out.M1stim;

freqs=out.freqs;
amps=out.amps;
durs=out.durs;
nreps=out.nreps;
numfreqs=out.numfreqs;
numamps=out.numamps;
numdurs=out.numdurs;
samprate=out.samprate; %in Hz
% M1=out.M1;
% scaledtrace=out.scaledtrace;
% traces_to_keep=out.traces_to_keep;
% mM1=out.mM1;
mM1ON=out.mM1ON;
mM1OFF=out.mM1OFF;
% M1OFF=out.M1OFF;
mM1ONLaser=out.mM1ONLaser;
% M1ONLaser=out.M1ONLaser;
mM1OFFLaser=out.mM1OFFLaser;
mM1ONStim=out.mM1ONStim;
mM1OFFStim=out.mM1OFFStim;
nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;
IL=out.IL; %whether there were interleaved laser trials or not
if isempty(xlimits) xlimits=out.xlimits; end
fprintf('\nusing xlimits [%d-%d]\n', xlimits(1), xlimits(2))

% %find optimal axis limits
if isempty(ylimits)
    ylimits=[0 0];
    for dindex=1:numdurs
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                trace1=squeeze(mM1OFF(findex, aindex, dindex, :));
                if flag.filt
                    %trace1=filtfilt(b,a,trace1);
                    trace1 = bandpass(trace1,[lo_pass_cutoff hi_pass_cutoff],samprate/2);
                end
                trace1=trace1-mean(trace1(1:100));
                if min([trace1])<ylimits(1); ylimits(1)=min(trace1);end
                if max([trace1])>ylimits(2); ylimits(2)=max(trace1);end
            end
        end
    end
end

% ylimits=round(ylimits*100)/100;

%plot the mean tuning curve BOTH
if IL
    for dindex=1:numdurs
        figure
        p=0;
        subplot1(numamps,numfreqs)
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                traces1=squeeze(mM1ON(:,findex, aindex, dindex, :));
                traces2=squeeze(mM1OFF(:,findex, aindex, dindex, :));

                traces1=traces1 -mean(traces1(:,1:100), 2);
                traces1=traces1(out.channelorder,:); %re-order to be in probe channel order
                inc=ylimits(2);
                offsets1=inc*[1:size(traces1, 1)]';
                traces2=traces2 -mean(traces2(:,1:100), 2);
                traces2=traces2(out.channelorder,:); %re-order to be in probe channel order
                inc=ylimits(2);
                offsets1=inc*[1:size(traces2, 1)]';

                if flag.filt
                    %trace1=filtfilt(b,a,trace1);
                    trace1 = bandpass(trace1,[lo_pass_cutoff hi_pass_cutoff],samprate/2);
                    %trace2=filtfilt(b,a,trace2);
                    trace2 = bandpass(trace2,[lo_pass_cutoff hi_pass_cutoff],samprate/2);
                end

                t=1:length(traces(1,:));
                t=1000*t/samprate; %convert to ms
                t=t+xlimits(1); %correct for xlim in original processing call
                line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                hold on;
                plot(t, traces1+offsets, 'b');
                plot(t, traces2+offsets, 'k');
                %   ylim(ylimits)
                xlim(xlimits)
                box off

            end
        end
        subplot1(1)
        h=title(sprintf('%s: %dms, nreps: %d-%d, ON&OFF, ch%d',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF))), channel));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')

        %label amps and freqs
        p=0;
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                if findex==1
                    text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    if freqs(findex)>0
                        text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                    elseif freqs(findex)==-1000
                        text(xlimits(1), vpos, 'WN')
                    elseif freqs(findex)==-2000
                        text(xlimits(1), vpos, 'SS')
                    end
                end
            end
        end
    end
end

%samprate=30e3;
reps_to_use=[];


%plot the mean tuning curve OFF
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            traces=squeeze(mM1OFF(:,findex, aindex, dindex, :));
            if flag.filt
                %trace1=filtfilt(b,a,trace1);
                trace1 = bandpass(trace1,[lo_pass_cutoff hi_pass_cutoff],samprate/2);
            end
            traces=traces -mean(traces(:,1:100), 2);
            traces=traces(out.channelorder,:); %re-order to be in probe channel order
            inc=ylimits(2);
            offsets=inc*[1:size(traces, 1)]';
            offsets(chans_per_shank)=inc*20;

            if ~isempty(mM1OFFLaser)
                Lasertrace=squeeze(mM1OFFLaser(findex, aindex, dindex, :));
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=.05*diff(ylimits)*Lasertrace;
            else
                Lasertrace = traces(1,:)*0;
            end
            if ~isempty(mM1OFFStim)
                Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.05*diff(ylimits)*Stimtrace;
            else
                Stimtrace = traces(1,:)*0;
            end
            t=1:length(traces(1,:));
            t=1000*t/samprate; %convert to ms
            t=t+xlimits(1); %correct for xlim in original processing call
            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            hold on; plot(t, traces+offsets, 'k');
            offset=ylimits(1)+.1*diff(ylimits);
            plot(t, Stimtrace+offset, 'm', t, Lasertrace+offset, 'c')
            try
                % ylim(ylimits)
            end
            xlim(xlimits)
            xlabel off
            ylabel off

            %label channels and shanks
            for i=1:128
                str{i}=sprintf('ch%d', i);
                if i<=64 sh{i}=1;
                elseif i>64 sh{i}=2;
                end
            end
            text(repmat(-30, 1, 128), offsets, str, 'fontsize', 6)
            text(repmat(-25, 1, 128), offsets, sh, 'fontsize', 6)

            %label amps and freqs

            subplot1(p)
            if findex==1
                text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
            end
            %mM1ONLaser=mM1ONLaser;
            axis off
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s: %dms, nreps: %d-%d',datadir,durs(dindex),min(nrepsOFF(:)),max(nrepsOFF(:))));
    %h=title(sprintf('OFF %s: %dms, nreps: %d-%d',datadir,durs(dindex), reps_to_use));
    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
    set(h,  'interpreter', 'none')

    %label amps and freqs
    p=0;
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            if findex==1
                text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
            end

            if aindex==1
                if mod(findex,2) %odd freq
                    vpos=ylimits(1)-mean(ylimits);
                else
                    vpos=ylimits(1)-mean(ylimits);
                end
                if freqs(findex)==-2000
                    text(xlimits(1), vpos, 'SS')
                elseif freqs(findex)==-1000
                    text(xlimits(1), vpos, 'WN')
                elseif freqs(findex)>0
                    text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                end
            else
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    if freqs(findex)>0
                        text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                    elseif freqs(findex)==-1000
                        text(xlimits(1), vpos, 'WN')
                    elseif freqs(findex)==-2000
                        text(xlimits(1), vpos, 'SS')
                    end
                end

                %                 if freqs(findex)>0
                %                     text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                %                 elseif freqs(findex)==-1000
                %                     text(xlimits(1), vpos, 'WN')
                %                 elseif freqs(findex)==-2000
                %                     text(xlimits(1), vpos, 'SS')
                %                 end
            end
            %             if findex==numfreqs && aindex==numamps
            %                 axis on
            %                 ylab=[ceil(ylimits(1)*10)/10 floor(ylimits(2)*10)/10];
            %                 set(gca,'ytick',ylab,'yticklabel',ylab,'YAxisLocation','right')
            %             end
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s \n%dms, nreps: %d-%d',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')
end
set(gcf, 'pos', [ 998         198         660        1100])



% plot on
if IL
    error ('plot laser on not finished yet')
    for dindex=1:numdurs
        figure
        p=0;
        subplot1(numamps,numfreqs)
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                axis off

                traces=squeeze(mM1ON(:,findex, aindex, dindex, :));
                if flag.filt
                    %trace1=filtfilt(b,a,trace1);
                    trace1 = bandpass(trace1,[lo_pass_cutoff hi_pass_cutoff],samprate/2);
                end
                traces=traces -mean(traces(:,1:100), 2);
                traces=traces(out.channelorder,:); %re-order to be in probe channel order
                inc=ylimits(2);
                offsets=inc*[1:size(traces, 1)]';

                Stimtrace=squeeze(mM1ONStim(findex, aindex, dindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.05*diff(ylimits)*Stimtrace;

                t=1:length(traces(1,:));
                t=1000*t/samprate; %convert to ms
                t=t+xlimits(1); %correct for xlim in original processing call
                line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                hold on; plot(t, traces+offsets, 'k');
                offset=ylimits(1)+.1*diff(ylimits);
                plot(t, Stimtrace+offset, 'm', t, Lasertrace+offset, 'c')



                for rep=1:nrepsON(findex, aindex, dindex)
                    Lasertrace=squeeze(M1ONLaser(findex, aindex, dindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=.05*diff(ylimits)*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end

                % ylim(ylimits)
                xlim(xlimits)
                axis off

            end
        end
        subplot1(1)
        h=title(sprintf('ON %s: %dms, nreps: %d-%d, ch%d',datadir,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON))), channel));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')

        %label amps and freqs
        p=0;
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                if findex==1
                    text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                end
                %             if findex==numfreqs && aindex==numamps
                %                 axis on
                %                 ylab=[ceil(ylimits(1)*10)/10 floor(ylimits(2)*10)/10];
                %                 set(gca,'ytick',ylab,'yticklabel',ylab,'YAxisLocation','right')
                %             end
            end
        end
    end
end



traces=squeeze(mM1OFF(:,findex, aindex, dindex, :));
%traces=traces -mean(traces(:,1:100), 2);
traces=traces(out.channelorder,:);
clear depthstr


%Ira smoothed the traces like this:
% smooth_traces=[];
% for j=2:size(traces,1)-1 %smoothing to get rid of variance, -2 channels
%     smooth_traces(j-1,:)=(1/4)*(traces(j+1,:)+2*traces(j,:)+traces(j-1,:));
% end
% for j=2:size(smooth_traces,1)-1 %cannot compute CSD for first and last channels, - 2 channels
%     CSD(j-1,:)=(smooth_traces(j-1,:)+smooth_traces(j+1,:)-2*smooth_traces(j,:))/distance^2;
% end

traces_by_shank(1, 1:chans_per_shank, :)=traces(1:chans_per_shank, :);
traces_by_shank(2, 1:chans_per_shank, :)=traces(chans_per_shank+1:2*chans_per_shank, :);
%putting into traces_by_shank because Ira's smoothing (and possibly the csd calculation) was wrapping
%around to inputs from channels on the next shank

% % %Ira's smoothing method to get rid of variance, -2 channels
% % smooth_traces=traces_by_shank;
% % for shank=1:2
% % for j=2:chans_per_shank-1 %smoothing to get rid of variance, -2 channels
% %     smooth_traces(shank,j-1,:)=(1/4)*(traces_by_shank(shank, j+1,:)+2*traces_by_shank(shank, j,:)+traces_by_shank(shank, j-1,:));
% % end
% % end
% % traces_by_shank=smooth_traces;

%Ira's smoothing method to get rid of variance, -2 channels, trying to
%include edge channels too
smooth_traces=traces_by_shank;
for shank=1:2
for j=2:chans_per_shank-1 %smoothing to get rid of variance, -2 channels
    smooth_traces(shank,j,:)=(1/4)*(traces_by_shank(shank, j+1,:)+2*traces_by_shank(shank, j,:)+traces_by_shank(shank, j-1,:));
end
end
traces_by_shank=smooth_traces;


for shank=1:2
    for j=2:chans_per_shank-1 %cannot compute CSD for first and last channels, - 2 channels
            CSD(shank, j-1,:)=(traces_by_shank(shank, j-1,:)+traces_by_shank(shank, j+1,:)-2*traces_by_shank(shank,j,:))/pitch^2;
            depth(shank, j-1)=pitch*j;
            depthstr{j-1}=sprintf('%.0f',depth(shank, j-1));
    end
end

% % old way
% for shank=1:2
%     for j=2:chans_per_shank-1 %cannot compute CSD for first and last channels, - 2 channels
%         if shank==1
%             CSD(shank, j-1,:)=(traces(j-1,:)+traces(j+1,:)-2*traces(j,:))/pitch^2;
%             depth(shank, j-1)=pitch*j;
%             depthstr{j-1}=sprintf('%.0f',depth(shank, j-1));
%         elseif shank==2
%             k=j+64;
%             CSD(shank, j-1,:)=(traces(k-1,:)+traces(k+1,:)-2*traces(k,:))/pitch^2;
%             depth(shank, j-1)=pitch*j;
%         end
%     end
% end
% depthstr=fliplr(depthstr);

%this puts a marker in the lower left corner to check that we're correctly
%plotting the deepest channel at the bottom
% CSD(1,62,1:100)=.05; %check depth
% CSD(2,62,1:100)=-.05; %check depth

numxticks=7;
numyticks=6;
figure
subplot(121)
imagesc(squeeze(CSD(1,:,:)));
caxis([min(CSD(:)) max(CSD(:))])
xticks(linspace(1, size(CSD,3), numxticks))
xticklabels(linspace(xlimits(1), xlimits(2), numxticks))
yticks(linspace(1, 60, numyticks))
yticklabels(depthstr(round(yticks)))
xlabel('time, ms')
ylabel('depth, um')
h=title(sprintf('CSD shank 1 %s \n %dms, nreps: %d',BonsaiFolder,durs(dindex),min(nrepsOFF(:))));
set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')

subplot(122)
imagesc(squeeze(CSD(2,:,:)));
caxis([min(CSD(:)) max(CSD(:))])
colorbar
colormap(parula)
xticks(linspace(1, size(CSD,3), numxticks))
xticklabels(linspace(xlimits(1), xlimits(2), numxticks))
yticklabels('')
xlabel('time, ms')
ylabel('depth, um')
title('shank 2')
set(gca, 'pos', [ 0.5300    0.1100    0.3347    0.8150])

if printtofile
    pdffilename=sprintf('%s-LFP-figs.pdf', BonsaiFolder);
    %print figures to postscript file
    f=findobj('type', 'figure');
    for idx=1:length(f)
        pause(.5)
        %figure(f(idx))
        % orient landscape
        % % print figs -dpsc2 -append -bestfit
        exportgraphics(f(idx),pdffilename,'Append',true)
        pause(.5)

        if closewindows
            close
        end
    end
end

