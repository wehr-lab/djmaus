function PlotTC_LFP2(varargin)

% plots continuous tuning curve data from djmaus
% detects bad channels based on impedance test or manually, and replaces with
% means of neighboring channels
% plots CSDs vs depth assuming 128-2 probe and channel map stored in outfile
%   note that channels 1 and 65 are at the bottom, and 64 and 128 at the
%   surface.
% asks you to select the L3/4 sink, and the L5/6 sink, and writes it to sink_chans.txt
% calculates corrected depths and writes to file depths.mat
%   corrected_depth is depth along the shank, corrected so that the L3/4 sink is at 400 µm
%   depths.mat also contains corrected_depth sink_chans
%
% using new OpenEphys and kilosort file formats and hierarchy
%
% usage: PlotTC_LFP([datapath], [xlimits],[ylimits])
% all inputs are optional
% defaults to datapath=pwd, xlimits [0 200], ylimits autoscaled
% Plots all continuous channels
% Processes data if outfile is not found;a

% set filtering flag and cutoffs (if isempty(lo_pass_cutoff), then hipass created)
flag.filt = 0;
hi_pass_cutoff=3000;
lo_pass_cutoff=300;
printtofile=1; %print figures to postscript file
closewindows=1; %close windows as soon as you print them
%write_depth_textfile=1; %create or edit depth.txt file
force_reprocess=0;
interactive=1; %1 asks user to confirm bad channels and impedance file, asks user to
%     select sinks, and saves bad channels, sinks, and depths files.
%if you set interactive=0, it will run without any user input, and save an
%     outfile, but will not save bad_channels or calculate corrected depth. This
%     is just to precompute the outfile to save time when later plotting it
%     interactively to select the sinks. It saves about 4-5 minutes.


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
    load dirs.mat
    cd(Bdirs{1}) %go to Bonsai Folder
catch
    try
        load bdirs.mat
        cd(Bdirs{1}) %go to Bonsai Folder
    catch
        ProcessTC_LFP2(datadir, xlimits, ylimits);
    end
end

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

tmp=split(macifypath(BonsaiPath), '/');
BonsaiFolder=tmp{end}; %remove absolute path
DataRoot=macifypath(DataRoot); %does nothing if you're on windows
cd(DataRoot)
cd(BonsaiFolder)
cd(EphysPath)

if printtofile
    pdffilename=sprintf('%s-LFP-figs.pdf', BonsaiFolder);
    delete(pdffilename)
end

d=dir(outfilename);
if ~isempty(d)
    tic
    fprintf('\nloading outfile %s (%.3f Gb) ...', outfilename, d.bytes/1e9);
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
            traces=traces -mean(traces(:,1:100), 2);
            traces=traces(out.channelorder,:); %re-order to be in probe channel order
            if flag.filt
                %trace1=filtfilt(b,a,trace1);
                trace1 = bandpass(trace1,[lo_pass_cutoff hi_pass_cutoff],samprate/2);
            end
            inc=1.5*ylimits(2);
            offsets=inc*[1:size(traces, 1)]';
            offsets(chans_per_shank+1:end)=offsets(chans_per_shank+1:end)+inc*20;

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
            text(repmat(-15, 1, 128), offsets, str, 'fontsize', 6)
            text(repmat(-19, 1, 128), offsets, sh, 'fontsize', 6)

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
    h=title(sprintf('OFF %s: %dms, nreps: %d-%d\nno channel replacement or smoothing',datadir,durs(dindex),min(nrepsOFF(:)),max(nrepsOFF(:))));
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
end
set(gcf, 'pos', [ 998         198         660        1100])

%plot heatmap to look for bad channels
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            traces=squeeze(mM1OFF(:,findex, aindex, dindex, :));
            traces=traces -mean(traces(:,1:100), 2);
            traces=traces(out.channelorder,:); %re-order to be in probe channel order

            pcolor(traces)
            shading interp
            %label channels and shanks
            for i=1:128
                str{i}=sprintf('ch%d', i);
                if i<=64 sh{i}=1;
                elseif i>64 sh{i}=2;
                end
            end
            text(repmat(-15, 1, 128), 1:128, str, 'fontsize', 6)
            text(repmat(-19, 1, 128), 1:128, sh, 'fontsize', 6)
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s: %dms, nreps: %d-%d\nno channel replacement or smoothing',datadir,durs(dindex),min(nrepsOFF(:)),max(nrepsOFF(:))));
    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
    set(h,  'interpreter', 'none')
    set(gcf, 'pos', [ 330         198         660        1100])
end

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

% check for bad channels and replace them with neighborhood mean

%plot raw voltage for each channel for entire recording duration
% offset=0;
% inc=1200;
% figure
% hold on
% tracelength=length(mM1OFF);
% for i=1:128
%     offset=offset+inc;
%     trace1=reshape(squeeze(out.M1OFF(out.channelorder(i),1,1,1,1:out.nrepsOFF,:))', 1, nreps*tracelength);
%     plot(trace1+offset)
% end

%I tried to detect bad channels as those that have std>thresh
%this doesn't work very well on most sessions, and is definitely not robust across sessions
%switching to manually entering bad channels into a text file
%based on impedance test, visual inspection, or other


%plot voltage stats for each channel to help inspect for bad channels
tracelength=length(mM1OFF);
for i=1:128
    trace1=reshape(squeeze(out.M1OFF(out.channelorder(i),1,1,1,1:out.nrepsOFF,:))', 1, nreps*tracelength);
    chan_std(i)=std(trace1);
    chan_max(i)=max(abs(trace1));
end

figure
hold on
sh1=1:64;
sh2=65:128;
all_fchan_std=[];
all_fchan_max=[];
for shank=1:2
    %I had to filter each shank separately because otherwise the 64/65
    %discontinuity throws a bad-channel-like spike
    if shank==1 shx=sh1; else shx=sh2; end
    plot(shx, chan_std(shx)./max(chan_std), 'b', 'linew', 2)
    plot(shx, chan_max(shx)./max(chan_max), 'r', 'linew', 2)
    [b,a]=butter(1, .5, 'high');
    fchan_std=filtfilt(b,a,chan_std(shx));
    fchan_max=filtfilt(b,a,chan_max(shx));
    all_fchan_std=[all_fchan_std fchan_std./max(chan_std)];
    all_fchan_max=[all_fchan_max fchan_max./max(chan_max)];
    plot(shx, fchan_std./max(chan_std), 'b')
    plot(shx, fchan_max./max(chan_max), 'r')
end
legend('std voltage', 'max voltage')
guess_bad_chans=unique([find(all_fchan_std>.1) find(all_fchan_max>.1)]);
fprintf('\nguess that these might be bad channels: ')
fprintf('%d ', guess_bad_chans)


bad_channels_from_xml=ReadImpedanceTestXML; %reads impedance test xml file generated by open ephys, prints and plots results
try
    load bad_channels.txt
    fprintf('\nfound and loaded local bad_channels.txt file')

catch
    if interactive
        ButtonName = questdlg('Could not find bad_channels.txt, do you want to locate the file, manually enter bad channels, or cancel? ', ...
            'bad_channels.txt', ...
            'choose file', 'enter manually', 'cancel', 'cancel');
        switch ButtonName
            case 'choose file'
                [filename, pathname] = uigetfile('bad_channels.txt', 'choose bad_channels.txt file appropriate for this recording session')
                if filename==0     error('user cancelled, no bad_channels.txt file '), end
                load(fullfile(pathname, filename));
                ButtonName = questdlg('save this bad_channels.txt to this session folder?', ...
                    int2str(bad_channels), 'yes', 'no', 'no');
                switch ButtonName
                    case 'yes'
                        bad_channels=bad_channels(:);
                        save bad_channels.txt bad_channels -ascii
                        fprintf('\nsaved bad channels to local bad_channels.txt file')
                end

            case 'enter manually'
                prompt={sprintf('Enter bad channels:\ndefault channels are populated from ReadImpedanceTestXML \n(separated by space,comma,return,etc.)')};
                name='bad channels input';
                numlines=10+ length(bad_channels_from_xml);

                str=inputdlg(prompt,name,numlines, {int2str(bad_channels_from_xml)});
                if isempty(str)   error('user cancelled, no bad_channels.txt file '), end
                bad_channels=str2num(str{:});
                ButtonName = questdlg('save these to bad_channels.txt file in this session folder?', ...
                    int2str(bad_channels), 'yes', 'no', 'no');
                switch ButtonName
                    case 'yes'
                        bad_channels=bad_channels(:);
                        save bad_channels.txt bad_channels -ascii
                        fprintf('\nsaved manually entered bad channels to bad_channels.txt file')
                end

            case 'cancel'
                abort = questdlg('Abort, or Skip bad channel entry and proceed as if there were no bad channels?', ...
                    'Abort or Skip?', 'Abort', 'Skip', 'Skip');
                switch abort
                    case 'Abort'
                        error('user cancelled ')
                    case 'Skip'
                        fprintf('\nskip bad channel entry, pretending there are no bad channels')
                        bad_channels=[];
                end
        end % switch
    else %if not interactive
        if ~isempty(bad_channels_from_xml)
            fprintf('\nnon-interactive mode, using bad channels from xml impedance file');
        end
        bad_channels=bad_channels_from_xml; %will be empty if no xml
    end %if interactive
end


% enforce bad_channels is nx1
bad_channels=bad_channels(:);


% bad_chan_thresh=.1;
% bad_chans=unique([find(all_fchan_std>bad_chan_thresh) find(all_fchan_max>bad_chan_thresh)]);
num_bad_chans=length(bad_channels);
text(bad_channels, 1+0*bad_channels, 'B')
fprintf('\nusing %d bad channels: ', num_bad_chans)
fprintf('%d ', bad_channels)



figure
hold on
bar(chan_std./max(chan_std))
bar(chan_max./max(chan_max))
text(bad_channels, 1+0*bad_channels, 'B')
xlabel('channel')
ylabel('normalized max/std signal')
title(['channel amplitudes. B=detected bad channels: ', sprintf('%d ', bad_channels)])

bad_channels_sh1=bad_channels(find(bad_channels<=64)); %separated by shank for plot labeling
bad_channels_sh2=bad_channels(find(bad_channels>64))-64;



% replace bad chans with neighborhood mean
%this should be able to handle multiple bad channels in a row, and bad channels 1,64,65,128.
for b=bad_channels'
    if b==1
        g=b;
        while ismember(g, bad_channels)
            g=g+1;
        end
        new_trace=traces(g,:);
        traces(b,:)=new_trace;
    elseif b==128
        g=b;
        while ismember(g, bad_channels)
            g=g-1;
        end
        new_trace=traces(g,:);
        traces(b,:)=new_trace;
    elseif b==64
        g=b;
        while ismember(g, bad_channels)
            g=g-1;
        end
        new_trace=traces(g,:);
        traces(b,:)=new_trace;
    elseif b==65
        g=b;
        while ismember(g, bad_channels)
            g=g+1;
        end
        new_trace=traces(g,:);
        traces(b,:)=new_trace;
    else
        if ~ismember([b-1, b+1], bad_channels) %simplest case
            new_trace=mean(traces([b-1, b+1], :));
            traces(b,:)=new_trace;
        else
            g1=b-1;
            while ismember(g1, bad_channels)
                g1=g1-1;
            end
            g2=b+1;
            while ismember(g2, bad_channels)
                g2=g2+1;
            end
            if g2>128 %oops there's no good channel above
                g2=g1;
            end
            new_trace=mean(traces([g1, g2], :));
            traces(b,:)=new_trace;
        end
    end
end

%traces now has bad channels replaced by neighborhood means
traces=traces -mean(traces(:,1:100), 2);
offsets=inc*[1:size(traces, 1)]';
offsets(chans_per_shank+1:end)=offsets(chans_per_shank+1:end)+inc*20;

%plot the mean tuning curve OFF with bad channels replaced
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            inc=ylimits(2);
            t=1:length(traces(1,:));
            t=1000*t/samprate; %convert to ms
            t=t+xlimits(1); %correct for xlim in original processing call
            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            hold on;
            plot(t, traces+offsets, 'k');

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
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s \n%dms, nreps: %d-%d\nbad chans replaced, no channel smoothing',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')
end
set(gcf, 'pos', [ 998         198         660        1100])


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

figure
offsets1=inc*[1:64]';
subplot(121)
plot(t, squeeze(traces_by_shank(1,:,:))+offsets1, 'k');
text(repmat(-30, 1, 64), offsets1, int2str([1:64]'), 'fontsize', 6)
title('smoothed traces by shank')
axis off
subplot(122)
plot(t, squeeze(traces_by_shank(2,:,:))+offsets1, 'k');
text(repmat(-30, 1, 64), offsets1, int2str([65:128]'), 'fontsize', 6)
set(gcf, 'pos', [ 998         198         660        1100])
axis off

for shank=1:2
    for j=2:chans_per_shank-1 %cannot compute CSD for first and last channels, - 2 channels
        CSD(shank, j-1,:)=(traces_by_shank(shank, j-1,:)+traces_by_shank(shank, j+1,:)-2*traces_by_shank(shank,j,:))/pitch^2;
        % depth(shank, j-1)=pitch*j;
        % depthstr{j-1}=sprintf('%.0f',depth(shank, j-1));
    end
end


offsets2=.1*max(abs(CSD(:)))*[1:size(CSD,2)]'; %positive inc so ch1 is at bottom and 64 is at top
%trim end channels (1, 65) from bad_channels
bad_channels_sh1_trim=bad_channels_sh1(find(bad_channels_sh1>1));
bad_channels_sh1_trim=bad_channels_sh1_trim(find(bad_channels_sh1_trim<64));
bad_channels_sh2_trim=bad_channels_sh2(find(bad_channels_sh2>1));
bad_channels_sh2_trim=bad_channels_sh2_trim(find(bad_channels_sh2_trim<64));

figure
subplot(121)
%shank1
plot(t, squeeze(CSD(1,:,:))+offsets2, 'k');
ylim([0 1.1*max(offsets2)])
title('CSD')
text(repmat(-30, 1, 62), offsets2, int2str([2:63]'), 'fontsize', 6)
text(0*bad_channels_sh1_trim-15, offsets2(bad_channels_sh1_trim), 'X', 'fontsize', 12, 'color', 'r')
axis off

subplot(122)
%shank2
plot(t, squeeze(CSD(2,:,:))+offsets2, 'k');
ylim([0 1.1*max(offsets2)])
text(repmat(-30, 1, 62), offsets2, int2str([66:127]'), 'fontsize', 6)
text(0*bad_channels_sh2_trim-15, offsets2(bad_channels_sh2_trim-1), 'X', 'fontsize', 12, 'color', 'r')
set(gcf, 'pos', [ 998         198         660        1100])
axis off


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

%caxis limits excluding bad chans and edges (combined across both shanks)
cmin=min( min(min(min(CSD(1,setdiff(1:62, bad_channels_sh1_trim-1),:)))), ...
    min(min(min(CSD(2,setdiff(1:62, bad_channels_sh2_trim-1),:)))));
cmax=max( max(max(max(CSD(1,setdiff(1:62, bad_channels_sh1_trim-1),:)))), ...
    max(max(max(CSD(2, setdiff(1:62, bad_channels_sh2_trim-1),:)))));

%caxis limits for each shank individually, if that looks better
cmin1=min( min(min(min(CSD(1,setdiff(1:62, bad_channels_sh1_trim-1),:)))));
cmax1=max( max(max(max(CSD(1,setdiff(1:62, bad_channels_sh1_trim-1),:)))));
cmin2= min(min(min(CSD(2,setdiff(1:62, bad_channels_sh2_trim-1),:))));
cmax2= max(max(max(CSD(2, setdiff(1:62, bad_channels_sh2_trim-1),:))));


numxticks=7;
numyticks=6;
figure
subplot(121)
imagesc(squeeze(CSD(1,:,:)));
set(gca, 'ydir', 'normal')
% caxis([min(CSD(:)) max(CSD(:))])
caxis(.5*[cmin cmax])
xticks(linspace(1, size(CSD,3), numxticks))
xticklabels(linspace(xlimits(1), xlimits(2), numxticks))
% yticks(linspace(1, 60, numyticks))
% yticklabels(depthstr(round(yticks)))
yticks(1:size(CSD, 2))
yticklabels(yticks+1)
xlabel('time, ms')
ylabel('channel')
h=title(sprintf('CSD shank 1 %s \n %dms, nreps: %d',BonsaiFolder,durs(dindex),min(nrepsOFF(:))));
set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')
text(0*bad_channels_sh1-20, bad_channels_sh1, 'X', 'fontsize', 12, 'color', 'r')

subplot(122)
imagesc(squeeze(CSD(2,:,:)));
set(gca, 'ydir', 'normal')
caxis(.5*[cmin cmax])
colorbar
colormap(parula)
xticks(linspace(1, size(CSD,3), numxticks))
xticklabels(linspace(xlimits(1), xlimits(2), numxticks))
yticks(1:size(CSD, 2))
yticklabels(yticks+65)
xlabel('time, ms')
ylabel('channel')
title('shank 2')
text(0*bad_channels_sh2-20, bad_channels_sh2, 'X', 'fontsize', 12, 'color', 'r')
set(gca, 'pos', [ 0.5300    0.1100    0.3347    0.8150])
set(gcf, 'pos', [93 141 895 1156])

try
    fid=fopen('sink_chans.txt', 'r');
    fgetl(fid);
    sink_chans.shank1L34=str2num(extractAfter(fgetl(fid),'shank1 L3/4:'));
    sink_chans.shank1L56=str2num(extractAfter(fgetl(fid),'shank1 L5/6:'));
    sink_chans.shank2L34=str2num(extractAfter(fgetl(fid),'shank2 L3/4:'));
    sink_chans.shank2L56=str2num(extractAfter(fgetl(fid),'shank2 L5/6:'));
    fclose(fid);
    fprintf('\nfound and loaded sink channels from local sink_chans.txt file')
    %expecting this txt file to have the following structure:
    % channel of shortest-latency positive-going CSD sink % header comment
    % shank1: 44 % (a number between 2-63)
    % shank2: 108 % (a number between 65-127)

catch
    if interactive
        ButtonName = questdlg('        did not find sink_chans.txt file, would you like to identify the channels with shortest-latency current sink? ', ...
            'locate sink channel?', ...
            'Sure', 'Not now', 'Not now');
        switch ButtonName
            case 'Sure'
                prompt={'Enter the channel number for shank 1 L3/4 sink',...
                    'Enter the channel number for shank 1 L5/6 sink',...
                    'Enter the channel number for shank 2 L3/4 sink',...
                    'Enter the channel number for shank 2 L5/6 sink',...
                    };
                name='Identify the channels with shortest-latency positive-going (yellow) current sink';
                numlines=[1 120];
                 opts.Resize='on';
                 def={'','','',''}
                sink_chans_str=inputdlg(prompt,name,numlines, def, opts);

                if isempty(sink_chans_str)   warning('user cancelled, no depth.txt file written'), end
                sink_chans.shank1L34=str2num(sink_chans_str{1});
                sink_chans.shank1L56=str2num(sink_chans_str{2});
                sink_chans.shank2L34=str2num(sink_chans_str{3});
                sink_chans.shank2L56=str2num(sink_chans_str{4});
                if sink_chans.shank1L34<2 | sink_chans.shank1L34>63 | sink_chans.shank1L56<2 | sink_chans.shank1L56>63 | ...
                   sink_chans.shank2L34<66 | sink_chans.shank2L34>127 | sink_chans.shank2L56<66 | sink_chans.shank2L56>127
                    warning('invalid channels entered, must be 2-63/65-127, no sink_chans.txt file written')
                    clear sink_chans
                else

                    %write out sink_chans.txt file
                    fid=fopen('sink_chans.txt', 'w');
                    fprintf(fid, 'channel of shortest-latency positive-going CSD sink');
                    fprintf(fid, '\nshank1 L3/4: %d', sink_chans.shank1L34);
                    fprintf(fid, '\nshank1 L5/6: %d', sink_chans.shank1L56);
                    fprintf(fid, '\nshank2 L3/4: %d', sink_chans.shank2L34);
                    fprintf(fid, '\nshank2 L5/6: %d', sink_chans.shank2L56);
                    fclose(fid);
                end
            case 'Not now'
                %do nothing
        end % switch
    else %if non-interactive
        %do nothing
        fprintf('\nnon-interactive mode, skipping user selection of sink channels')
    end
end


try
    fid=fopen('penetration_angles.txt', 'r');
    fgetl(fid);
    penetration_angles.shank1=str2num(extractAfter(fgetl(fid),'shank1:'));
    penetration_angles.shank2=str2num(extractAfter(fgetl(fid),'shank2:'));
    fclose(fid);
    fprintf('\nfound and loaded penetration angles from local penetration_angles.txt file')
    %expecting this txt file to have the following structure:
    % angles (in °) of penetration with respect to orthogonal to cortical surface
    % shank1: 8 % (a number between 0-90)
    % shank2: 13 % (a number between 0-90)
    % note that these angles are reconstructed from di-I or di-O histology and recorded in
    % Projects/5XFAD/Rig3Phys/Summary of Recording Sites.xlsx among other places

catch
    if interactive
        ButtonName = questdlg('        did not find penetration_angles.txt file, would you like to enter them? ', ...
            'enter penetration angles?', ...
            'Sure', 'Not now', 'Not now');
        switch ButtonName
            case 'Sure'
                prompt={'Enter the penetration angle for shank 1',...
                    'Enter the penetration angle for shank 2';...
                    };
                name=sprintf('Enter the penetration angles for %s', BonsaiFolder);
                numlines=[1 120; 1 120];
                 opts.Resize='on';

                penetration_angles_str=inputdlg(prompt,name,numlines,{'',''}, opts);

                if isempty(penetration_angles_str)   warning('user cancelled, no penetration_angles.txt file written'), end
                penetration_angles.shank1=str2num(penetration_angles_str{1});
                penetration_angles.shank2=str2num(penetration_angles_str{2});
                if penetration_angles.shank1<0 | penetration_angles.shank1>90 | penetration_angles.shank2<0 | penetration_angles.shank2>90 
                    warning('invalid angles entered, must be 0-90, no penetration_angles.txt file written')
                    clear penetration_angles
                else

                    %write out penetration_angles.txt file
                    fid=fopen('penetration_angles.txt', 'w');
                    fprintf(fid, 'angles (in °) of penetration with respect to orthogonal to cortical surface');
                    fprintf(fid, '\nshank1: %d', penetration_angles.shank1);
                    fprintf(fid, '\nshank2: %d', penetration_angles.shank2);
                    fclose(fid);
                end
            case 'Not now'
                %do nothing
        end % switch
    else %if non-interactive
        %do nothing
        fprintf('\nnon-interactive mode, skipping user entry of penetration angles ')
    end
end



if exist('sink_chans', 'var')
    for j=1:64
        corrected_depth(j)=-pitch*j+pitch*sink_chans.shank1L34+400;
    end
    for j=65:128
        corrected_depth(j)=-pitch*j+pitch*sink_chans.shank2L34+400;
    end

    if exist('penetration_angles', 'var')
        for j=1:64
            angle_corrected_depth(j)=(-pitch*j+pitch*sink_chans.shank1L34)*cosd(penetration_angles.shank1)+400;
        end
        for j=65:128
            angle_corrected_depth(j)=(-pitch*j+pitch*sink_chans.shank2L34)*cosd(penetration_angles.shank2)+400;
        end
    else
        angle_corrected_depth=[];
    end

    generated_on=datestr(now);
    generated_by=mfilename;
    readme='corrected_depth = subtracted so that L3/4 sink is at 400µm, ignoring L5/6 sink. angle_corrected_depth = subtracted so that L3/4 sink is at 400µm and corrected for penetration angle, ignoring L5/6 sink';
    save depths.mat readme corrected_depth angle_corrected_depth sink_chans BonsaiFolder generated_on generated_by

   %sanity check. comment this out if it gets tedious
    fprintf('\nangle correction sanity check:\n\n')
    for j=1:128
        fprintf('\nchannel %d corrected_depth %.0f angle_corrected_depth %.0f', j, corrected_depth(j), angle_corrected_depth(j))
    end


    %replot with corrected depth labels
    figure
    subplot(121)
    %shank1
    plot(t, squeeze(CSD(1,:,:))+offsets2, 'k');
    ylim([0 1.1*max(offsets2)])
    title('CSD')
    if isempty(angle_corrected_depth)
        text(repmat(-30, 1, 62), offsets2, int2str(corrected_depth([2:63])'), 'fontsize', 6)
        h=title(sprintf('CSD shank 1 %s \nwith corrected depths',BonsaiFolder));
    else
        text(repmat(-30, 1, 62), offsets2, int2str(angle_corrected_depth([2:63])'), 'fontsize', 6)
        h=title(sprintf('CSD shank 1 %s \nwith angle-corrected depths',BonsaiFolder));
    end
    text(0*sink_chans.shank1L34-10, offsets2(sink_chans.shank1L34-1), '*', 'fontsize', 18, 'color', 'g')
    text(0*bad_channels_sh1_trim-15, offsets2(bad_channels_sh1_trim), 'X', 'fontsize', 12, 'color', 'r')
    axis off
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')

    subplot(122)
    %shank2
    plot(t, squeeze(CSD(2,:,:))+offsets2, 'k');
    ylim([0 1.1*max(offsets2)])
    if isempty(angle_corrected_depth)
        text(repmat(-30, 1, 62), offsets2, int2str(corrected_depth([66:127])'), 'fontsize', 6)
    else
        text(repmat(-30, 1, 62), offsets2, int2str(angle_corrected_depth([66:127])'), 'fontsize', 6)
    end
    text(0*sink_chans.shank2L34-10, offsets2(sink_chans.shank2L34-65), '*', 'fontsize', 18, 'color', 'g')
    text(0*bad_channels_sh2_trim-15, offsets2(bad_channels_sh2_trim-1), 'X', 'fontsize', 12, 'color', 'r')
    set(gcf, 'pos', [ 998         198         660        1100])
    axis off
    title('shank 2')


    for shank=1:2
        for j=2:chans_per_shank-1 %cannot compute CSD for first and last channels, - 2 channels
            depth(shank, j-1)=pitch*j;
            depthstr{j-1}=sprintf('%.0f',depth(shank, j-1));
        end
    end

    figure
    subplot(121)
    %shank1
    imagesc(squeeze(CSD(1,:,:)));
    set(gca, 'ydir', 'normal')
    caxis(.5*[cmin1 cmax1])
    xticks(linspace(1, size(CSD,3), numxticks))
    xticklabels(linspace(xlimits(1), xlimits(2), numxticks))
    yt=[(sink_chans.shank1L34-50-1):10:(sink_chans.shank1L34+50)];
    yt=yt(find(yt>1 & yt<64));
    yticks(yt)
    if isempty(angle_corrected_depth)
        yticklabels(round(corrected_depth(round(yticks)+1)))
        h=title(sprintf('CSD shank 1 %s \nwith corrected depths',BonsaiFolder));
    else
        yticklabels((angle_corrected_depth(round(yticks)+1)))
        h=title(sprintf('CSD shank 1 %s \nwith angle-corrected depths',BonsaiFolder));
    end
    xlabel('time, ms')
    ylabel('depth, um')
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none')
    text(0*bad_channels_sh1-20, bad_channels_sh1, 'X', 'fontsize', 12, 'color', 'r')
    text(0*sink_chans.shank1L34-10, sink_chans.shank1L34-1, '*', 'fontsize', 18, 'color', 'g')
    text(0*sink_chans.shank1L56-10, sink_chans.shank1L56-1, '*', 'fontsize', 18, 'color', 'w')


    subplot(122)
    %shank2
    imagesc(squeeze(CSD(2,:,:)));
    set(gca, 'ydir', 'normal')
    caxis(.5*[cmin2 cmax2])
    colorbar
    colormap(parula)
    xticks(linspace(1, size(CSD,3), numxticks))
    xticklabels(linspace(xlimits(1), xlimits(2), numxticks))
    yt=[(sink_chans.shank2L34-50-0):10:(sink_chans.shank2L34+50)];
    yt=yt(find(yt>65 & yt<128));
    yticks(yt-65)
    if isempty(angle_corrected_depth)
        yticklabels(round(corrected_depth(yt)))
    else
        yticklabels(round(angle_corrected_depth(yt)))
    end
    xlabel('time, ms')
    ylabel('depth, um')
    title('shank 2')
    text(0*bad_channels_sh2-20, bad_channels_sh2, 'X', 'fontsize', 12, 'color', 'r')
    text(0*sink_chans.shank2L34-10, sink_chans.shank2L34-65, '*', 'fontsize', 18, 'color', 'g')
    text(0*sink_chans.shank2L56-10, sink_chans.shank2L56-65, '*', 'fontsize', 18, 'color', 'w')
    set(gca, 'pos', [ 0.5300    0.1100    0.3347    0.8150])


end %if exist sink chans

if printtofile
    fprintf('\nprinting figs to pdf...')
    cd(Bdirs{1}) %go to Bonsai Folder
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
        fprintf('\t done')
end

