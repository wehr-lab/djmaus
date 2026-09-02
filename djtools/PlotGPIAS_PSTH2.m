function PlotGPIAS_PSTH2(varargin)

%plots clustered spiking GPIAS data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%
% usage: PlotGPIAS_PSTH2([datapath], [cells], [xlimits],[ylimits], [binwidth])
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
    if isempty(binwidth), binwidth=20;end

catch
    binwidth=20;
end

outfilename='outPSTH.mat';
cd(datadir)

%single decision point for (re)processing -- avoids ever calling
%ProcessGPIAS_PSTH2 more than once per invocation, and avoids a
%redundant save-then-reload round trip through outfilename when we do
%process (ProcessGPIAS_PSTH2 now returns out directly instead of only
%saving it). Directory resolution (BonsaiDir vs EphysDir) is unchanged
%from before. -mike/claude 9.2.26
if force_reprocess
    need_process=true;
elseif exist(outfilename,'file')
    need_process=false;
elseif exist('Bdirs.mat', 'file') %we're in BonsaiDir, look for outfile in EphysDir
    load('Bdirs.mat')
    cd(dirs{1})
    need_process=~exist(outfilename,'file');
    if need_process
        cd .. %cd(Bdirs{1}) -- back to BonsaiDir before (re)processing, below
    end
else
    need_process=true;
end

if need_process
    if force_reprocess
        fprintf('\nForce re-process\n')
        warning('Force re-process')
    else
        fprintf('\ncould not find outfile, calling ProcessSession...')
    end
    if exist('./notebook.mat')==2
        %we're in Ephys folder, so go up one
        cd ..
    end
    ProcessSession
    if ~exist('SortedUnitsFile') fprintf('\nNo Sorted Units File, probably because there is no kilosort data. \nPlotGPIAS_PSTH2 will fail.'); end
    load(SortedUnitsFile)
    out=ProcessGPIAS_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits);
else
    load(outfilename)
    fprintf('\nloaded outfile.')
end

if isempty(ylimits)
    autoscale_ylimits=1;
else
    autoscale_ylimits=0;
end

%if xlimits are requested but don't match those in outfile, reprocess once more
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%.1f %.1f] but xlimits in outfile are [%.1f %.1f], re-processing...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        out=ProcessGPIAS_PSTH2(out.SortedUnits, out.BonsaiPath, out.EphysPath, out.EphysPath_KS, xlimits, ylimits);
    end
else %xlimits is empty
    xlimits=out.xlimits;
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
    %pdf file goes in the Bonsai folder, ensured by using absolute path in filename
    pdffilename=fullfile(macifypath(BonsaiPath), sprintf('%s-figs.pdf', BonsaiFolder));
    delete(pdffilename)
    close all
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
if isfield(out, 'LaserRecorded')
    LaserRecorded=out.LaserRecorded;
    M1ONLaser=out.M1ONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
    mM1ONLaser=out.mM1ONLaser;
    M1OFFLaser=out.M1OFFLaser;
    mM1OFFLaser=out.mM1OFFLaser;
else
    LaserRecorded=0;
end
if isfield(out, 'StimRecorded')
    StimRecorded=out.StimRecorded;
    M1ONStim=out.M1ONStim;
    mM1ONStim=out.mM1ONStim;
    M1OFFStim=out.M1OFFStim;
    mM1OFFStim=out.mM1OFFStim;
else
    StimRecorded=0;
end

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
fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

%hardcoding for probe P128-2
distance=20;
chans_per_shank=64;
% look for depths.mat in local directory or master directory
try
    corrected_depths_from_file=load('depths.mat');
    fprintf('\nfound and loaded depths.mat')
catch
    %if this was batch kilosorted, there might be a depths.mat in the
    %LFP directory which is somewhere in dirs{:}
        fprintf('\nno depths.mat in this directory, searching dirs/bdirs...')

    if exist('Bdirs.mat', 'file') load('Bdirs.mat')
    elseif exist('dirs.mat', 'file') load('dirs.mat')
    end
    try
        for dirs_idx=1:length(dirs)
            if exist(fullfile(dirs{dirs_idx}, 'depths.mat'))
                corrected_depths_from_file=load(fullfile(dirs{dirs_idx}, 'depths.mat'));
                fprintf('\nfound and loaded depths.mat from directory %s', fullfile(dirs{dirs_idx}, 'depths.mat'))
            end
        end
    catch
        warning('could not find depths.mat file in this directory or in any of the directories in dirs or bdirs')
    end
end
try
    corrected_depths=corrected_depths_from_file.corrected_depth;
    angle_corrected_depths=corrected_depths_from_file.angle_corrected_depth;
    fprintf('\nfound and loaded corrected depths file')
catch
    corrected_depths=[];
    angle_corrected_depths=[];
    warning('\ncould not find corrected depths file, using raw depth')
end

for cellnum=cells
    tic
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

    chan=out.SortedUnits(cellnum).channel;
    if chan<=chans_per_shank
        shank=1;
        depth=distance*chan;
    else
        shank=2;
        depth=distance*chan-chans_per_shank;
    end
    if ~isempty(corrected_depths)
        corrected_depth=corrected_depths(chan);
    else
        corrected_depth=nan;
    end
    if ~isempty(angle_corrected_depths)
        angle_corrected_depth=angle_corrected_depths(chan);
    else
        angle_corrected_depth=nan;
    end
    fprintf('\ncell %d, chan %d, shank %d, raw depth %d, corrected depth %.0f, angle_corrected_depth %.0f', cellnum, chan, shank, depth, corrected_depth, angle_corrected_depth)

    % %find optimal axis limits
    if autoscale_ylimits
        ymax=0;
        for paindex=1:numpulseamps
            for gdindex=1:numgapdurs
                st=mM1OFF(cellnum,gdindex,paindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=histcounts(st, X);
                N=N./nrepsOFF(cellnum, gdindex,paindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));

                if IL
                    st=mM1ON(cellnum,gdindex,paindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=histcounts(st, X);
                    N=N./nrepsON(cellnum, gdindex,paindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));
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

    %plot the mean tuning curve OFF


    figure('position',windowpos, 'PaperOrientation', 'landscape')
    p=0;
    subplot1(numpulseamps, numgapdurs, 'Max', [.95 .9])
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            hold on
            spiketimes1=mM1OFF(cellnum,gdindex,paindex).spiketimes; %spiketimes are in ms relative to gap termination
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=histcounts(spiketimes1, X);
            N=N./nrepsOFF(cellnum,gdindex,paindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            offset=0;
            yl=ylimits;
            inc=(yl(2))/max(nrepsOFF(:));
            if rasters==1
                for n=1:nrepsOFF(cellnum,gdindex,paindex)
                    spiketimes2=M1OFF(cellnum,gdindex,paindex, n).spiketimes;
                    offset=offset+inc;
                    h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                end
            end
            bar(X(1:end-1), N,1,'facecolor','k','edgecolor','k');
            line(xlimits, [0 0], 'color', 'k')
            ylimits2(2)=ylimits(2)*2.2;
            ylimits2(1)=-2;
            ylim(ylimits2)
            xlim(xlimits)
            set(gca, 'fontsize', fs)

            if StimRecorded %plot actual recorded stimulus trace in magenta
                Stimtrace=squeeze(mM1OFFStim(gdindex,paindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace);
                height=.1*diff(ylimits2); %reasonable height
                %height=.5*diff(ylimits2); %magnified height for scrutinizing stimulus
                Stimtrace=height*Stimtrace;
                t=1:length(Stimtrace);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                offset=-range(Stimtrace);
                ylimits2(1)=ylimits2(1)+2*offset;
                %  offset=ylimits(1)+.1*diff(ylimits);
                plot(t, Stimtrace+offset, 'm')
                ylimits2(1)=offset;
                ylim(ylimits2)
            end

            %draw line to indicate gaps
            if gapdurs(gdindex)>0
                line([0 0],[ylim],'color','m')
                line(-[(gapdurs(gdindex)) (gapdurs(gdindex))],[ylim],'color','m')
            end
            line(xlimits, [0 0], 'color', 'k')
            ylimits2(2)=ylimits(2)*2.2;
            ylim(ylimits2)


            if LaserRecorded & IL %plot laser square pulses
                height=.05*diff(ylimits2);
                offset=-height;
                for rep=1:nrepsOFF(cellnum,gdindex,paindex)
                    Lasertrace=squeeze(M1OFFLaser(gdindex,paindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=height*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end
                ylimits2(1)=ylimits2(1)-2*range(Lasertrace);
                ylim(ylimits2)
            end
        end
    end

    ymin=inf;
    ymax=0;
    p=0;
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            yl=ylim;
            if yl(1)<ymin ymin=yl(1);end
            if yl(2)>ymax ymax=yl(2);end
        end
    end
    p=0;
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            ylim([ymin ymax])
        end
    end


    subplot1(1)
    h=title(sprintf('GPIAS OFF %s: \ncell %d, chan %d, shank %d, raw depth %d um, corrected depth %.0f um, angle-corrected depth %.0f um, %d spikes, nreps: %d-%d, OFF ',datadir,cellnum, chan, shank, depth, corrected_depth, angle_corrected_depth, length(out.SortedUnits(cellnum).spiketimes),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

    %label amps and freqs
    p=0;
    for paindex=1:numpulseamps
        for gdindex=1:numgapdurs
            p=p+1;
            subplot1(p)
            vpos=.8*ylimits(2);
            text(xlimits(1)+100, vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r', 'fontsiz', 14)
            set(gca, 'yticklabel', '');
        end
    end

    %turn on ytick for left-most plot
    subplot1(1)
    set(gca, 'yticklabelmode', 'auto');

    if IL
        %plot the mean tuning curve ON
        if show_plots == 0
            figure('Visible','off')
        else
            figure
        end
        p=0;
        subplot1(numpulseamps, numgapdurs, 'Max', [.95 .9])
        for paindex=1:numpulseamps
            for gdindex=1:numgapdurs
                p=p+1;
                subplot1(p)
                hold on
                spiketimes1=mM1ON(cellnum,gdindex,paindex).spiketimes; %spiketimes are in ms relative to gap termination
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=histcounts(spiketimes1, X);
                N=N./nrepsON(cellnum,gdindex,paindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                offset=0;
                yl=ylimits;
                inc=(yl(2))/max(nrepsON(:));
                if rasters==1
                    for n=1:nrepsON(cellnum,gdindex,paindex)
                        spiketimes2=M1ON(cellnum,gdindex,paindex, n).spiketimes;
                        offset=offset+inc;
                        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    end
                end
                bar(X(1:end-1), N,1,'facecolor','g','edgecolor','k');
                line(xlimits, [0 0], 'color', 'k')
                ylimits2(2)=ylimits(2)*3;
                ylimits2(1)=-2;
                ylim(ylimits2)

                xlim(xlimits)
                set(gca, 'fontsize', fs)

                if StimRecorded
                    Stimtrace=squeeze(mM1ONStim(gdindex,paindex, :));
                    Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                    Stimtrace=.05*diff(ylimits)*Stimtrace;
                    t=1:length(Stimtrace);
                    t=1000*t/out.samprate; %convert to ms
                    t=t+out.xlimits(1); %correct for xlim in original processing call
                    offset=ylimits(1)+.1*diff(ylimits);
                    plot(t, Stimtrace+offset, 'm')
                else
                    %do nothing
                end

                %draw line to indicate gaps
                if gapdurs(gdindex)>0
                    line([0 0],[ylim],'color','m')
                    line(-[(gapdurs(gdindex)) (gapdurs(gdindex))],[ylim],'color','m')
                end
                line(xlimits, [0 0], 'color', 'k')
                ylimits2(2)=ylimits(2)*2.2;
                ylimits2(1)=-2;
                ylim(ylimits2)


                if LaserRecorded
                    for rep=1:nrepsOFF(cellnum,gdindex,paindex)
                        Lasertrace=squeeze(M1ONLaser(gdindex,paindex,rep, :));
                        Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                        Lasertrace=.05*diff(ylimits)*Lasertrace;
                        plot( t, Lasertrace+offset, 'c')
                    end
                end

            end %for gd
        end %for pa

        subplot1(1)
        h=title(sprintf('GPIAS ON %s: \nchannel%d cell%d, nreps: %d-%d, ON',datadir,channel,out.cluster,min(nrepsON(:)),max(nrepsON(:))));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

        %label amps and freqs
        p=0;
        for paindex=1:numpulseamps
            for gdindex=1:numgapdurs
                p=p+1;
                subplot1(p)
                vpos=.8*ylimits(2);
                text(xlimits(1)+100, vpos, sprintf('%d', gapdurs(gdindex)), 'color', 'r', 'fontsiz', 14)
                set(gca, 'yticklabel', '');
            end
        end

        %turn on ytick for left-most plot
        subplot1(1)
        set(gca, 'yticklabelmode', 'auto');


    end %            %if IL plot the mean tuning curve ON


    if printtofile
        f=findobj('type', 'figure');
        for idx=1:length(f)
            exportgraphics(f(idx),pdffilename,'Append',true)
        end
        if closewindows %doing this in a separate step to ensure matlab prints all figs first
            for idx=1:length(f)
                close(f(idx))
            end
        end
    end
    fprintf('\telapsed time %.0fs', toc)

end %for cellnum








