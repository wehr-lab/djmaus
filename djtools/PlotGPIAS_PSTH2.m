function PlotGPIAS_PSTH2(varargin)
%PlotGPIAS_PSTH2 - plot PSTHs and rasters of clustered spiking GPIAS (gap-in-noise) data,
%new OpenEphys/kilosort file formats and hierarchy. Processes the data (ProcessSession,
%ProcessGPIAS_PSTH2) if outPSTH.mat is not found.
%
%usage:
%   PlotGPIAS_PSTH2()                                  %datapath = pwd, all cells
%   PlotGPIAS_PSTH2(datapath)                          %Bonsai or Ephys folder
%   PlotGPIAS_PSTH2(datapath, 'cells', 1:20, 'xlimits', [-300 500], 'binwidth', 10)
%name-value options (defaults):
%   'cells'          []      cell indices into SortedUnits (default all)
%   'xlimits'        []      ms re gap onset (default: from outfile; reprocesses if wider)
%   'ylimits'        []      firing-rate axis (default autoscaled)
%   'binwidth'       20      PSTH bin width, ms
%   'rasters'        true    overlay spike rasters
%   'pdfMode'        'jpeg'  'jpeg': one JPEG per page via Booklet (small, quick to browse over the
%                            network); 'vector': the old exportgraphics vector booklet (very large)
%   'pdfResolution'  150     dpi for 'jpeg' pages
%   'pdfQuality'     85      JPEG quality for 'jpeg' pages
%   'forceReprocess' false   re-run ProcessGPIAS_PSTH2 even if outPSTH.mat exists
%   'printToFile'    true    append every figure to <BonsaiFolder>-figs.pdf in the Bonsai folder
%   'showWindows'    false   show figures on screen (false: invisible, pdf output only)
%   'openPdf'        true    open the pdf booklet when done (if printToFile)
%   'maxWindows'     20      ask before continuing beyond this many open figures
%   'fontSize'       12
%   'windowPos'      figure position (default depends on platform)
%The old positional form PlotGPIAS_PSTH2(datapath, cells, xlimits, ylimits, binwidth)
%still works (with a warning).
%-mike/claude 10.2.26: inputParser name-value interface

%legacy positional call: PlotGPIAS_PSTH2(datapath, cells, xlimits, ylimits, binwidth)
if numel(varargin)>=2 && ~(ischar(varargin{2}) || isstring(varargin{2}))
    warning('PlotGPIAS_PSTH2: positional arguments are deprecated -- use name-value pairs, e.g. PlotGPIAS_PSTH2(path, ''cells'', 1:20)');
    legacy={'cells','xlimits','ylimits','binwidth'};
    nv={};
    for k=2:min(numel(varargin), 5)
        if ~isempty(varargin{k}), nv=[nv, legacy(k-1), varargin(k)]; end %#ok<AGROW>
    end
    varargin=[varargin(1), nv];
end

if ismac, defaultpos=[-1415 479 1381 680]; else, defaultpos=[200 100 1381 680]; end %mike's pref on mac
p=inputParser;
p.FunctionName=mfilename;
addOptional(p, 'datapath', pwd, @(x) ischar(x) || isstring(x));
addParameter(p, 'cells', [], @(x) isempty(x) || (isnumeric(x) && all(x>0)));
addParameter(p, 'xlimits', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'ylimits', [], @(x) isempty(x) || (isnumeric(x) && numel(x)==2));
addParameter(p, 'binwidth', 20, @(x) isnumeric(x) && isscalar(x) && x>0);
addParameter(p, 'rasters', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'pdfMode', 'jpeg', @(x) any(strcmpi(x, {'jpeg','vector'})));
addParameter(p, 'pdfResolution', 150, @isscalar);
addParameter(p, 'pdfQuality', 85, @isscalar);
addParameter(p, 'forceReprocess', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'printToFile', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'showWindows', false, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'openPdf', true, @(x) islogical(x) || isnumeric(x));
addParameter(p, 'maxWindows', 20, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'fontSize', 12, @(x) isnumeric(x) && isscalar(x));
addParameter(p, 'windowPos', defaultpos, @(x) isnumeric(x) && numel(x)==4);
parse(p, varargin{:});
P=p.Results;

datadir=char(P.datapath); if isempty(datadir), datadir=pwd; end
cells=P.cells(:)';
xlimits=P.xlimits;
ylimits=P.ylimits;
binwidth=P.binwidth;
rasters=logical(P.rasters);
force_reprocess=logical(P.forceReprocess);
printtofile=logical(P.printToFile);
showwindows=logical(P.showWindows);
openpdf=logical(P.openPdf) && logical(P.printToFile);
maxwindows=P.maxWindows;
fs=P.fontSize;
windowpos=P.windowPos;
if showwindows || ~printtofile, figvis='on'; else, figvis='off'; end   %figures stay visible if not printing
already_printed=gobjects(0);   %figures already appended to the pdf (when showwindows keeps them open)

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
    sess_idx=find(strcmpi(Bdirs, datadir)); %strcmpi: macOS paths are case-insensitive, and Bdirs may differ in case from datadir (e.g. 5xfad/rig3phys vs 5XFAD/Rig3Phys)
    cd(dirs{sess_idx})
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
    pdffilename=fullfile(macifypath(out.BonsaiPath), sprintf('%s-figs.pdf', out.BonsaiFolder));
    if isfile(pdffilename), delete(pdffilename), end
    if strcmpi(P.pdfMode, 'jpeg'), Booklet('start', pdffilename); end
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
    %chan = depth-ordered index 1..128 (1-64 shank 1 from the tip up, 65-128 shank 2).
    %raw depth = distance below the top recording site, same convention as PlotTC_PSTH2
    %(fixed mw 12.17.24). Previously distance*chan / distance*chan-chans_per_shank. -mike/claude 10.2.26
    if chan<=chans_per_shank
        shank=1;
        depth=distance*(chans_per_shank-chan);
    else
        shank=2;
        depth=distance*(2*chans_per_shank-chan);
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


    figure('position',windowpos, 'PaperOrientation', 'landscape', 'Visible',figvis)
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
            if any(f(idx)==already_printed), continue, end   %shown figures from earlier cells: already in the pdf
            set(f(idx), 'Visible', 'off')   %export from an invisible figure (avoids the macOS HiDPI cropping bug)
            if strcmpi(P.pdfMode, 'jpeg')   %10.8.26: JPEG pages via Booklet instead of a vector booklet
                Booklet('add', pdffilename, f(idx), 'Resolution', P.pdfResolution, 'Quality', P.pdfQuality);
            else
                exportgraphics(f(idx), pdffilename, 'Append', true)
            end
            if showwindows
                set(f(idx), 'Visible', 'on'); already_printed(end+1)=f(idx);
            else
                close(f(idx))
            end
        end
    end
    fprintf('\telapsed time %.0fs', toc)

end %for cellnum

if printtofile && strcmpi(P.pdfMode, 'jpeg')
    np=Booklet('finish', pdffilename); d=dir(pdffilename);
    if ~isempty(d), fprintf('\nwrote %s (%d pages, %.1f MB)', pdffilename, np, d.bytes/1e6); end
end
if openpdf
    if ismac
        system(sprintf('open "%s"', pdffilename));
    elseif ispc
        winopen(char(pdffilename));
    end
end

