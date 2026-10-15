function PlotTC_PSTH2(varargin)
%PlotTC_PSTH2 - plot PSTHs and rasters of clustered spiking tuning-curve data,
%new OpenEphys/kilosort file formats and hierarchy. Processes the data (ProcessSession and
%the matching Process function) if outPSTH.mat is not found.
%
%usage:
%   PlotTC_PSTH2()                                  %datapath = pwd, all cells
%   PlotTC_PSTH2(datapath)                          %Bonsai or Ephys folder
%   PlotTC_PSTH2(datapath, 'cells', 1:20, 'xlimits', [-100 300], 'binwidth', 10)
%name-value options (defaults):
%   'cells'          []      cell indices into SortedUnits (default all)
%   'xlimits'        []      ms (default: from outfile; reprocesses if wider)
%   'ylimits'        []      firing-rate axis (default autoscaled)
%   'binwidth'       5       PSTH bin width, ms
%   'rasters'        true    overlay spike rasters
%   'pdfMode'        'jpeg'  'jpeg': one JPEG per page via Booklet (~0.2-0.4 MB/page); 'vector': the old
%                            exportgraphics booklet (~20 MB/page with rasters; 'image' was ignored)
%   'pdfResolution'  150     dpi for 'jpeg' pages
%   'pdfQuality'     85      JPEG quality for 'jpeg' pages
%   'forceReprocess' false   re-process even if outPSTH.mat exists
%   'printToFile'    true    append every figure to <BonsaiFolder>-figs.pdf in the Bonsai folder
%   'showWindows'    false   show figures on screen (false: invisible, pdf output only)
%   'openPdf'        true    open the pdf booklet when done (if printToFile)
%   'maxWindows'     20      ask before continuing beyond this many open figures
%   'fontSize'       12
%   'windowPos'      figure position (default depends on platform)
%The old positional form PlotTC_PSTH2(datapath, cells, xlimits, ylimits, binwidth)
%still works (with a warning).
%-mike/claude 10.2.26: inputParser name-value interface

%legacy positional call: PlotTC_PSTH2(datapath, cells, xlimits, ylimits, binwidth)
if numel(varargin)>=2 && ~(ischar(varargin{2}) || isstring(varargin{2}))
    warning('PlotTC_PSTH2: positional arguments are deprecated -- use name-value pairs, e.g. PlotTC_PSTH2(path, ''cells'', 1:20)');
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
addParameter(p, 'binwidth', 5, @(x) isnumeric(x) && isscalar(x) && x>0);
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

if force_reprocess
    fprintf('\nForce re-process\n')
    warning('Force re-process')
    if exist('./notebook.mat')==2
        %we're in Ephys folder, so go up one
        cd ..
    end
    ProcessSession
    if ~exist('SortedUnitsFile') fprintf('\nNo Sorted Units File, probably because there is no kilosort data. \nPlotTC_PSTH2 will fail.'); end
    load(SortedUnitsFile)
    ProcessTC_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits);
    load(outfilename);
end

if exist(outfilename,'file')
    load(outfilename)
    fprintf('\nloaded outfile.')
elseif exist('Bdirs.mat', 'file') %we're in BonsaiDir, look for outfile in EphysDir
    load('Bdirs.mat')
    sess_idx=find(strcmpi(Bdirs, datadir)); %strcmpi: macOS paths are case-insensitive, and Bdirs may differ in case from datadir (e.g. 5xfad/rig3phys vs 5XFAD/Rig3Phys)
    cd(dirs{sess_idx})
    try
        load(outfilename)
        fprintf('\nloaded outfile.')
    catch
        cd .. %cd(Bdirs{1})
        ProcessSession
        if ~exist('SortedUnitsFile') fprintf('\nNo Sorted Units File, probably because there is no kilosort data. PlotTC_PSTH2 will fail.'); end
        load(SortedUnitsFile)
        ProcessTC_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits);
        load(outfilename);
    end

else
    fprintf('\ncould not find outfile, calling ProcessSession...')
    if exist('./notebook.mat')==2
        %we're in Ephys folder, so go up one
        cd ..
    end
    ProcessSession
    if ~exist('SortedUnitsFile') fprintf('\nNo Sorted Units File, probably because there is no kilosort data. PlotTC_PSTH2 will fail.'); end
    load(SortedUnitsFile)
    ProcessTC_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, xlimits,ylimits);
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
    %delete figs.ps
    %delete figs.pdf
    pdffilename=fullfile(macifypath(out.BonsaiPath), sprintf('%s-figs.pdf', out.BonsaiFolder)); %Bonsai folder, absolute path
    if isfile(pdffilename), delete(pdffilename), end
    if strcmpi(P.pdfMode, 'jpeg'), Booklet('start', pdffilename); end
    close all
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
fprintf('\nusing xlimits [%d %d]', xlimits)

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

%hardcoding for probe P128-2
pitch=20; %distance between recording sites in µm
chans_per_shank=64;
% look for depths.mat in local directory or master directory
try
    corrected_depths_from_file=load('depths.mat');
catch
    %if this was batch kilosorted, there might be a depths.mat in the
    %LFP directory which is dirs{1}
    if exist('Bdirs.mat', 'file') load('Bdirs.mat')
    elseif exist('dirs.mat', 'file') load('dirs.mat')
    end
    try
        corrected_depths_from_file=load(fullfile(dirs{1}, 'depths.mat'));
    end
end
try
    corrected_depths=corrected_depths_from_file.corrected_depth;
    angle_corrected_depths=corrected_depths_from_file.angle_corrected_depth;
    fprintf('\nfound and loaded corrected depths file')
catch
    corrected_depths=[];
    angle_corrected_depths=[];
    warning('\ncould not find corrected depths file')
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
        depth=pitch*(chans_per_shank-chan);
    else
        shank=2;
        depth=pitch*(2*chans_per_shank-chan); %fixed mw 12.17.24
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
        figure('position',windowpos, 'PaperOrientation', 'landscape', 'Visible',figvis)
        
        p=0;
        subplot1(numy,numx, 'Max', [.95 .9])
        set(gcf, 'visible', 'off')
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
                [N, x]=histcounts(spiketimes1, X);
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
                b=bar(X(1:end-1),N,1);
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
                    h=title(sprintf('%s: \ncell %d, chan %d, shank %d, raw depth %d um, corrected depth %.0f um, angle-corrected depth %.0f um, %d spikes, %dms, nreps: %d-%d, %s, OFF ',datadir,cellnum, chan, shank, depth, corrected_depth, angle_corrected_depth, length(out.SortedUnits(cellnum).spiketimes), durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF))), out.SortedUnits(cellnum).rating));
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
                    [N, x]=histcounts(spiketimes1, X);
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
                    b=bar(X(1:end-1),N,1);
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
            figure('position',windowpos, 'PaperOrientation', 'landscape', 'Visible',figvis)

            p=0;
            subplot1(numy,numx, 'Max', [.95 .9])
            set(gcf, 'visible', 'off')
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
                    [N, x]=histcounts(spiketimes1, X);
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
                    b=bar(X(1:end-1), N,1);
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
                        for rep=1:nrepsON(cellnum, findex, aindex, dindex)
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
                    h=title(sprintf('%s: \ncell %d, chan %d, shank %d, raw depth %d um, corrected depth %.0f um, angle_corrected depth %.0f um, %d spikes, %dms, nreps: %d-%d, %s, ON ',out.BonsaiFolder,cellnum, chan, shank, depth, corrected_depth, angle_corrected_depth, length(out.SortedUnits(cellnum).spiketimes), durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON))), SortedUnits(cellnum).rating));
                    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

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
                        % [N, x]=hist(spiketimes1, X);
                        [N, x]=histcounts(spiketimes1, X);
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
                        b=bar(X(1:end-1), N,1);
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

        %plot OFF & ON overlayed, curves w/ no rasters
        for windex=1:numw
            figure('position',windowpos, 'PaperOrientation', 'landscape', 'Visible',figvis)

            p=0;
            subplot1(numy,numx, 'Max', [.95 .9])
            set(gcf, 'visible', 'off')
            for yindex=ystart:ystep:yend
                for findex=1:numfreqs
                    p=p+1;
                    subplot1(p)
                    hold on
                    if numamps>=numdurs, aindex=yindex; dindex=windex;
                    else aindex=windex; dindex=yindex;
                    end
                    spiketimesON=mM1ON(cellnum, findex, aindex, dindex).spiketimes;
                    spiketimesOFF=mM1OFF(cellnum, findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [Non, xon]=histcounts(spiketimesON, X);
                    [Noff, xoff]=histcounts(spiketimesOFF, X);
                    Non=Non./nrepsON(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    Non=1000*Non./binwidth; %normalize to spike rate in Hz
                    Noff=Noff./nrepsOFF(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    Noff=1000*Noff./binwidth; %normalize to spike rate in Hz
                    offset=ylimits(2);
                    yl=ylimits;
                    inc=(yl(2))/max(nrepsON(:));

                    Hon=plot(X(1:end-1), Non);
                    set(Hon, 'color', 'c', 'linew', 2);
                    hold on
                    Hoff=plot(X(1:end-1), Noff);
                    set(Hoff, 'color', 'k', 'linew', 2);
                    line(xlimits, [0 0], 'color', 'k')
                    ylimits2(2)=ylimits(2)+offset;
                    if StimRecorded & unique(Stimtrace)
                        Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
                        Stimtrace=squeeze(mM1ONStim(findex, aindex, dindex, :));
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

                    X=[xlimits(1), LaserStart, LaserStart, LaserStart+LaserWidth, LaserStart+LaserWidth, xlimits(2)];
                    height=diff(ylimits)/10;
                    Y=[0 0 height height 0 0];
                    line(X,Y, 'color', 'b', 'linewidth', 2)




                    xlim(xlimits)
                    %xlim([-50 100])
                    set(gca, 'fontsize', fs)
                    % set(gca, 'xticklabel', '')
                    % set(gca, 'yticklabel', '')

                    subplot1(1)
                    h=title(sprintf('%s: \ncell %d, chan %d, shank %d, raw depth %d, corrected depth %d %d spikes, %dms, nreps: %d-%d, %s, ON & OFF ',out.BonsaiFolder,cellnum, chan, shank, depth, corrected_depth, length(out.SortedUnits(cellnum).spiketimes), durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF))), SortedUnits(cellnum).rating));
                    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

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

                        b=plot(x,N);
                        set(b, 'color', 'c');
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

        %plot OFF & ON overlayed, curves WITH rasters
        for windex=1:numw
            figure('position',windowpos, 'PaperOrientation', 'landscape', 'Visible',figvis)

            p=0;
            subplot1(numy,numx, 'Max', [.95 .9])
            set(gcf, 'visible', 'off')
            for yindex=ystart:ystep:yend
                for findex=1:numfreqs
                    p=p+1;
                    subplot1(p)
                    hold on
                    if numamps>=numdurs, aindex=yindex; dindex=windex;
                    else aindex=windex; dindex=yindex;
                    end
                    spiketimesON=mM1ON(cellnum, findex, aindex, dindex).spiketimes;
                    spiketimesOFF=mM1OFF(cellnum, findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [Non, xon]=histcounts(spiketimesON, X);
                    [Noff, xoff]=histcounts(spiketimesOFF, X);
                    Non=Non./nrepsON(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    Non=1000*Non./binwidth; %normalize to spike rate in Hz
                    Noff=Noff./nrepsOFF(cellnum, findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    Noff=1000*Noff./binwidth; %normalize to spike rate in Hz
                    offset=ylimits(2);
                    yl=ylimits;
                    inc=(yl(2))/max(nrepsON(:));

                    Hon=plot(X(1:end-1), Non);
                    set(Hon, 'color', 'c', 'linew', 2);
                    hold on
                    Hoff=plot(X(1:end-1), Noff);
                    set(Hoff, 'color', 'k', 'linew', 2);
                    line(xlimits, [0 0], 'color', 'k')
                    ylimits2(2)=ylimits(2)+offset;

                    f=1.5; %raster offset fudge factor
                    if rasters==1
                        for n=1:nrepsON(cellnum, findex, aindex, dindex)
                            spiketimesON2=M1ON(cellnum, findex, aindex, dindex, n).spiketimes;
                            offset=offset+inc;
                            h=plot(spiketimesON2, f*yl(2)+ones(size(spiketimesON2))+offset, '.c');
                        end
                        for n=1:nrepsOFF(cellnum, findex, aindex, dindex)
                            spiketimesOFF2=M1OFF(cellnum, findex, aindex, dindex, n).spiketimes;
                            offset=offset+inc;
                            h=plot(spiketimesOFF2, f*yl(2)+ones(size(spiketimesOFF2))+offset, '.k');
                        end
                    end

                    if StimRecorded & unique(Stimtrace)
                        Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
                        Stimtrace=squeeze(mM1ONStim(findex, aindex, dindex, :));
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

                    ylimits2(2)=ylimits(2)+2+ offset;

                    %ylimits2(2)=2*ylimits(2);
                    try
                        ylimits(2) = max([ylimits(2) 1]);
                        ylim(ylimits2)
                    end

                    X=[xlimits(1), LaserStart, LaserStart, LaserStart+LaserWidth, LaserStart+LaserWidth, xlimits(2)];
                    height=diff(ylimits)/10;
                    Y=[0 0 height height 0 0];
                    line(X,Y, 'color', 'c', 'linewidth', 2)




                    xlim(xlimits)
                    %xlim([-50 100])
                    set(gca, 'fontsize', fs)
                    % set(gca, 'xticklabel', '')
                    % set(gca, 'yticklabel', '')

                    subplot1(1)
                    h=title(sprintf('%s: \ncell %d, chan %d, shank %d, raw depth %d, corrected depth %.0f, angle-corrected depth %.0f, %d spikes, %dms, nreps: %d-%d, %s, ON & OFF ',out.BonsaiFolder,cellnum, chan, shank, depth, corrected_depth, angle_corrected_depth, length(out.SortedUnits(cellnum).spiketimes), durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF))), SortedUnits(cellnum).rating));
                    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

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
                end
            end
        end


    end

    if printtofile
        %print figures to postscript file
        f=findobj('type', 'figure');
        %drawnow;
        for idx=1:length(f)
            %figure(f(idx))
            % orient landscape
            % % print figs -dpsc2 -append -bestfit
            if ~any(f(idx)==already_printed) %if printtofile=1 but showwindows=1, this keeps track so as to not re-print
               set(f(idx), 'visible', 'off')
                if strcmpi(P.pdfMode, 'jpeg')   %10.8.26: 'image' below was silently ignored -> 400 MB vector booklets
                    Booklet('add', pdffilename, f(idx), 'Resolution', P.pdfResolution, 'Quality', P.pdfQuality);
                else
                    exportgraphics(f(idx),pdffilename,'Append',true, 'ContentType', 'image', 'Resolution', 200)
                end
                % exportgraphics(f(idx),pdffilename, 'ContentType', 'vector', 'Append',true) 
                % 'vector' works, the details are more clear, but takes longer and file
                % size is more than 150x bigger, and you'll need to reduce
                %  the fontsize. 
                % 'image' seems to have crop issues on a mac 

                if showwindows
                    already_printed(end+1)=f(idx);
                    set(f(idx), 'visible', 'on')
                else
                    close(f(idx))   %was close (current figure), which need not be f(idx)
                end
            end
        end
    end
    fprintf('\telapsed time %.0fs', toc)

end
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

