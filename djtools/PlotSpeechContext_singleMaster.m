function PlotSpeechContext_single(varargin)

%plots a single file of clustered spiking tuning curve data from djmaus
%
% usage: PlotSpeechContext_single(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% or     PlotSpeechContext_single(datapath, outfilename, [xlimits],[ylimits], [binwidth])
% (xlimits, ylimits, binwidth are optional)
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};
t_filename=varargin{2};

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
catch
    binwidth=5;
end

if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessSoundfile_singleMaster(datadir,  t_filename, xlimits, ylimits);
end

if ischar(varargin{2}) & strfind(varargin{2}, 'out')
    outfilename=varargin{2};
else
    %Nick addition 8/31/18 - accomodates kilosort input
    if ischar(t_filename)
        [p,f,ext]=fileparts(t_filename);
        split=strsplit(f, '_');
        ch=strsplit(split{1}, 'ch');
        channel=str2num(ch{2});
        clust=str2num(split{end});
    else %reads kilosort input, which is [clust, channel, cellnum]
        channel=t_filename(1,2);
        clust=t_filename(1,1);
    end
    %end of Nick addition 8/31/18.
    
    outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
    fprintf('\nchannel %d, cluster %d', channel, clust)
end

fprintf('\n%s', t_filename)
fprintf('\n%s', outfilename)

cd(datadir)

fprintf('\nloading outfile')
if exist(outfilename,'file')
    d=dir(outfilename); fprintf(' (%.1f MB) ...', d.bytes/(1024^2))
    load(outfilename)
    fprintf('\nloaded outfile.')
else
    fprintf('\noutfile not found, calling ProcessSoundfile')
    ProcessSoundfile_singleMaster(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%d %d] but xlimits in outfile are [%d %d], re-processing...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        ProcessSoundfile_singleMaster(datadir,  t_filename, xlimits, ylimits);
        load(outfilename);
    end
end

IL=out.IL; %whether there are any interleaved laser trials
sourcefiles=out.sourcefiles;
numsourcefiles=out.numsourcefiles;
amps=out.amps;
durs=out.durs;
nreps=out.nreps;
numamps=out.numamps;
numdurs=out.numdurs;
samprate=out.samprate; %in Hz
if isfield(out, 'mM1ON')
    mM1ON=out.mM1ON;
end
if isfield(out, 'M1ON')
    M1ON=out.M1ON;
end
mM1OFF=out.mM1OFF;
M1OFF=out.M1OFF;
nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;
if isempty(xlimits) xlimits=out.xlimits;end

LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
try
    M_LaserStart=out.M_LaserStart;
    M_LaserWidth=out.M_LaserWidth;
    M_LaserNumPulses=out.M_LaserNumPulses;
    M_LaserISI=out.M_LaserISI;
    LaserStart=out.LaserStart;
    LaserWidth=out.LaserWidth;
    LaserNumPulses=out.LaserNumPulses;
    LaserISI=out.LaserISI;
    LaserStart=LaserStart(~isnan(LaserStart));
    LaserWidth=LaserWidth(~isnan(LaserWidth));
end
M1ONStim=out.M1ONStim;
M1ONLaser=out.M1ONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
% mM1ONStim=out.mM1ONStim;
% mM1ONLaser=out.mM1ONLaser;
M1OFFStim=out.M1OFFStim;
M1OFFLaser=out.M1OFFLaser;
mM1OFFStim=out.mM1OFFStim;
mM1OFFLaser=out.mM1OFFLaser;
fs=10; %fontsize

% %find optimal axis limits
if isempty(ylimits)
    ymax=0;
    for aindex=[numamps:-1:1]
        for sourcefileindex=1:numsourcefiles
            for dindex=1:numdurs
                st=mM1OFF(sourcefileindex, aindex, dindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(st, X);
                N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));
                
                if IL
                    st=mM1ON(sourcefileindex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(st, X);
                    N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));
                end
            end
        end
    end
    ylimits=[-.3 ymax];
end

aindex=2; %hard-coding to avoid aindex1 (which is silent sound amplitude=-1000)
dindex=2;

%plot the mean tuning curve OFF
% figure('position',[200 100 600 600])
p=0;
% subplot1(numsourcefiles,1, 'Max', [.95 .9])
for sourcefileindex=1:numsourcefiles
    p=p+1;
    figure
    % subplot1(p)
    hold on
    spiketimes1=mM1OFF(sourcefileindex, aindex, dindex).spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    offset=0;
    yl=ylimits;
    xlabel('time, ms')
    inc=(yl(2))/max(max(max(nreps)));
    if rasters==1
        for n=1:nrepsOFF(sourcefileindex, aindex, dindex)
            spiketimes2=M1OFF(sourcefileindex, aindex, dindex, n).spiketimes;
            offset=offset+inc;
            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
        end
    end
    b=bar(x,N,1);
    set(b, 'facecolor', 'k');
    
    offsetS=ylimits(1)+.05*diff(ylimits);
    if StimRecorded
        Stimtrace=squeeze(mM1OFFStim(sourcefileindex, aindex, dindex, :));
        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
        Stimtrace=1.25*diff(ylimits)*Stimtrace;
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
    %             if LaserRecorded
    %                 for rep=1:nrepsOFF(findex, aindex, dindex)
    %                     Lasertrace=squeeze(M1OFFLaser(findex, aindex, dindex,rep, :));
    %                     Lasertrace=Lasertrace -mean(Lasertrace(1:100));
    %                     Lasertrace=.05*diff(ylimits)*Lasertrace;
    %                     plot( t, Lasertrace+offsetS, 'c')
    %                 end
    %             end
    %                 line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
    %                 line(xlimits, [0 0], 'color', 'k')
    ylimits2(2)=ylimits(2)+offset;
    ylimits2(2)=2*ylimits(2);
    %             ylim([-2 1.1*(yl(2)+offset)])
    try
        ylim(ylimits2)
    catch
        ylim auto
    end
    
    xlim(xlimits)
    set(gca, 'fontsize', fs)
    %set(gca, 'xticklabel', '')
    %set(gca, 'yticklabel', '')
    if p==1
        h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, OFF',datadir,out.channel,out.cluster,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
        set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    end
    vpos=mean(ylimits);
    text(xlimits(1), vpos, sprintf('%s', sourcefiles{sourcefileindex}), 'interpreter', 'none')
    
end

% %label amps and freqs
% p=0;
% for sourcefileindex=1:numsourcefiles
%     p=p+1;
%     %     subplot1(p)
%     vpos=mean(ylimits);
%     text(xlimits(1), vpos, sprintf('%s', sourcefiles{sourcefileindex}), 'interpreter', 'none')
% end

if IL
    %plot the mean tuning curve ON
    figure('position',[1250 100 600 700])
    p=0;
    subplot1(numsourcefiles,1, 'Max', [.95 .9])
    for sourcefileindex=1:numsourcefiles
        p=p+1;
        subplot1(p)
        hold on
        spiketimes1=mM1ON(sourcefileindex, aindex, dindex).spiketimes;
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(spiketimes1, X);
        N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        offset=0;
        yl=ylimits;
        inc=(yl(2))/max(max(max(nreps)));
        if rasters==1
            for n=1:nrepsON(sourcefileindex, aindex, dindex)
                spiketimes2=M1ON(sourcefileindex, aindex, dindex, n).spiketimes;
                offset=offset+inc;
                %this should plot a cyan line for every trial among the rasters
                %it should accomodate trial-by-trial changes to
                %Laser params
                MLaserStart=M_LaserStart(sourcefileindex,aindex,dindex, n);
                MLaserWidth=M_LaserWidth(sourcefileindex,aindex,dindex, n);
                MLaserNumPulses=M_LaserNumPulses(sourcefileindex,aindex,dindex, n);
                MLaserISI=M_LaserISI(sourcefileindex,aindex,dindex, n);
                for np=1:MLaserNumPulses
                    plot([MLaserStart+(np-1)*(MLaserWidth+MLaserISI) MLaserStart+(np-1)*(MLaserWidth+MLaserISI)+MLaserWidth], [1 1]+yl(2)+offset, 'c')
                end
                h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
            end
        end
        b=bar(x, N,1);
        set(b, 'facecolor', 'k');
        line(xlimits, [0 0], 'color', 'k')
        
        if StimRecorded
            Stimtrace=squeeze(mM1OFFStim(sourcefileindex, aindex, dindex, :));
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
        
        ylimits2(2)=ylimits(2)+offset;
        ylimits2(2)=2*ylimits(2);
        try
            ylim(ylimits2)
        catch
            ylim auto
        end
        %                 ylim([-2 1.1*(yl(2)+offset)])
        
        if 0*LaserRecorded
            for rep=1:nrepsON(sourcefileindex, aindex, dindex)
                Lasertrace=squeeze(M1ONLaser(sourcefileindex, aindex, dindex,rep, :));
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
        %this should plot a cyan line at the unique Laser
        %params - not sure what will happen if not scalar
        %                 for np=1:LaserNumPulses
        %                     plot([LaserStart+(np-1)*(LaserWidth+LaserISI) LaserStart+(np-1)*(LaserWidth+LaserISI)+LaserWidth], [-2 -2], 'c', 'linewidth', 2)
        %                 end
        %                 line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
        
        xlim(xlimits)
        set(gca, 'fontsize', fs)
        % set(gca, 'xticklabel', '')
        % set(gca, 'yticklabel', '')
        
    end
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, ON',datadir,out.channel,out.cluster,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %label amps and freqs
    p=0;
    for sourcefileindex=1:numsourcefiles
        p=p+1;
        subplot1(p)
        vpos=mean(ylimits);
        text(xlimits(1), vpos, sprintf('%s', sourcefiles{sourcefileindex}), 'interpreter', 'none')
    end
    
end %if IL


fprintf('\nlist of sourcefiles:')
for sourcefileindex=1:numsourcefiles
    fprintf('\n%s',    sourcefiles{sourcefileindex})
end
% 1:10: ba-da
% 11:20: iba-ida
% 21:30: uba-uda
% note that 10 is out of order
peakidx=70:90;%for one cell, bin 76 is a response peak
order=[1 3 4 5 6 7 8 9 10 2];
cmap=jet(10);
%plot overlayed curves for ba-da continua
figure
c=0;
for sourcefileindex=order
    c=c+1;
    fprintf('\n%s',    sourcefiles{sourcefileindex})
    hold on
    spiketimes1=mM1OFF(sourcefileindex, aindex, dindex).spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    b=plot(x,N);
    bada(c)=max(N(peakidx));
    set(b, 'color', cmap(c,:), 'linewid', 2);
end
legend ('ba', '', '','','','','','','','da')
title('ba-da')

figure
c=0;
for sourcefileindex=order+10
    c=c+1;
    hold on
    spiketimes1=mM1OFF(sourcefileindex, aindex, dindex).spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    b=plot(x,N);
    ibaida(c)=max(N(peakidx));
    set(b, 'color', cmap(c,:), 'linewid', 2);
end
legend ('ba', '', '','','','','','','','da')
title('iba-ida')

figure
c=0;
for sourcefileindex=order+20
    c=c+1;
    hold on
    spiketimes1=mM1OFF(sourcefileindex, aindex, dindex).spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nreps(sourcefileindex, aindex, dindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    b=plot(x,N);
    ubauda(c)=max(N(peakidx));
    set(b, 'color', cmap(c,:), 'linewid', 2);
end
xlabel('time, ms')
title('uba-uda')
legend ('ba', '', '','','','','','','','da')
figure
hold on
x=1:10;
p=plot(x, bada, x, ibaida, x, ubauda);
set(p, 'linewidth', 2)
legend('ba-da', 'iba-ida', 'uba-uda')
ylabel('firing rate')
xlabel('stimulus (ba-da continuum)')

% keyboard




