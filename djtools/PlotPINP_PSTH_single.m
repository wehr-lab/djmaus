function PlotPINP_PSTH_single(varargin)

%plots a single file of clustered spiking tuning curve data from djmaus
%this function plots silent sounds with laser on, that's it.
%
% usage: PlotPINP_PSTH_single(datapath, t_filename, [xlimits],[ylimits], [binwidth])
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
    ProcessTC_PSTH_single(datadir,  t_filename, xlimits, ylimits);
end

[p,f,ext]=fileparts(t_filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\n%s', t_filename)
fprintf('\n%s', outfilename)

cd(datadir)

if exist(outfilename,'file')
    load(outfilename)
else
    ProcessTC_PSTH_single(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        ProcessTC_PSTH_single(datadir,  t_filename, xlimits, ylimits);
        load(outfilename);
    end
end

try
    IL=out.IL; %whether there are any interleaved laser trials
    nreps_ssON=out.nreps_ssON;
    nreps_ssOFF=out.nreps_ssOFF;
    SilentSoundON=out.SilentSoundON;
    SilentSoundOFF=out.SilentSoundOFF;
    mSilentSoundON=out.mSilentSoundON;
    mSilentSoundOFF=out.mSilentSoundOFF;
    SilentSoundONspikecount=out.SilentSoundONspikecount;
    SilentSoundOFFspikecount=out.SilentSoundOFFspikecount;
    SilentSoundONStim=out.SilentSoundONStim;
    SilentSoundOFFStim=out.SilentSoundOFFStim;
    SilentSoundONLaser=out.SilentSoundONLaser;
    SilentSoundOFFLaser=out.SilentSoundOFFLaser;
    M_LaserStart=out.M_LaserStart;
    LaserStart=out.LaserStart;
catch
    fprintf('\nre-processing %s', outfilename)
    ProcessTC_PSTH_single(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
    IL=out.IL; %whether there are any interleaved laser trials
    nreps_ssON=out.nreps_ssON;
    nreps_ssOFF=out.nreps_ssOFF;
    SilentSoundON=out.SilentSoundON;
    SilentSoundOFF=out.SilentSoundOFF;
    mSilentSoundON=out.mSilentSoundON;
    mSilentSoundOFF=out.mSilentSoundOFF;
    SilentSoundONspikecount=out.SilentSoundONspikecount;
    SilentSoundOFFspikecount=out.SilentSoundOFFspikecount;
    SilentSoundONStim=out.SilentSoundONStim;
    SilentSoundOFFStim=out.SilentSoundOFFStim;
    SilentSoundONLaser=out.SilentSoundONLaser;
    SilentSoundOFFLaser=out.SilentSoundOFFLaser;
    M_LaserStart=out.M_LaserStart;
    LaserStart=out.LaserStart;
end
LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
samprate=out.samprate; %in Hz
if isempty(xlimits) xlimits=out.xlimits;end

%plot the mean tuning curve OFF
figure
fs=10;
hold on
spiketimes1=mSilentSoundOFF.spiketimes;
X=xlimits(1):binwidth:xlimits(2); %specify bin centers
[N, x]=hist(spiketimes1, X);
N=N./nreps_ssOFF; %normalize to spike rate (averaged across trials)
N=1000*N./binwidth; %normalize to spike rate in Hz
bar(x,N,1);
offset=0;
yl=ylim;
inc=(yl(2))/max(nreps_ssOFF);
if rasters==1
    for n=1:nreps_ssOFF
        spiketimes2=SilentSoundOFF(n).spiketimes;
        offset=offset+inc;
        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
    end
end

if StimRecorded
    yl(1)=yl(1)-.1*diff(yl);
    offsetS=yl(1);
    for rep=1:nreps_ssOFF
        Stimtrace=SilentSoundOFFStim(rep,:);
        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
        Stimtrace=.25*diff(yl)*Stimtrace;
        t=1:length(Stimtrace);
        t=1000*t/out.samprate; %convert to ms
        t=t+out.xlimits(1); %correct for xlim in original processing call
        plot(t, Stimtrace+offsetS, 'm')
    end
else
    line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
end
if LaserRecorded
    yl(1)=yl(1)-.1*diff(yl);
    offsetL=yl(1);
    for rep=1:nreps_ssOFF
        Lasertrace=SilentSoundOFFLaser(rep, :);
        Lasertrace=Lasertrace -mean(Lasertrace(1:100));
        Lasertrace=.05*diff(yl)*Lasertrace;
        plot( t, Lasertrace+offsetL, 'c')
    end
end
if ~isempty(ylimits) ylim(ylimits); end

xlim(xlimits)
set(gca, 'fontsize', fs)
h=title(sprintf('%s: \nSilent Sound tetrode%d cell%d, nreps: %d-%d, OFF',datadir,channel,out.cluster,min(nreps_ssON),max(nreps_ssON)));
set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')

if IL
    %plot the mean tuning curve ON
    figure
    hold on
    spiketimes1=mSilentSoundON.spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nreps_ssON; %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    bar(x,N,1);
    offset=0;
    yl=ylim;
    inc=(yl(2))/max(nreps_ssON);
    if rasters==1
        for n=1:nreps_ssON
            spiketimes2=SilentSoundON(n).spiketimes;
            offset=offset+inc;
            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
        end
    end
    
    if StimRecorded
        yl(1)=yl(1)-.1*diff(yl);
        offsetS=yl(1);
        for rep=1:nreps_ssON
            Stimtrace=SilentSoundONStim(rep,:);
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=.25*diff(yl)*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            plot(t, Stimtrace+offsetS, 'm')
        end
    else
        line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
    end
    if LaserRecorded
        yl(1)=yl(1)-.1*diff(yl);
        offsetL=yl(1);
        for rep=1:nreps_ssON
            Lasertrace=SilentSoundONLaser(rep, :);
            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
            Lasertrace=.05*diff(yl)*Lasertrace;
            plot( t, Lasertrace+offsetL, 'c')
        end
    end
    if ~isempty(ylimits) ylim(ylimits); end
    
    xlim(xlimits)
    set(gca, 'fontsize', fs)
    
    h=title(sprintf('%s: \nSilent Sound tetrode%d cell%d, nreps: %d-%d, ON',datadir,channel,out.cluster,min(nreps_ssON),max(nreps_ssON)));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
end

%How about a t-test for effect of Laser in first 25 ms

for rep=1:nreps_ssON
    stop=LaserStart+25;
    st=SilentSoundON(n).spiketimes;
    spiketimes=st(st>LaserStart & st<stop); % spiketimes in region
    ON(rep)=length(spiketimes);
end
for rep=1:nreps_ssOFF
    stop=LaserStart+25;
    st=SilentSoundOFF(n).spiketimes;
    spiketimes=st(st>LaserStart & st<stop); % spiketimes in region
    OFF(rep)=length(spiketimes);
end
[h,p]=ttest2(ON, OFF, 'tail', 'right')
yl=ylim;
text(xlimits(1)+25, .95*yl(2), sprintf('h=%d, p=%.4f effect of laser (t-test)', h,p), 'fontsize', 14)








