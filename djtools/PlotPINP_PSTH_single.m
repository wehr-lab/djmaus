function PlotPINP_PSTH_single(varargin)

%plots a single file of clustered spiking tuning curve data from djmaus
%this function plots silent sounds with laser on, that's it.
%
% usage: PlotPINP_PSTH_single(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% (xlimits, ylimits, binwidth are optional)
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=1;

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
    ProcessPINP_PSTH_single(datadir,  t_filename, xlimits, ylimits);
end

[p,f,ext]=fileparts(t_filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\n%s', t_filename)
fprintf('\nlooking for %s', outfilename)

cd(datadir)

if exist(outfilename,'file')
    load(outfilename)
    fprintf('\nloaded outfile')
else
    fprintf('\ndid not find outfile')
    ProcessPINP_PSTH_single(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        ProcessPINP_PSTH_single(datadir,  t_filename, xlimits, ylimits);
        load(outfilename);
    end
end

tetrode=out.channel;
mMPulse=out.mMPulse;
mMSilentSoundOFF=out.mMSilentSoundOFF;
mMTrain=out.mMTrain;
mMPulseLasertrace=out.mMPulseLasertrace;
mMSilentSoundOFFLasertrace=out.mMSilentSoundOFFLasertrace;
mMTrainLasertrace=out.mMTrainLasertrace;
mMPulseStimtrace=out.mMPulseStimtrace;
mMSilentSoundOFFStimtrace=out.mMSilentSoundOFFStimtrace;
mMTrainStimtrace=out.mMTrainStimtrace;

silentsounddurs=out.silentsounddurs;
trainnumpulses=out.trainnumpulses;
pulsewidths=out.pulsewidths;
trainisis=out.trainisis;
trainpulsewidths=out.trainpulsewidths;
numsilentsounddurs=out.numsilentsounddurs;
numtrainnumpulses=out.numtrainnumpulses;
numpulsewidths=out.numpulsewidths;
numtrainisis=out.numtrainisis;
numtrainpulsewidths=out.numtrainpulsewidths;

nrepsOFF=out.nrepsOFF;
nrepsPulse=out.nrepsPulse;
nrepsTrain=out.nrepsTrain;

spiketimes=out.spiketimes;
numnexts=out.numnexts;
nexts=out.nexts;
next=min(nexts);
laserstarts=out.laserstarts;
numlaserstarts=out.numlaserstarts;


LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
samprate=out.samprate; %in Hz
if isempty(xlimits) xlimits=out.xlimits;end

% %find optimal axis limits
if isempty(ylimits)
    ymax=0;
    for dindex=1:numsilentsounddurs
        st=mMSilentSoundOFF(dindex).spiketimes;
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(st, X);
        N=N./nrepsOFF(dindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        ymax= max(ymax,max(N));
    end
    for pwindex=1:numpulsewidths
        st=mMPulse(pwindex).spiketimes;
        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
        [N, x]=hist(st, X);
        N=N./nrepsPulse(dindex); %normalize to spike rate (averaged across trials)
        N=1000*N./binwidth; %normalize to spike rate in Hz
        ymax= max(ymax,max(N));
    end
    for tnpindex=1:numtrainnumpulses
        for tpwindex=1:numtrainpulsewidths
            for tiindex=1:numtrainisis
                st=mMTrain(tnpindex,tpwindex,tiindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(st, X);
                N=N./nrepsOFF(dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));
            end
        end
    end
    
    
    ylimits=[-.3 ymax*2.5];
end


%plot the mean tuning curve silent sound no laser
figure
subplot1(numsilentsounddurs,1, 'Max', [.95 .9])
fs=10;
p=0;
for dindex=1:numsilentsounddurs
    p=p+1;
    subplot1(p)
    
    hold on
    spiketimes1=mMSilentSoundOFF(dindex).spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nrepsOFF(dindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    bar(x,N,1);
    %offset=0;
    %yl=ylim;
    yl=ylimits;
    offset=-yl(2)/2;
    inc=(yl(2)/2.5)/max(nrepsOFF);
    if rasters==1
        for n=1:nrepsOFF
            spiketimes2=MSilentSoundOFF(dindex, n).spiketimes;
            offset=offset+inc;
            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
        end
    end
    
    if StimRecorded
        yl(1)=yl(1)-.05*diff(yl);
        offsetS=yl(1);
        for rep=1:nrepsOFF
            Stimtrace=MSilentSoundOFFStimtrace(dindex, rep, :);
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=.25*diff(yl)*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            plot(t, Stimtrace+offsetS, [.5 .5 .5])
        end
    else
        line([0 0+silentsounddurs(dindex)], ylimits(1)+[0 0], 'color', [.5 .5 .5], 'linewidth', 5)
    end
    if LaserRecorded
        yl(1)=yl(1)-.05*diff(yl);
        offsetL=yl(1);
        for rep=1:nrepsOFF
            Lasertrace=MSilentSoundOFFLasertrace(dindex, rep, :);
            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
            Lasertrace=.05*diff(yl)*Lasertrace;
            plot( t, Lasertrace+offsetL, 'c')
        end
    end
    if ~isempty(yl) ylim(yl); end
    
    %     xlim(xlimits) %(makes less sense for silent sound)
    xlim([-next/2 silentsounddurs(dindex)+next/2]);
    set(gca, 'fontsize', fs)
    h=title(sprintf('%s: \nSilent Sound no laser, tetrode%d cell%d, nreps: %d-%d',datadir,channel,out.cluster,min(nrepsOFF),max(nrepsOFF)));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
end

%plot the psth for silent sound with single laser pulse
figure
subplot1(numpulsewidths,1, 'Max', [.95 .9])
fs=10;
p=0;
if numlaserstarts>1 warning('\nusing only minimum laserstart for t-test');end
laserstart=min(laserstarts);
if numsilentsounddurs>1 warning('\nusing only minimum silentsound duration for t-test');end
ssdindex=1;

for pwindex=1:numpulsewidths
    p=p+1;
    subplot1(p)
    
    hold on
    spiketimes1=mMPulse(pwindex).spiketimes;
    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
    [N, x]=hist(spiketimes1, X);
    N=N./nrepsPulse(pwindex); %normalize to spike rate (averaged across trials)
    N=1000*N./binwidth; %normalize to spike rate in Hz
    bar(x,N,1);
    %offset=0;
    %yl=ylim;
    yl=ylimits;
    offset=-yl(2)/2;
    inc=(yl(2)/2.5)/max(nrepsPulse);
    if rasters==1
        for n=1:nrepsPulse(pwindex)
            spiketimes2=MPulse(pwindex, n).spiketimes;
            offset=offset+inc;
            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
        end
    end
    
    if StimRecorded
        yl(1)=yl(1)-.05*diff(yl);
        offsetS=yl(1);
        for rep=1:nrepsPulse(pwindex)
            Stimtrace=MPulseStimtrace(pwindex, rep, :);
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=.25*diff(yl)*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            plot(t, Stimtrace+offsetS, [.5 .5 .5])
        end
    else
        % do nothing
    end
    if LaserRecorded
        yl(1)=yl(1)-.05*diff(yl);
        offsetL=yl(1);
        for rep=1:nrepsPulse(pwindex)
            Lasertrace=MPulseLasertrace(pwindex, rep, :);
            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
            Lasertrace=.05*diff(yl)*Lasertrace;
            plot( t, Lasertrace+offsetL, 'c')
        end
    else
        line([laserstart laserstart+pulsewidths(pwindex)], ylimits(1)+[0 0], 'color', 'c', 'linewidth', 5)
    end
    if ~isempty(yl) ylim(yl); end
    
    %xlim(xlimits)
        xlim([-200 pulsewidths(pwindex)+200]);
    set(gca, 'fontsize', fs)
    h=title(sprintf('%s: \nSilent Sound with single laser pulse, tetrode%d cell%d, nreps: %d-%d',datadir,channel,out.cluster,min(nrepsOFF),max(nrepsOFF)));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    
    %How about a t-test for effect of Laser in first 25 ms
    for rep=1:nrepsPulse(pwindex)
        stop=laserstart+25;
        st=MPulse(pwindex,rep).spiketimes;
        spiketimes=st(st>laserstart & st<stop); % spiketimes in region
        ON(rep)=length(spiketimes);
    end
    for rep=1:nrepsOFF(ssdindex)
        stop=laserstart+25;
        st=MSilentSoundOFF(rep).spiketimes;
        spiketimes=st(st>laserstart & st<stop); % spiketimes in region
        OFF(rep)=length(spiketimes);
    end
    [h,p]=ttest2(ON, OFF, 'tail', 'right');
    yl=ylim;
    xl=xlim;
    text(xl(1)+25, .95*yl(2), sprintf('h=%d, p=%.4f effect of laser (t-test)', h,p), 'fontsize', 14)
    fprintf('\nch%d cell %d: h=%d, p=%.4f effect of laser (1-tailed t-test)',channel, clust, h,p)
end

%plot the psth for silent sound with laser train
for tpwindex=1:numtrainpulsewidths
    
    figure
    subplot1(numtrainnumpulses,numtrainisis, 'Max', [.95 .9])
    
    fs=10;
    p=0;
    for tnpindex=1:numtrainnumpulses
        for tiindex=1:numtrainisis
            p=p+1;
            subplot1(p)
            
            hold on
            spiketimes1=mMTrain(tnpindex,tpwindex,tiindex).spiketimes;
            X=xlimits(1):binwidth:xlimits(2); %specify bin centers
            [N, x]=hist(spiketimes1, X);
            N=N./nrepsTrain(tnpindex,tpwindex,tiindex); %normalize to spike rate (averaged across trials)
            N=1000*N./binwidth; %normalize to spike rate in Hz
            bar(x,N,1);
            %offset=0;
            %yl=ylim;
            yl=ylimits;
            offset=-yl(2)/2;
            inc=(yl(2)/2.5)/max(nrepsTrain(tnpindex,tpwindex,tiindex));
            if rasters==1
                for rep=1:nrepsTrain(tnpindex,tpwindex,tiindex)
                    spiketimes2=MTrain(tnpindex,tpwindex,tiindex, rep).spiketimes;
                    offset=offset+inc;
                    h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                end
            end
            
            if StimRecorded
                yl(1)=yl(1)-.05*diff(yl);
                offsetS=yl(1);
                for rep=1:nrepsTrain(tnpindex,tpwindex,tiindex)
                    Stimtrace=MTrainStimtrace(tnpindex,tpwindex,tiindex, rep, :);
                    Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                    Stimtrace=.25*diff(yl)*Stimtrace;
                    t=1:length(Stimtrace);
                    t=1000*t/out.samprate; %convert to ms
                    t=t+out.xlimits(1); %correct for xlim in original processing call
                    plot(t, Stimtrace+offsetS, [.5 .5 .5])
                end
            else
                %do nothing
                %line([0 0+silentsounddurs(dindex)], ylimits(1)+[0 0], 'color', [.5 .5 .5], 'linewidth', 5)
            end
            if LaserRecorded
                yl(1)=yl(1)-.05*diff(yl);
                offsetL=yl(1);
                for rep=1:nrepsTrain(tnpindex,tpwindex,tiindex)
                    Lasertrace=MTrainLasertrace(tnpindex,tpwindex,tiindex, rep, :);
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=.05*diff(yl)*Lasertrace;
                    plot( t, Lasertrace+offsetL, 'c')
                end
            else
                trainisi=trainisis(tiindex);
                trainpulsewidth=trainpulsewidths(tpwindex);
                for pnum=1:trainnumpulses(tnpindex)
%                     line([laserstart laserstart+trainpulsewidth], ylimits(1)+[0 0], 'color', 'c', 'linewidth', 5)
                    line([laserstart+trainisi*(pnum-1) laserstart+trainisi*(pnum-1)+trainpulsewidth], ylimits(1)+[0 0], 'color', 'c', 'linewidth', 5)
                end
            end
            if ~isempty(yl) ylim(yl); end
            
            xlim(xlimits)
            %xlim([-next/2 silentsounddurs(dindex)+next/2]);
            set(gca, 'fontsize', fs)
            h=title(sprintf('%s: \nSilent Sound with laser train, tetrode%d cell%d, nreps: %d-%d',datadir,channel,out.cluster,min(nrepsOFF),max(nrepsOFF)));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
        end
    end
end



