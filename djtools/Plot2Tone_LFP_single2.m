function Plot2Tone_LFP_single2(varargin)

%plots a single file of clustered spiking 2 tone data from djmaus
%
% usage: Plot2Tone_PSTH_single(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% (xlimits, ylimits, binwidth are optional)
%
%Processes data if outfile is not found;

rasters=1;
force_reprocess=0;

datadir=varargin{1};

try
    channel=varargin{2};
catch
    prompt=('please enter channel number: ');
    channel=input(prompt);
end
if strcmp('char',class(channel))
    channel=str2num(channel);
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

high_pass_cutoff=400;
low_pass_cutoff=300;
[a,b]=butter(1, high_pass_cutoff/(30e3/2), 'high');


if force_reprocess
    fprintf('\nForce re-process\n')
    Process2Tone_LFP_single2(datadir,  channel, xlimits, ylimits);
end
cd(datadir)
outfilename=sprintf('outLFP_ch%d.mat',channel);
d=dir(outfilename);
if ~isempty(d)
    load(outfilename)
else
    Process2Tone_LFP_single2(datadir,  channel, xlimits, ylimits);
    load(outfilename);
end

SOAs=out.SOAs; % SOAs and LaserISIs are the same
%if xlimits are requested but don't match those in outfile, force preprocess
% xlimits=[-50 max(SOAs)+150];
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%d %d] but xlimits in outfile are [%d %d], re-processing...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        Process2Tone_LFP_single2(datadir,  channel, xlimits, ylimits);
        load(outfilename);
    end
end

IL=out.IL; %whether there are any interleaved laser trials
freqs=out.freqs; % only one -1 WN
probefreqs=out.probefreqs; % probe freq -1 WN or 0 silent sound
durs=out.durs; % 25 ms one dur for WN and silent sound

nreps=out.nreps;
numfreqs=out.numfreqs;
numprobefreqs=out.numprobefreqs;
numamps=out.numamps;
numdurs=out.numdurs;
numSOAs=out.numSOAs;
samprate=out.samprate; %in Hz
M1=out.M1;
mM1=out.mM1;

numLaserNumPulses=out.numLaserNumPulses;
numLaserISIs=out.numLaserISIs;
LaserNumPulses=out.LaserNumPulses;
LaserISIs=out.LaserISIs;

if isempty(xlimits) xlimits=out.xlimits;end
if numdurs>1 error('cannot handle multiple durations'), end
dindex=1;
LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
M1Laser=out.M1Laser;
mM1Laser=out.mM1Laser;
M1Stim=out.M1Stim;
mM1Stim=out.mM1Stim;
if isempty(xlimits);
    xlimits=out.xlimits;
end
fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))
fs=10; %fontsize
% %find optimal axis limits
if isempty(ylimits)
    ylimits=[0 0];
    for findex=1:numfreqs
        for pfindex=1:numprobefreqs
            for SOAindex=1:numSOAs
                trace1=squeeze(mM1(findex, pfindex, SOAindex,:));
                trace1=trace1-mean(trace1(1:100));
                if min([trace1])<ylimits(1); ylimits(1)=min([trace1]);end
                if max([trace1])>ylimits(2); ylimits(2)=max([trace1]);end
                
            end
        end
    end
end

ylimits=round(ylimits*100)/100;

%plot the mean tuning curve OFF
figure('position',[650 100 600 600])
p=0;

%since SOA and laserISIs are the same, those are the two paramaters varying
%in this stimulus, i will plot them as one


% Plot OFF WN and ON silent sound laser pulses on the same plot, 2 tone
% only
subplot1(numSOAs-1,1, 'Max', [.95 .9])
ylimits2=ylimits;
for findex=[numfreqs:-1:1]
    pfindex=1;
    p=0;
    for SOAindex=2:numSOAs
        p=p+1;
        subplot1(p)
        hold on
        trace1=squeeze(mM1(findex, pfindex, SOAindex,:));
        [b,a]=butter(1, low_pass_cutoff/(samprate/2), 'low');
        trace1=filtfilt(b,a,trace1);
        trace1=trace1-mean(trace1(1:10));
        
        t=1:length(trace1);
        t=1000*t/out.samprate; %convert to ms
        t=t+out.xlimits(1); %correct for xlim in original processing call
        if findex==1
            plot(t, trace1, 'k');
        else
            plot(t, trace1, 'g');
        end
        vpos=ylimits(1);
        text(0, vpos*2, sprintf('%d ms', SOAs(SOAindex)))
       
        offsetS=ylimits(1)-.1*diff(ylimits);
        if StimRecorded
            Stimtrace=squeeze(mM1Stim(findex, pfindex, SOAindex, :));
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=1.25*diff(ylimits)*Stimtrace;
            t=1:length(Stimtrace);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            stimh=plot(t, Stimtrace+offsetS, 'm');
            uistack(stimh, 'bottom')
            ylimits2(1)=2*offsetS;
        else
            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            ylimits2(1)=-2;
        end
        if LaserRecorded
            for rep=1:nreps(findex, pfindex, SOAindex)
                Lasertrace=squeeze(M1Laser(findex, pfindex, SOAindex,rep, :));
                Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                Lasertrace=.05*diff(ylimits)*Lasertrace;
                plot( t, Lasertrace+offsetS, 'c')
            end
        end
        line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
        line(xlimits, [0 0], 'color', 'k')
        ylim(ylimits2)
        
        xlim(xlimits)
        set(gca, 'fontsize', fs)
        %set(gca, 'xticklabel', '')
        %set(gca, 'yticklabel', '')
        if p==1
            h=title(sprintf('%s: \nchannel%d %dms, nreps: %d-%d',datadir,channel,SOAs(SOAindex),min(min(min(nreps))),max(max(max(nreps)))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            
        end
        if p==1
        ylabel('FR Hz')
        end
%         xl = xlim; yl = ylim;
%         h=text(0, range(yl),sprintf('%d ms',SOAs(SOAindex)), 'color', 'r');
%         set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','position',[0 range(yl)*.8 0])
    end
end

%plot single pulse or single WN only
figure;
ylimits2=ylimits;
for findex=[numfreqs:-1:1] %plot laser first
    pfindex=2; %single stimulus
    SOAindex=1; % 0 SOA
    
    hold on
    trace1=squeeze(mM1(findex, pfindex, SOAindex,:));
    [b,a]=butter(1, low_pass_cutoff/(samprate/2), 'low');
    trace1=filtfilt(b,a,trace1);
    trace1=trace1-mean(trace1(1:10));
    
    t=1:length(trace1);
    t=1000*t/out.samprate; %convert to ms
    t=t+out.xlimits(1); %correct for xlim in original processing call
    if findex==1
    plot(t, trace1, 'k');
    else
        plot(t, trace1, 'g');
    end
        
    offsetS=ylimits(1)-.1*diff(ylimits);
    
    if StimRecorded
        Stimtrace=squeeze(mM1Stim(findex, pfindex, SOAindex, :));
        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
        Stimtrace=1.25*diff(ylimits)*Stimtrace;
        t=1:length(Stimtrace);
        t=1000*t/out.samprate; %convert to ms
        t=t+out.xlimits(1); %correct for xlim in original processing call
        stimh=plot(t, Stimtrace+offsetS, 'm');
        uistack(stimh, 'bottom')
        ylimits2(1)=2*offsetS;
    else
        line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
        ylimits2(1)=-2;
    end
    if LaserRecorded
        for rep=1:nreps(findex, pfindex, SOAindex)
            Lasertrace=squeeze(M1Laser(findex, pfindex, SOAindex,rep, :));
            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
            Lasertrace=.05*diff(ylimits)*Lasertrace;
            plot( t, Lasertrace+offsetS, 'c')
        end
    end
    line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
    line(xlimits, [0 0], 'color', 'k')
    ylim(ylimits2)
    
    xlim(xlimits)
    set(gca, 'fontsize', fs)
        h=title('WN burst and laser pulse alone');
        set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    ylabel('FR Hz')
    xl = xlim; yl = ylim;
    h=text(0, range(yl),sprintf('%d ms',SOAs(SOAindex)), 'color', 'r');
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','position',[0 range(yl)*.8 0])
    
end
