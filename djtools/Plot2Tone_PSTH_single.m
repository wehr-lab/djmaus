function Plot2Tone_PSTH_single(varargin)

%plots a single file of clustered spiking 2 tone data from djmaus
%
% usage: Plot2Tone_PSTH_single(datapath, t_filename, [xlimits],[ylimits], [binwidth])
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
    Process2Tone_PSTH_single(datadir,  t_filename, xlimits, ylimits);
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
    fprintf('\nloaded outfile.')
else
    Process2Tone_PSTH_single(datadir,  t_filename, xlimits, ylimits);
    load(outfilename);
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        fprintf('\nPlot called with xlimits [%d %d] but xlimits in outfile are [%d %d], re-processing...', xlimits(1), xlimits(2), out.xlimits(1), out.xlimits(2))
        Process2Tone_PSTH_single(datadir,  t_filename, xlimits, ylimits);
        load(outfilename);
    end
end

IL=out.IL; %whether there are any interleaved laser trials
freqs=out.freqs;
probefreqs=out.probefreqs;
amps=out.amps;
durs=out.durs;
SOAs=out.SOAs;
nreps=out.nreps;
numfreqs=out.numfreqs;
numprobefreqs=out.numprobefreqs;
numamps=out.numamps;
numdurs=out.numdurs;
numSOAs=out.numSOAs;
samprate=out.samprate; %in Hz
mM1ON=out.mM1ON;
mM1OFF=out.mM1OFF;
M1ON=out.M1ON;
M1OFF=out.M1OFF;
nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;
if isempty(xlimits) xlimits=out.xlimits;end
if numdurs>1 error('cannot handle multiple durations'), end
dindex=1;

LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
M_LaserStart=out.M_LaserStart;
M_LaserWidth=out.M_LaserWidth;
M_LaserNumPulses=out.M_LaserNumPulses;
M_LaserISI=out.M_LaserISI;
LaserStart=out.LaserStart;
LaserWidth=out.LaserWidth;
LaserNumPulses=out.LaserNumPulses;
LaserISI=out.LaserISI;
M1ONStim=out.M1ONStim;
M1ONLaser=out.M1ONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
mM1ONStim=out.mM1ONStim;
mM1ONLaser=out.mM1ONLaser;
M1OFFStim=out.M1OFFStim;
M1OFFLaser=out.M1OFFLaser;
mM1OFFStim=out.mM1OFFStim;
mM1OFFLaser=out.mM1OFFLaser;
fs=10; %fontsize

% %find optimal axis limits
if isempty(ylimits)
    ymax=0;
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            for pfindex=1:numprobefreqs
                for SOAindex=1:numSOAs
                    st=mM1OFF(findex, pfindex, aindex, SOAindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(st, X);
                    N=N./nrepsOFF(findex, pfindex, aindex, SOAindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    ymax= max(ymax,max(N));
                    
                    if IL
                        st=mM1ON(findex, pfindex, aindex, SOAindex).spiketimes;
                        X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                        [N, x]=hist(st, X);
                        N=N./nrepsON(findex, pfindex, aindex, SOAindex); %normalize to spike rate (averaged across trials)
                        N=1000*N./binwidth; %normalize to spike rate in Hz
                        ymax= max(ymax,max(N));
                    end
                end
            end
        end
    end
    ylimits=[-.3 ymax];
end

if freqs(1)==-1
    numy=1+(numfreqs-1)*(numprobefreqs-1);
else %no WN
    numy=numfreqs*numprobefreqs;
end


%plot the mean tuning curve OFF
figure('position',[650 100 600 600])
p=0;

subplot1(numy,numSOAs, 'Max', [.95 .9])
for findex=1:numfreqs
    for pfindex=1:numprobefreqs
        for aindex=1:numamps
            for SOAindex=1:numSOAs
                if xor(freqs(findex)==-1 , probefreqs(pfindex)==-1)
                    break
                end
                
                p=p+1;
                subplot1(p)
                hold on
                spiketimes1=mM1OFF(findex, pfindex, aindex, SOAindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(spiketimes1, X);
                N=N./nrepsOFF(findex, pfindex, aindex, SOAindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                offset=0;
                yl=ylimits;
                inc=(yl(2))/max(max(max(nrepsOFF)));
                if rasters==1
                    for n=1:nrepsOFF(findex, pfindex, aindex, SOAindex)
                        spiketimes2=M1OFF(findex, pfindex, aindex, SOAindex, n).spiketimes;
                        offset=offset+inc;
                        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    end
                end
                bar(x,N,1);
                %vpos=mean(ylimits);
                vpos=ylimits(2)-2*inc;

                if freqs(findex)==-1
                    text(0, vpos, 'WN', 'color', 'k')
                else
                    text(0, vpos, sprintf('%.1f kHz', freqs(findex)/1000), 'color', 'k')
                end
                if freqs(findex)==-1
                    text(SOAs(SOAindex), vpos-2*inc, 'WN', 'color', 'r')
                else
                    
                    text(SOAs(SOAindex), vpos-2*inc, sprintf('%.1f kHz', probefreqs(pfindex)/1000), 'color', 'r')
                end
                offsetS=ylimits(1)-.1*diff(ylimits);
                if StimRecorded
                    Stimtrace=squeeze(mM1OFFStim(findex, pfindex, aindex, SOAindex, :));
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
                ylim(ylimits2)
                
                xlim(xlimits)
                set(gca, 'fontsize', fs)
                %set(gca, 'xticklabel', '')
                %set(gca, 'yticklabel', '')
                if p==1
                    h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, OFF',datadir,channel,out.cluster,SOAs(SOAindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
                    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
                end
                xl = xlim; yl = ylim;
                h=text(0, range(yl),sprintf('%d ms',SOAs(SOAindex)), 'color', 'r');
                set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','position',[0 range(yl)*.8 0])
            end
        end
    end
end


if IL
    %plot the mean tuning curve ON
    figure('position',[1250 100 600 700])
    p=0;
    subplot1(numy,numSOAs, 'Max', [.95 .9])
    for findex=1:numfreqs
        for pfindex=1:numprobefreqs
            for aindex=1:numamps
                for SOAindex=1:numSOAs
                    if xor(freqs(findex)==-1 , probefreqs(pfindex)==-1)
                        break
                    end
                    p=p+1;
                    subplot1(p)
                    hold on
                    spiketimes1=mM1ON(findex, pfindex, aindex, SOAindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nrepsON(findex, pfindex, aindex, SOAindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=0;
                    yl=ylimits;
                    inc=(yl(2))/max(max(max(nrepsON)));
                    if rasters==1
                        for n=1:nrepsON(findex, pfindex, aindex, SOAindex)
                            spiketimes2=M1ON(findex, pfindex, aindex, SOAindex, n).spiketimes;
                            offset=offset+inc;
                            %this should plot a cyan line for every trial among the rasters
                            %it should accomodate trial-by-trial changes to
                            %Laser params
                            MLaserStart=M_LaserStart(findex, pfindex,aindex,SOAindex, n);
                            MLaserWidth=M_LaserWidth(findex, pfindex,aindex,SOAindex, n);
                            MLaserNumPulses=M_LaserNumPulses(findex, pfindex,aindex,SOAindex, n);
                            MLaserISI=M_LaserISI(findex, pfindex,aindex,SOAindex, n);
                            for np=1:MLaserNumPulses
                                plot([MLaserStart+(np-1)*(MLaserWidth+MLaserISI) MLaserStart+(np-1)*(MLaserWidth+MLaserISI)+MLaserWidth], [1 1]+yl(2)+offset, 'c')
                            end
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    bar(x, N,1);
                    line(xlimits, [0 0], 'color', 'k')
                    %vpos=mean(ylimits);
                    vpos=ylimits(2);
                    text(0, vpos, sprintf('%.1f kHz', freq/1000), 'color', 'k')
                    text(SOAs(SOAindex), vpos, sprintf('%.1f kHz', probefreq/1000),'color', 'r')
                    
                    if StimRecorded
                        Stimtrace=squeeze(mM1OFFStim(findex, pfindex, aindex, SOAindex, :));
                        Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                        Stimtrace=.25*diff(ylimits)*Stimtrace;
                        t=1:length(Stimtrace);
                        t=1000*t/out.samprate; %convert to ms
                        t=t+out.xlimits(1); %correct for xlim in original processing call
                        offset=ylimits(1)-.1*diff(ylimits);
                        stimh=plot(t, Stimtrace+offset, 'm');
                        uistack(stimh, 'bottom')

                    else
                        line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                    end
                    
                    ylimits2(2)=ylimits(2)+offset;
                    ylimits2(2)=2*ylimits(2);
                    ylimits2(1)=-2;
                    ylim(ylimits2)
                    %                 ylim([-2 1.1*(yl(2)+offset)])
                    
                    if LaserRecorded
                        for rep=1:nrepsON(findex, pfindex, aindex, SOAindex)
                            Lasertrace=squeeze(M1ONLaser(findex, pfindex, aindex, SOAindex,rep, :));
                            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                            Lasertrace=.05*diff(ylimits)*Lasertrace;
                            plot( t, Lasertrace+offset, 'c')
                        end
                    end
                    %this should plot a cyan line at the unique Laser
                    %params - not sure what will happen if not scalar
                    for np=1:LaserNumPulses
                        plot([LaserStart+(np-1)*(LaserWidth+LaserISI) LaserStart+(np-1)*(LaserWidth+LaserISI)+LaserWidth], [-2 -2], 'c', 'linewidth', 2)
                    end
                    line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
                    
                    xlim(xlimits)
                    set(gca, 'fontsize', fs)
                    set(gca, 'xticklabel', '')
                    set(gca, 'yticklabel', '')
                    
                end
            end
        end
    end
    subplot1(1)
    h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, ON',datadir,channel,out.cluster,SOAs(SOAindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
    set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
    
    %label amps and freqs
    p=0;
    for yindex=ystart:ystep:yend
        for findex=1:numfreqs
            for pfindex=1:numprobefreqs
                p=p+1;
                subplot1(p)
                
                if findex==1
                    if numamps>=numSOAs
                        text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(yindex)))
                    else
                        text(xlimits(1)+diff(xlimits)/20, mean(ylimits), int2str(SOAs(yindex)))
                    end
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                end
            end
        end
    end
end