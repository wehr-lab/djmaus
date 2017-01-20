function out=PlotArchPVRev3(varargin)

%like PlotArchPVRev3 but only plots smoothed firing rate curves with
%different colors for laser start times
%
%plots clustered spiking tuning curve data from djmaus
%same as PlotArchPVRev1 but takes a single filename as input (instead of channel &
%clust)
% usage: PlotArchPVRev1(datapath, t_filename, [xlimits],[ylimits], [binwidth])
% ( xlimits, ylimits are optional)
% xlimits default to [0 200]
%
% if no outfile is found, does nothing


rasters=1;

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};
t_filename=varargin{2};
try
    xlimits=varargin{3};
catch
    xlimits=[-100 300];
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

[p,f,ext]=fileparts(t_filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\n%s', t_filename)
fprintf('\n%s', outfilename)


fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

djPrefs;
global pref
cd(datadir)


d=dir(outfilename);
if isempty(d)
    ProcessArchPVRev2(datadir, t_filename, xlimits)
end

load(outfilename)

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

M_LaserStart=out.M_LaserStart;
M_LaserWidth=out.M_LaserWidth;
M_LaserNumPulses=out.M_LaserNumPulses;
M_LaserISI=out.M_LaserISI;
LaserStart=out.LaserStart;
LaserWidth=out.LaserWidth;
LaserNumPulses=out.LaserNumPulses;
LaserISI=out.LaserISI;
fs=10; %fontsize

% %find optimal axis limits
if isempty(ylimits)
    ymax=0;
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            for dindex=1:numdurs
                st=mM1OFF(findex, aindex, dindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(st, X);
                N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));
                
                st=mM1ON(findex, aindex, dindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                [N, x]=hist(st, X);
                N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                ymax= max(ymax,max(N));
            end
        end
    end
    ylimits=[-.3 ymax];
end

% if numamps>=numdurs
%     %plot freqs x amps, with each dur in a separate window
%     numy=numamps;
%     ystep=-1;
%     ystart=numy;
%     yend=1;
%     numw=numdurs;
% elseif numdurs>numamps
%     %plot freqs x durs, with each amp in a separate window
%     numw=numamps;
%     numy=numdurs;
%     ystart=1;
%     ystep=1;
%     yend=numy;
% end

ystart=1;
yend=numamps;
numy=numamps;
ystep=1;
numw=1;

%plot the mean tuning curve OFF
if(0)
    for windex=1:numw
        figure
        p=0;
        subplot1(numy,numfreqs, 'Max', [.95 .9])
        for yindex=ystart:ystep:yend
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                hold on
                if numamps>=numdurs, aindex=yindex; dindex=windex;
                else aindex=windex; dindex=yindex;
                end
                spiketimes1=mM1OFF(findex, aindex, dindex).spiketimes;
                X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                %                 [N, x]=hist(spiketimes1, X);
                %                 N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                %                 N=1000*N./binwidth; %normalize to spike rate in Hz
                %                 offset=0;
                %                 yl=ylimits;
                %                 inc=(yl(2))/max(max(max(nreps)));
                %                 if rasters==1
                %                     for n=1:nrepsOFF(findex, aindex, dindex)
                %                         spiketimes2=M1OFF(findex, aindex, dindex, n).spiketimes;
                %                         offset=offset+inc;
                %                         h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                %                     end
                %                 end
                %                 ph=plot(x, N);
                set(ph, 'color', 'k')
                line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
                line(xlimits, [0 0], 'color', 'k')
                ylimits2(2)=ylimits(2)*3;
                ylimits2(1)=-2;
                ylim(ylimits2)
                
                xlim(xlimits)
                set(gca, 'fontsize', fs)
                %set(gca, 'xticklabel', '')
                %set(gca, 'yticklabel', '')
                
            end
        end
        subplot1(1)
        h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, OFF',datadir,channel,out.cluster,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
        set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
        
        %label amps and freqs
        p=0;
        for yindex=ystart:ystep:yend
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                if findex==1
                    if numamps>=numdurs
                        text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(yindex)))
                    else
                        text(xlimits(1)+diff(xlimits)/20, mean(ylimits), int2str(durs(yindex)))
                    end
                    
                end
                if yindex==1
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

%do a little test to see if the WN response is significant
for rep=1:median(out.nrepsON) %was min, but there is a zero rep condition that's irksome
    spiketimes=out.M1ON(1,:,1,rep).spiketimes;
    WNresponse(rep)=length(find(spiketimes>0 & spiketimes<100));
    spont(rep)=length(find(spiketimes>-100 & spiketimes<0));
end
[ptest,htest]=ranksum(WNresponse, spont, 'tail', 'right');
fprintf('\nIs there a significant WN response? h=%d, p=%.4f', htest, ptest)

htest=1; %HACK

% if ~htest
%     Durn=nan;
%     tFR=nan;
%     return
% end

%do a little test to see if the laser has a significant effect on th eWN response
%for this I compare FR in 0-100 ms for laserstart=-50 (ON) to OFF
for rep=1:median(out.nrepsON)
    spiketimesON=out.M1ON(1,1,1,rep).spiketimes;
    WNresponseON(rep)=length(find(spiketimesON>0 & spiketimesON<100));
end
for rep=1:median(out.nrepsOFF(find(out.nrepsOFF))) %was min, but there is are many zero rep conditions that're irksome
    condition=find(out.nrepsOFF); %find the "laserstart" or "aindex" that is actually populated in M1OFF
    spiketimesOFF=out.M1OFF(1,condition,1,rep).spiketimes;
    WNresponseOFF(rep)=length(find(spiketimesOFF>0 & spiketimesOFF<100));
end

[ptestArcheffect,htestArcheffect]=ranksum(WNresponseON, WNresponseOFF, 'tail', 'right');
fprintf('\nDoes laser have an effect on WN response? h=%d, p=%.4f', htestArcheffect, ptestArcheffect)
fprintf('\n')

% htestArcheffect=1; %HACK
% 
% if ~htestArcheffect
%     Durn=nan;
%     tFR=nan;
%     return
% end

% cmap=colormap(jet(numy));
cmap=colormap(cool(numy));
PH2=[];
if IL
    %     get a max firing rate for normalization purposes
    
    maxFR=0;
    for y=1:numy
        spiketimes1=mM1ON(1, y, 1).spiketimes;
        sigma=10;
        [t, fr]=GaussSmooth(spiketimes1, sigma, xlimits);
        fr=fr./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
        if max(fr)>maxFR
            maxFR=max(fr);
        end
    end
    
    spiketimes2=[];
    for y=1:numy
        spiketimes2=[spiketimes2 mM1ON(1, y, 1).spiketimes];
    end
    [t, fr]=GaussSmooth(spiketimes2, sigma, xlimits);
    ti_peakFR=find(fr==max(fr), 1); %time index at which fr hits peak
        
    if 1        %plot rasters for ON
        for windex=1:numw
            figure
            p=0;
            subplot1(numy,numfreqs, 'Max', [.95 .9], 'Gap', [.01 0])
            for yindex=ystart:ystep:yend
                for findex=1:numfreqs
                    p=p+1;
                    subplot1(p)
                    hold on
                    if numamps>=numdurs, aindex=yindex; dindex=windex;
                    else aindex=windex; dindex=yindex;
                    end
                    spiketimes1=mM1ON(findex, aindex, dindex).spiketimes;
                    X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    [N, x]=hist(spiketimes1, X);
                    N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    N=1000*N./binwidth; %normalize to spike rate in Hz
                    offset=0;
                    yl=ylimits;
                    inc=(yl(2))/max(max(max(nreps)));
                    if rasters==1
                        for n=1:nrepsON(findex, aindex, dindex)
                            spiketimes2=M1ON(findex, aindex, dindex, n).spiketimes;
                            offset=offset+inc;
                            %this should plot a cyan line for every trial among the rasters
                            %it should accomodate trial-by-trial changes to
                            %Laser params
                            MLaserStart=M_LaserStart(findex,aindex,dindex, n);
                            MLaserWidth=M_LaserWidth(findex,aindex,dindex, n);
                            MLaserNumPulses=M_LaserNumPulses(findex,aindex,dindex, n);
                            MLaserISI=M_LaserISI(findex,aindex,dindex, n);
                            for np=1:MLaserNumPulses
                                ph=plot([MLaserStart+(np-1)*(MLaserWidth+MLaserISI) MLaserStart+(np-1)*(MLaserWidth+MLaserISI)+MLaserWidth], [1 1]+yl(2)+offset);
                                set(ph, 'color',  [0 1 1 .15])
                            end
                            h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                        end
                    end
                    bar(x, N,1);
                    line(xlimits, [0 0], 'color', 'k')
                    
                    ylimits2(2)=ylimits(2)*3;
                    ylimits2(1)=-2;
                    ylim(ylimits2(:))
                    
                    %this should plot a cyan line at the unique Laser
                    %params - not sure what will happen if not scalar
                    for np=1:LaserNumPulses
                        %    plot([LaserStart+(np-1)*(LaserWidth+LaserISI) LaserStart+(np-1)*(LaserWidth+LaserISI)+LaserWidth], [-2 -2], 'c', 'linewidth', 2)
                    end
                    line([0 0+durs(dindex)], [-.2 -.2], 'color', 'm', 'linewidth', 4)
                    
                    xlim(xlimits)
                    set(gca, 'fontsize', fs)
                    set(gca, 'xticklabel', '')
                    set(gca, 'yticklabel', '')
                    
                end
            end
            subplot1(1)
            h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, ON',datadir,channel,out.cluster,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            orient tall
        end %if plot rasters
        
        %plot the mean tuning curve ON
        for windex=1:numw
            figure
            p=0;
            %             subplot1(numy,numfreqs, 'Max', [.95 .9], 'Gap', [.01 0])
            for yindex=ystart:ystep:yend
                for findex=1:numfreqs
                    p=p+1;
                    %                     subplot1(p)
                    hold on
                    if numamps>=numdurs, aindex=yindex; dindex=windex;
                    else aindex=windex; dindex=yindex;
                    end
                    spiketimes1=mM1ON(findex, aindex, dindex).spiketimes;
                    %                     X=xlimits(1):binwidth:xlimits(2); %specify bin centers
                    %                     [N, x]=hist(spiketimes1, X);
                    %                     N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    %                     N=1000*N./binwidth; %normalize to spike rate in Hz
                    %                     offset=0;
                    %                     yl=ylimits;
                    %                     inc=(yl(2))/max(max(max(nreps)));
                    %                     if rasters==1
                    %                         for n=1:nrepsON(findex, aindex, dindex)
                    %                             spiketimes2=M1ON(findex, aindex, dindex, n).spiketimes;
                    %                             offset=offset+inc;
                    %                             %this should plot a cyan line for every trial among the rasters
                    %                             %it should accomodate trial-by-trial changes to
                    %                             %Laser params
                    %                             MLaserStart=M_LaserStart(findex,aindex,dindex, n);
                    %                             MLaserWidth=M_LaserWidth(findex,aindex,dindex, n);
                    %                             MLaserNumPulses=M_LaserNumPulses(findex,aindex,dindex, n);
                    %                             MLaserISI=M_LaserISI(findex,aindex,dindex, n);
                    %                             for np=1:MLaserNumPulses
                    %                                 ph=plot([MLaserStart+(np-1)*(MLaserWidth+MLaserISI) MLaserStart+(np-1)*(MLaserWidth+MLaserISI)+MLaserWidth], [1 1]+yl(2)+offset);
                    %                             set(ph, 'color',  [0 1 1 .15])
                    %                             end
                    %                             h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    %                         end
                    %                     end
                    
                    
                    
                    sigma=10;
                    [t, fr]=GaussSmooth(spiketimes1, sigma, xlimits);
                    fr=fr./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                    ph=plot(t(150:end), fr(150:end)); %crop first 15ms
                    c=cmap(yindex,:);
                    set(ph, 'color', c)
                    fig=gcf;
                    
                    %extract some estimate of duration or integrated FR
                    Intg(p)=sum(fr(1000:2000));
                    
                    %                 nfr=fr./max(fr);
                    nfr=fr./maxFR; %global normalization
                    t_poststim=find(t>0);
                    THRESH=.25;
                    above=find(nfr(t_poststim)>THRESH);
                    if ~isempty(above)
                        upcross=min(above); %index within t_poststim
                        d_above=diff(above);
                        if all(d_above==1) %unimodal
                            downcross=max(above); %index within t_poststim
                        else
                            downcross=upcross+find(d_above>1, 1, 'first');
                        end
                        
                        Durn(p)=(downcross-upcross)/10;
                    else
                        Durn(p)=nan;
                    end
                    
                    %thought from run: I could plot FR at a specific latency, as a different sort of proxy for response duration.
                    %looks like a good time point = 35ms
                    %t is sampled at 0.1 ms (from gauss)
                    ti=find(t>20, 1, 'first');
                    tFR20(p)=nfr(ti); %FR at t=20 ms
                    tFR25(p)=nfr(ti+50); %FR at t
                    tFR30(p)=nfr(ti+100); %FR at t
                    tFR35(p)=nfr(ti+150); %FR at t
                    tFR40(p)=nfr(ti+200); %FR at t
                    tFR45(p)=nfr(ti+250); %FR at t
                    tFR50(p)=nfr(ti+300); %FR at t
                    tFR55(p)=nfr(ti+350); %FR at t
                    tFR60(p)=nfr(ti+400); %FR at t
                    
                    tFRp0(p)=nfr(ti_peakFR); %FR at t=peak FR
                    tFRp5(p)=nfr(ti_peakFR+50); %FR at t=5 ms after peak FR 
                    tFRp10(p)=nfr(ti_peakFR+100); %FR at t=10 ms after peak FR
                    tFRp15(p)=nfr(ti_peakFR+150); %FR at t=
                    tFRp20(p)=nfr(ti_peakFR+200); %FR at t=
                    tFRp25(p)=nfr(ti_peakFR+250); %FR at t=
                    tFRp30(p)=nfr(ti_peakFR+300); %FR at t=
                    tFRp35(p)=nfr(ti_peakFR+350); %FR at t=
                    tFRp40(p)=nfr(ti_peakFR+400); %FR at t=
                    
                    
                    
                    %                 figure(100)
                    %                 hold on
                    %                 plot(t, nfr, 'k')
                    %                 plot(t(t_poststim), nfr(t_poststim), 'r')
                    %                 plot(t(t_poststim(upcross)), nfr(t_poststim(upcross)), 'o')
                    %                 plot(t(t_poststim(downcross)), nfr(t_poststim(downcross)), 'o')
                    %                 figure(fig)
                    
                    if aindex<=size(M_LaserStart, 2)
                        MLaserStart=M_LaserStart(findex,aindex,dindex, 1);
                        yl=ylim;
                        vpos=yl(1)+.05*diff(yl);
                        ph2=plot(MLaserStart, vpos, '^');
                        set(ph2, 'color', c)
                        PH2=[PH2 ph2];
                        line(xlimits, [0 0], 'color', 'k')
                        
                        %                     ylimits2(2)=ylimits(2)*3;
                        %                     ylimits2(1)=-2;
                        %                     ylim(ylimits)
                        
                        %this should plot a cyan line at the unique Laser
                        %params - not sure what will happen if not scalar
                        for np=1:LaserNumPulses
                            %    plot([LaserStart+(np-1)*(LaserWidth+LaserISI) LaserStart+(np-1)*(LaserWidth+LaserISI)+LaserWidth], [-2 -2], 'c', 'linewidth', 2)
                        end
                        %set(gca, 'xticklabel', '')
                        %set(gca, 'yticklabel', '')
                    end
                end
            end
            vpos=-.05*diff(yl);
            line([0 0+durs(dindex)], vpos*[1 1], 'color', 'm', 'linewidth', 4)
            vpos=-.025*diff(yl);
            set(PH2, 'ydata', vpos)
            xlim(xlimits)
            set(gca, 'fontsize', fs)
            
            subplot1(1)
            h=title(sprintf('%s: \ntetrode%d cell%d %dms, nreps: %d-%d, ON',datadir,channel,out.cluster,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
            set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal')
            
            
            %label amps and freqs
            p=0;
            %             for yindex=ystart:ystep:yend
            %                 for findex=1:numfreqs
            %                     p=p+1;
            %                     subplot1(p)
            %
            %                     if findex==1
            %                         if numamps>=numdurs
            %                             text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(yindex)))
            %                         else
            %                             text(xlimits(1)+diff(xlimits)/20, mean(ylimits), int2str(durs(yindex)))
            %                         end
            %                     end
            %                     if aindex==1
            %                         if mod(findex,2) %odd freq
            %                             vpos=ylimits(1)-mean(ylimits);
            %                         else
            %                             vpos=ylimits(1)-mean(ylimits);
            %                         end
            %                         text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
            %                     end
            %                 end
            %             end
        end
    end
end

laserstarts=squeeze(M_LaserStart(1,:,1,1));
if length(laserstarts) ~= length(Intg)
    %working through some special cases here
    if length(laserstarts) == length(Intg) -1 & Intg(end)==0
        Intg=Intg(1:end-1);
        Durn=Durn(1:end-1);
        tFR35=tFR35(1:end-1);
    end
end

figure
plot(laserstarts, Intg)
xlabel('laser start relative to tone onset')
ylabel('integrated FR 0-100 ms')

figure
plot(laserstarts, Durn)
xlabel('laser start relative to tone onset')
ylabel('duration, ms')
title(sprintf('%s %s', datadir, outfilename), 'interpreter', 'none')

figure
plot(laserstarts, tFR35)
xlabel('laser start relative to tone onset')
ylabel(sprintf('nFR at %.1f ms', t(ti)))
title(sprintf('%s %s', datadir, outfilename), 'interpreter', 'none')

out.Durn=Durn;
out.Intg=Intg;
out.tFR20=tFR20;
out.tFR25=tFR25
out.tFR30=tFR30;
out.tFR35=tFR35;
out.tFR40=tFR40;
out.tFR45=tFR45;
out.tFR50=tFR50;
out.tFR55=tFR55;
out.tFR60=tFR60;
out.tFRp0=tFRp0;
out.tFRp5=tFRp5;
out.tFRp10=tFRp10;
out.tFRp15=tFRp15;
out.tFRp20=tFRp20;
out.tFRp25=tFRp25;
out.tFRp30=tFRp30;
out.tFRp35=tFRp35;
out.tFRp40=tFRp40;
out.ptest=ptest;
out.htest=htest;
out.htestArcheffect=htestArcheffect;
out.ptestArcheffect=ptestArcheffect;
out.laserstarts=laserstarts;
                    