function PlotArchPVRev2(varargin)

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
    
    %plot the mean tuning curve OFF
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
                [N, x]=hist(spiketimes1, X);
                N=N./nreps(findex, aindex, dindex); %normalize to spike rate (averaged across trials)
                N=1000*N./binwidth; %normalize to spike rate in Hz
                offset=0;
                yl=ylimits;
                inc=(yl(2))/max(max(max(nreps)));
                if rasters==1
                    for n=1:nrepsOFF(findex, aindex, dindex)
                        spiketimes2=M1OFF(findex, aindex, dindex, n).spiketimes;
                        offset=offset+inc;
                        h=plot(spiketimes2, yl(2)+ones(size(spiketimes2))+offset, '.k');
                    end
                end
                bar(x, N,1);
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
    
    if IL
        %plot the mean tuning curve ON
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
    
end %main cluster loop
