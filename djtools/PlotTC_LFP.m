function PlotTC_LFP(varargin)

% plots continuous tuning curve data from djmaus
%
% usage: PlotTC_LFP(datapath, [channel], [xlimits],[ylimits])
% (xlimits & ylimits are optional)
% xlimits default to [0 200]
% channel number should be an integer
%
% Processes data if outfile is not found;

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
    xlimits=[0 200];
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end

high_pass_cutoff=400;
[a,b]=butter(1, high_pass_cutoff/(30e3/2), 'high');
%[a,b]=butter(1, [ lo hi]/(30e3/2));

fprintf('\nusing xlimits [%d-%d]', xlimits(1), xlimits(2))

force_reprocess=1;
if force_reprocess
    fprintf('\nForce Re-process')
    fprintf('\ncalling ProcessTC_LFP')
    ProcessTC_LFP(datadir,  channel, xlimits, ylimits);
end

cd(datadir)
outfilename=sprintf('outLFP_ch%d.mat',channel);
d=dir(outfilename);
if ~isempty(d)
    load(outfilename)
else
    fprintf('\ncalling ProcessTC_LFP')
    ProcessTC_LFP(datadir,  channel, xlimits, ylimits);
    load(outfilename);
end



M1stim=out.M1stim;
freqs=out.freqs;
amps=out.amps;
durs=out.durs;
nreps=out.nreps;
numfreqs=out.numfreqs;
numamps=out.numamps;
numdurs=out.numdurs;
samprate=out.samprate; %in Hz
M1=out.M1;
scaledtrace=out.scaledtrace;
traces_to_keep=out.traces_to_keep;
mM1=out.mM1;
mM1ON=out.mM1ON;
mM1OFF=out.mM1OFF;
mM1ONLaser=out.mM1ONLaser;
M1ONLaser=out.M1ONLaser;
mM1OFFLaser=out.mM1OFFLaser;
mM1ONStim=out.mM1ONStim;
mM1OFFStim=out.mM1OFFStim;
nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;
IL=out.IL; %whether there were interleaved laser trials or not

% %find optimal axis limits
if isempty(ylimits)
    ylimits=[0 0];
    for dindex=1:numdurs
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                trace1=squeeze(mM1(findex, aindex, dindex, :));
                %                 trace1=filtfilt(b,a,trace1);
                trace1=trace1-mean(trace1(1:100));
                if min([trace1])<ylimits(1); ylimits(1)=min([trace1]);end
                if max([trace1])>ylimits(2); ylimits(2)=max([trace1]);end
            end
        end
    end
end

ylimits=round(ylimits*100)/100;

%plot the mean tuning curve BOTH
if IL
    for dindex=1:numdurs
        figure
        p=0;
        subplot1(numamps,numfreqs)
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                %
                %             for i=1:length(nrepsON)
                %             endmM1ONLaser=out.mM1ONLaser;
                
                %             trace1=squeeze(squeeze(meanONbl(findex, aindex, dindex, :)));
                %             trace2=(squeeze(meanOFFbl(findex, aindex, dindex, :)));
                trace1=squeeze(squeeze(out.mM1ON(findex, aindex, dindex, :)));
                trace2=(squeeze(out.mM1OFF(findex, aindex, dindex, :)));
                
                trace1=trace1 -mean(trace1(1:10));
                trace2=trace2-mean(trace2(1:10));
                %             trace1=filtfilt(b,a,trace1);
                %             trace2=filtfilt(b,a,trace2);
                t=1:length(trace1);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                line([0 0+durs(dindex)], [0 0], 'color', 'm', 'linewidth', 5)
                plot(t, trace1, 'b');
                hold on; plot(t, trace2, 'k');
                %ylim([-5000 5000])
                ylim(ylimits)
                xlim(xlimits)
                box off
                
            end
        end
        subplot1(1)
        h=title(sprintf('%s: %dms, nreps: %d-%d, ON&OFF',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
        
        %label amps and freqs
        p=0;
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                if findex==1
                    text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    if freqs(findex)>0
                        text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                    elseif freqs(findex)==-1000
                        text(xlimits(1), vpos, 'WN')
                    elseif freqs(findex)==-2000
                        text(xlimits(1), vpos, 'SS')
                    end
                end
            end
        end
    end
end

%plot the mean tuning curve OFF
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            trace1=squeeze(mM1OFF(findex, aindex, dindex, :));
            %             trace1=filtfilt(b,a,trace1);
            trace1=trace1 -mean(trace1(1:100));
            
            Lasertrace=squeeze(mM1OFFLaser(findex, aindex, dindex, :));
            Lasertrace=Lasertrace -mean(Lasertrace(1:100));
            Lasertrace=.05*diff(ylimits)*Lasertrace;
            Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
            Stimtrace=Stimtrace -mean(Stimtrace(1:100));
            Stimtrace=.05*diff(ylimits)*Stimtrace;
            
            
            t=1:length(trace1);
            t=1000*t/out.samprate; %convert to ms
            t=t+out.xlimits(1); %correct for xlim in original processing call
            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            hold on; plot(t, trace1, 'k');
            offset=ylimits(1)+.1*diff(ylimits);
            plot(t, Stimtrace+offset, 'm', t, Lasertrace+offset, 'c')
            ylim(ylimits)
            xlim(xlimits)
            xlabel off
            ylabel off
            axis off
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s: %dms, nreps: %d-%d',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF)))));
    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
    
    %label amps and freqs
    p=0;
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            if findex==1
                text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                endmM1ONLaser=out.mM1ONLaser;
                
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    if freqs(findex)>0
                        text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                    elseif freqs(findex)==-1000
                        text(xlimits(1), vpos, 'WN')
                    elseif freqs(findex)==-2000
                        text(xlimits(1), vpos, 'SS')
                    end
                end
                %             if findex==numfreqs && aindex==numamps
                %                 axis on
                %                 ylab=[ceil(ylimits(1)*10)/10 floor(ylimits(2)*10)/10];
                %                 set(gca,'ytick',ylab,'yticklabel',ylab,'YAxisLocation','right')
                %             end
            end
        end
    end
end

% plot on
if IL
    for dindex=1:numdurs
        figure
        p=0;
        subplot1(numamps,numfreqs)
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                axis off
                
                trace1=squeeze(mM1ON(findex, aindex, dindex, :));
                %             trace1=filtfilt(b,a,trace1);
                trace1=trace1 -mean(trace1(1:100));
                
                
                Stimtrace=squeeze(mM1ONStim(findex, aindex, dindex, :));
                Stimtrace=Stimtrace -mean(Stimtrace(1:100));
                Stimtrace=.05*diff(ylimits)*Stimtrace;
                
                t=1:length(trace1);
                t=1000*t/out.samprate; %convert to ms
                t=t+out.xlimits(1); %correct for xlim in original processing call
                line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
                hold on; plot(t, trace1, 'b');
                offset=ylimits(1)+.1*diff(ylimits);
                plot(t, Stimtrace+offset, 'm')
                
                for rep=1:nrepsON(findex, aindex, dindex)
                    Lasertrace=squeeze(M1ONLaser(findex, aindex, dindex,rep, :));
                    Lasertrace=Lasertrace -mean(Lasertrace(1:100));
                    Lasertrace=.05*diff(ylimits)*Lasertrace;
                    plot( t, Lasertrace+offset, 'c')
                end
                
                ylim(ylimits)
                xlim(xlimits)
                axis off
                
            end
        end
        subplot1(1)
        h=title(sprintf('ON %s: %dms, nreps: %d-%d',datadir,durs(dindex),min(min(min(nrepsON))),max(max(max(nrepsON)))));
        set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
        
        %label amps and freqs
        p=0;
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                if findex==1
                    text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                end
                if aindex==1
                    if mod(findex,2) %odd freq
                        vpos=ylimits(1)-mean(ylimits);
                    else
                        vpos=ylimits(1)-mean(ylimits);
                    end
                    text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                end
                %             if findex==numfreqs && aindex==numamps
                %                 axis on
                %                 ylab=[ceil(ylimits(1)*10)/10 floor(ylimits(2)*10)/10];
                %                 set(gca,'ytick',ylab,'yticklabel',ylab,'YAxisLocation','right')
                %             end
            end
        end
    end
end
