function PlotTC_2P(varargin)




% plots 2P mesoscope tuning curve data from djmaus
%
% usage: PlotTC_LFP(datapath, [channel], [xlimits],[ylimits])
% (xlimits & ylimits are optional)
% channel number should be an integer
%
% Processes data if outfile is not found;

force_reprocess=0;

if nargin==0
    %select outfile with GUI
    [outfilename, datadir]=uigetfile('o*.mat', 'select outfile to plot');
    
    ylimits=[];
    xlimits=[];
else
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
    outfilename=sprintf('outLFP_ch%d.mat',channel);
    try
        xlimits=varargin{3};
    catch
        xlimits=[-20 20];
    end
    try
        ylimits=varargin{4};
    catch
        ylimits=[];
    end
end



if force_reprocess
    fprintf('\nForce Re-process')
    fprintf('\ncalling ProcessTC_2P')
    ProcessTC_2P(datadir,  channel, xlimits, ylimits);
end

cd(datadir)
d=dir(outfilename);
if ~isempty(d)
    fprintf('\nloading outfile %s ...', outfilename);
    load(outfilename)
    fprintf('\tdone');
else
    fprintf('\ncalling ProcessTC_2P')
    ProcessTC_2P(datadir,  channel, xlimits, ylimits);
    load(outfilename);
end

% M1stim=out.M1stim;

freqs=out.freqs;
amps=out.amps;
durs=out.durs;
nreps=out.nreps;
numfreqs=out.numfreqs;
numamps=out.numamps;
numdurs=out.numdurs;
samprate=out.samprate; %in Hz
% M1=out.M1;
% scaledtrace=out.scaledtrace;
% traces_to_keep=out.traces_to_keep;
% mM1=out.mM1;
mM1OFF=out.mM1OFF;
M1OFF=out.M1OFF;
% mM1ONStim=out.mM1ONStim;
% mM1OFFStim=out.mM1OFFStim;
% nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;
% IL=out.IL; %whether there were interleaved laser trials or not
if isempty(xlimits) xlimits=out.xlimits; end
fprintf('\nusing xlimits [%d-%d]\n', xlimits(1), xlimits(2))

% %find optimal axis limits
if isempty(ylimits)
    ylimits=[0 0];
    for dindex=1:numdurs
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                trace1=squeeze(mM1OFF(findex, aindex, dindex, :));
                if min([trace1])<ylimits(1); ylimits(1)=min(trace1);end
                if max([trace1])>ylimits(2); ylimits(2)=max(trace1);end
            end
        end
    end
end

% ylimits=round(ylimits*100)/100;



%samprate=30e3;
reps_to_use=[];


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
            %            [b,a]=butter(1, low_pass_cutoff/(samprate/2), 'low');
            %             reps_to_use=1000;
            %             trace1=squeeze(mean(M1OFF(findex, aindex, dindex, 1:reps_to_use,:), 4));
           
%             trace1=trace1 -mean(trace1(1:100));
           
%             if ~isempty(mM1OFFStim)
%                 Stimtrace=squeeze(mM1OFFStim(findex, aindex, dindex, :));
%                 Stimtrace=Stimtrace -mean(Stimtrace(1:100));
%                 Stimtrace=.05*diff(ylimits)*Stimtrace;
%             else
%                 Stimtrace = trace1*0;
%             end
            t=1:length(trace1);
            t=1000*t/samprate; %convert to ms
            t=t+xlimits(1); %correct for xlim in original processing call
            line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            hold on; plot(t, trace1, 'k');
            offset=ylimits(1)+.1*diff(ylimits);
%             plot(t, Stimtrace+offset, 'm', t, Lasertrace+offset, 'c')
            try
                ylim(ylimits)
            end
            xlim(xlimits)
            xlabel off
            ylabel off
            
            %label amps and freqs
            
            subplot1(p)
            if findex==1
                text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
            end
            %mM1ONLaser=mM1ONLaser;
            axis off
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s: %dms, nreps: %d-%d',datadir,durs(dindex),min(nrepsOFF(:)),max(nrepsOFF(:))));
    %h=title(sprintf('OFF %s: %dms, nreps: %d-%d',datadir,durs(dindex), reps_to_use));
    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
    set(h,  'interpreter', 'none')
    
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
                if freqs(findex)==-2000
                    text(xlimits(1), vpos, 'SS')
                elseif freqs(findex)==-1000
                    text(xlimits(1), vpos, 'WN')
                elseif freqs(findex)>0
                    text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                end
            else
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
                
                %                 if freqs(findex)>0
                %                     text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
                %                 elseif freqs(findex)==-1000
                %                     text(xlimits(1), vpos, 'WN')
                %                 elseif freqs(findex)==-2000
                %                     text(xlimits(1), vpos, 'SS')
                %                 end
            end
            %             if findex==numfreqs && aindex==numamps
            %                 axis on
            %                 ylab=[ceil(ylimits(1)*10)/10 floor(ylimits(2)*10)/10];
            %                 set(gca,'ytick',ylab,'yticklabel',ylab,'YAxisLocation','right')
            %             end
        end
    end
    subplot1(1)
    h=title(sprintf('OFF %s: %dms, nreps: %d-%d, ch%d',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF))), channel));
    set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
end


