function PlotTC_2P(varargin)




% plots 2P mesoscope tuning curve data from djmaus
%plots all cells
%
% usage: PlotTC_2P(datapath, [xlimits],[ylimits])
% (xlimits & ylimits are optional)
% channel number should be an integer
%
% Processes data if outfile is not found;

force_reprocess=0;

if nargin==0
    if exist('./out2P.mat')==2
        datadir=pwd;
        outfilename='out2P.mat';
    else
        %select outfile with GUI
        [outfilename, datadir]=uigetfile('o*.mat', 'select outfile to plot');
    end
    ylimits=[];
    xlimits=[];
else
    datadir=varargin{1};


    outfilename=sprintf('out2P.mat');
    try
        xlimits=varargin{2};
    catch
        xlimits=[];
    end
    try
        ylimits=varargin{3};
    catch
        ylimits=[];
    end
end



if force_reprocess
    fprintf('\nForce Re-process')
    fprintf('\ncalling ProcessTC_2P')
    ProcessTC_2P(datadir, xlimits, ylimits);
end

cd(datadir)
d=dir(outfilename);
if ~isempty(d)
    fprintf('\nloading outfile %s ...', outfilename);
    load(outfilename)
    fprintf('\tdone');
else
    fprintf('\ncalling ProcessTC_2P')
    ProcessTC_2P(datadir, xlimits, ylimits);
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
numframes=out.numframes;
xlimits=out.xlimits;
% M1=out.M1;
% scaledtrace=out.scaledtrace;
% traces_to_keep=out.traces_to_keep;
% mM1=out.mM1;
mM1=out.mM1;
M1=out.M1;
% mM1ONStim=out.mM1ONStim;
% mM1OFFStim=out.mM1OFFStim;
% nrepsON=out.nrepsON;
nreps=out.nreps;
% IL=out.IL; %whether there were interleaved laser trials or not
if isempty(xlimits) xlimits=out.xlimits; end
fprintf('\nusing xlimits [%d-%d]\n', xlimits(1), xlimits(2))

% %find optimal axis limits
% if isempty(ylimits)
%     ylimits=[0 0];
%     for dindex=1:numdurs
%         for aindex=numamps:-1:1
%             for findex=1:numfreqs
%                 trace1=squeeze(mM1(findex, aindex, dindex, :));
%                 if min([trace1])<ylimits(1); ylimits(1)=min(trace1);end
%                 if max([trace1])>ylimits(2); ylimits(2)=max(trace1);end
%             end
%         end
%     end
% end

% ylimits=round(ylimits*100)/100;



%samprate=30e3;
reps_to_use=[];


%plot the mean tuning curve
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            trace1=squeeze(mM1(findex, aindex, dindex, :, :));
            if unique(trace1(:))==0
                %empty matrix, do nothing
                axis off
            else
                imagesc(trace1)
              

xticks([0:samprate:numframes])
xticklabels([0:1:numframes/samprate])
if aindex==1 xlabel('time, s');end
if findex==1 ylabel('cell #');end

                %label amps and freqs

                if freqs(findex)==-2000 & aindex==1
                    title( 'SS')
                elseif freqs(findex)==-1000
                    title(sprintf('WN %ddB', amps(aindex)))
                elseif freqs(findex)>0
                    title( sprintf('%.1fkHz %ddB', freqs(findex)/1000, amps(aindex)))
                end

                %             subplot1(p)
                %             if findex==1
                %                 %text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
                %                 text(xlimits(1), mean(ylimits), int2str(amps(aindex)))
                %             end
                %mM1ONLaser=mM1ONLaser;
                %axis off
            end
        end
    end
end
[~,fname, ~]=fileparts(datadir);
h=suptitle(sprintf('%s\n%dms, nreps: %d-%d',fname,durs(dindex),min(nreps(:)),max(nreps(:))));
%h=title(sprintf('OFF %s: %dms, nreps: %d-%d',datadir,durs(dindex), reps_to_use));
set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
set(h,  'interpreter', 'none')

%label amps and freqs
%     p=0;
%     for aindex=numamps:-1:1
%         for findex=1:numfreqs
%             p=p+1;
%             subplot1(p)
%             if findex==1
%                 text(xlimits(1)-diff(xlimits)/2, mean(ylimits), int2str(amps(aindex)))
%             end
%
%             if aindex==1
%                 if mod(findex,2) %odd freq
%                     vpos=ylimits(1)-.2*mean(ylimits);
%                 else
%                     vpos=ylimits(1)-.2*mean(ylimits);
%                 end
%                 if freqs(findex)==-2000
%                     text(xlimits(1), vpos, 'SS')
%                 elseif freqs(findex)==-1000
%                     text(xlimits(1), vpos, 'WN')
%                 elseif freqs(findex)>0
%                     text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
%                 end
%             else
%                 if aindex==1
%                     if mod(findex,2) %odd freq
%                         vpos=ylimits(1)-mean(ylimits);
%                     else
%                         vpos=ylimits(1)-mean(ylimits);
%                     end
%                     if freqs(findex)>0
%                         text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
%                     elseif freqs(findex)==-1000
%                         text(xlimits(1), vpos, 'WN')
%                     elseif freqs(findex)==-2000
%                         text(xlimits(1), vpos, 'SS')
%                     end
%                 end
%
%                 %                 if freqs(findex)>0
%                 %                     text(xlimits(1), vpos, sprintf('%.1f', freqs(findex)/1000))
%                 %                 elseif freqs(findex)==-1000
%                 %                     text(xlimits(1), vpos, 'WN')
%                 %                 elseif freqs(findex)==-2000
%                 %                     text(xlimits(1), vpos, 'SS')
%                 %                 end
%             end
%             %             if findex==numfreqs && aindex==numamps
%             %                 axis on
%             %                 ylab=[ceil(ylimits(1)*10)/10 floor(ylimits(2)*10)/10];
%             %                 set(gca,'ytick',ylab,'yticklabel',ylab,'YAxisLocation','right')
%             %             end
%         end
%     end
%     subplot1(1)
%     h=title(sprintf('OFF %s: %dms, nreps: %d-%d, ch%d',datadir,durs(dindex),min(min(min(nrepsOFF))),max(max(max(nrepsOFF))), channel));
%     set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
% end


