function PlotTC_2P(varargin)

% plots 2P mesoscope tuning curve data from djmaus
% plots all cells
%
% usage: PlotTC_2P(datapath, [xlimits],[ylimits])
% (xlimits & ylimits are optional)
% datapath should be the path to an outfile (out2P.mat)
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
numisis=out.numisis;
isi=out.isis;
if numisis>1 warning('multiple isis, unsupported case');end
samprate=out.samprate; %in Hz
numframes=out.numframes;
xlimits=out.xlimits;
mM1=out.mM1;
M1=out.M1;
nreps=out.nreps;
if isempty(xlimits) xlimits=out.xlimits; end
fprintf('\nusing xlimits [%d-%d]ms', xlimits(1), xlimits(2))
reps_to_use=[];
fprintf('\ncells sorted by pca separately for each stimulus\n')

%plot the mean dF/F for all cells for each freq/amp 
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            X=squeeze(mM1(findex, aindex, dindex, :, :));
            [pcs, score, latent]=pca(X);
            [~, I]=sort(score(:,1));
            X=X(I,:); %sort by pc1
            %note that each freq/amp combo is pca-sorted separately
            %this means the cell order is different for each combo
            %and the pc1 it is sorted by is also different
            if unique(X(:))==0
                %empty matrix, do nothing
                axis off
            else
                imagesc(X)

                xticks([(xlimits(1)/1000):samprate:numframes])
                xticklabels([(xlimits(1)/1000):1:numframes/samprate])
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


            end
        end
    end
end
[~,fname, ~]=fileparts(datadir);
h=suptitle(sprintf('%s\n%dms, nreps: %d-%d',fname,durs(dindex),min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')

%plot the population mean dF/F (averaged across all cells) for each freq/amp 
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
                t=1:length(nanmean(trace1));
                t=t/samprate;
                t=t+xlimits(1)/1000;
                plot(t, nanmean(trace1))

                xticks([(xlimits(1)/1000):1:xlimits(2)])
                
                line([0 0], ylim, 'linestyle', '--')
                line((isi/1000)+[0 0], ylim, 'linestyle', '--')
                line((2*isi/1000)+[0 0], ylim, 'linestyle', '--')


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


            end
        end
    end
end
[~,fname, ~]=fileparts(datadir);
h=suptitle(sprintf('population mean (all cells)%s\n%dms, nreps: %d-%d',fname,durs(dindex),min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')

