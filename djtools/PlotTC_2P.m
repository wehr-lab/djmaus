function PlotTC_2P(varargin)

% plots 2P mesoscope tuning curve data from djmaus
% plots all cells
%
% usage: PlotTC_2P([datapath], [xlimits],[ylimits])
% (xlimits & ylimits are optional)
% datapath should be the path to an outfile (out2P.mat)
%
% Processes data if outfile is not found;

force_reprocess=0;

if nargin==0
    if exist('./out2P.mat')==2
        datadir=pwd;
        outfilename='out2P.mat';
    elseif exist('./Fall.mat')==2 & exist('../../out2P.mat')==2
        %we're in suite2p/plane0, go up
        cd ../..
        datadir=pwd;
        outfilename='out2P.mat';
    else
        datadir=pwd;
        fprintf('\nno outfile found...')
        outfilename=sprintf('out2P.mat');
        % %select outfile with GUI
        % [outfilename, datadir]=uigetfile('o*.mat', 'select outfile to plot');
        % if outfilename==0
        %    fprintf('\nuser pressed cancel')
        %    return
        % end
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
mM1dff=out.mM1dff;
mM1f=out.mM1f;
M1f=out.M1f;
M1dff=out.M1dff;
nreps=out.nreps;
if isempty(xlimits) xlimits=out.xlimits; end
fprintf('\nusing xlimits [%d-%d]ms', xlimits(1), xlimits(2))
reps_to_use=[];
fprintf('\ncells sorted by pca separately for each stimulus\n')

%sort by WN pca loading
X=squeeze(mM1dff(1, 2, 1, :, :)); 
            %f0=median(X(:,1:15), 2); %compute f0, assuming mM1 is raw F
            %X=(X-f0)./f0;
             [pcs, score, latent]=pca(X);
                [~, Iwn]=sort(score(:,1));
                %X=X(I,:); %sort by pc1


%plot the mean dF/F for all cells for each freq/amp 
for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            X=squeeze(mM1f(findex, aindex, dindex, :, :));
            f0=median(X(:,1:15), 2); %compute f0, assuming mM1 is raw F
            X=(X-f0)./f0;

          
            if 0
                %sort by pca
                [pcs, score, latent]=pca(X);
                [~, I]=sort(score(:,1));
                X=X(I,:); %sort by pc1
                %note that each freq/amp combo is pca-sorted separately
                %this means the cell order is different for each combo
                %and the pc1 it is sorted by is also different
            else
                %sort by tone response
                [~, I]=sort(mean(X(:, 18:22), 2)); %tone response seems to be around frames 18-22
                X=X(I,:);
            end

            if unique(X(:))==0
                %empty matrix, do nothing
                axis off
            else
                imagesc(X)
                cl=clim;
                caxis([0 .25*cl(2)]) %turn up gain on colormap

                xticks([(xlimits(1)/1000):samprate:numframes])
                xticklabels([(xlimits(1)/1000):1:numframes/samprate])
                if aindex==1 xlabel('time, s');end
                if findex==1 ylabel('cell # (sorted separately in each panel)');end

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
h=suptitle(sprintf('%s\n%dms, nreps: %d-%d',fname,1000*durs(dindex),min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
pos=get(gcf, 'pos'); pos(3)=1200; pos(4)=1000;
set(gcf, 'pos', pos)

%plot dF/F tuning curve for each cell 
%blue is louder, green is quiet
%mean reponse for all freqs are concatenated
dindex=1;
numcells=length(Iwn);
rootnumcells=ceil(sqrt(numcells));
figure
subplot1(rootnumcells,rootnumcells)
for p=1:numcells
    subplot1(p)
    hold on
    aindex=1;
    x1=squeeze(mM1f(:, aindex, dindex, Iwn(p), 1:45));
    x1=[x1 nan(6,15)];
    x1row=reshape(x1', 1, prod(size(x1)));
    plot(x1row, 'b')
    aindex=2;
    x2=squeeze(mM1f(:, aindex, dindex, Iwn(p), 1:45));
    x2=[x2 nan(6,15)];
    x2row=reshape(x2', 1, prod(size(x2)));
    plot(x2row, 'color', [0 .65 0])
    axis off
    ylim([0 max(mM1f(:))])
    text(1, 1, int2str(Iwn(p)))
end
for p=numcells+1:rootnumcells^2
    subplot1(p)
    axis off
end

[~,fname, ~]=fileparts(datadir);
h=suptitle(sprintf('%s\n%dms, nreps: %d-%d',fname,1000*durs(dindex),min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
pos=get(gcf, 'pos'); pos(3)=1200; pos(4)=1000;
set(gcf, 'pos', pos)

%plot the population mean dF/F (averaged across all cells) for each freq/amp 

%get ylimits
if isempty(ylimits)
    ylimits=[0 0];
    for dindex=1:numdurs
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                % trace1=nanmean(squeeze(mM1(findex, aindex, dindex, :, :)));
                X=squeeze(mM1f(findex, aindex, dindex, :, :));
                f0=mean(X(:,1:15), 2); %compute f0, assuming mM1 is raw F
                X=(X-f0)./f0;
                trace1=nanmean(X);

                if min(trace1)<ylimits(1) ylimits(1)=min(trace1); end
                if max(trace1)>ylimits(2) ylimits(2)=max(trace1); end
            end
        end
    end
end


for dindex=1:numdurs
    figure
    p=0;
    subplot1(numamps,numfreqs)
    for aindex=numamps:-1:1
        for findex=1:numfreqs
            p=p+1;
            subplot1(p)
            X=squeeze(mM1f(findex, aindex, dindex, :, :));
            f0=mean(X(:,1:15), 2); %compute f0, assuming mM1 is raw F
            X=(X-f0)./f0;
            trace1=nanmean(X);
            if unique(trace1(:))==0
                %empty matrix, do nothing
                axis off
            else
                t=1:length(trace1);
                t=t/samprate;
                t=t+xlimits(1)/1000;
                plot(t, trace1)

                xticks([(xlimits(1)/1000):1:xlimits(2)])
                ylim(ylimits)

                line([0 0], ylim, 'linestyle', '--')
                line((isi/1000)+[0 0], ylim, 'linestyle', '--')
                line((2*isi/1000)+[0 0], ylim, 'linestyle', '--')


                if aindex==1 xlabel('time, s');end
                if findex==1 ylabel('mean dff');end


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
h=suptitle(sprintf('population mean (all cells)%s\n%dms, nreps: %d-%d',fname,1000*durs(dindex),min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'left', 'interpreter', 'none')
pos=get(gcf, 'pos'); pos(3)=2000;
set(gcf, 'pos', pos)

% plot tonotopy
% assume tones of a single duration and amplitude (should check this)
aindex=1;
dindex=1;
f0=squeeze(mM1f(:, aindex, dindex, :, 1:15));
f0=mean(f0, 3);
f=squeeze(mM1(:, aindex, dindex, :, 18:22));
f=mean(f, 3);
dff=(f-f0)./f0;


[mag, CFidx]=(max(dff));
mag=mag-min(mag);
mag=mag./max(mag);
mag=mag.^.5;

for i=1:length(out.iscell)
    fprintf('\n CF %d %.1fkHz x %d y %d', CFidx(i), freqs(CFidx(i))/1000, out.stat{i}.xpix(1), out.stat{i}.ypix(1))
end
cmap=jet(numfreqs);
figure
hold on
for i=1:length(out.iscell)
    h=plot(out.stat{i}.xpix(1), -out.stat{i}.ypix(1), '.');
    set(h, 'Color', mag(i)*cmap(CFidx(i),:), 'MarkerSize', 20)
end
xlabel('x position, in pixels')
ylabel('y position, in pixels')
title(sprintf('%s tonotopic map with color-coded best frequency of cells', fname))

          
figure('position',[1561  817 222  420])
hold on
for i=1:length(freqs)
    h=plot(1,i, '.');
    set(h, 'Color', cmap(i,:), 'MarkerSize', 20)
    text(1.2, i, sprintf('%.1f kHz',freqs(i)/1000))
end
title('frequency color code ')







