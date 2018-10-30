function PlotGPIAS_Behavior_varlaser(datadir)

%plots GPIAS behavioral data from djmaus and open-ephys accelerometers
%
%usage: PlotGPIAS_Behavior_varlaser(datadir)
%
%modified from PlotGPIAS_Behavior_kip mw 10.29.2018
%this version is for behavioral data with variable laser pulse trains or
%pulse timing
%
%Processes data if outfile is not found;
% adds on GTR and burst data from all cells with outfiles in pwd
% _kip version adds on GTR and burst data from all cells with outfiles in pwd
% should be cleaned up

if nargin==0 datadir=pwd;end

flag.plot = 0;  % plot all startles for each gapdur?

PreStartleWindowms=[-100 0]; % in ms relative to onset of startle pulse
PostStartleWindowms=[0 100]; % in ms relative to onset of startle-pulse
ISIWindowms=[0 60]; % in ms relative to onset of pre-pulse    %added by APW 3_31_14

% ind = strfind(datadir,'\');
% FNtemp = datadir(ind(end)+1:end);
[~,FNtemp,~]=fileparts(datadir);
force_reprocess=1;
if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessGPIAS_BehaviorTiltvarlaser(datadir,4);
end

cd(datadir)
outfilename=sprintf('outGPIAS_Behavior.mat');
d=dir(outfilename);
if ~isempty(d)
    load(outfilename)
else
    ProcessGPIAS_BehaviorTiltvarlaser(datadir)
    load(outfilename);
end



IL=out.IL; %whether there are any interleaved laser trials
numpulseamps=out.numpulseamps;
numgapdurs=out.numgapdurs;
pulseamps=out.pulseamps;
gapdurs=out.gapdurs;
gapdelay=out.gapdelay;
samprate=out.samprate; %in Hz
mM1=out.mM1;
M1=out.M1;
nreps=out.nreps;
soa=out.soa;
isi=out.isi;
soaflag=out.soaflag;
PeakTilt=out.Peak;
mPeak=out.mPeak;
semPeak=out.semPeak;
percentGPIAS=out.percentGPIAS;
pTilt = out.p;
samprate=out.samprate;
M1stim=out.M1stim;
mM1stim=out.mM1stim;
xlimits=out.xlimits;

% %find optimal axis limits
ylimits(1)=min(mM1(:));
ylimits(2)=max(mM1(:));

fs=10;

%plot the tuning curves

% Plot the actual trace with mean trace overlayed
% Separated into 2 figures with laser OFF/ON
if true
    
    %plot the mean response
    figure;
    p=0;
    subplot1(numgapdurs*numpulsewidths*numlaserstarts,1)
    for pwindex=1:numpulsewidths
        for lsindex=1:numlaserstarts
            for gdindex=1:numgapdurs
                p=p+1;
                subplot1(p)
                hold on
                
                offset=10*std(M1(:));
                
                % add the stimulus in magenta
                stimtrace=squeeze(mM1stim(gdindex, pwindex,lsindex,:));
                stimtrace=stimtrace-median(stimtrace(1:1000));
                t=1:length(stimtrace);
                t=1e3*t/samprate;
                t=t+xlimits(1);
                plot(t, stimtrace-offset, 'm')
                
                % plot each trial in blue
                for i=1:nreps(gdindex, pwindex, lsindex)
                    trace1=squeeze(M1(gdindex, pwindex, lsindex,i,:));
                    trace1=trace1-median(trace1(1:1000));
                    plot(t, trace1+i*offset, 'b');
                end
                
                % plot the mean trace in red
                trace1=squeeze(mM1(gdindex, pwindex, lsindex,:));
                trace1=trace1-median(trace1(1:1000));
                
                plot(t, trace1, 'r')
                
                %ylim([ylimits(1)-3*offset ylimits(2)])
                xlim(xlimits)
                ylabel(sprintf('%d ms',gapdurs(gdindex)));
                text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d',nreps(gdindex, pwindex, lsindex)))
                %axis off
                
            end
        end
    end
end
subplot1(1)
h=title(sprintf('%s:\nnreps: %d-%d',FNtemp,min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','interpreter','none')

subplot1(numgapdurs)
xlabel('Time (ms)');


%plot the percent GPIAS
Hfig_percentGPIAS = figure;
hold on
gd=1:numgapdurs;
cmap=copper(numpulsewidths*numlaserstarts);
j=0;
% percentGPIAS(numgapdurs, numpulsewidths, numlaserstarts)
legstr={};
for pwindex=1:numpulsewidths
    for lsindex=1:numlaserstarts
        j=j+1;
        h=plot(gd, squeeze(percentGPIAS(:,pwindex, lsindex)), 'k-o');
        set(h, 'color', cmap(j,:))
        legstr{j}=sprintf('pwd%.1f ls%d');
    end
end
set(gca, 'xtick', 1:numgapdurs)
set(gca, 'xticklabel', gapdurs)
h = title ([FNtemp ' percent GPIAS']);
set(h,'interpreter','none')
xlabel('gap duration')
ylabel('percent GPIAS')
legend(legstr)

%plot the mean peak rectified startle
figure
hold on
j=0;
legstr={};
for pwindex=1:numpulsewidths
    for lsindex=1:numlaserstarts
        j=j+1;
        e=  errorbar(gd, mPeak(:,pwindex, lsindex),semPeak(:,pwindex, lsindex), 'k-o');
        set(e, 'color', cmap(j,:))
        legstr{j}=sprintf('pwd%.1f ls%d');
    end
end
set(gca, 'xtick', 1:numgapdurs)
set(gca, 'xticklabel', gapdurs)
h = title ([FNtemp ' startle']);
set(h,'interpreter','none')
xlabel('gap duration')
ylabel('startle response +- sem')


%plot all trials peak rectified startle
if flag.plot
    figure
    hold on
    j=0;
    legstr={};
    for pwindex=1:numpulsewidths
        for lsindex=1:numlaserstarts
            j=j+1;
            
            h1=errorbar(gd, mPeak(:,pwindex, lsindex),semPeak(:,pwindex, lsindex), 'k-o');
            h2=plot(gd, squeeze(PeakTilt(:,pwindex, lsindex,:)), 'ko')
            set([h1 h2], 'color', cmap(j,:))
            legstr{j}=sprintf('pwd%.1f ls%d');
            
            for igap = 2:numgapdurs
                str = {['GPIAS ' num2str(round(percentGPIAS(igap),1)) ' %']; ['p ~ ' num2str(round(pTilt(igap),3))]};
                text(igap-.25,mPeak(igap,pwindex, lsindex)*.95,str)
            end
            
        end
    end
    set(gca, 'xtick', 1:numgapdurs)
    set(gca, 'xticklabel', gapdurs)
    h = title ([FNtemp ' startle']);
    set(h,'interpreter','none')
    xlabel('gap duration')
    ylabel('startle responses all trials')
    
    
    
    
    
    
