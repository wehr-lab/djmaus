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
force_reprocess=0;
if force_reprocess
    fprintf('\nForce re-process\n')
    ProcessGPIAS_BehaviorTiltvarlaser(datadir,4);
end

cd(datadir)
outfilename=sprintf('outGPIAS_Behavior.mat');
d=dir(outfilename);
if ~isempty(d)
    fprintf('\nloading outfile...')
    load(outfilename)
else
    fprintf('\noutfile not found, calling ProcessGPIAS_BehaviorTiltvarlaser\n')
    ProcessGPIAS_BehaviorTiltvarlaser(datadir)
    load(outfilename);
end



numpulseamps=out.numpulseamps;
numpulsewidths=out.numpulsewidths;
numlaserstarts=out.numlaserstarts;
numgapdurs=out.numgapdurs;
pulseamps=out.pulseamps;
pulsewidths=out.pulsewidths;
laserstarts=out.laserstarts;
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

fprintf('\n%d gapdurs: ', numgapdurs)
fprintf('%d ', gapdurs)

fprintf('\nLaser conditions: ')
if out.LaserNumPulses==1
    fprintf('\nlaser is a single pulse, not a train')
else
    fprintf('\nlaser flashtrain')
end
fprintf('\n%d laser pulse widths (pw): ', numpulsewidths)
fprintf('%g ', pulsewidths)
fprintf('\n%d laser start times (ls): ', numlaserstarts)
fprintf('%d ', laserstarts)



%plot the tuning curves

% Plot the actual trace with mean trace overlayed
if true
    
    %plot the mean response
    figure;
    p=0;
    numconditions=length(find(nreps));
    %numconditions is less than numgapdurs*numpulsewidths*numlaserstarts
    %because there are some empty conditions
    subplot1(numconditions,1)
    for pwindex=1:numpulsewidths
        for lsindex=1:numlaserstarts
            for gdindex=1:numgapdurs
                
                if nreps(gdindex, pwindex,lsindex)>0 %skip empty conditions
                    
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
                text(xlimits(1)+10,ylimits(2)/2,sprintf('n=%d pw=%g ls=%d',nreps(gdindex, pwindex, lsindex), pulsewidths(pwindex), laserstarts(lsindex)), 'color', 'r', 'fontsize', 14)
                %axis off
                end
            end
        end
    end
end
subplot1(1)
h=title(sprintf('%s:\nnreps: %d-%d',FNtemp,min(nreps(:)),max(nreps(:))));
set(h, 'HorizontalAlignment', 'center', 'interpreter', 'none', 'fontsize', fs, 'fontw', 'normal','interpreter','none')

subplot1(numgapdurs)
xlabel('Time (ms)');

%depending on whether we're varying the flashtrain pulsewidth, or the laser
%start time, we will get nans in one of the columns.
% There isn't a full M x N matrix

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
        pg=squeeze(percentGPIAS(:,pwindex, lsindex));
        if ~all(isnan(pg))
            j=j+1;
            h=plot(gd, pg, 'k-o');
            set(h, 'color', cmap(j,:))
            legstr{j}=sprintf('pw%g ls%d', pulsewidths(pwindex), laserstarts(lsindex));
        end
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
        mp=mPeak(:,pwindex, lsindex);
        semp=semPeak(:,pwindex, lsindex);
        if all(~isnan(mp)) 
            j=j+1;
            e=  errorbar(gd, mp,semp, 'k-o');
            set(e, 'color', cmap(j,:))
            legstr{j}=sprintf('pw%.1f ls%d', pulsewidths(pwindex), laserstarts(lsindex));
        end
    end
end
set(gca, 'xtick', 1:numgapdurs)
set(gca, 'xticklabel', gapdurs)
h = title ([FNtemp ' startle']);
set(h,'interpreter','none')
xlabel('gap duration')
ylabel('startle response +- sem')
legend(legstr)


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
        legstr{j}=sprintf('pw%.1f ls%d', pwindex, lsindex);
            
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
    legend(legstr)
end
    
    
    
    
