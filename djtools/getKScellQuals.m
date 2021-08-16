function [channel, meanSpikeAmp, stdSpikeAmp, SNR]= getKScellQuals(cellnum);


% read cluster ID and group id. 0= noise, 1= MUA, 2= good (single cell)

sp = loadKSdir(pwd);
cell_ID=sp.cids(cellnum); %get Kilosort id

[spikeAmps, spikeDepths, templateYpos, tempAmps, tempsUnW, tempDur, tempPeakWF] = ...
    templatePositionsAmplitudes(sp.temps, sp.winv, sp.ycoords, sp.spikeTemplates(sp.clu==cell_ID), sp.tempScalingAmps(sp.clu==cell_ID));



% create gwfparams to be able to plot waveforms for one cell
gwfparams.dataDir = pwd;
gwfparams.fileName = '64OE.dat';         % .dat file containing the raw data
gwfparams.dataType = 'int16';            % Data type of .dat file (this should be BP filtered)
gwfparams.nCh = 64;                      % Number of channels that were streamed to disk in .dat file
gwfparams.wfWin = [-40 41];              % Number of samples before and after spiketime to include in waveform
gwfparams.nWf = 2000;                    % Number of waveforms per unit to pull out
gwfparams.spikeTimes =   ceil(sp.st(sp.clu==cell_ID)*30000); % Vector of cluster spike times (in samples) same length as .spikeClusters
gwfparams.spikeClusters = sp.clu(sp.clu==cell_ID); % Vector of cluster IDs (Phy nomenclature)   same length as .spikeTimes

% get waveforms for the cell
wf = getWaveForms(gwfparams);
%OUTPUT:
% wf.unitIDs                               % [nClu,1]            List of cluster IDs; defines order used in all wf.* variables
% wf.spikeTimeKeeps                        % [nClu,nWf]          Which spike times were used for the waveforms
% wf.waveForms                             % [nClu,nWf,nCh,nSWf] Individual waveforms
% wf.waveFormsMean                         % [nClu,nCh,nSWf]     Average of all waveforms (per channel)
% nClu: number of different clusters in .spikeClusters
% nSWf: number of samples per waveform

%find the channel with the largest spike
[r,c]=find(min(min(squeeze(wf.waveFormsMean)))==squeeze(wf.waveFormsMean));
channel= r;

% plot mean waveform
figure; 
imagesc(squeeze(wf.waveFormsMean))
set(gca, 'YDir', 'normal'); xlabel('time (samples)'); ylabel('channel number'); 
colormap(colormap_BlueWhiteRed); caxis([-1 1]*max(abs(caxis()))/2); box off; yticks(1:size(wf.waveForms,3));
title(sprintf(' Cell # %d, KS ID %d, peak spike on Ch%d ', cellnum, cell_ID, channel))
set(gcf, 'Position', [1190 72 562 904])

% plot all traces on peak channel, calculate quality measures
figure; hold on;
for i =1:size(wf.waveForms,2)
    trace=squeeze(wf.waveForms(:,i,channel,:));
    trace=trace-mean(trace(1:10));
    plot(trace,'k')
end
xlabel('time (samples)');
ylabel('Spike Amplitude (\muV)')
meanSpikeAmp=mean(spikeAmps);
stdSpikeAmp=std(spikeAmps);
traces=mean(squeeze(wf.waveForms(:,:,channel,:)),4);
SNR= abs(mean(traces)/std(traces));
title(sprintf('Mean spike amplitude = %.1f, st dev = %.1f, STN ratio = %.3f', meanSpikeAmp, stdSpikeAmp, SNR))



