function ProcessPINP_PSTH2(varargin)

%processes kilosorted spiking data from a djmaus PINP protocol
%using new OpenEphys and kilosort file formats and hierarchy
%creates only a single outfile per session
%(automatically processes all cells)
%-mike 04.28.2024
%
% usage: ProcessPINP_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, [xlimits],[ylimits])
%   SortedUnits is created by ProcessSpikes
% BonsaiPath, EphysDir, EphysDir_KS are folders, they are stored in OEinfo-xxx.mat in the Bonsai directory
% (xlimits, ylimits are optional)
% xlimits default to [0 200]
%
% saves to outfile
%
%customized to process a PINP protocol. It will only look for  (1) silent sound, (2) silent
%sound with a single laser pulse, and (3) silent sound with laser flash trains

if nargin==0
    fprintf('\nno input');
    return;
end
SortedUnits=varargin{1};
BonsaiPath=varargin{2};
EphysPath=varargin{3};
EphysPath_KS=varargin{4};

xlimits=[];
try
    xlimits=varargin{5};
end
if isempty(xlimits)
    xlimits=[-100 200];
    s=GetStimParams(fullfile(BonsaiPath, EphysPath));
    if isempty(s) %in new OE format, notebook is one up from kilosorted data, so fileparts looks in ..
        s=GetStimParams(fileparts(datadir));
    end
    durs=s.durs;
    dur=max(durs);
    xlimits=[-.5*dur 1.5*dur]; %default x limits for axis
end
try
    ylimits=varargin{6};
catch
    ylimits=[];
end
fprintf('\nprocessing with xlimits [%.1f-%.1f]', xlimits(1), xlimits(2))


try
    cd(BonsaiPath)
    cd(EphysPath)
    load notebook.mat
    %check if nb and stimlog are actually there
    if ~exist('stimlog','var')
        stimlog=[];
        fprintf('\nfound notebook file but there was no stimlog in it!!!');
        fprintf('\n in principle it should be possible to reconstruct all stimuli\nfrom network events, but that would take a lot of coding, does this problem happen often enough to be worth the effort???');

        error('found notebook file but there was no stimlog in it!!!');
    end
    if ~exist('nb','var')
        nb=[];
        warning('found notebook file but there was no nb notebook information in it!!!');
    end
catch
    warning('could not find notebook file')
end

%read messages
cd(BonsaiPath)
behavior_filename=dir('Behavior_*.mat');
load(behavior_filename.name);


%get Events and soundcard trigger timestamps
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps

if exist('Events.mat')
    load('Events.mat')
    fprintf('loaded Events file \n')
else
    [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps2(Sky);
    save('Events.mat','Events')
    save('StartAcquisitionSec.mat','StartAcquisitionSec')
end
if exist('StartAcquisitionSec.mat')
    load('StartAcquisitionSec.mat')
else
    [~, StartAcquisitionSec] = GetEventsAndSCT_Timestamps2(Sky);
    save('StartAcquisitionSec.mat','StartAcquisitionSec')
end

try
    fprintf('\nNumber of logged stimuli in notebook: %d', length(stimlog));
catch
    fprintf('\nCould not find stimlog, no logged stimuli in notebook!!');
end


%check if this is an appropriate stimulus protocol
switch (GetPlottingFunction(fullfile(BonsaiPath, EphysPath)))
    case 'PlotPINP_PSTH'
        %fine
    otherwise
        error('This stimulus protcol does not appear to be a PINP protocol.')
        %you could add an exception here
end

fprintf('\ncomputing tuning curve...');

samprate=Sky.OEsamplerate;

%get freqs/amps
j=0;
alltrainnumpulses=[];
alltrainisis=[];
alltrainpulsewidths=[];
allpulsewidths=[];
allsilentsounddurs=[];
if exist('Events.mat')
    load('Events.mat')
else
    save('Events.mat','Events')
end

%apparently we sometimes call VarLaser params VL params instead, so copy
%them to old names

for i=1:length(Events)
    if  strcmp(Events(i).type, 'silentsound')
        if isfield(Events(i), 'VL')
            Events(i).VarLaser=Events(i).VL;
            Events(i).VarLaserstart=Events(i).VLstart;
            Events(i).VarLasernumpulses=Events(i).VLnumpulses;
            Events(i).VarLaserpulsewidth=Events(i).VLpulsewidth;
            Events(i).VarLaserisi=Events(i).VLisi;
        end
    end
end

for i=1:length(Events)
    if  strcmp(Events(i).type, 'silentsound')
        j=j+1;
        alldurs(j)=Events(i).duration;
        allnexts(j)=Events(i).next;
        %I should probably check that these fields exist before capturing them
        alllasers(j)=Events(i).laser;
        allVarLasers(j)=Events(i).VarLaser;
        allVarLaserstarts(j)=Events(i).VarLaserstart;
        if Events(i).VarLaser ==0 & Events(i).laser ==0 %silent sound no laser
            allsilentsounddurs=[allsilentsounddurs Events(i).duration];
        end
        if Events(i).VarLasernumpulses ==1 %single pulse
            allpulsewidths=[ allpulsewidths Events(i).VarLaserpulsewidth];
        elseif Events(i).VarLasernumpulses>1 %train
            alltrainpulsewidths=[alltrainpulsewidths Events(i).VarLaserpulsewidth];
            alltrainisis=[alltrainisis Events(i).VarLaserisi]; %isi in train
            alltrainnumpulses=[alltrainnumpulses Events(i).VarLasernumpulses]; %isi in train
        end

    end
end

lasers=unique(alllasers);
laserstarts=unique(allVarLaserstarts);
trainnumpulses=unique(alltrainnumpulses);
trainpulsewidths=unique(alltrainpulsewidths);
trainisis=unique(alltrainisis);
pulsewidths=unique(allpulsewidths);
silentsounddurs=unique(allsilentsounddurs);

VarLaserstarts=unique(allVarLaserstarts);
durs=unique(alldurs);
nexts=unique(allnexts);
numnexts=length(nexts);

numlaserstarts=length(laserstarts);
numtrainnumpulses=length(trainnumpulses);
numtrainpulsewidths=length(trainpulsewidths);
numtrainisis=length(trainisis);
numpulsewidths=length(pulsewidths);
numdurs=length(durs);
numsilentsounddurs=length(silentsounddurs);

%try to load laser and stimulus monitor
warning('still need to figure out how to load laser and stimulus traces')
if 1%isempty(getLaserfile('.'))
    LaserRecorded=0;
    warning('Laser monitor channel not recorded')
else
    LaserRecorded=1;
end
if 1%isempty(getStimfile('.'))
    StimRecorded=0;
    warning('Stimulus monitor channel not recorded')
else
    StimRecorded=1;
end

if LaserRecorded
    try
        [Lasertrace, Lasertimestamps, Laserinfo] =load_open_ephys_data(getLaserfile('.'));
        Lasertimestamps=Lasertimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
        Lasertrace=Lasertrace./max(abs(Lasertrace));
        fprintf('\nsuccessfully loaded laser trace')
    catch
        fprintf('\nfound laser file %s but could not load laser trace', getLaserfile('.'))
    end
else
    fprintf('\nLaser trace not recorded')
end
if StimRecorded
    try
        [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(getStimfile('.'));
        Stimtimestamps=Stimtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
        Stimtrace=Stimtrace./max(abs(Stimtrace));
        fprintf('\nsuccessfully loaded stim trace')
    catch
        fprintf('\nfound stim file %s but could not load stim trace', getStimfile('.'))
    end
else
    fprintf('\nSound stimulus trace not recorded')
end

numcells=length(SortedUnits);

nrepsPulse=zeros(numcells, numpulsewidths);
nrepsTrain=zeros(numcells, numtrainnumpulses, numtrainpulsewidths, numtrainisis);
nrepsOFF=zeros(numcells, numsilentsounddurs);

%extract the traces into a big matrix M
j=0;
inRange=0;

MPulse=[];
MPulseLasertrace=[];
MPulseStimtrace=[];
MTrain=[];
MTrainLasertrace=[];
MTrainStimtrace=[];
MSilentSoundOFF=[];
MSilentSoundOFFLasertrace=[];
MSilentSoundOFFStimtrace=[];
mMSilentSoundOFFLasertrace=[];
mMSilentSoundOFFStimtrace=[];
mMPulseLasertrace=[];
mMPulseStimtrace=[];
mMTrainLasertrace=[];
mMTrainStimtrace=[];
mMSilentSoundOFF=[];
mMPulse=[];
mMTrain=[];

for cellnum=1:length(SortedUnits);
    spiketimes=[];
    fprintf('\nreading SortedUnits cell %d', cellnum)
    spiketimes= SortedUnits(cellnum).spiketimes;

    totalnumspikes=length(spiketimes);
    fprintf('\nsuccessfully loaded spike data with %d spikes\n',totalnumspikes)


    for i=1:length(Events)
        if strcmp(Events(i).type, 'silentsound')
            if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
                pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
            else
                error('???')
            end
            start=(pos+xlimits(1)*1e-3); %in seconds
            stop=(pos+xlimits(2)*1e-3);
            region=round(start*samprate)+1:round(stop*samprate);
            if start>0 %(disallow negative or zero start times)


                st=spiketimes; %are in seconds
                spiketimes1=st(st>start & st<stop); % spiketimes in region
                spikecount=length(spiketimes1); % No. of spikes fired in response to this rep of this stim.
                inRange=inRange+ spikecount; %accumulate total spikecount in region
                spiketimes1=(spiketimes1-pos)*1000;%covert to ms after tone onset
                spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window

                if Events(i).VarLaser ==0 & Events(i).laser ==0 %silent sound no laser
                    silentsounddur = Events(i).duration;
                    dindex= find(silentsounddur==silentsounddurs);
                    nrepsOFF(cellnum, dindex)=nrepsOFF(cellnum, dindex)+1;
                    MSilentSoundOFF(cellnum, dindex, nrepsOFF(cellnum, dindex)).spiketimes=spiketimes1;
                    if LaserRecorded
                        MSilentSoundOFFLasertrace(dindex, nrepsOFF(dindex),:)=Lasertrace(region);
                    else
                        MSilentSoundOFFLasertrace=[];
                    end
                    if StimRecorded
                        MSilentSoundOFFStimtrace(dindex, nrepsOFF(dindex),:)=Stimtrace(region);
                    else
                        MSilentSoundOFFStimtrace=[];
                    end
                end
                if Events(i).VarLasernumpulses ==1 %single pulse
                    pulsewidth= Events(i).VarLaserpulsewidth;
                    pwindex= find(pulsewidth==pulsewidths);
                    nrepsPulse(cellnum, pwindex)=nrepsPulse(cellnum, pwindex)+1;
                    MPulse(cellnum, pwindex, nrepsPulse(cellnum, pwindex)).spiketimes=spiketimes1;
                    if LaserRecorded
                        MPulseLasertrace(pwindex,nrepsPulse(pwindex),:)=Lasertrace(region);
                    else
                        MPulseLasertrace=[];
                    end
                    if StimRecorded
                        MPulseStimtrace(pwindex,nrepsPulse(pwindex),:)=Stimtrace(region);
                    else
                        MPulseStimtrace=[];
                    end


                elseif Events(i).VarLasernumpulses>1 %train
                    thistrainnumpulses= Events(i).VarLasernumpulses; %numpulses in train
                    trainpulsewidth= Events(i).VarLaserpulsewidth; %pulse width in train
                    trainisi= Events(i).VarLaserisi; %isi in train
                    tnpindex= find(thistrainnumpulses==trainnumpulses);
                    tpwindex= find(trainpulsewidth==trainpulsewidths);
                    tiindex= find(trainisi==trainisis);
                    nrepsTrain(cellnum, tnpindex, tpwindex, tiindex)=nrepsTrain(cellnum, tnpindex,tpwindex, tiindex)+1;
                    MTrain(cellnum, tnpindex, tpwindex, tiindex, nrepsTrain(cellnum, tnpindex, tpwindex, tiindex)).spiketimes=spiketimes1;
                    if LaserRecorded
                        MTrainLasertrace(tnpindex,tpwindex,tiindex, nrepsTrain(tnpindex, tpwindex, tiindex),:)=Lasertrace(region);
                    else
                        MTrainLasertrace=[];
                    end
                    if StimRecorded
                        MTrainStimtrace(tnpindex,tpwindex,tiindex, nrepsTrain(tnpindex, tpwindex, tiindex),:)=Stimtrace(region);
                    else
                        MTrainStimtrace=[];
                    end
                end
            end
        end
    end


  


    % Accumulate spiketimes across trials, for psth...
    % ... and average laser and stimulus monitor M matrices across trials, if recorded


    %silent sound no laser
    for dindex=1:numsilentsounddurs
        spiketimesOFF=[];
        for rep=1:nrepsOFF(cellnum, dindex)
            spiketimesOFF=[spiketimesOFF MSilentSoundOFF(cellnum, dindex, rep).spiketimes];
        end
        mMSilentSoundOFF(cellnum, dindex).spiketimes=spiketimesOFF;
        if LaserRecorded
            mMSilentSoundOFFLasertrace(dindex,:)=mean(MSilentSoundOFFLasertrace(dindex, 1:nrepsOFF(cellnum, dindex),:), 2);
        end
        if StimRecorded
            mMSilentSoundOFFStimtrace(dindex,:)=mean(MSilentSoundOFFStimtrace(dindex, 1:nrepsOFF(cellnum, dindex),:), 2);
        end
    end

    %single laser pulse
    for pwindex=1:numpulsewidths
        spiketimesPulse=[];
        for rep=1:nrepsPulse(cellnum, pwindex)
            spiketimesPulse=[spiketimesPulse MPulse(cellnum, pwindex, rep).spiketimes];
        end
        mMPulse(cellnum, pwindex).spiketimes=spiketimesPulse;
        if LaserRecorded
            mMPulseLasertrace(pwindex,:)=mean(MPulseLasertrace(pwindex, 1:nrepsPulse(cellnum, pwindex),:), 2);
        end
        if StimRecorded
            mMPulseStimtrace(pwindex,:)=mean(MPulseStimtrace(pwindex, 1:nrepsPulse(cellnum, pwindex),:), 2);
        end
    end

    %laser train
    for tnpindex=1:numtrainnumpulses
        for tpwindex=1:numtrainpulsewidths
            for tiindex=1:numtrainisis
                spiketimesTrain=[];
                for rep=1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex)
                    spiketimesTrain=[spiketimesTrain MTrain(cellnum, tnpindex,tpwindex,tiindex, rep).spiketimes];
                end
                mMTrain(cellnum, tnpindex,tpwindex,tiindex).spiketimes=spiketimesTrain;
                if LaserRecorded
                    mMTrainLasertrace(tnpindex,tpwindex,tiindex,:)=mean(MTrainLasertrace(tnpindex,tpwindex,tiindex, 1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex),:), 4);
                end
                if StimRecorded
                    mMTrainStimtrace(tnpindex,tpwindex,tiindex,:)=mean(MTrainStimtrace(tnpindex,tpwindex,tiindex, 1:nrepsTrain(cellnum, tnpindex,tpwindex,tiindex),:), 4);
                end
            end
        end
    end
end % for cellnum=1:length(SortedUnits);

  fprintf('\nnreps silent sound no laser: %d - %d', min(nrepsOFF(:)), max(nrepsOFF(:)))
    fprintf('\nnreps single laser pulse: %d - %d',min(nrepsPulse(:)), max(nrepsPulse(:)))
    fprintf('\nnreps laser train: %d - %d',min(nrepsTrain(:)), max(nrepsTrain(:)))
    fprintf('\ntotal num spikes: %d', length(spiketimes))
    fprintf('\nIn range: %d', inRange)
%save to outfile


out.MPulse=MPulse;
out.MSilentSoundOFF=MSilentSoundOFF;
out.MTrain=MTrain;
out.MPulseLasertrace=MPulseLasertrace;
out.MSilentSoundOFFLasertrace=MSilentSoundOFFLasertrace;
out.MTrainLasertrace=MTrainLasertrace;
out.MPulseStimtrace=MPulseStimtrace;
out.MSilentSoundOFFStimtrace=MSilentSoundOFFStimtrace;
out.MTrainStimtrace=MTrainStimtrace;

out.mMPulse=mMPulse;
out.mMSilentSoundOFF=mMSilentSoundOFF;
out.mMTrain=mMTrain;
out.mMPulseLasertrace=mMPulseLasertrace;
out.mMSilentSoundOFFLasertrace=mMSilentSoundOFFLasertrace;
out.mMTrainLasertrace=mMTrainLasertrace;
out.mMPulseStimtrace=mMPulseStimtrace;
out.mMSilentSoundOFFStimtrace=mMSilentSoundOFFStimtrace;
out.mMTrainStimtrace=mMTrainStimtrace;

out.silentsounddurs=silentsounddurs;
out.trainnumpulses=trainnumpulses;
out.pulsewidths=pulsewidths;
out.trainisis=trainisis;
out.trainpulsewidths=trainpulsewidths;
out.numsilentsounddurs=numsilentsounddurs;
out.numtrainnumpulses=numtrainnumpulses;
out.numpulsewidths=numpulsewidths;
out.numtrainisis=numtrainisis;
out.numtrainpulsewidths=numtrainpulsewidths;

out.nrepsOFF=nrepsOFF;
out.nrepsPulse=nrepsPulse;
out.nrepsTrain=nrepsTrain;

out.xlimits=xlimits;
out.nexts=nexts;
out.numnexts=numnexts;

out.samprate=samprate;

out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel
out.laserstarts=laserstarts;
out.numlaserstarts=numlaserstarts;

[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysPath=EphysPath;
out.EphysPath_KS=EphysPath_KS;
out.SortedUnits=SortedUnits;

out.generated_by=mfilename;
out.generated_on=datestr(now);

try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
end
outfilename='outPSTH.mat';
save (outfilename, 'out')

