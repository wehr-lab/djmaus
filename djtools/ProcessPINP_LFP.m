function ProcessPINP_LFP(varargin)
if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};
channel=varargin{2};
xlimits=varargin{3};

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

%check if this is an appropriate stimulus protocol
switch (GetPlottingFunction(datadir))
    case 'PlotPINP_PSTH'
        %fine
    otherwise
        error('This stimulus protcol does not appear to be a PINP protocol.')
        %you could add an exception here
end

try
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
filename=getContinuousFilename('.', channel);
[scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(filename); %data
%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);
if exist(messagesfilename, 'file')~=2
    error('Could not find messages.events file. Is this the right data directory?')
end

%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

%get Events and soundcard trigger timestamps
[Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
%there are some general notes on the format of Events and network messages in help GetEventsAn
filename=getContinuousFilename('.', channel);

if exist('Events.mat')
    load('Events.mat')
    fprintf('\nloading fixed Events files with correct number of SCTs\n')
end

try
    fprintf('\nNumber of logged stimuli in notebook: %d', length(stimlog));
catch
    fprintf('\nCould not find stimlog, no logged stimuli in notebook!!');
end


j=0;
alltrainnumpulses=[];
alltrainisis=[];
alltrainpulsewidths=[];
allpulsewidths=[];
allsilentsounddurs=[];
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
if isempty(getLaserfile('.'))
    LaserRecorded=0;
    warning('Laser monitor channel not recorded')
else
    LaserRecorded=1;
end
if isempty(getStimfile('.'))
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

nrepsPulse=zeros(numpulsewidths);
nrepsTrain=zeros(numtrainnumpulses, numtrainpulsewidths, numtrainisis);
nrepsOFF=zeros(numsilentsounddurs);

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
samprate=sampleRate;
for i=1:length(Events)
    if strcmp(Events(i).type, 'silentsound')
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate;
        else
            error('???')
        end
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        if isempty(find(region<1)) %(disallow negative or zero start times)
            if Events(i).VarLaser ==0 & Events(i).laser ==0 %silent sound no laser
                silentsounddur = Events(i).duration;
                dindex= find(silentsounddur==silentsounddurs);
                nrepsOFF(dindex)=nrepsOFF(dindex)+1;
                MSilentSoundOFF(dindex, nrepsOFF(dindex),:)=scaledtrace(region);
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
                nrepsPulse(pwindex)=nrepsPulse(pwindex)+1;
                MPulse(pwindex, nrepsPulse(pwindex),:)=scaledtrace(region);
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
                nrepsTrain(tnpindex, tpwindex, tiindex)=nrepsTrain(tnpindex,tpwindex, tiindex)+1;
                MTrain(tnpindex, tpwindex, tiindex, nrepsTrain(tnpindex, tpwindex, tiindex),:)=scaledtrace(region);
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


%silent sound no laser
mMSilentSoundOFF=[];
for dindex=1:numsilentsounddurs
     mMSilentSoundOFF(dindex,:)=mean(MSilentSoundOFF(dindex, 1:nrepsOFF(dindex),:),2);
    if LaserRecorded
        mMSilentSoundOFFLasertrace(dindex,:)=mean(MSilentSoundOFFLasertrace(dindex, 1:nrepsOFF(dindex),:), 2);
    end
    if StimRecorded
        mMSilentSoundOFFStimtrace(dindex,:)=mean(MSilentSoundOFFStimtrace(dindex, 1:nrepsOFF(dindex),:), 2);
    end
end

%single laser pulse
mMPulse=[];
for pwindex=1:numpulsewidths
    mMPulse(pwindex,:)=mean(MPulse(pwindex,1:nrepsPulse(pwindex),:),2);
    if LaserRecorded
        mMPulseLasertrace(pwindex,:)=mean(MPulseLasertrace(pwindex, 1:nrepsPulse(pwindex),:), 2);
    end
    if StimRecorded
        mMPulseStimtrace(pwindex,:)=mean(MPulseStimtrace(pwindex, 1:nrepsPulse(pwindex),:), 2);
    end
end

%laser train
    mMTrain=[];
    for tnpindex=1:numtrainnumpulses
        for tpwindex=1:numtrainpulsewidths
            for tiindex=1:numtrainisis
                mMTrain(tnpindex,tpwindex,tiindex,:)=mean(MTrain(tnpindex,tpwindex,tiindex, 1:nrepsTrain(tnpindex,tpwindex,tiindex),:), 4);
                
                if LaserRecorded
                    mMTrainLasertrace(tnpindex,tpwindex,tiindex,:)=mean(MTrainLasertrace(tnpindex,tpwindex,tiindex, 1:nrepsTrain(tnpindex,tpwindex,tiindex),:), 4);
                end
                if StimRecorded
                    mMTrainStimtrace(tnpindex,tpwindex,tiindex,:)=mean(MTrainStimtrace(tnpindex,tpwindex,tiindex, 1:nrepsTrain(tnpindex,tpwindex,tiindex),:), 4);
                end
            end
        end
    end
out.channel=channel;
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
out.samprate=sampleRate;
out.datadir=datadir;
out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel
out.laserstarts=laserstarts;
out.numlaserstarts=numlaserstarts;

try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
end
out.filename=filename;
outfilename=sprintf('outLFP_ch%d.mat',channel);
save (outfilename, 'out')
