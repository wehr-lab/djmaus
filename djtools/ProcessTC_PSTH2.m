function out = ProcessTC_PSTH2(varargin)

%processes clustered spiking tuning curve data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%creates only a single outfile per session
%(automatically processes all cells)
%-mike 04.28.2024
%
% usage: ProcessTC_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, [xlimits],[ylimits])
%   SortedUnits is created by ProcessSpikes
% BonsaiPath, EphysDir, EphysDir_KS are folders, they are stored in OEinfo-xxx.mat in the Bonsai directory
% (xlimits, ylimits are optional)
% xlimits default to [0 200]
%
% saves to outfile


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

try
    cd(BonsaiPath)
    cd(EphysPath)
    sessionfilename=['session-',EphysPath];
    load(sessionfilename)
    fprintf('\nloaded session object')
catch
    warning('could not load session object')
end

%read messages
cd(BonsaiPath)
behavior_filename=dir('Behavior_*.mat');
load(behavior_filename.name);


%get Events and soundcard trigger timestamps
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps

if exist('Events.mat')
    load('Events.mat')
    fprintf('\nloaded Events file \n')
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
    case 'PlotTC_PSTH'
        %that's fine
    case 'PlotMGB_LGN'
        %also fine, for now
    case 'PlotPINP_PSTH'
        %fine
    otherwise
        error('This stimulus protcol does not appear to have any tones or whitenoise.')
end

%read  Kilosort data


%uncomment this to run some sanity checks
%SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)


fprintf('\ncomputing tuning curve...');

samprate=Sky.OEsamplerate;

%get freqs/amps
j=0;
allfreqs=[];
allamps=[];
alltrainnumpulses=[];
alltrainisis=[];
alltrainpulsewidths=[];
allpulsewidths=[];
allsilentsounddurs=[];
allVarLasers=[];
allVarLaserstarts=[];
alllasers=[];

for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') ...
            |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise') ...
            | strcmp(Events(i).type, 'grating') |strcmp(Events(i).type, 'silentsound')
        j=j+1;
        alldurs(j)=Events(i).duration;
        if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, '2tone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).frequency;
        elseif strcmp(Events(i).type, 'whitenoise')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=-1;
        elseif strcmp(Events(i).type, 'fmtone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).carrier_frequency;
        elseif strcmp(Events(i).type, 'grating')
            allfreqs(j)=Events(i).angle*1000;
            allamps(j)=Events(i).spatialfrequency;
        elseif strcmp(Events(i).type, 'silentsound')
            allfreqs(j)=-1;
            allamps(j)=-1000;
            alldurs(j)=Events(i).duration;
            %PINP protocol params
            alllasers(j)=Events(i).laser;
            if isfield(Events(i), 'VarLaser')
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
            else
                %allVarLasers(j)=[];
                %allVarLaserstarts(j)=[];
                allsilentsounddurs=[allsilentsounddurs Events(i).duration];
            end
        end
    end
end
freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
numfreqs=length(freqs);
numamps=length(amps);
numdurs=length(durs);

lasers=unique(alllasers);
laserstarts=unique(allVarLaserstarts);
trainnumpulses=unique(alltrainnumpulses);
trainpulsewidths=unique(alltrainpulsewidths);
trainisis=unique(alltrainisis);
pulsewidths=unique(allpulsewidths);
silentsounddurs=unique(allsilentsounddurs);
VarLaserstarts=unique(allVarLaserstarts);
numlaserstarts=length(laserstarts);
numtrainnumpulses=length(trainnumpulses);
numtrainpulsewidths=length(trainpulsewidths);
numtrainisis=length(trainisis);
numpulsewidths=length(pulsewidths);
numsilentsounddurs=length(silentsounddurs);

%check for laser in Eventsdjtools/ProcessTC_PSTH_single.m
for i=1:length(Events)
    if isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
        LaserScheduled(i)=Events(i).laser; %whether the stim protocol scheduled a laser for this stim
        LaserOnOffButton(i)=Events(i).LaserOnOff; %whether the laser button was turned on
        LaserTrials(i)=LaserScheduled(i) & LaserOnOffButton(i);
        if isempty(stimlog(i).LaserStart)
            LaserStart(i)=nan;
            LaserWidth(i)=nan;
            LaserNumPulses(i)=nan;
            LaserISI(i)=nan;
        else
            %these laser params are from the djmaus gui, they would be
            %overwritten if there are also varlaser params in the protocol
            LaserStart(i)=stimlog(i).LaserStart; 
            LaserWidth(i)=stimlog(i).LaserWidth;
            LaserNumPulses(i)=stimlog(i).LaserNumPulses;
            LaserISI(i)=stimlog(i).LaserISI;
        end
        if isfield(Events(i), 'VL') %this is a varlaser protocol
            VL(i)=Events(i).VL;
            VLstart(i)=Events(i).VLstart;
            VLpulsewidth(i)= Events(i).VLpulsewidth;
            VLnumpulses(i)= Events(i).VLnumpulses;
            VLisi(i)= Events(i).VLisi;
        end

    elseif isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %Not sure about this one. Assume no laser for now, but investigate.
        warning('ProcessTC_PSTH2: Cannot tell if laser button was turned on in djmaus GUI');
        LaserTrials(i)=0;
        Events(i).laser=0;
    elseif ~isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %if neither of the right fields are there, assume no laser
        LaserTrials(i)=0;
        Events(i).laser=0;
    else
        error('wtf?')
    end
end
fprintf('\n%d laser pulses in this Events file', sum(LaserTrials))
try
    if sum(LaserOnOffButton)==0
        fprintf('\nLaser On/Off button remained off for entire file.')
    end
end
if sum(LaserTrials)>0
    IL=1;
else
    IL=0;
end
%if lasers were used, we'll un-interleave them and save ON and OFF data


%try to load laser and stimulus monitor
if exist('lasertrace')==1
    %we loaded lasertrace from open ephys Session object
        LaserRecorded=1;
elseif isempty(getLaserfile('.'))
    LaserRecorded=0;
    warning('Laser monitor channel not recorded')
    warning('ProcessTC_PSTH_single2: still need to figure out how to find/plot Laser monitor channel ')
else
    LaserRecorded=1;
end
if exist('stimtrace')==1
    %we loaded stimtrace from open ephys Session object
        StimRecorded=1;
elseif isempty(getStimfile('.'))
    StimRecorded=0;
    warning('Stimulus monitor channel not recorded')
    warning('ProcessTC_PSTH_single2: still need to figure out how to find/plot Stimulus monitor channel ')
else
    StimRecorded=1;
end

if LaserRecorded
    if exist('lasertrace')==1
        Lasertimestamps=timestamps-timestamps(1);
        Lasertrace=lasertrace./max(abs(lasertrace));
    else
        try
            [Lasertrace, Lasertimestamps, Laserinfo] =load_open_ephys_data(getLaserfile('.'));
            %Lasertimestamps=Lasertimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
            Lasertimestamps=Lasertimestamps-Lasertimestamps(1);
            Lasertrace=Lasertrace./max(abs(Lasertrace));
            fprintf('\nsuccessfully loaded laser trace')
        catch
            [Lasertrace, Lasertimestamps, Laserinfo] =load_open_ephys_data('105_ADC2.continuous');
            %fprintf('\nfound laser file %s but could not load laser trace',getStimfile('.'))
        end
    end
else
    fprintf('\nLaser trace not recorded')
end
if StimRecorded
     if exist('stimtrace')==1
        Stimtimestamps=timestamps-timestamps(1);
        Stimtrace=stimtrace./max(abs(stimtrace));
    else
    try
        [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(getStimfile('.'));
        %Stimtimestamps=Stimtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
        Stimtimestamps=Stimtimestamps-Stimtimestamps(1);
        Stimtrace=Stimtrace./max(abs(Stimtrace));
        fprintf('\nsuccessfully loaded stim trace')
    catch
        [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data('105_ADC1.continuous');
        %fprintf('\nfound stim file %s but could not load stim trace', getStimfile('.'))
    end
     end
else
    fprintf('\nSound stimulus trace not recorded')
end

numcells=length(SortedUnits);
M1=[];M1ON=[];M1OFF=[];
mM1ON=[];
mM1OFF=[];
nreps=zeros(numcells, numfreqs, numamps, numdurs);
nrepsON=zeros(numcells, numfreqs, numamps, numdurs);
nrepsOFF=zeros(numcells, numfreqs, numamps, numdurs);
nreps_ssON=zeros(numcells);
nreps_ssOFF=zeros(numcells);
%extract the traces into a big matrix M
j=0;
inRange=0;

M1ONLaser=[];
M1ONStim=[];
SilentSoundON=[];
SilentSoundOFF=[];
SilentSoundONspikecount=[];
SilentSoundOFFspikecount=[];

SilentSoundONStim=[];
SilentSoundOFFStim=[];
SilentSoundONLaser=[];
SilentSoundOFFLaser=[];

for cellnum=1:length(SortedUnits);
    spiketimes=[];
    fprintf('\nreading SortedUnits cell %d', cellnum)
    spiketimes= SortedUnits(cellnum).spiketimes;

    totalnumspikes=length(spiketimes);
    fprintf('\nsuccessfully loaded spike data with %d spikes\n',totalnumspikes)

    for i=1:length(Events)
        if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') | ...
                strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating') | ...
                strcmp(Events(i).type, 'silentsound')
            if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
                pos=Events(i).soundcard_trigger_timestamp_sec;
            else
                error('???')
                pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
            end
            laser=LaserTrials(i);
            start=(pos+xlimits(1)*1e-3); %in seconds
            stop=(pos+xlimits(2)*1e-3);
            region=round(start*samprate)+1:round(stop*samprate);
            %if this wierd case where start was clicked twice
            %region=region+abs(StartAcquisitionSamples-StartAcquisitionSamples_wrong);
            if start>0 %(disallow negative or zero start times)
                %             if stop>lostat
                %                 fprintf('\ndiscarding trace (after lostat)')
                %             elseif start<gotat
                %                fprintf('\ndiscarding trace (before gotat)')
                %                %commented out by ira 09-05-2013
                %             else
                switch Events(i).type
                    case {'tone', '2tone'}
                        freq=Events(i).frequency;
                        amp=Events(i).amplitude;
                    case 'fmtone'
                        freq=Events(i).carrier_frequency;
                        amp=Events(i).amplitude;
                    case 'whitenoise'
                        freq=-1;
                        amp=Events(i).amplitude;
                    case 'grating'
                        amp=Events(i).spatialfrequency;
                        freq=Events(i).angle*1000;
                end

                st=spiketimes; %are in seconds
                spiketimes1=st(st>start & st<stop); % spiketimes in region
                spikecount=length(spiketimes1); % No. of spikes fired in response to this rep of this stim.
                inRange=inRange+ spikecount; %accumulate total spikecount in region
                spiketimes1=(spiketimes1-pos)*1000;%convert to ms after tone onset
                spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window

                if strcmp(Events(i).type, 'silentsound')
                    if laser
                        nreps_ssON(cellnum)=nreps_ssON(cellnum)+1;
                        SilentSoundON(cellnum, nreps_ssON(cellnum)).spiketimes=spiketimes1;
                        SilentSoundONspikecount(cellnum, nreps_ssON(cellnum))=spikecount;

                        if LaserRecorded
                            SilentSoundONLaser(nreps_ssON(cellnum),:)=Lasertrace(region);
                        else
                            SilentSoundONLaser=[];
                        end
                        if StimRecorded
                            SilentSoundONStim(nreps_ssON(cellnum),:)=Stimtrace(region);
                        else
                            SilentSoundONStim=[];
                        end
                    else
                        nreps_ssOFF(cellnum)=nreps_ssOFF(cellnum)+1;
                        SilentSoundOFF(cellnum, nreps_ssOFF(cellnum)).spiketimes=spiketimes1;
                        SilentSoundOFFspikecount(cellnum, nreps_ssOFF(cellnum))=spikecount;
                        if LaserRecorded
                            SilentSoundOFFLaser(nreps_ssOFF(cellnum),:)=Lasertrace(region);
                        else
                            SilentSoundOFFLaser=[];
                        end
                        if StimRecorded
                            SilentSoundOFFStim(nreps_ssOFF(cellnum),:)=Stimtrace(region);
                        else
                            SilentSoundOFFStim=[];
                        end
                    end
                else
                    dur=Events(i).duration;
                    findex= find(freqs==freq);
                    aindex= find(amps==amp);
                    dindex= find(durs==dur);
                    nreps(cellnum, findex, aindex, dindex)=nreps(cellnum, findex, aindex, dindex)+1;


                    if laser
                        nrepsON(cellnum,findex, aindex, dindex)=nrepsON(cellnum,findex, aindex, dindex)+1;
                        M1ON(cellnum, findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex)).spiketimes=spiketimes1; % Spike times
                        M1ONspikecounts(cellnum, findex,aindex,dindex,nrepsON(cellnum, findex, aindex, dindex))=spikecount; % No. of spikes
                        M1spontON(cellnum, findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                        M_LaserStart(findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex))=LaserStart(i);
                        M_LaserWidth(findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex))= LaserWidth(i);
                        M_LaserNumPulses(findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex))= LaserNumPulses(i);
                        M_LaserISI(findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex))= LaserISI(i);
                        if LaserRecorded
                            M1ONLaser(findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex),:)=Lasertrace(region);
                        end
                        if StimRecorded
                            M1ONStim(findex,aindex,dindex, nrepsON(cellnum, findex, aindex, dindex),:)=Stimtrace(region);
                        end
                    else
                        nrepsOFF(cellnum,findex, aindex, dindex)=nrepsOFF(cellnum,findex, aindex, dindex)+1;
                        M1OFF(cellnum, findex,aindex,dindex, nrepsOFF(cellnum, findex, aindex, dindex)).spiketimes=spiketimes1;
                        M1OFFspikecounts(cellnum, findex,aindex,dindex,nrepsOFF(cellnum, findex, aindex, dindex))=spikecount;
                        M1spontOFF(cellnum, findex,aindex,dindex, nrepsOFF(cellnum, findex, aindex, dindex))=spont_spikecount;
                        if LaserRecorded
                            M1OFFLaser(findex,aindex,dindex, nrepsOFF(cellnum, findex, aindex, dindex),:)=Lasertrace(region);
                        end
                        if StimRecorded
                            M1OFFStim(findex,aindex,dindex, nrepsOFF(cellnum, findex, aindex, dindex),:)=Stimtrace(region);
                        end
                    end
                end
            end
        end
    end


    fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
    fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))
    fprintf('\ncell %d, total num spikes: %d', cellnum, length(spiketimes))
    fprintf('\nIn range: %d', inRange)

    % Accumulate spiketimes across trials, for psth...
    for dindex=1:length(durs); % Hardcoded.
        for aindex=[numamps:-1:1]
            for findex=1:numfreqs

                % on
                spiketimesON=[];
                for rep=1:nrepsON(cellnum, findex, aindex, dindex)
                    spiketimesON=[spiketimesON M1ON(cellnum, findex, aindex, dindex, rep).spiketimes];
                end

                % All spiketimes for a given f/a/d combo, for psth:
                mM1ON(cellnum, findex, aindex, dindex).spiketimes=spiketimesON;

                % off
                spiketimesOFF=[];
                for rep=1:nrepsOFF(cellnum, findex, aindex, dindex)
                    spiketimesOFF=[spiketimesOFF M1OFF(cellnum, findex, aindex, dindex, rep).spiketimes];
                end
                mM1OFF(cellnum, findex, aindex, dindex).spiketimes=spiketimesOFF;
            end
        end
    end

    spiketimesON=[];
    for rep=1:nreps_ssON(cellnum)
        spiketimesON=[spiketimesON SilentSoundON(cellnum, rep).spiketimes];
    end
    mSilentSoundON(cellnum).spiketimes=spiketimesON;
    spiketimesOFF=[];
    for rep=1:nreps_ssOFF(cellnum)
        spiketimesOFF=[spiketimesOFF SilentSoundOFF(cellnum, rep).spiketimes];
    end
    mSilentSoundOFF(cellnum).spiketimes=spiketimesOFF;


    %average laser and stimulus monitor M matrices across trials
    if LaserRecorded
        for aindex=1:numamps
            for findex=1:numfreqs
                for dindex=1:numdurs
                    if nrepsON(cellnum, findex, aindex, dindex)>0
                        mM1ONLaser(findex, aindex, dindex,:)=mean(M1ONLaser(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                    else %no reps for this stim, since rep=0
                        mM1ONLaser(findex, aindex, dindex,:)=zeros(size(region));
                    end
                    if nrepsOFF(cellnum, findex, aindex, dindex)>0
                        mM1OFFLaser(findex, aindex, dindex,:)=mean(M1OFFLaser(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    else %no reps for this stim, since rep=0
                        mM1OFFLaser(findex, aindex, dindex,:)=zeros(size(region));
                    end
                end
            end
        end
    end
    if StimRecorded
        for aindex=1:numamps
            for findex=1:numfreqs
                for dindex=1:numdurs
                    if nrepsON(cellnum, findex, aindex, dindex)>0
                        mM1ONStim(findex, aindex, dindex,:)=mean(M1ONStim(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                    else %no reps for this stim, since rep=0
                        mM1ONStim(findex, aindex, dindex,:)=zeros(size(region));
                    end
                    if nrepsOFF(cellnum, findex, aindex, dindex)>0
                        mM1OFFStim(findex, aindex, dindex,:)=mean(M1OFFStim(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    else %no reps for this stim, since rep=0
                        mM1OFFStim(findex, aindex, dindex,:)=zeros(size(region));
                    end
                end
            end
        end
    end

    if ~IL %no laser pulses in this file
        mM1ONspikecount=[];
        sM1ONspikecount=[];
        semM1ONspikecount=[];
        M1ONspikecounts=[];
    else
        if exist('M1ONspikecounts') & sum(nrepsON(:)) 
            %for a PINP protocol with only silent sounds, sum(nrepsON(:)
            %will be zero
            mM1ONspikecount(cellnum,:,:,:)=mean(M1ONspikecounts(cellnum,:,:,:,:),5); % Mean spike count
            sM1ONspikecount(cellnum, :,:,:)=std(M1ONspikecounts(cellnum,:,:,:,:),[],5); % Std of the above
            semM1ONspikecount(cellnum,:,:,:)=sM1ONspikecount(cellnum,:,:,:,:)./sqrt(max(nrepsON(:))); % Sem of the above

            % Spont
            mM1spontON(cellnum,:,:,:)=mean(M1spontON(cellnum,:,:,:,:),5);
            sM1spontON(cellnum,:,:,:)=std(M1spontON(cellnum,:,:,:,:),[],5);
            semM1spontON(cellnum,:,:,:)=sM1spontON(cellnum,:,:,:,:)./sqrt(max(nrepsON(:)));
        else
            mM1ONspikecount=[];
            sM1ONspikecount=[];
            semM1ONspikecount=[];
            M1ONspikecounts=[];
            mM1spontON=[];
            sM1spontON=[];
            semM1spontON=[];
        end
    end
    if isempty(mM1OFF) %no laser pulses in this file
        mM1OFFspikecount=[];
        sM1OFFspikecount=[];
        semM1OFFspikecount=[];
        M1OFFspikecounts=[];
        mM1spontOFF=[];
        sM1spontOFF=[];
        semM1spontOFF=[];
    else
        if exist('M1OFFspikecounts') & sum(nrepsOFF(:)) 
            mM1OFFspikecount(cellnum,:,:,:)=mean(M1OFFspikecounts(cellnum,:,:,:,:),5); % Mean spike count
            sM1OFFspikecount(cellnum,:,:,:)=std(M1OFFspikecounts(cellnum,:,:,:,:),[],5); % Std of the above
            semM1OFFspikecount(cellnum,:,:,:)=sM1OFFspikecount(cellnum,:,:,:,:)./sqrt(max(nrepsOFF(:))); % Sem of the above
            % Spont
            mM1spontOFF(cellnum,:,:,:)=mean(M1spontOFF(cellnum,:,:,:,:),5);
            sM1spontOFF(cellnum,:,:,:)=std(M1spontOFF(cellnum,:,:,:,:),[],5);
            semM1spontOFF(cellnum,:,:,:)=sM1spontOFF(cellnum,:,:,:,:)./sqrt(max(nrepsOFF(:)));
        else
            mM1OFFspikecount=[];
            sM1OFFspikecount=[];
            semM1OFFspikecount=[];
            M1OFFspikecounts=[];
            mM1spontOFF=[];
            sM1spontOFF=[];
            semM1spontOFF=[];
        end
    end

end % for cellnum=1:length(SortedUnits);


%save to outfiles
%one outfile for all cells
% previously existing outfiles for this session will be overwritten

%after squeezing cluster, saves with the following dimensions:
% M1ON(cellnum, findex,aindex,dindex, nrepsON).spiketimes
% mM1ON(cellnum, findex,aindex,dindex).spiketimes
% mM1ONspikecount(cellnum, findex,aindex,dindex)

out.IL=IL;
out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
out.LaserWidth=unique(LaserWidth);
out.LaserNumPulses=unique(LaserNumPulses);
out.LaserISI=unique(LaserISI);

if IL
    out.M1ON=M1ON;
    out.mM1ON=mM1ON; %Accumulated spike times for *all* presentations of each laser/f/a/d combo.
    out.mM1ONspikecount=mM1ONspikecount;
    out.sM1ONspikecount=sM1ONspikecount;
    out.semM1ONspikecount=semM1ONspikecount;
    out.mM1spontON=mM1spontON;
    out.sM1spontON=sM1spontON;
    out.semM1spontON=semM1spontON;

    if exist('M_LaserStart')
        out.M_LaserStart=M_LaserStart;
        out.M_LaserWidth=M_LaserWidth;
        out.M_LaserNumPulses=M_LaserNumPulses;
        out.M_LaserISI=M_LaserISI;
    end

else
    out.M1ON=[];
    out.mM1ONspikecount=[];
    out.sM1ONspikecount=[];
    out.semM1ONspikecount=[];
    out.mM1ON=[];
    out.mM1spontON=[];
    out.sM1spontON=[];
    out.semM1spontON=[];
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
    out.M_LaserStart=[];
    out.M_LaserWidth=[];
    out.M_LaserNumPulses=[];
    out.M_LaserISI=[];

end
out.M1OFF=M1OFF;
out.mM1OFF=mM1OFF;
out.mM1OFFspikecount=mM1OFFspikecount;
out.sM1OFFspikecount=sM1OFFspikecount;
out.semM1OFFspikecount=semM1OFFspikecount;
out.mM1spontOFF=mM1spontOFF;
out.sM1spontOFF=sM1spontOFF;
out.semM1spontOFF=semM1spontOFF;
out.amps=amps;
out.freqs=freqs;
out.durs=durs;
out.numamps=numamps;
out.numfreqs=numfreqs;
out.numdurs=numdurs;
out.nreps=nreps;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.xlimits=xlimits;
out.samprate=samprate;

[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysPath=EphysPath;
out.EphysPath_KS=EphysPath_KS;
out.SortedUnits=SortedUnits;
out.spiketimes=spiketimes;
out.nreps_ssON=nreps_ssON;
out.nreps_ssOFF=nreps_ssOFF;
out.SilentSoundON=SilentSoundON;
out.SilentSoundOFF=SilentSoundOFF;
out.SilentSoundONspikecount=SilentSoundONspikecount;
out.SilentSoundOFFspikecount=SilentSoundOFFspikecount;
out.mSilentSoundON=mSilentSoundON;
out.mSilentSoundOFF=mSilentSoundOFF;

out.SilentSoundONStim=SilentSoundONStim;
out.SilentSoundOFFStim=SilentSoundOFFStim;
out.SilentSoundONLaser=SilentSoundONLaser;
out.SilentSoundOFFLaser=SilentSoundOFFLaser;

out.generated_by=mfilename;
out.generated_on=datestr(now);


out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel
if LaserRecorded
    out.M1ONLaser=M1ONLaser;
    if exist('mM1ONLaser')
        out.mM1ONLaser=mM1ONLaser;
    end
    if exist('M1OFFLaser')
        out.M1OFFLaser=M1OFFLaser;
    end
    if exist('mM1OFFLaser')
        out.mM1OFFLaser=mM1OFFLaser;
    end
    if exist('mM1ONLaser')
        out.mM1ONLaser=mM1ONLaser;
    end
else
    out.M1ONLaser=[];
    out.mM1ONLaser=[];
    out.M1OFFLaser=[];
    out.mM1OFFLaser=[];
end

if exist('VL')
    out.VL=VL;
    out.VLstart=VLstart;
    out.VLpulsewidth=VLpulsewidth;
    out.VLnumpulses=VLnumpulses;
    out.VLisi=VLisi;
end

if StimRecorded
    out.M1ONStim=M1ONStim;
    if exist('mM1ONStim')
        out.mM1ONStim=mM1ONStim;
    end
    if exist('M1OFFStim')
        out.M1OFFStim=M1OFFStim;
    end
    if exist('mM1OFFStim')
        out.mM1OFFStim=mM1OFFStim;
    end
else
    out.M1ONStim=[];
    out.mM1ONStim=[];
    out.M1OFFStim=[];
    out.mM1OFFStim=[];
end

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

