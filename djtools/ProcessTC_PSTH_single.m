function ProcessTC_PSTH_single(varargin)

%processes a single .t file of clustered spiking tuning curve data from djmaus
%
% usage: ProcessTC_PSTH(datadir, t_filename, [xlimits],[ylimits])
% (xlimits, ylimits are optional)
% xlimits default to [0 200]
% channel (tetrode) number should be an integer, if omitted you will be prompted for it
% automatically processes dafunction ProcessTC_PSTH_single(varargin)

%processes a single .t file of clustered spiking tuning curve data from djmaus
%
% usage: ProcessTC_PSTH(datadir, t_filename, [xlimits],[ylimits])
% (xlimits, ylimits are optional)
% xlimits default to [0 200]
% channel (tetrode) number should be an integer, if omitted you will be prompted for it
% automatically processes data for all cells clustered for that tetrode
% saves to outfile


if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

xlimits=[];
try
    xlimits=varargin{3};
end
if isempty(xlimits)
    xlimits=[-100 200];
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end

filename=varargin{2};
[p,f,ext]=fileparts(filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

fprintf('\nchannel %d, cluster %d', channel, clust)
fprintf('\nprocessing with xlimits [%d-%d]', xlimits(1), xlimits(2))

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

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
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps

try
    fprintf('\nNumber of logged stimuli in notebook: %d', length(stimlog));
catch
    fprintf('\nCould not find stimlog, no logged stimuli in notebook!!');
end

%check if this is an appropriate stimulus protocol
if ~strcmp(GetPlottingFunction(datadir), 'PlotTC_PSTH')
    error('This stimulus protcol does not appear to have any tones or whitenoise.')
end

%read MClust .t file
fprintf('\nreading MClust output file %s', filename)
spiketimes=read_MClust_output(filename)'/10000; %spiketimes now in seconds
%correct for OE start time, so that time starts at 0
spiketimes=spiketimes-StartAcquisitionSec;
totalnumspikes=length(spiketimes);
fprintf('\nsuccessfully loaded MClust spike data')
Nclusters=1;

%uncomment this to run some sanity checks
%SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)


fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') ...
            |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise')| strcmp(Events(i).type, 'grating')
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
        end
    end
end
freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
numfreqs=length(freqs);
numamps=length(amps);
numdurs=length(durs);

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
            LaserStart(i)=stimlog(i).LaserStart;
            LaserWidth(i)=stimlog(i).LaserWidth;
            LaserNumPulses(i)=stimlog(i).LaserNumPulses;
            LaserISI(i)=stimlog(i).LaserISI;
        end
        
    elseif isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %Not sure about this one. Assume no laser for now, but investigate.
        warning('ProcessTC_LFP: Cannot tell if laser button was turned on in djmaus GUI');
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

M1=[];M1ON=[];M1OFF=[];
nreps=zeros(numfreqs, numamps, numdurs);
nrepsON=zeros(numfreqs, numamps, numdurs);
nrepsOFF=zeros(numfreqs, numamps, numdurs);
nreps_ssON=0;
nreps_ssOFF=0;
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
            spiketimes1=(spiketimes1-pos)*1000;%covert to ms after tone onset
            spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window
            
            if strcmp(Events(i).type, 'silentsound')
                if laser
                    nreps_ssON=nreps_ssON+1;
                    SilentSoundON(nreps_ssON).spiketimes=spiketimes1;
                    SilentSoundONspikecount(nreps_ssON)=spikecount;
                    if LaserRecorded
                        SilentSoundONLaser(nreps_ssON,:)=Lasertrace(region);
                    else
                        SilentSoundONLaser=[];
                    end
                    if StimRecorded
                        SilentSoundONStim(nreps_ssON,:)=Stimtrace(region);
                    else
                        SilentSoundONStim=[];
                    end
                else
                    nreps_ssOFF=nreps_ssOFF+1;
                    SilentSoundOFF(nreps_ssOFF).spiketimes=spiketimes1;
                    SilentSoundOFFspikecount(nreps_ssOFF)=spikecount;
                    if LaserRecorded
                        SilentSoundOFFLaser(nreps_ssOFF,:)=Lasertrace(region);
                    else
                        SilentSoundOFFLaser=[];
                    end
                    if StimRecorded
                        SilentSoundOFFStim(nreps_ssOFF,:)=Stimtrace(region);
                    else
                        SilentSoundOFFStim=[];
                    end
                end
            else
                dur=Events(i).duration;
                findex= find(freqs==freq);
                aindex= find(amps==amp);
                dindex= find(durs==dur);
                nreps(findex, aindex, dindex)=nreps(findex, aindex, dindex)+1;
                
                
                if laser
                    nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
                    M1ON(findex,aindex,dindex, nrepsON(findex, aindex, dindex)).spiketimes=spiketimes1; % Spike times
                    M1ONspikecounts(findex,aindex,dindex,nrepsON(findex, aindex, dindex))=spikecount; % No. of spikes
                    M1spontON(findex,aindex,dindex, nrepsON(findex, aindex, dindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                    M_LaserStart(findex,aindex,dindex, nrepsON(findex, aindex, dindex))=LaserStart(i);
                    M_LaserWidth(findex,aindex,dindex, nrepsON(findex, aindex, dindex))= LaserWidth(i);
                    M_LaserNumPulses(findex,aindex,dindex, nrepsON(findex, aindex, dindex))= LaserNumPulses(i);
                    M_LaserISI(findex,aindex,dindex, nrepsON(findex, aindex, dindex))= LaserISI(i);
                    if LaserRecorded
                        M1ONLaser(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        M1ONStim(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=Stimtrace(region);
                    end
                else
                    nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
                    M1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex)).spiketimes=spiketimes1;
                    M1OFFspikecounts(findex,aindex,dindex,nrepsOFF(findex, aindex, dindex))=spikecount;
                    M1spontOFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex))=spont_spikecount;
                    if LaserRecorded
                        M1OFFLaser(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        M1OFFStim(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Stimtrace(region);
                    end
                end
            end
        end
    end
end


fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))
fprintf('\ntotal num spikes: %d', length(spiketimes))
fprintf('\nIn range: %d', inRange)


mM1ON=[];
mM1OFF=[];

% Accumulate spiketimes across trials, for psth...
for dindex=1:length(durs); % Hardcoded.
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            
            % on
            spiketimesON=[];
            spikecountsON=[];
            for rep=1:nrepsON(findex, aindex, dindex)
                spiketimesON=[spiketimesON M1ON(findex, aindex, dindex, rep).spiketimes];
            end
            
            % All spiketimes for a given f/a/d combo, for psth:
            mM1ON(findex, aindex, dindex).spiketimes=spiketimesON;
            
            % off
            spiketimesOFF=[];
            spikecountsOFF=[];
            for rep=1:nrepsOFF(findex, aindex, dindex)
                spiketimesOFF=[spiketimesOFF M1OFF(findex, aindex, dindex, rep).spiketimes];
            end
            mM1OFF(findex, aindex, dindex).spiketimes=spiketimesOFF;
        end
    end
end

spiketimesON=[];
for rep=1:nreps_ssON
    spiketimesON=[spiketimesON SilentSoundON(rep).spiketimes];
end
mSilentSoundON.spiketimes=spiketimesON;
spiketimesOFF=[];
for rep=1:nreps_ssOFF
    spiketimesOFF=[spiketimesOFF SilentSoundOFF(rep).spiketimes];
end
mSilentSoundOFF.spiketimes=spiketimesOFF;


%average laser and stimulus monitor M matrices across trials
if LaserRecorded
    for aindex=1:numamps
        for findex=1:numfreqs
            for dindex=1:numdurs
                if nrepsON(findex, aindex, dindex)>0
                    mM1ONLaser(findex, aindex, dindex,:)=mean(M1ONLaser(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1ONLaser(findex, aindex, dindex,:)=zeros(size(region));
                end
                if nrepsOFF(findex, aindex, dindex)>0
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
                if nrepsON(findex, aindex, dindex)>0
                    mM1ONStim(findex, aindex, dindex,:)=mean(M1ONStim(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1ONStim(findex, aindex, dindex,:)=zeros(size(region));
                end
                if nrepsOFF(findex, aindex, dindex)>0
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
    mM1ONspikecount=mean(M1ONspikecounts,4); % Mean spike count
    sM1ONspikecount=std(M1ONspikecounts,[],4); % Std of the above
    semM1ONspikecount=sM1ONspikecount./sqrt(max(nrepsON(:))); % Sem of the above
    
    % Spont
    mM1spontON=mean(M1spontON,4);
    sM1spontON=std(M1spontON,[],4);
    semM1spontON(:,:)=sM1spontON./sqrt(max(nrepsON(:)));
end
if isempty(mM1OFF) %no laser pulses in this file
    mM1OFFspikecount=[];
    sM1OFFspikecount=[];
    semM1OFFspikecount=[];
    M1OFFspikecounts=[];
else
    mM1OFFspikecount=mean(M1OFFspikecounts,4); % Mean spike count
    sM1OFFspikecount=std(M1OFFspikecounts,[],4); % Std of the above
    semM1OFFspikecount=sM1OFFspikecount./sqrt(max(nrepsOFF(:))); % Sem of the above
    % Spont
    mM1spontOFF=mean(M1spontOFF,4);
    sM1spontOFF=std(M1spontOFF,[],4);
    semM1spontOFF(:,:)=sM1spontOFF./sqrt(max(nrepsOFF(:)));
end

%save to outfiles
%one outfile for each cell
% previously existing outfiles for this cell will be overwritten

%after squeezing cluster, saves with the following dimensions:
% M1ON(findex,aindex,dindex, nrepsON).spiketimes
% mM1ON(findex,aindex,dindex).spiketimes
% mM1ONspikecount(findex,aindex,dindex)

out.IL=IL;
out.Nclusters=Nclusters;
out.tetrode=channel;
out.channel=channel;
out.cluster=clust; %there are some redundant names here
out.cell=clust;
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
    
    out.M_LaserStart=M_LaserStart;
    out.M_LaserWidth=M_LaserWidth;
    out.M_LaserNumPulses=M_LaserNumPulses;
    out.M_LaserISI=M_LaserISI;
    
    
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
out.datadir=datadir;
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


out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel
if LaserRecorded
    out.M1ONLaser=M1ONLaser;
    out.mM1ONLaser=mM1ONLaser;
    out.M1OFFLaser=M1OFFLaser;
    out.mM1OFFLaser=mM1OFFLaser;
else
    out.M1ONLaser=[];
    out.mM1ONLaser=[];
    out.M1OFFLaser=[];
    out.mM1OFFLaser=[];
end
if StimRecorded
    out.M1ONStim=M1ONStim;
    out.mM1ONStim=mM1ONStim;
    out.M1OFFStim=M1OFFStim;
    out.mM1OFFStim=mM1OFFStim;
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
out.t_filename=filename;
outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
save (outfilename, 'out')

