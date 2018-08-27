function Process2tone_PSTH_single2(varargin)
% tailored to plot Yashar's stimuli, comparing tones vs laser pulses
%processes a single .t file of clustered spiking 2 tone curve data from djmaus
%
% usage: Process2tone_PSTH_single(datadir, t_filename, [xlimits],[ylimits])
% (xlimits, ylimits are optional)
% xlimits default to [-.5*dur SOA+1.5*dur]
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
dbstop if error

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
    %     xlimits=[-100 200];
    s=GetStimParams(datadir);
    durs=s.durs;
    dur=max(durs);
    SOAs=s.SOAs;
    SOA=max(SOAs);
    xlimits=[-.5*dur SOA+1.5*dur]; %default x limits for axis
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
switch (GetPlottingFunction(datadir))
    case 'Plot2Tone_PSTH2'
      case 'Plot2Tone_PSTH_single2'
    otherwise
        warning('GetPlottingFunction does not check out.')
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
VarLaser=0;
%get freqs/amps
j=0;
allfreqs=[];
allamps=[];
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone')
        j=j+1;
        if Events(i).frequency==0
            keyboard
        end
        alldurs(j)=Events(i).duration;
        allamps(j)=Events(i).amplitude;
        allfreqs(j)=Events(i).frequency;
        allprobefreqs(j)=Events(i).probefreq;
        allprobeamps(j)=Events(i).probeamp;
        allprobedurs(j)=Events(i).duration;
        allSOAs(j)=Events(i).SOA;
    elseif strcmp(Events(i).type, 'whitenoise')
        j=j+1;
        if Events(i).frequency==0
            keyboard
        end
        alldurs(j)=Events(i).duration;
        allamps(j)=Events(i).amplitude;
        allfreqs(j)=Events(i).frequency;
    elseif strcmp(Events(i).type, 'silentsound')
        j=j+1;
        alldurs(j)=Events(i).duration;
        if Events(i).VarLaser==1
            allVarLaserstart(j)=Events(i).VarLaserstart;
            allVarLaserpulsewidth(j)=Events(i).VarLaserpulsewidth;
            allVarLasernumpulses(j)=Events(i).VarLasernumpulses;
            allVarLaserisi(j)=Events(i).VarLaserisi;
        end
    end
end
freqs=nonzeros(unique(allfreqs));
amps=nonzeros(unique(allamps));
durs=unique(alldurs);
SOAs=unique(allSOAs);
probefreqs=unique(allprobefreqs);
probeamps=unique(allprobeamps);
probedurs=unique(allprobedurs);
numfreqs=length(freqs);
numprobefreqs=length(probefreqs);
numamps=length(amps);
numdurs=length(durs);
numSOAs=length(SOAs);
numprobeamps=length(probeamps);
numprobedurs=length(probedurs);

%check for laser in Eventsdjtools/ProcessTC_PSTH_single.m
for i=1:length(Events)
    if Events(i).VarLaser==1
        LaserTrials(i)=1;
        allLaserStart(i)=stimlog(i).LaserStart;
        allLaserWidth(i)=stimlog(i).LaserWidth;
        allLaserNumPulses(i)=stimlog(i).LaserNumPulses;
        allLaserISI(i)=stimlog(i).LaserISI;
        VarLaser=1;
    else
        LaserTrials(i)=0;
    end
end

if VarLaser==1
    LaserStarts=unique(allVarLaserstart);
    LaserWidths=unique(allVarLaserpulsewidth);
    LaserWidths=nonzeros(LaserWidths);
    LaserNumPulses=unique(allVarLasernumpulses);
    LaserNumPulses=nonzeros(LaserNumPulses);
    LaserISIs=unique(allVarLaserisi);
    numLaserStarts=length(LaserStarts);
    numLaserWidths=length(LaserWidths);
    numLaserNumPulses=length(LaserNumPulses);
    numVarLaserISIs=length(LaserISIs);
    
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
nreps=zeros(numfreqs, numprobefreqs, numamps, numSOAs);
nrepsON=[];
nrepsOFF=zeros(numfreqs, numprobefreqs, numamps, numSOAs);
nreps_ssON=0;
nreps_ssOFF=0;
if VarLaser==1
    
    nreps_ssON= zeros( numLaserStarts, numLaserWidths, numLaserNumPulses, numVarLaserISIs); %nreps for ss with varying laser
end
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
    if strcmp(Events(i).type, '2tone') | ...
            strcmp(Events(i).type, 'silentsound') | strcmp(Events(i).type, 'whitenoise')%%%I should pull out silent sound processing into separate stanza
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
        if start>0 %(disallow negative or zero start times)
            
            
            st=spiketimes; %are in seconds
            spiketimes1=st(st>start & st<stop); % spiketimes in region
            spikecount=length(spiketimes1); % No. of spikes fired in response to this rep of this stim.
            inRange=inRange+ spikecount; %accumulate total spikecount in region
            spiketimes1=(spiketimes1-pos)*1000;%covert to ms after tone onset
            spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window
            if strcmp(Events(i).type, '2tone')
                freq=Events(i).frequency;
                probefreq=Events(i).probefreq;
                dur=Events(i).duration;
                amp=Events(i).amplitude;
                findex= find(freqs==freq);
                aindex= find(amps==amp);
                dindex= find(durs==dur);
                pfindex= find(probefreqs==probefreq);
                SOAindex= find(SOAs==SOA);
                nreps(findex, pfindex, aindex, SOAindex)+1;
                
                nrepsOFF(findex, pfindex, aindex, SOAindex)=nrepsOFF(findex, pfindex, aindex, SOAindex)+1;
                M1OFF(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex)).spiketimes=spiketimes1;
                M1OFFspikecounts(findex, pfindex, aindex, SOAindex,nrepsOFF(findex, pfindex, aindex, SOAindex))=spikecount;
                M1spontOFF(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex))=spont_spikecount;
                if LaserRecorded
                    M1OFFLaser(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex),:)=Lasertrace(region);
                end
                if StimRecorded
                    M1OFFStim(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex),:)=Stimtrace(region);
                end
            elseif strcmp(Events(i).type, 'silentsound')
                if laser && VarLaser==1 %this should only be the case with silent sound and var laser
                    LaserStart=Events(i).VarLaserstart;
                    LaserWidth=Events(i).VarLaserpulsewidth;
                    LaserNumPulse=Events(i).VarLasernumpulses;
                    LaserISI=Events(i).VarLaserisi;
                    
                    lsindex= find (LaserStarts== LaserStart);
                    lwindex= find( LaserWidths== LaserWidth);
                    lnpindex= find(LaserNumPulses== LaserNumPulse) ;
                    liindex=  find(LaserISIs == LaserISI);
                    
                    SS=1;
                    nreps_ssON(lsindex, lwindex, lnpindex, liindex)=nreps_ssON(lsindex, lwindex, lnpindex, liindex)+1;
                    SilentSoundON(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex)).spiketimes=spiketimes1; % Spike times
                    SilentSoundONspikecount(lsindex, lwindex, lnpindex, liindex,nreps_ssON(lsindex, lwindex, lnpindex, liindex))=spikecount; % No. of spikes
                    M1spontONSS(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                    M_LaserStartSS(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex))=LaserStart;
                    M_LaserWidthSS(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex))= LaserWidth;
                    M_LaserNumPulsesSS(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex))= LaserNumPulse;
                    M_LaserISISS(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex))= LaserISI;
                    if LaserRecorded
                        SilentSoundONLaser(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        SilentSoundONStim(lsindex, lwindex, lnpindex, liindex, nreps_ssON(lsindex, lwindex, lnpindex, liindex),:)=Stimtrace(region);
                    end
                end
            else
                SOA=Events(i).SOA;
                if isempty(SOA)
                    SOA=0;
                    probefreq=0;
                    probeamp=0;
                    probedur=0;
                end
                freq=Events(i).frequency;
                dur=Events(i).duration;
                amp=Events(i).amplitude;
                findex= find(freqs==freq);
                aindex= find(amps==amp);
                dindex= find(durs==dur);
                pfindex= find(probefreqs==probefreq);
                SOAindex= find(SOAs==SOA);
                nreps(findex, pfindex, aindex, SOAindex)+1;
                
                nrepsOFF(findex, pfindex, aindex, SOAindex)=nrepsOFF(findex, pfindex, aindex, SOAindex)+1;
                M1OFF(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex)).spiketimes=spiketimes1;
                M1OFFspikecounts(findex, pfindex, aindex, SOAindex,nrepsOFF(findex, pfindex, aindex, SOAindex))=spikecount;
                M1spontOFF(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex))=spont_spikecount;
                if LaserRecorded
                    M1OFFLaser(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex),:)=Lasertrace(region);
                end
                if StimRecorded
                    M1OFFStim(findex, pfindex, aindex, SOAindex, nrepsOFF(findex, pfindex, aindex, SOAindex),:)=Stimtrace(region);
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
for SOAindex=1:length(SOAs); % Hardcoded.
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            for pfindex=1:numprobefreqs
                % off
                spiketimesOFF=[];
                for rep=1:nrepsOFF(findex, pfindex, aindex, SOAindex)
                    spiketimesOFF=[spiketimesOFF M1OFF(findex, pfindex, aindex, SOAindex, rep).spiketimes];
                end
                mM1OFF(findex, pfindex, aindex, SOAindex).spiketimes=spiketimesOFF;
            end
        end
    end
end
mSilentSoundON=[];

spiketimesON=[];
lsindex=1;
lwindex=1;
for lnpindex=1:numLaserNumPulses
    for liindex=1:numVarLaserISIs
        for rep=1:nreps_ssON(lsindex, lwindex, lnpindex, liindex)
            spiketimesON=[spiketimesON SilentSoundON(lsindex, lwindex, lnpindex, liindex, rep).spiketimes];
        end
        mSilentSoundON(lnpindex, liindex).spiketimes=spiketimesON;
    end
end



%average laser and stimulus monitor M matrices across trials
if LaserRecorded
    for aindex=1:numamps
        for findex=1:numfreqs
            for pfindex=1:numprobefreqs
                for SOAindex=1:numSOAs
                    if nrepsOFF(findex, pfindex, aindex, SOAindex)>0
                        mM1OFFLaser(findex, pfindex, aindex, SOAindex,:)=mean(M1OFFLaser(findex, pfindex, aindex, SOAindex, 1:nrepsOFF(findex, pfindex, aindex, SOAindex),:), 5);
                    else %no reps for this stim, since rep=0
                        mM1OFFLaser(findex, pfindex, aindex, SOAindex,:)=zeros(size(region));
                    end
                end
            end
        end
    end
    lsindex=1;
    lwindex=1;
    for lnpindex=1:numLaserNumPulses
        for liindex=1:numVarLaserISIs
            for rep=1:nreps_ssON(lsindex, lwindex, lnpindex, liindex)
                mSilentSoundONLaser(lsindex, lwindex, lnpindex, liindex,:)=mean(SilentSoundONLaser(lsindex, lwindex, lnpindex, liindex, 1:nreps_ssON(lsindex, lwindex, lnpindex, liindex),:), 5);
            end
        end
    end
    
end

if StimRecorded 
    for aindex=1:numamps
        for findex=1:numfreqs
            for pfindex=1:numprobefreqs
                for SOAindex=1:numSOAs
                    if nrepsOFF(findex, pfindex, aindex, SOAindex)>0
                        mM1OFFStim(findex, pfindex, aindex, SOAindex,:)=mean(M1OFFStim(findex, pfindex, aindex, SOAindex, 1:nrepsOFF(findex, pfindex, aindex, SOAindex),:), 5);
                    else %no reps for this stim, since rep=0
                        mM1OFFStim(findex, pfindex, aindex, SOAindex,:)=zeros(size(region));
                    end
                end
            end
        end
    end
    lsindex=1;
    lwindex=1;
    for lnpindex=1:numLaserNumPulses
        for liindex=1:numVarLaserISIs
            for rep=1:nreps_ssON(lsindex, lwindex, lnpindex, liindex)
                mSilentSoundONStim(lsindex, lwindex, lnpindex, liindex,:)=mean(SilentSoundONStim(lsindex, lwindex, lnpindex, liindex, 1:nreps_ssON(lsindex, lwindex, lnpindex, liindex),:), 5);
            end
        end
    end
end


if exist('M1ONspikecounts')
    mM1ONspikecount=mean(M1ONspikecounts,4); % Mean spike count
    sM1ONspikecount=std(M1ONspikecounts,[],4); % Std of the above
    semM1ONspikecount=sM1ONspikecount./sqrt(max(nrepsON(:))); % Sem of the above
    
    % Spont
    mM1spontON=mean(M1spontON,4);
    sM1spontON=std(M1spontON,[],4);
    semM1spontON(:,:)=sM1spontON./sqrt(max(nrepsON(:)));
else
    mM1ONspikecount=[];
    sM1ONspikecount=[];
    semM1ONspikecount=[];
    M1ONspikecounts=[];
    mM1spontON=[];
    sM1spontON=[];
    semM1spontON=[];
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
    if exist('M1OFFspikecounts')
        mM1OFFspikecount=mean(M1OFFspikecounts,5); % Mean spike count
        sM1OFFspikecount=std(M1OFFspikecounts,[],5); % Std of the above
        semM1OFFspikecount=sM1OFFspikecount./sqrt(max(nrepsOFF(:))); % Sem of the above
        % Spont
        mM1spontOFF=mean(M1spontOFF,5);
        sM1spontOFF=std(M1spontOFF,[],5);
        semM1spontOFF=sM1spontOFF./sqrt(max(nrepsOFF(:)));
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
%save to outfiles
%one outfile for each cell
% previously existing outfiles for this cell will be overwritten

%after squeezing cluster, saves with the following dimensions:
% M1ON(findex, pfindex,aindex,SOAindex, nrepsON).spiketimes
% mM1ON(findex, pfindex,aindex,SOAindex).spiketimes
% mM1ONspikecount(findex, pfindex,aindex,SOAindex)

out.IL=IL;
out.Nclusters=Nclusters;
out.tetrode=channel;
out.channel=channel;
out.cluster=clust; %there are some redundant names here
out.cell=clust;
out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
out.LaserWidth=unique(LaserWidth);
out.LaserNumPulses=unique(LaserNumPulses);
out.LaserISIs=unique(LaserISIs);

if IL
    if exist('M_LaserStartSS')
        out.M_LaserStartSS=M_LaserStartSS;
        out.M_LaserWidthSS=M_LaserWidthSS;
        out.M_LaserNumPulsesSS=M_LaserNumPulsesSS;
        out.M_LaserISISS=M_LaserISISS;
    end
    
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
out.probefreqs=probefreqs;
out.SOAs=SOAs;
out.durs=durs;
out.numamps=numamps;
out.numfreqs=numfreqs;
out.numprobefreqs=numprobefreqs;
out.numdurs=numdurs;
out.numSOAs=numSOAs;
out.nreps=nreps;
out.nrepsOFF=nrepsOFF;
out.xlimits=xlimits;
out.samprate=samprate;
out.datadir=datadir;
out.spiketimes=spiketimes;
out.nreps_ssON=nreps_ssON;
out.SilentSoundON=SilentSoundON;
out.SilentSoundONspikecount=SilentSoundONspikecount;
out.mSilentSoundON=mSilentSoundON;

out.SilentSoundONStim=SilentSoundONStim;
out.mSilentSoundONStim=mSilentSoundONStim;
out.SilentSoundONLaser=SilentSoundONLaser;
out.mSilentSoundONLaser=mSilentSoundONLaser;
out.M1spontONSS=M1spontONSS;
out.M_LaserStartSS=M_LaserStartSS;
out.M_LaserWidthSS=M_LaserWidthSS;
out.M_LaserNumPulsesSS=M_LaserNumPulsesSS;
out.M_LaserISISS=M_LaserISISS;


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
out.t_filename=filename;
outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
save (outfilename, 'out', '-v7.3')

