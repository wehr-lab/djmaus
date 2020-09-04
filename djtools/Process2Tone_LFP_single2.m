function Process2Tone_LFP_single2(varargin)
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
djPrefs;
global pref

if nargin==0
    fprintf('\nno input\n')
    return
else
    datadir=varargin{1};
end
if nargin==1
    xlimits=[-300 1100]; %x limits for axis
    ylimits=[-.1 .2];
    prompt=('Please enter channel number: ');
    channel=input(prompt) ;
elseif nargin==2
    channel=varargin{2};
    xlimits=[-300 1100]; %x limits for axis
    ylimits=[-.1 .2];
elseif nargin==3
    channel=varargin{2};
    xlimits=varargin{3};
    ylimits=[-.1 .2];
elseif nargin==4
    channel=varargin{2};
    xlimits=varargin{3};
    ylimits=varargin{4};
else
    error('wrong number of arguments');
end
if ischar(channel); channel=str2num(channel);end

cd(pref.datapath)
cd(datadir)
filename=getContinuousFilename('.', channel);
if exist(filename, 'file')~=2 %couldn't find it
    filename=sprintf('114_CH%d.continuous', channel);
end
if exist(filename, 'file')~=2 %couldn't find it
    error(sprintf('could not find data file %s in datadir %s', filename, datadir))
end

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
if exist('Events.mat')
    load('Events.mat')
    try 
        load('StartAcquisitionSec.mat')
    catch
    [~, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
    end
else
[Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
save('Events.mat', 'Events')
save('StartAcquisitionSec.mat', 'StartAcquisitionSec')
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps
end
try
    fprintf('\nNumber of logged stimuli in notebook: %d', length(stimlog));
catch
    fprintf('\nCould not find stimlog, no logged stimuli in notebook!!');
end

[scaledtrace, datatimestamps, datainfo] =load_open_ephys_data(filename);
datatimestamps=datatimestamps-StartAcquisitionSec;
%should record a copy of stim
stim=0*scaledtrace;



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
    if strcmp(Events(i).type, '2tone') %get 2 tone paramaters: masker freq, probe freq, dur, masker amp and probe amp, SOAs
        j=j+1;
        alldurs(j)=Events(i).duration;
        allamps(j)=Events(i).amplitude;
        allfreqs(j)=Events(i).frequency;
        allprobefreqs(j)=Events(i).probefreq;
        allprobeamps(j)=Events(i).probeamp;
        allSOAs(j)=Events(i).SOA;
        LaserTrials(i)=0;
    elseif strcmp(Events(i).type, 'whitenoise')
        j=j+1;
        alldurs(j)=Events(i).duration;
        allamps(j)=Events(i).amplitude;
        allfreqs(j)=Events(i).frequency;
        LaserTrials(i)=0;
    elseif strcmp(Events(i).type, 'silentsound')
        j=j+1;
        alldurs(j)=Events(i).duration;
        allamps(j)=0;
        allfreqs(j)=0;
        allprobefreqs(j)=0;
        allprobeamps(j)=0;
        LaserTrials(i)=1;
        allVarLaserstart(j)=Events(i).VarLaserstart;
        allVarLaserpulsewidth(j)=Events(i).VarLaserpulsewidth;
        allVarLasernumpulses(j)=Events(i).VarLasernumpulses;
        allVarLaserisi(j)=Events(i).VarLaserisi;
        VarLaser=1;
    end
    

end

freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
SOAs=unique(allSOAs);
probefreqs=unique(allprobefreqs);
probeamps=unique(allprobeamps);

numfreqs=length(freqs);
numprobefreqs=length(probefreqs);
numamps=length(amps);
numdurs=length(durs);
numSOAs=length(SOAs);
numprobeamps=length(probeamps);

%check for laser
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
    numLaserISIs=length(LaserISIs);
end


fprintf('\n%d silent sound laser trials (single and double) in this Events file', sum(LaserTrials))
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

%load Laser
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


%Load stimulus
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


M1=[];M1ON=[];M1=[];

nreps=zeros(numfreqs, numLaserNumPulses, numLaserISIs);
nreps2=zeros(numfreqs, numprobefreqs, numSOAs);
if size(nreps)~=size(nreps2)
    error('laser and wn parameters dont match')
end


%extract the traces into a big matrix M
j=0;
inRange=0;
M1ONLaser=[];
M1ONStim=[];
if exist('Events.mat')
    load('Events.mat')
end
t=0; w=0; m=0; n=0;
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') | ...
            strcmp(Events(i).type, 'silentsound') | strcmp(Events(i).type, 'whitenoise')%%%I should pull out silent sound processing into separate stanza
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate; %pos is in samples
        end
        laser=LaserTrials(i);
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        if start>0 %(disallow negative or zero start times)
          
            if strcmp(Events(i).type, '2tone') || strcmp(Events(i).type, 'whitenoise')
                freq=-1;
                probefreq=Events(i).probefreq;
                SOA=Events(i).SOA;
                t=t+1;
                if strcmp(Events(i).type, 'whitenoise')
                    probefreq=0;
                    SOA=0;
                    w=w+1;
                end
                findex= find(freqs==freq);
                pfindex= find(probefreqs==probefreq);
                SOAindex= find(SOAs==SOA);
                
                
                nreps(findex, pfindex, SOAindex)=nreps(findex, pfindex, SOAindex)+1;
                M1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex),:)=scaledtrace(region);;
               
                if LaserRecorded
                    M1Laser(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex),:)=Lasertrace(region);
                end
                if StimRecorded
                    %figure(1); hold on;plot(Stimtrace(region));
                    M1Stim(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex),:)=Stimtrace(region);
                end
                
            elseif strcmp(Events(i).type, 'silentsound')
                LaserStart=Events(i).VarLaserstart;
                LaserWidth=Events(i).VarLaserpulsewidth;
                LaserNumPulse=Events(i).VarLasernumpulses;
                LaserISI=Events(i).VarLaserisi;
                freq=0;
                findex= find(freqs==freq);
                
                
                %reverse N of laser pulses to match pfindex 1=2tone, 2-single
                %should be lnindex=1 = 2tone if lnindex=2 single
                lnpindex= find(LaserNumPulses== LaserNumPulse);
                if lnpindex==1
                    lnpindex=2;
                    LaserISI=0;
                    n=n+1;
                elseif lnpindex==2
                    lnpindex=1;
                    m=m+1;
                end
                liindex=  find(LaserISIs == LaserISI);
                
                nreps(findex, lnpindex, liindex)=nreps(findex, lnpindex, liindex)+1;
                M1(findex,lnpindex, liindex, nreps(findex, lnpindex, liindex),:)=scaledtrace(region);
                M_LaserStart(findex, lnpindex, liindex, nreps(findex, lnpindex, liindex))=LaserStart;
                M_LaserWidth(findex, lnpindex, liindex, nreps(findex, lnpindex, liindex))= LaserWidth;
                M_LaserNumPulses(findex, lnpindex, liindex, nreps(findex, lnpindex, liindex))= LaserNumPulse;
                M_LaserISI(findex, lnpindex, liindex, nreps(findex, lnpindex, liindex))= LaserISI;
                if LaserRecorded
                    %figure(2); hold on; plot(Lasertrace(region));
                    M1Laser(findex,lnpindex, liindex, nreps(findex, lnpindex, liindex),:)=Lasertrace(region);
                end
                if StimRecorded
                    M1Stim(findex,lnpindex, liindex, nreps(findex, lnpindex, liindex),:)=Stimtrace(region);
                end
                
            end
            
        end
    end
end
fprintf('\nmin reps: %d\nmax reps: %d',min(nreps(:)), max(nreps(:)))

mM1=[];

% Accumulate spiketimes across trials, for psth...
for SOAindex=1:numSOAs; % or numLaserISIs
    for findex=1:numfreqs %-1 or 0
        for pfindex=1:numprobefreqs % or numLaserPulses, for 2tone 1 means one WN burst, for laser 1 mean one laser pulse
            mM1(findex, pfindex, SOAindex,:)=mean(M1(findex, pfindex, SOAindex, 1:nreps(findex, pfindex, SOAindex),:),4);
        end
    end
    
end

%average laser and stimulus monitor M matrices across trials
if LaserRecorded
    for findex=1:numfreqs
        for pfindex=1:numprobefreqs
            for SOAindex=1:numSOAs
                if nreps(findex, pfindex, SOAindex)>0
                    mM1Laser(findex, pfindex, SOAindex,:)=mean(M1Laser(findex, pfindex, SOAindex, 1:nreps(findex, pfindex, SOAindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1Laser(findex, pfindex, SOAindex,:)=zeros(size(region));
                end
            end
        end
    end
end


if StimRecorded
    for findex=1:numfreqs
        for pfindex=1:numprobefreqs
            for SOAindex=1:numSOAs
                if nreps(findex, pfindex, SOAindex)>0
                    mM1Stim(findex, pfindex, SOAindex,:)=mean(M1Stim(findex, pfindex, SOAindex, 1:nreps(findex, pfindex, SOAindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1Stim(findex, pfindex, SOAindex,:)=zeros(size(region));
                end
            end
        end
    end
end



%save to outfiles

out.IL=IL;
out.channel=channel;

%spikes
out.M1=M1;
out.mM1=mM1;


out.amps=amps;
out.freqs=freqs;
out.probefreqs=probefreqs;
out.probeamps=probeamps;
out.durs=durs;
out.SOAs=SOAs;

out.numamps=numamps;
out.numfreqs=numfreqs;
out.numprobefreqs=numprobefreqs;
out.numprobeamps=numprobeamps; %not used to parse data but just in case
out.numdurs=numdurs;
out.numSOAs=numSOAs;

out.nreps=nreps;
out.xlimits=xlimits;
out.samprate=samprate;
out.datadir=datadir;

out.LaserStarts=LaserStarts;
out.LaserWidths=LaserWidths;
out.LaserNumPulses=LaserNumPulses;
out.LaserISIs=LaserISIs;
out.numLaserStarts=numLaserStarts;
out.numLaserWidths=numLaserWidths;
out.numLaserNumPulses=numLaserNumPulses;
out.numLaserISIs=numLaserISIs;

out.M_LaserStart=M_LaserStart;
out.M_LaserWidth=M_LaserWidth;
out.M_LaserNumPulses=M_LaserNumPulses;
out.M_LaserISI=M_LaserISI;

out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel
if LaserRecorded
    if exist('M1Laser')
        out.M1Laser=M1Laser;
    end
    if exist('mM1Laser')
        out.mM1Laser=mM1Laser;
    end
else
    out.M1Laser=[];
    out.mM1Laser=[];
end
if StimRecorded
    if exist('M1Stim')
        out.M1Stim=M1Stim;
    end
    if exist('mM1Stim')
        out.mM1Stim=mM1Stim;
    end
else
    out.M1Stim=[];
    out.mM1Stim=[];
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
outfilename=sprintf('outLFP_ch%d.mat',channel);
save(outfilename, 'out')
fprintf('\n saved to %s', outfilename)

