function ProcessGPIASflashtrain_PSTH_single(varargin)

%processes a single .t file of clustered spiking GPIAS data from djmaus
%
% usage: ProcessGPIAS_PSTH_single(datadir, t_filename, [xlimits],[ylimits])
% (xlimits, ylimits are optional)
% xlimits default to [-1.5*max(gapdurs) 2*soa]
% saves to outfile

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};

try
    xlimits=varargin{3};
catch
    xlimits=[];
end
if isempty(xlimits)
    xlimits=[];
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
if exist(messagesfilename, 'file')~=2
    error('Could not find messages.events file. Is this the right data directory?')
end
[messages] = GetNetworkEvents(messagesfilename);


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


%read MClust .t file or Kilosort

if exist('params.py','file') || exist('dirs.mat','file')
    fprintf('\nreading KiloSort output cell %d', clust)
    spiketimes=readKiloSortOutput(clust, sampleRate);
else
    fprintf('\nreading MClust output file %s', filename)
    spiketimes=read_MClust_output(filename)'/10000; %spiketimes now in seconds
    %correct for OE start time, so that time starts at 0
    spiketimes=spiketimes-StartAcquisitionSec;
    fprintf('\nsuccessfully loaded MClust spike data')
end
totalnumspikes=length(spiketimes);

Nclusters=1;

%%%uncomment this to run some sanity checks
%SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)

fprintf('\ncomputing tuning curve...');
samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS')
        j=j+1;
        allsoas(j)=Events(i).soa;
        allgapdurs(j)=Events(i).gapdur;
        allgapdelays(j)=Events(i).gapdelay;
        allpulseamps(j)=Events(i).pulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allnoiseamps(j)=Events(i).amplitude;
    elseif strcmp(Events(i).type, 'gapinnoise')
        allsoas(j)=Events(i).soa;
        allgapdurs(j)=Events(i).gapdur;
        allgapdelays(j)=Events(i).gapdelay;
        allpulseamps(j)=Events(i).pulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allnoiseamps(j)=Events(i).amplitude;
    end
    
end
gapdurs=unique(allgapdurs);
pulsedurs=unique(allpulsedurs);
soas=unique(allsoas);
gapdelays=unique(allgapdelays);
pulseamps=unique(allpulseamps);
pulsedurs=unique(allpulsedurs);
noiseamps=unique(allnoiseamps);
numgapdurs=length(gapdurs);
numpulseamps=length(pulseamps);

if length(noiseamps)~=1
    error('not able to handle multiple noiseamps')
end
if length(gapdelays)~=1
    error('not able to handle multiple gapdelays')
end
if length(pulsedurs)~=1
    error('not able to handle multiple pulsedurs')
end
if length(soas)~=1
    error('not able to handle multiple soas')
end
noiseamp=noiseamps;
soa=soas;
pulsedur=pulsedurs;
gapdelay=gapdelays;

%check for laser in Events
LaserScheduled = zeros(1,length(Events));
LaserOnOffButton = zeros(1,length(Events));
LaserTrials = zeros(1,length(Events));

allVarLaserstarts=[];
alltrainnumpulses=[];
alltrainisis=[];
alltrainpulsewidths=[];
allpulsewidths=[];
for i=1:length(Events)
    VarLaser(i)=Events(i).VarLaser;
    if Events(i).VarLaser
        %        Events(i).VarLaserstart
        %        Events(i).VarLaserpulsewidth
        %        Events(i).VarLaserisi
        %        Events(i).VarLasernumpulses
        allVarLaserstarts=[allVarLaserstarts Events(i).VarLaserstart];
        allpulsewidths=[ allpulsewidths Events(i).VarLaserpulsewidth];
        alltrainisis=[alltrainisis Events(i).VarLaserisi]; %isi in train
        alltrainnumpulses=[alltrainnumpulses Events(i).VarLasernumpulses]; %isi in train
        LaserPulsewidth(i)=Events(i).VarLaserpulsewidth;
    else
        LaserPulsewidth(i)=0;
        allpulsewidths=[ allpulsewidths 0];

    end
end
laserstarts=unique(allVarLaserstarts);
trainnumpulses=unique(alltrainnumpulses);
trainisis=unique(alltrainisis);
pulsewidths=unique(allpulsewidths);
numpulsewidths=length(pulsewidths);

% for i=1:length(Events)
%     if isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
%         if ~isempty(Events(i).laser)
%             LaserScheduled(i)=Events(i).laser; %whether the stim protocol scheduled a laser for this stim
%         end
%         if ~isempty(Events(i).LaserOnOff)
%             LaserOnOffButton(i)=Events(i).LaserOnOff; %whether the laser button was turned on
%         end
%         LaserTrials(i)=LaserScheduled(i) & LaserOnOffButton(i);
%         if isempty(stimlog(i).LaserStart)
%             LaserStart(i)=nan;
%             LaserWidth(i)=nan;
%             LaserNumPulses(i)=nan;
%             LaserISI(i)=nan;
%         else
%             LaserStart(i)=stimlog(i).LaserStart;
%             LaserWidth(i)=stimlog(i).LaserWidth;
%             LaserNumPulses(i)=stimlog(i).LaserNumPulses;
%             LaserISI(i)=stimlog(i).LaserISI;
%         end
%         
%     elseif isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
%         %Not sure about this one. Assume no laser for now, but investigate.
%         warning('ProcessGPIAS_PSTH_single: Cannot tell if laser button was turned on in djmaus GUI');
%         LaserTrials(i)=0;
%         Events(i).laser=0;
%     elseif ~isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
%         %if neither of the right fields are there, assume no laser
%         LaserTrials(i)=0;
%         Events(i).laser=0;
%     else
%         error('wtf?')
%     end
% end
fprintf('\n%d laser trials in this Events file', sum(VarLaser))
% try
%     if sum(LaserOnOffButton)==0
%         fprintf('\nLaser On/Off button remained off for entire file.')
%     end
% end
% if sum(LaserTrials)>0
     IL=1;
% else
%     IL=0;
% end
%if lasers were used, we'll un-interleave them and save ON and OFF data
%try to load laser and stimulus monitor files
if isempty(getLaserfile('.'))
    LaserRecorded=0;
else
    LaserRecorded=1;
end
if isempty(getStimfile('.'))
    StimRecorded=0;
else
    StimRecorded=1;
end

if LaserRecorded
    try
        [Lasertrace, Lasertimestamps, Laserinfo] =load_open_ephys_data(getLaserfile('.'));
        Lasertimestamps=Lasertimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
        Lasertrace=Lasertrace./max(abs(Lasertrace));
        fprintf('\nsuccessfully loaded laser trace\n')
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

%there is only 1 matrix, containing both laser on and off
%pulsewidth 0 is where we put the laser-off trials
M1=[];
nreps=zeros(numgapdurs, numpulsewidths); %I don't think we need to add 1

% %find optimal axis limits
if isempty(xlimits)
    xlimits(1)=-1.5*max(gapdurs);
    xlimits(2)=2*soa;
end
fprintf('\nprocessing with xlimits [%d-%d]', xlimits(1), xlimits(2))

%extract the traces into a big matrix M
j=0;
inRange=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') | strcmp(Events(i).type, 'gapinnoise')
        
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=VarLaser(i);
        start=pos + gapdelay/1000 +xlimits(1)/1000; %start is in seconds
        stop=pos+ gapdelay/1000 + xlimits(2)/1000; %stop is in seconds
        region=round(start*samprate)+1:round(stop*samprate);
        if start>0 %(disallow negative or zero start times)
            gapdur=Events(i).gapdur;
            gdindex= find(gapdur==gapdurs);
            pulsewidth=LaserPulsewidth(i);
            pwindex= find(pulsewidth==pulsewidths);
            st=spiketimes; %are in seconds
            st_inrange=st(st>start & st<stop); % spiketimes in region, in seconds relative to start of acquisition
            spikecount=length(st_inrange); % No. of spikes fired in response to this rep of this stim.
            inRange=inRange+ spikecount; %accumulate total spikecount in region
            spiketimes1=st_inrange*1000 - pos*1000 - gapdelay;%covert to ms after gap termination
            spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window
            nreps(gdindex,pwindex)=nreps(gdindex,pwindex)+1;
            M1(gdindex,pwindex, nreps(gdindex,pwindex)).spiketimes=spiketimes1; % Spike times
            M1spikecounts(gdindex,pwindex,nreps(gdindex,pwindex))=spikecount; % No. of spikes
            M1spont(gdindex,pwindex, nreps(gdindex,pwindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
           % M_LaserStart(gdindex,pwindex, nreps(gdindex,pwindex))=Events(i).VarLaserstart;
           % M_LaserPulsewidth(gdindex,pwindex, nreps(gdindex,pwindex))= Events(i).VarLaserpulsewidth;
           % M_LaserNumPulses(gdindex,pwindex, nreps(gdindex,pwindex))= Events(i).VarLasernumpulses;
           % M_LaserISI(gdindex,pwindex, nreps(gdindex,pwindex))= Events(i).VarLaserisi;
            if LaserRecorded
                M1Laser(gdindex, pwindex, nreps(gdindex,pwindex),:)=Lasertrace(region);
            end
            if StimRecorded
                M1Stim(gdindex,pwindex, nreps(gdindex,pwindex),:)=Stimtrace(region);
            end
        end
    end
end

fprintf('\nmin num reps: %d\nmax num reps: %d', min(nreps(:)), max(nreps(:)))
fprintf('\ntotal num spikes: %d', length(spiketimes))
fprintf('\nIn range: %d', inRange)

% Accumulate spiketimes across trials, for psth...
for gdindex=1:numgapdurs; % Hardcoded.
    for pwindex=1:numpulsewidths
        % on
        spiketimes=[];
        for rep=1:nreps(gdindex,pwindex)
            spiketimes=[spiketimes M1(gdindex,pwindex, rep).spiketimes];
        end
        
        % All spiketimes for a given f/a/d combo, for psth:
        mM1(gdindex,pwindex).spiketimes=spiketimes;
        
    end
end



mM1spikecount=mean(M1spikecounts,3); % Mean spike count
sM1spikecount=std(M1spikecounts,[],3); % Std of the above
semM1spikecount=sM1spikecount./sqrt(max(nreps(:))); % Sem of the above
% Spont
mM1spont=mean(M1spont,3);
sM1spont=std(M1spont,[],3);
semM1spont=sM1spont./sqrt(max(nreps(:)));




%average laser and stimulus monitor M matrices across trials
if LaserRecorded
    for gdindex=1:numgapdurs
        for pwindex=1:numpulsewidths
            if nreps(gdindex,pwindex)>0
                mM1Laser(gdindex, pwindex, :)=mean(M1Laser(gdindex,pwindex, 1:nreps(gdindex,pwindex),:), 3);
            else %no reps for this stim, since rep=0
                mM1Laser(gdindex,pwindex,:)=zeros(size(region));
            end
        end
    end
end
if StimRecorded
    for gdindex=1:numgapdurs
        for pwindex=1:numpulsewidths
            if nreps(gdindex,pwindex)>0
                mM1Stim(gdindex,pwindex,:)=mean(M1Stim(gdindex,pwindex, 1:nreps(gdindex,pwindex),:), 3);
            else %no reps for this stim, since rep=0
                mM1Stim(gdindex,pwindex,:)=zeros(size(region));
            end
        end
    end
end

%sanity check - are the stimuli where we think they are?
if StimRecorded
    figure %
    hold on
    offset2=0;
    offset=1.5*range(M1Stim(:));
    for gdindex=1:numgapdurs
        for pwindex =1:numpulsewidths
            for r=1:nreps(gdindex,pwindex)
                stim=squeeze(M1Stim(gdindex,pwindex,r,:));
                t=1:length(stim);t=1000*t/samprate; %in ms
                t=t+xlimits(1); %correct for xlim
                offset2=offset2+offset;
                plot(t, stim+offset2, 'm')
            end
        end
    end
    title(' stimulus monitor')

end

%save to outfiles
%one outfile for each cell

%saves with the following dimensions:
% M1(numgapdurs, numpulsewidths, nreps).spiketimes
% mM1(numgapdurs, numpulsewidths).spiketimes
% mM1spikecount(numgapdurs, numpulsewidths)
%note that numpulsewidths includes 0, the laser-off condition
%thus there is no M1 OFF

out.IL=IL;
out.Nclusters=Nclusters;
out.tetrode=channel;
out.channel=channel;
out.cluster=clust; %there are some redundant names here
out.cell=clust;
out.M1=M1; %isn't this so much easier?
out.mM1=mM1;
out.mM1spikecount=mM1spikecount;
out.sM1spikecount=sM1spikecount;
out.semM1spikecount=semM1spikecount;
out.mM1spont=mM1spont;
out.sM1spont=sM1spont;
out.semM1spont=semM1spont;
% out.M_LaserStart=M_LaserStart;
% out.M_LaserPulsewidth=M_LaserPulsewidth;
% out.M_LaserNumPulses=M_LaserNumPulses;
% out.M_LaserISI=M_LaserISI;

out.laserstarts=laserstarts;
out.trainnumpulses=trainnumpulses;
out.trainisis=trainisis;
out.pulsewidths=pulsewidths;
out.numpulsewidths = numpulsewidths;

out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;
out.nreps=nreps;
out.xlimits=xlimits;
out.samprate=samprate;
out.datadir=datadir;
out.spiketimes=spiketimes;

out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel

if LaserRecorded
    if exist('M1Laser')
        out.M1Laser=M1Laser;
        out.mM1Laser=mM1Laser;
    else
        out.M1Laser=[];
        out.mM1Laser=[];
    end
else
    out.M1Laser=[];
    out.mM1Laser=[];
end
if StimRecorded
    if exist('M1Stim')
        out.M1Stim=M1Stim;
        out.mM1Stim=mM1Stim;
    else
        out.M1Stim=[];
        out.mM1Stim=[];
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
outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
save (outfilename, 'out')
fprintf('\nsaved %s', outfilename)





