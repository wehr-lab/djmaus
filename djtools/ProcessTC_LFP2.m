function ProcessTC_LFP2(varargin)
%processes continuous tuning curve data from djmaus.
%sorts Open Ephys continuous data into a big response matrix.
%using new OpenEphys file formats and hierarchy
%
%use for LFP, WC, or other continuous data (i.e. does not extract spikes)
%
%usage: ProcessTC_LFP2([datapath], [xlimits], [ylimits])
%all inputs are optional, defaults to current directory
%saves output in an outfile
%
% notes: whitenoise plotted as freq=-1 kHz, silent sound as -2 kHz
%       data are in microvolts


if nargin==0
    datadir=pwd;
    xlimits=[-20 100]; %x limits for axis
    ylimits=[-.1 .2];
elseif nargin==1
    datadir=varargin{1};
    xlimits=[-20 100]; %x limits for axis
    ylimits=[-.1 .2];
elseif nargin==2
    datadir=varargin{1};
    xlimits=varargin{2};
    ylimits=[-.1 .2];
elseif nargin==3
    datadir=varargin{1};
    xlimits=varargin{2};
    ylimits=varargin{3};
else
    error('wrong number of arguments');
end
if isempty(xlimits)
    xlimits=[-20 100]; %default xlimits
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cd(datadir)
try
    load dirs.mat
catch
    load bdirs.mat
end
BonsaiPath=Bdirs{1};
EphysPath=dirs{1};
tmp=split(macifypath(BonsaiPath), '/');
BonsaiFolder=tmp{end}; %remove absolute path
DataRoot=macifypath(DataRoot); %does nothing if you're on windows
try
    cd(DataRoot)
    cd(BonsaiFolder)
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
    cd(DataRoot)
    cd(BonsaiFolder)
    cd(EphysPath)
    sessionfilename=['session-',EphysPath];
    load(sessionfilename)
    fprintf('\nloaded session object')
catch
    warning('could not load session object')
end

%read messages
cd(DataRoot)
cd(BonsaiFolder)
behavior_filename=dir('Behavior_*.mat');
if isempty(behavior_filename)
    fprintf('\nno behavior file found, calling ProcessSession ')

    ProcessSession
    behavior_filename=dir('Behavior_*.mat');
end
load(behavior_filename.name);
fprintf('\nloaded behavior file. ')

if exist('Events.mat')
    load('Events.mat')
    fprintf('\nloaded Events file \n')
else
    [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps2(BonsaiPath);
    save('Events.mat','Events')
    save('StartAcquisitionSec.mat','StartAcquisitionSec')
end
if exist('StartAcquisitionSec.mat')
    load('StartAcquisitionSec.mat')
else
    [~, StartAcquisitionSec] = GetEventsAndSCT_Timestamps2(BonsaiPath);
    save('StartAcquisitionSec.mat','StartAcquisitionSec')
end

try
    fprintf('\nNumber of logged stimuli in notebook: %d', length(stimlog));
catch
    fprintf('\nCould not find stimlog, no logged stimuli in notebook!!');
end

%check if this is an appropriate stimulus protocol
% switch GetPlottingFunction(datadir)
%     case {'PlotTC_LFP', 'PlotTC_PSTH'}
%     otherwise
%         error('This does not appear to be a tuning curve stimulus protcol');
% end



%load OpenEphys session information
sessionfilename=['session-',EphysPath];
samplesfilename=['samples-',EphysPath];
cd(DataRoot)
cd(BonsaiFolder)
cd(EphysPath)
fprintf('\nloading Open Ephys data ...')
try
    tic
    load(sessionfilename)
    load(samplesfilename)
    fprintf(' done.  ')
    toc
catch
    fprintf('\ndid not find saved Open Ephys session object for this directory ...')
    fprintf('\nloading Open Ephys data (this will take a couple minutes) ...\n')
    tic
    session = Session(pwd); %open-ephys-matlab-tools function
    num_channels=session.recordNodes{1}.recordings{1}.info.continuous.num_channels;
    for ch=1:num_channels
        bit_volts(ch)=session.recordNodes{1}.recordings{1}.info.continuous.channels(ch).bit_volts;
    end
    keys=session.recordNodes{1}.recordings{1}.continuous.keys();
    key=keys{1};
    timestamps=session.recordNodes{1}.recordings{1}.continuous(key).timestamps;
    samples=session.recordNodes{1}.recordings{1}.continuous(key).samples(:,:);
    if num_channels==78 %for 64 channel neuronexus probe
        stimtracech=71;
        soundcardtriggerch=72;
        lasertracech=73;
    elseif num_channels==142 %for 128 channel diagnostic biochips probe
        stimtracech=135;
        soundcardtriggerch=136;
        lasertracech=137;
    elseif num_channels==136 %for 128 channel diagnostic biochips probe where we forgot to record the 6 AUX channels
        stimtracech=129;
        soundcardtriggerch=130;
        lasertracech=131;
    end
    stimtrace=session.recordNodes{1}.recordings{1}.continuous(key).samples(stimtracech,:);
    soundcardtrigger=session.recordNodes{1}.recordings{1}.continuous(key).samples(soundcardtriggerch,:);
    lasertrace=session.recordNodes{1}.recordings{1}.continuous(key).samples(lasertracech,:);
    session.recordNodes{1}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.spikes=[];
    fprintf('\nsaving Open Ephys session object in OpenEphys folder...')
    save(sessionfilename, 'session', 'stimtrace','soundcardtrigger','lasertrace','timestamps','num_channels', 'bit_volts')
    save(samplesfilename, 'samples', '-v7.3')
    fprintf(' done. ')
    toc
end

%try to load laser and stimulus monitor
if exist('lasertrace')==1
    %we loaded lasertrace from open ephys Session object
    LaserRecorded=1;
else    LaserRecorded=0;
    warning('Laser monitor channel not recorded')
end
if exist('stimtrace')==1
    %we loaded stimtrace from open ephys Session object
    StimRecorded=1;
else
    StimRecorded=0;
    warning('Stimulus monitor channel not recorded')
end

if LaserRecorded
    Lasertimestamps=timestamps-timestamps(1);
    Lasertrace=lasertrace./max(abs(lasertrace));
else
    fprintf('\nLaser trace not recorded')
end
if StimRecorded
    Stimtimestamps=timestamps-timestamps(1);
    Stimtrace=stimtrace./max(abs(stimtrace));
else
    fprintf('\nSound stimulus trace not recorded')
end




fprintf('\ncomputing tuning curve...');

samprate=Sky.OEsamplerate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'silentsound') ...
            |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise')| strcmp(Events(i).type, 'grating')
        j=j+1;
        alldurs(j)=Events(i).duration;
        if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, '2tone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).frequency;
        elseif strcmp(Events(i).type, 'whitenoise')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=-1000;
        elseif strcmp(Events(i).type, 'silentsound')
            allfreqs(j)=-2000;
            allamps(j)=nan; %flagging silent sound amp as nan
        elseif strcmp(Events(i).type, 'fmtone')
            allamps(j)=Events(i).amplitude;
            allfreqs(j)=Events(i).carrier_frequency;
        elseif strcmp(Events(i).type, 'grating')
            allfreqs(j)=Events(i).angle*1000;
            allamps(j)=Events(i).spatialfrequency;
        end
    end
end
allamps=allamps(~isnan(allamps)); %strip out nans from silent sound
freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
numfreqs=length(freqs);
numamps=length(amps);
numdurs=length(durs);

%check for laser in Events
%for now, setting AOPulseOn to 0 for all events
for i=1:length(Events)
    if isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
        LaserScheduled(i)=Events(i).laser; %whether the stim protocol scheduled a laser for this stim
        LaserOnOffButton(i)=Events(i).LaserOnOff; %whether the laser button was turned on
        LaserTrials(i)=LaserScheduled(i) & LaserOnOffButton(i);
    elseif isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %Not sure about this one. Assume no laser for now, but investigate.
        warning('ProcessTC_LFP: Cannot tell if laser button was turned on in djmaus GUI');
        LaserTrials(i)=0;
        Events(i).laser=0;%

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

M1=[];
M1ON=[];M1OFF=[];
M1ONLaser=[]; M1OFFLaser=[];
M1ONStim=[];M1OFFStim=[];

nreps=zeros(numfreqs, numamps, numdurs);
nrepsON=zeros(numfreqs, numamps, numdurs);
nrepsOFF=zeros(numfreqs, numamps, numdurs);

fprintf('\nextracting traces into BIG matrix...');

samprate=Sky.OEsamplerate;

channelorder=1:128; %default
probe='';
if num_channels==78 %for 64 channel neuronexus probe
    % ch 1-64 = headstage channels
    % ch  = A1_AUX channels 1-3
    % ch  = A2_AUX channels 1-3
    % ch  =  ADC channels ADC1-ADC8
elseif num_channels==142 | num_channels==136 %for 128 channel diagnostic biochips probe
    probe='P128-2';
    % ch 1-128 = headstage channels
    % ch 129-131 = A1_AUX channels 1-3
    % ch 132-134 = A2_AUX channels 1-3
    % ch 135-142 =  ADC channels ADC1-ADC8
    try
        load(fullfile(DataRoot, 'chanMap128.mat'))
        chans=1:128;
        channelorder=chanMap;
        out.channelmap_file=fullfile(DataRoot, 'chanMap128.mat');
        fprintf('\nusing channelmap %s', out.channelmap_file)
    catch
        warning('could not load channel map, defaulting to wrong order 1-128')
        chans=1:128;
        channelorder=1:128;
        out.channelmap_warning='could not load channel map, defaulting to wrong order 1-128';
    end
else
    warning(['unrecognized number of channels', int2str(num_channels)])
    warning('could not load channel map, defaulting to wrong order 1-128')
    out.channelmap_warning='could not load channel map, defaulting to wrong order 1-128';

end

%extract the traces into a big matrix M
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') |  strcmp(Events(i).type, 'silentsound') | ...
            strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating')
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

        if isempty(find(region<1)) %(disallow negative or zero start times)
            switch Events(i).type
                case {'tone', '2tone'}
                    freq=Events(i).frequency;
                    amp=Events(i).amplitude;
                case 'fmtone'
                    freq=Events(i).carrier_frequency;%
                    amp=Events(i).amplitude;
                case 'whitenoise'
                    freq=-1000;
                    amp=Events(i).amplitude;
                case 'silentsound'
                    freq=-2000;
                    amp=min(amps); %put silentsound in it's own column (freq=-2) in the lowest row
                case 'grating'
                    amp=Events(i).spatialfrequency;
                    freq=Events(i).angle*1000;
            end

            dur=Events(i).duration;
            findex= find(freqs==freq);
            aindex= find(amps==amp);
            dindex= find(durs==dur);
            nreps(findex, aindex, dindex)=nreps(findex, aindex, dindex)+1;
            M1(chans, findex,aindex,dindex, nreps(findex, aindex, dindex),:)=double(samples(chans, region))*bit_volts(1);
            M1stim(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=Stimtrace(region);
            if laser
                nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
                M1ON(chans, findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=double(samples(chans, region))*bit_volts(1);                M1ONLaser(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=Lasertrace(region);
                M1ONStim(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=Stimtrace(region);
            else
                nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
                M1OFF(chans, findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=double(samples(chans, region))*bit_volts(1);                M1OFFLaser(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Lasertrace(region);
                M1OFFStim(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Stimtrace(region);
            end
        end
    end
end

%region=length(M1OFF);
traces_to_keep=[];
if ~isempty(traces_to_keep)
    fprintf('\n using only traces %d to %d, discarding others\n', traces_to_keep(1),traces_to_keep(end));
    mM1=mean(M1(:,:,:,traces_to_keep,:), 4);
    if ~isempty(M1ON)
        mM1ON=mean(M1ON(:,:,:,traces_to_keep,:), 4);
    else
        mM1ON = [];
        mM1stim = [];
        mM1ONStim = [];
        mM1ONLaser = [];
        mM1OFFStim = [];
        mM1OFFLaser = [];
    end
    mM1OFF=mean(M1OFF(:,:,:,traces_to_keep,:), 4);
else
    for aindex=1:numamps
        for findex=1:numfreqs
            for dindex=1:numdurs
                if nreps(findex, aindex, dindex)>0
                    mM1(chans, findex, aindex, dindex,:)=mean(M1(chans, findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 5);
                    mM1stim(findex, aindex, dindex,:)=mean(M1stim(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1(findex, aindex, dindex,:)=zeros(size(region));
                    mM1stim(findex, aindex, dindex,:)=zeros(size(region));
                end
                if nrepsON(findex, aindex, dindex)>0
                    mM1ON(chans, findex, aindex, dindex,:)=mean(M1ON(chans, findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 5);
                    mM1ONStim(findex, aindex, dindex,:)=mean(M1ONStim(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                    mM1ONLaser(findex, aindex, dindex,:)=mean(M1ONLaser(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1ON(findex, aindex, dindex,:)=zeros(size(region));
                    mM1ONStim(findex, aindex, dindex,:)=zeros(size(region));
                    mM1ONLaser(findex, aindex, dindex,:)=zeros(size(region));
                end
                if nrepsOFF(findex, aindex, dindex)>0
                    mM1OFF(chans, findex, aindex, dindex,:)=mean(M1OFF(chans, findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 5);
                    mM1OFFStim(findex, aindex, dindex,:)=mean(M1OFFStim(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    mM1OFFLaser(findex, aindex, dindex,:)=mean(M1OFFLaser(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1OFF(findex, aindex, dindex,:)=zeros(size(region));
                    mM1OFFStim(findex, aindex, dindex,:)=zeros(size(region));
                    mM1OFFLaser(findex, aindex, dindex,:)=zeros(size(region));
                end
            end
        end
    end
end


%find optimal axis ylimits
if ylimits<0
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            trace=mM1(findex, aindex, dindex,:);
            ylimits(1)=min(ylimits(1), min(trace));
            ylimits(2)=max(ylimits(2), max(trace));
        end
    end
end

fprintf('\nsaving out file...');
tic
%assign outputs
out.M1=M1;
out.M1ON=M1ON;
out.M1ONStim=M1ONStim;
out.M1ONLaser=M1ONLaser;
out.mM1ON=mM1ON;
out.mM1ONStim=mM1ONStim;
out.mM1ONLaser=mM1ONLaser;
out.M1OFF=M1OFF;
out.M1OFFStim=M1OFFStim;
out.M1OFFLaser=M1OFFLaser;
out.mM1OFF=mM1OFF;
out.mM1OFFStim=mM1OFFStim;
out.mM1OFFLaser=mM1OFFLaser;
out.M1stim=M1stim;
out.mM1stim=mM1stim;
out.mM1=mM1;
out.datadir=datadir;
out.freqs=freqs;
out.amps=amps;
out.durs=durs;
out.nreps=nreps;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.numfreqs=numfreqs;
out.numamps=numamps;
out.numdurs=numdurs;
out.traces_to_keep=traces_to_keep;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.IL=IL;
out.xlimits=xlimits;
out.ylimits=ylimits;
out.samprate=samprate;
% out.nstd=nstd;
out.chans=chans;
out.channelorder=channelorder;
out.num_channels=num_channels;
out.probe=probe;

[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysPath=EphysPath;

%should probably save header info and stuff like that



outfilename=sprintf('outLFP.mat');
% save(outfilename, 'out')
save(outfilename, 'out', '-v7.3' )
fprintf('\n saved to %s', outfilename)
toc
