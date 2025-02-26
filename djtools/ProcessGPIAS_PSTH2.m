function ProcessGPIAS_PSTH2(varargin)

%processes clustered spiking tuning curve data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%creates only a single outfile per session
%(automatically processes all cells)
%-mike 07.19.2024
%
% usage: ProcessGPIAS_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, [xlimits],[ylimits])
%   SortedUnits is created by ProcessSpikes
% BonsaiPath, EphysDir, EphysDir_KS are folders, they are stored in OEinfo-xxx.mat in the Bonsai directory
% (xlimits, ylimits are optional)
% xlimits default to [0 200]
%
% saves to outfile
%


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
    s=GetStimParams(fullfile(BonsaiPath, EphysPath));
    if isempty(s) %in new OE format, notebook is one up from kilosorted data, so fileparts looks in ..
        s=GetStimParams(fileparts(datadir));
    end
    xlimits(1)=-3*max(s.gapdurs);
    xlimits(2)=800;
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
    case 'PlotGPIAS_PSTH'
        %that's fine
    otherwise
        error('This stimulus protcol does not appear to be GPIAS protocol.')
end

%read  Kilosort data
fprintf('\ncomputing tuning curve...');
samprate=Sky.OEsamplerate;

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
                j=j+1;
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
nrepsON=zeros( numgapdurs, numpulseamps);
nrepsOFF=zeros( numgapdurs, numpulseamps);

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
for i=1:length(Events)
    if isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
        if isempty(Events(i).laser)
            Events(i).laser=0;
        end
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
        warning('ProcessGPIAS_PSTH2: Cannot tell if laser button was turned on in djmaus GUI');
        LaserTrials(i)=0;
        Events(i).laser=0;
    elseif ~isfield(Events(i), 'laser') & ~isfield(Events(i), 'LaserOnOff')
        %if neither of the right fields are there, assume no laser
        LaserTrials(i)=0;
        Events(i).laser=0;
    elseif ~isfield(Events(i), 'laser') & isfield(Events(i), 'LaserOnOff')
        %if laser field is not there, assume no laser
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

%try to load laser and stimulus monitor
if exist('lasertrace')==1
    %we loaded lasertrace from open ephys Session object
        LaserRecorded=1;
elseif isempty(getLaserfile('.'))
    LaserRecorded=0;
    warning('Laser monitor channel not recorded')
    warning('ProcessGPIAS_PSTH2: still need to figure out how to find/plot Laser monitor channel ')
else
    LaserRecorded=1;
end
if exist('stimtrace')==1
    %we loaded stimtrace from open ephys Session object
        StimRecorded=1;
elseif isempty(getStimfile('.'))
    StimRecorded=0;
    warning('Stimulus monitor channel not recorded')
    warning('ProcessGPIAS_PSTH2: still need to figure out how to find/plot Stimulus monitor channel ')
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


%sizes:
% M1ON(numcells, numgapdurs, numpulseamps, nrepsON).spiketimes
% mM1ON(numcells, numgapdurs, numpulseamps).spiketimes
% mM1ONspikecount(numcells, numgapdurs, numpulseamps)
% Mstim: stimulus matrix in same format as M1

%first sort trains into matrix M1

M1ON=[];M1OFF=[];M1ONStim=[];M1OFFStim=[];

numcells=length(SortedUnits);
M1=[];M1ON=[];M1OFF=[];
mM1ON=[];
mM1OFF=[];
nreps=zeros(numcells, numgapdurs, numpulseamps);
nrepsON=zeros(numcells, numgapdurs, numpulseamps);
nrepsOFF=zeros(numcells, numgapdurs, numpulseamps);
nreps_ssON=zeros(numcells);
nreps_ssOFF=zeros(numcells);
%extract the traces into a big matrix M
j=0;
inRange=0;



%extract the traces into a big matrix M
for cellnum=1:numcells;
    spiketimes=[];
    %fprintf('\nreading SortedUnits cell %d/%d', cellnum, numcells)
    spiketimes= SortedUnits(cellnum).spiketimes;

    totalnumspikes=length(spiketimes);
    % fprintf('\nsuccessfully loaded spike data with %d spikes\n',totalnumspikes)

    for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') | strcmp(Events(i).type, 'gapinnoise')

            pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
            laser=LaserTrials(i);
        start=pos + gapdelay/1000 +xlimits(1)/1000; %start is in seconds
        stop=pos+ gapdelay/1000 + xlimits(2)/1000; %stop is in seconds
            region=round(start*samprate)+1:round(stop*samprate);
            if start>0 %(disallow negative or zero start times)
            gapdur=Events(i).gapdur;
            gdindex= find(gapdur==gapdurs);
            pulseamp=Events(i).pulseamp;
            paindex= find(pulseamp==pulseamps);
                st=spiketimes; %are in seconds
                spiketimes1=st(st>start & st<stop); % spiketimes in region
                spikecount=length(spiketimes1); % No. of spikes fired in response to this rep of this stim.
                inRange=inRange+ spikecount; %accumulate total spikecount in region
                spiketimes1=spiketimes1*1000 - pos*1000 - gapdelay;%convert to ms after tone onset
                spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window


                if laser
                    nrepsON(cellnum, gdindex,paindex)=nrepsON(cellnum, gdindex,paindex)+1;
                    M1ON(cellnum, gdindex,paindex, nrepsON(cellnum, gdindex,paindex)).spiketimes=spiketimes1; % Spike times
                    M1ONspikecount(cellnum, gdindex,paindex,nrepsON(cellnum, gdindex,paindex))=spikecount; % No. of spikes
                    M1spontON(cellnum, gdindex,paindex, nrepsON(cellnum, gdindex,paindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                    if LaserRecorded
                        M1ONLaser(gdindex,paindex, nrepsON(cellnum, gdindex,paindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        M1ONStim(gdindex,paindex, nrepsON(cellnum, gdindex,paindex),:)=Stimtrace(region);
                    end
                else
                    nrepsOFF(cellnum, gdindex,paindex)=nrepsOFF(cellnum, gdindex,paindex)+1;

                    M1OFF(cellnum, gdindex,paindex, nrepsOFF(cellnum, gdindex,paindex)).spiketimes=spiketimes1;
                    M1OFFspikecount(cellnum, gdindex,paindex,nrepsOFF(cellnum, gdindex,paindex))=spikecount;
                    M1spontOFF(cellnum, gdindex,paindex, nrepsOFF(cellnum, gdindex,paindex))=spont_spikecount;
                    if LaserRecorded
                        M1OFFLaser(gdindex,paindex, nrepsOFF(cellnum, gdindex,paindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        M1OFFStim(gdindex,paindex, nrepsOFF(cellnum, gdindex,paindex),:)=Stimtrace(region);
                    end

                                 end
            end
        end
    end

    fprintf('\ncell %d/%d, total num spikes: %d', cellnum, numcells, length(spiketimes))
    fprintf('\nIn range: %d', inRange)


    %accumulate across trials
for gdindex=1:numgapdurs; % Hardcoded.
    for paindex=1:numpulseamps
        spiketimesON=[];
        for rep=1:nrepsON(cellnum, gdindex,paindex)
            spiketimesON=[spiketimesON M1ON(cellnum, gdindex,paindex, rep).spiketimes];
        end
        mM1ON(cellnum, gdindex,paindex).spiketimes=spiketimesON;
        spiketimesOFF=[];
        for rep=1:nrepsOFF(cellnum, gdindex,paindex)
            spiketimesOFF=[spiketimesOFF M1OFF(cellnum, gdindex,paindex, rep).spiketimes];
        end
        mM1OFF(cellnum, gdindex,paindex).spiketimes=spiketimesOFF;
    end
end
for gdindex=1:numgapdurs; 
    for paindex=1:numpulseamps
        if IL
            %we don't need to save a copy for every cell, there's only one trace
        
            mM1ONStim(gdindex,paindex,:)=mean(M1ONStim(gdindex,paindex, 1:nrepsON(cellnum, gdindex,paindex),:), 3);
            mM1ONLaser(gdindex,paindex,:)=mean(M1ONLaser(gdindex,paindex, 1:nrepsON(cellnum, gdindex,paindex),:), 3);
        else
            mM1ONStim=[];
            mM1ONLaser=[];
        end
        mM1OFFStim(gdindex,paindex,:)=mean(M1OFFStim(gdindex,paindex, 1:nrepsOFF(cellnum, gdindex,paindex),:), 3);
        mM1OFFLaser(gdindex,paindex,:)=mean(M1OFFLaser(gdindex,paindex, 1:nrepsOFF(cellnum, gdindex,paindex),:), 3);
    end
end

    % %spikecount matrices
    if IL
        mM1ONspikecount=mean(M1ONspikecount, 4);
        sM1ONspikecount=std(M1ONspikecount, 0, 4);
        semM1ONspikecount=sM1ONspikecount./nrepsON(cellnum, gdindex,paindex);
        mM1spontON=mean(M1spontON, 4);
        sM1spontON=std(M1spontON, 0,4);
        semM1spontON=sM1spontON./nrepsON(cellnum, gdindex,paindex);
    else
        M1ONspikecount=[];
        mM1ONspikecount=[];
        sM1ONspikecount=[];
        semM1ONspikecount=[];
        mM1spontON=[];
        sM1spontON=[];
        semM1spontON=[];
    end
    mM1OFFspikecount=mean(M1OFFspikecount, 4);
    sM1OFFspikecount=std(M1OFFspikecount, 0, 4);
    semM1OFFspikecount=sM1OFFspikecount./nrepsOFF(cellnum, gdindex,paindex);
    mM1spontOFF=mean(M1spontOFF, 4);
    sM1spontOFF=std(M1spontOFF, 0,4);
    semM1spontOFF=sM1spontOFF./nrepsOFF(cellnum, gdindex,paindex);

   
  


end % for cellnum=1:length(SortedUnits);

%save to outfiles

%sizes:
% M1ON(numcells, numgapdurs, numpulseamps, nrepsON).spiketimes
% mM1ON(numcells, numgapdurs, numpulseamps).spiketimes
% mM1ONspikecount(numcells, numgapdurs, numpulseamps)

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

out.IL=IL;
out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;

if IL
    out.M1ON=M1ON;
    out.mM1ON=mM1ON;
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
    out.LaserStart=uniquelo(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);


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
if isempty(M1OFF)
    out.M1OFF=[];
    out.mM1OFF=[];
    out.mM1OFFspikecount=[];
    out.sM1OFFspikecount=[];
    out.semM1OFFspikecount=[];
    out.mM1spontOFF=[];
    out.sM1spontOFF=[];
    out.semM1spontOFF=[];
else
    out.M1OFF=M1OFF;
    out.mM1OFF=mM1OFF;
     out.mM1OFFspikecount=mM1OFFspikecount;
     out.sM1OFFspikecount=sM1OFFspikecount;
     out.semM1OFFspikecount=semM1OFFspikecount;
     out.mM1spontOFF=mM1spontOFF;
     out.sM1spontOFF=sM1spontOFF;
     out.semM1spontOFF=semM1spontOFF;
end
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.xlimits=xlimits;
out.samprate=samprate;

out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel

if LaserRecorded
    if exist('M1ONLaser')
        out.M1ONLaser=M1ONLaser;
        out.mM1ONLaser=mM1ONLaser;
    else
        out.M1ONLaser=[];
        out.mM1ONLaser=[];
    end
    if exist('M1OFFLaser')
        out.M1OFFLaser=M1OFFLaser;
        out.mM1OFFLaser=mM1OFFLaser;
    else
        out.M1OFFLaser=[];
        out.mM1OFFLaser=[];
    end
else
    out.M1ONLaser=[];
    out.mM1ONLaser=[];
    out.M1OFFLaser=[];
    out.mM1OFFLaser=[];
end
if StimRecorded
    if exist('M1ONStim')
        out.M1ONStim=M1ONStim;
        out.mM1ONStim=mM1ONStim;
    else
        out.M1ONStim=[];
        out.mM1ONStim=[];
    end
    if exist('M1OFFStim')
        out.M1OFFStim=M1OFFStim;
        out.mM1OFFStim=mM1OFFStim;
    else
        out.M1OFFStim=[];
        out.mM1OFFStim=[];
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


