function ProcessGPIAS_Behavior4Mice(BonsaiPath, EphysFolder, flag_accel)

%processes gap detection behavioral data from djmaus
%
% usage: ProcessGPIAS_Behavior4Mice2(BonsaiPath, EphysFolder, flag_accel))
% saves to separate outfiles for each mouse
% Use this code for 4 mouse simultaneous gap detection on rig 2
%
% using new OpenEphys and kilosort file formats and hierarchy -mike 08.21.25
% I deleted or commented out ACC code because I don't think we're using the
% headstage accelerometer channels to compute startle responses at all
% anymore. I also changed it to write all 4 outfiles as
% outGPIAS_BehaviorMouse#.mat whereas before the 1 was omitted for mouse1 

if nargin < 2 || isempty(flag_accel)
    flag.accel = 4;
else
    flag.accel = flag_accel;
end

flag.includeALL = 1;     % 0) eliminate traces with too much activity before trial
flag.unwrap = 0;
flag.plot = 1;

if nargin==0
    try
        load dirs.mat
    catch
        try
            load bdirs.mat
        catch
            ProcessSession
            load bdirs.mat
        end
    end

    BonsaiPath=Bdirs{1};
    EphysFolder=dirs{1};
end

% djPrefs;
% global pref
% cd (pref.datapath);
% cd(datadir)

try
    cd(BonsaiPath)
    cd(EphysFolder)
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

if flag.includeALL
    fprintf('using all traces\n')
else
    fprintf('eliminating traces with too much activity before test\n')
end

try
    cd(BonsaiPath)
    cd(EphysFolder)
    sessionfilename=['session-',EphysFolder];
    load(sessionfilename)
    fprintf('\nloaded session object')
catch
    warning('could not load session object')
end

%read messages
cd(BonsaiPath)
behavior_filename=dir('Behavior_*.mat');
try 
    load(behavior_filename.name);
catch
    warning('could not load behavior file, not sure if we need something from that file')
end

%get Events and soundcard trigger timestamps
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps

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
if ~strcmp(GetPlottingFunction(BonsaiPath), 'PlotGPIAS_PSTH')
    error('This does not appear to be a GPIAS stimulus protcol');
end

[~,BonsaiFolder,~]=fileparts(BonsaiPath);
OEinfofilename=sprintf('OEinfo-%s', BonsaiFolder);
load(fullfile(BonsaiPath, OEinfofilename))

samprate=OEsamplerate;


%load OpenEphys session information
sessionfilename=['session-',EphysFolder, '.mat'];
samplesfilename=['samples-',EphysFolder, '.mat'];
cd(DataRoot)
cd(BonsaiFolder)
cd(EphysFolder)
fprintf('\nloading Open Ephys data ...')
try
    tic
    load(sessionfilename)
    fprintf('\nfound and loaded Open Ephys session object.')

    load(samplesfilename)
    fprintf('\nfound and loaded Open Ephys samples file.')
    fprintf(' done.  ')
    toc
catch
    if ~exist(sessionfilename)
        fprintf('\ndid not find saved Open Ephys session object for this directory ...')
    end
    if ~exist(samplesfilename)
        fprintf('\ndid not find saved Open Ephys samples file for this directory ...')
    end
    fprintf('\nloading Open Ephys data (this will take a couple minutes) ...\n')
    tic
    session = Session(pwd); %open-ephys-matlab-tools function
    num_channels=session.recordNodes{1}.recordings{1}.info.continuous(1).num_channels;
    for ch=1:num_channels
        bit_volts(ch)=session.recordNodes{1}.recordings{1}.info.continuous(1).channels(ch).bit_volts;
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
    elseif num_channels==11 %Rig2 gap detection behavior, one config
        stimtracech=4;
        soundcardtriggerch=5;
        lasertracech=6; %I think, but haven't confirmed
        %piezo data from mice 1-4 are on chans 8,9,10,11. Chans 1-3 are empty
        mouse1ch=8;
        mouse2ch=9;
        mouse3ch=10;
        mouse4ch=11;
    elseif num_channels==27 %Rig2 gap detection behavior, another config
        stimtracech=20;
        soundcardtriggerch=21;
        lasertracech=22; %I think, but haven't confirmed
        %piezo data from mice 1-4 should be on chans 24,25,26,27. Chans 1-3 are empty
        mouse1ch=24;
        mouse2ch=25;
        mouse3ch=26;
        mouse4ch=27;
    end
    stimtrace=session.recordNodes{1}.recordings{1}.continuous(key).samples(stimtracech,:);
    soundcardtrigger=session.recordNodes{1}.recordings{1}.continuous(key).samples(soundcardtriggerch,:);
    lasertrace=session.recordNodes{1}.recordings{1}.continuous(key).samples(lasertracech,:);
    session.recordNodes{1}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.continuous=[];
    session.recordNodes{2}.recordings{1}.spikes=[];
    fprintf('\nsaving Open Ephys session object and samples in OpenEphys folder...')
    save(sessionfilename, 'session', 'stimtrace','soundcardtrigger','lasertrace',...
        'timestamps','num_channels', 'bit_volts', '*ch')
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
    stimtrace=stimtrace./max(abs(stimtrace));
else
    fprintf('\nSound stimulus trace not recorded')
end

mouse1scaledtrace=double(samples(mouse1ch, :))*bit_volts(mouse1ch);
mouse2scaledtrace=double(samples(mouse2ch, :))*bit_volts(mouse1ch);
mouse3scaledtrace=double(samples(mouse3ch, :))*bit_volts(mouse1ch);
mouse4scaledtrace=double(samples(mouse4ch, :))*bit_volts(mouse1ch);


%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') | strcmp(Events(i).type, 'toneGPIAS')
        j=j+1;
        allsoas(j)=Events(i).soa;
        allsoaflags{j}=Events(i).soaflag;
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
soaflags=unique(allsoaflags);
gapdelays=unique(allgapdelays);
pulseamps=unique(allpulseamps);
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
if length(soaflags)~=1
    error('not able to handle multiple soaflags')
end
noiseamp=noiseamps;
soa=soas;
pulsedur=pulsedurs;
gapdelay=gapdelays;
soaflag=soaflags{:};


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
        warning('ProcessGPIAS_Behavior: Cannot tell if laser button was turned on in djmaus GUI');
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

M1ON=[];M1OFF=[];
M1ONACC=[];M1OFFACC=[];
M1ONstim=[];M1OFFstim=[];
nrepsON=zeros(numgapdurs, numpulseamps);
nrepsOFF=zeros(numgapdurs, numpulseamps);

xlimits=[-350 350]; %xlimits for storing traces
startle_window=[0 100]; %hard coded integration region for startle response (ms)
bl_window=[-100 0];%hard coded integration region for pre-startle baseline (ms)
fprintf('\nprocessing with xlimits [%d - %d]', xlimits(1), xlimits(2))
fprintf('\nprocessing with startle integration window [%d - %d]', startle_window(1), startle_window(2))

%extract the traces into a big matrix M
j=0;
if flag.plot
    for idur = 1:length(gapdurs)
        Hfig(idur) = figure;
        hold on;
        title(['raw traces for gapdur = ' num2str(gapdurs(idur)) ' ms'])
    end
end

for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') |  strcmp(Events(i).type, 'toneGPIAS')
        %note: gapdelay is the time from the soundcardtrigger (pos) to the
        %gap termination. The time to startle onset should be
        %(gapdelay + soa) after pos, if soaflag=soa
        gapdur=Events(i).gapdur;
        switch soaflag
            case 'isi'
                isi=Events(i).soa;
                soa=isi+gapdur;
            case 'soa'
                soa=Events(i).soa;
                isi=soa-gapdur;
        end
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=LaserTrials(i);
        %since this is all about quantifying startle response, we want a trace
        %locked to startle pulse (not gap)
        
        startle_onset=pos+gapdelay/1000+isi/1000;
        start=startle_onset + xlimits(1)/1000; %start is in seconds, should be at xlimits relative to startle onset
        stop=startle_onset + xlimits(2)/1000; %stop is in seconds
        if start>0 %(disallow negative or zero start times)
            gdindex= find(gapdur==gapdurs);
            pulseamp=Events(i).pulseamp;
            paindex= find(pulseamp==pulseamps);
            %start=round(pos+xlimits(1)*1e-3*samprate);
            %stop=round(pos+xlimits(2)*1e-3*samprate)-1;
            region=round(start*samprate)+1:round(stop*samprate);
            if isempty(find(region<1))
                if laser
                    nrepsON(gdindex,paindex)=nrepsON(gdindex,paindex)+1;
                    M1ON(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse1scaledtrace(region);
                    %M1ONACC(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse1scaledtraceACC(region);
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stimtrace(region);
                else
                    if flag.unwrap
                        temp = mouse1scaledtrace(region);
                        tempD = [0; diff(temp)];
                        indUp = find(tempD>2);
                        indDown = find(tempD<-2)-1;
                        njumps = min1([length(indUp) length(indDown)]);
                        for ijump = 1:njumps
                            temp(indUp(ijump):indDown(ijump)) = temp(indUp(ijump):indDown(ijump)) - (10-.4096);
                        end
                        mouse1scaledtrace(region) = temp;
                    end
                    
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse1scaledtrace(region);
                    % M1OFFACC(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse1scaledtraceACC(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stimtrace(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, mouse1scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
                    end
                    if flag.plot
                        %some sanity check to see if the stimulus played
                        %and there's sensor data
                        figure(8)
                        hold on
                        offset=i*.1;
                        t=region;t=t/samprate; t=t-t(1);
                        plot(t, stimtrace(region)-stimtrace(1)+offset, 'm', t, mouse1scaledtrace(region)-mouse1scaledtrace(1)+offset, 'b')
                        gap_termination=pos+gapdelay/1000;
                        gap_onset=pos+gapdelay/1000-gapdur/1000;
%                         plot(gap_onset, 0, '^', gap_termination,0, 'v')
%                         plot(startle_onset, 0, 'bo')
%                         plot(pos, 0, 'r*')
                    %                     keyboard
                    end
                end
            end
        end
    end
end
if flag.plot
    str = cell(1);
    Hfig(length(gapdurs)+1) = figure; hold on
    for ifig = 1:length(gapdurs)
        figure(Hfig(ifig))
        t=region; t = t-t(1);t=t/samprate;
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,.5]),'linewidth',2)
        figure(Hfig(length(gapdurs)+1))
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,1]),'linewidth',2)
        str{ifig} = num2str(gapdurs(ifig));
    end
    legend(str)
    [~, FN, ~] = fileparts(pwd);
    
    title(FN,'interpreter','none')
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))

PeakON=[];
PeakOFF=[];
SumOFF=[];   %%%Here's what I added
mSumOFF=[];
bl_SumOFF=[];




mM1OFF=mean(M1OFF, 3);
mM1ON=mean(M1ON, 3);
mM1OFFACC=mean(M1OFFACC, 3);
mM1ONACC=mean(M1ONACC, 3);

mM1OFFstim=mean(M1OFFstim, 3);
mM1ONstim=mean(M1ONstim, 3);


% Accumulate startle response across trials using peak (or summed) rectified signal in region
%%%% Added throwing out oddballs 6/27/2018
bl_start=(bl_window(1)-xlimits(1))*samprate/1000;
bl_stop=bl_start+diff(bl_window)*samprate/1000;
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;

SumON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
SumONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

bl_SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

PeakON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
PeakONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

fprintf('\n')
for paindex=1:numpulseamps
    for gdindex=1:numgapdurs; % Hardcoded.
        for k=1:nrepsON(gdindex, paindex);
            %traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
            %PeakON(gdindex, paindex, k) = max(abs(traceON));
            temp = squeeze(M1ON(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakON(gdindex, paindex, k) = nan;
                SumON(gdindex, paindex, k) = nan;
            end
            temp = squeeze(M1ONACC(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                PeakONACC(gdindex, paindex, k) = nan;
                SumONACC(gdindex, paindex, k) = nan;
            end
        end
        for k=1:nrepsOFF(gdindex, paindex);
            temp = squeeze(M1OFF(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFF(gdindex, paindex, k) = nan;
                bl_SumOFF(gdindex, paindex, k) = nan;
                SumOFF(gdindex, paindex, k) = nan;
            end
            % temp = squeeze(M1OFFACC(gdindex,paindex, k, 1:10000));
            % if abs(mean(temp)) <.01 & std(temp) < .1
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % elseif flag.includeALL
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % else
            %     fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
            %     PeakOFFACC (gdindex, paindex, k) = nan;
            %     SumOFFACC (gdindex, paindex, k) = nan;
            % end
        end
        if isempty(PeakON)
            mPeakON=[];
            semPeakON=[];
            mPeakONACC=[];
            semPeakONACC=[];
        else
            mPeakON(gdindex, paindex)=median(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakON(gdindex, paindex)=std(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mPeakONACC(gdindex, paindex)=median(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakONACC(gdindex, paindex)=std(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(SumON)
            mSumON=[];
            semSumON=[];
            mSumONACC=[];
            semSumONACC=[];
        else
            mSumON(gdindex, paindex)=median(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumON(gdindex, paindex)=std(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mSumONACC(gdindex, paindex)=median(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumONACC(gdindex, paindex)=std(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(PeakOFF)
            mPeakOFF=[];
            semPeakOFF=[];
            mPeakOFFACC=[];
            semPeakOFFACC=[];
        else
            mPeakOFF(gdindex, paindex)=median(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFF(gdindex, paindex)=std(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mPeakOFFACC(gdindex, paindex)=median(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFFACC(gdindex, paindex)=std(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
        if isempty(SumOFF)
            mSumOFF=[];
            semSumOFF=[];
            mSumOFFACC=[];
            semSumOFFACC=[];
        else
            mSumOFF(gdindex, paindex)=median(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFF(gdindex, paindex)=std(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mSumOFFACC(gdindex, paindex)=median(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFFACC(gdindex, paindex)=std(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
    end
    
    %sanity check that first gapdur is 0 (i.e. control condition)
    if gapdurs(1)~=0
        error('first gapdur is not 0, what is wrong?')
    end
    
    fprintf('\nusing MEDIAN of peak(abs(trace)) or of sum(abs(trace))  responses\n')
    
    %only makes sense for numgapdurs >= 2
    if isempty(PeakON)
        percentGPIAS_ON=[];
        percentGPIAS_ONACC=[];
        percentGPIAS_ONsum=[];
        percentGPIAS_ONACCsum=[];
        pON=[];
    else
        percentGPIAS_ON(1)=nan;
        pON(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser ON  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakON(1, paindex);
                m2=mPeakON(p, paindex);
                percentGPIAS_ON(p)=((m1-m2)/m1)*100;
                A=PeakON(1,paindex, 1:nrepsON(1, paindex));
                B=PeakON(p,paindex, 1:nrepsON(p, paindex));
                [H,pON(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ON(p),H,pON(p));
                
                m1=mSumON(1, paindex);
                m2=mSumON(p, paindex);
                percentGPIAS_ONsum(p)=((m1-m2)/m1)*100;
                A=SumON(1,paindex, 1:nrepsON(1, paindex));
                B=SumON(p,paindex, 1:nrepsON(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONsum(p),H,temp);
            end
            
            m1=mPeakONACC(1, paindex);
            m2=mPeakONACC(p, paindex);
            percentGPIAS_ONACC(p)=((m1-m2)/m1)*100;
            A=PeakONACC(1,paindex, 1:nrepsON(1, paindex));
            B=PeakONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACC(p),H,temp);
            
            m1=mSumONACC(1, paindex);
            m2=mSumONACC(p, paindex);
            percentGPIAS_ONACCsum(p)=((m1-m2)/m1)*100;
            A=SumONACC(1,paindex, 1:nrepsON(1, paindex));
            B=SumONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACCsum(p),H,temp);
        end
    end
    if isempty(PeakOFF)
        percentGPIAS_OFF=[];
        percentGPIAS_OFFACC=[];
        percentGPIAS_OFFsum=[];
        percentGPIAS_OFFACCsum=[];
        pOFF=[];
    else
        percentGPIAS_OFF(1)=nan;
        pOFF(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser OFF  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakOFF(1, paindex);
                m2=mPeakOFF(p, paindex);
                percentGPIAS_OFF(p)=((m1-m2)/m1)*100;
                A=PeakOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=PeakOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,pOFF(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFF(p),H,pOFF(p));
                
                m1=mSumOFF(1, paindex);
                m2=mSumOFF(p, paindex);
                percentGPIAS_OFFsum(p)=((m1-m2)/m1)*100;
                A=SumOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=SumOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFsum(p),H,temp);
            end
            m1=mPeakOFFACC(1, paindex);
            m2=mPeakOFFACC(p, paindex);
            percentGPIAS_OFFACC(p)=((m1-m2)/m1)*100;
            A=PeakOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=PeakOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACC(p),H,temp);
            
            m1=mSumOFFACC(1, paindex);
            m2=mSumOFFACC(p, paindex);
            percentGPIAS_OFFACCsum(p)=((m1-m2)/m1)*100;
            A=SumOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=SumOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACCsum(p),H,temp);
        end
    end
    
end



%save to outfiles
[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysFolder=EphysFolder;
out.generated_by=mfilename;
out.generated_on=datestr(now);

out.IL=IL;
out.flag=flag;

out.M1ON=M1ON;    % mouse1scaledtrace (depends on flag.accel)
out.M1OFF=M1OFF;    % mouse1scaledtrace
out.mM1ON=mM1ON;    % mouse1scaledtrace
out.mM1OFF=mM1OFF;    % mouse1scaledtrace
% out.M1ONACC=M1ONACC;
% out.M1OFFACC=M1OFFACC;
% out.mM1ONACC=mM1ONACC;
% out.mM1OFFACC=mM1OFFACC;
out.M1ONstim=M1ONstim;
out.M1OFFstim=M1OFFstim;
out.mM1ONstim=mM1ONstim;
out.mM1OFFstim=mM1OFFstim;
out.PeakON=PeakON;
out.PeakOFF=PeakOFF;
out.mPeakON=mPeakON;
out.mPeakOFF=mPeakOFF;
out.semPeakON=semPeakON;
out.semPeakOFF=semPeakOFF;
out.SumON=SumON;
out.SumOFF=SumOFF;
out.bl_SumOFF=bl_SumOFF;
out.mSumON=mSumON;
out.mSumOFF=mSumOFF;
out.semSumON=semSumON;
out.semSumOFF=semSumOFF;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.samprate=samprate;

out.percentGPIAS_OFF=percentGPIAS_OFF;
out.pOFF=pOFF;
out.percentGPIAS_ON=percentGPIAS_ON;
out.pON=pON;

if IL
    out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);
else
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
end


out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;
out.isi=isi;
out.soaflag=soaflag;
out.xlimits=xlimits;
out.startle_window=startle_window;
out.samprate=samprate;
out.mouse_number=1;
try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
    out.mouseID=nb.mouseID;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
    out.mouseID='unknown';
end
outfilename=sprintf('outGPIAS_BehaviorMouse1.mat');
out.outfilename=outfilename;

save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)


close all



%%%%Here is the non-looping code for mouse 2.

fprintf('\ncomputing mouse 2...');

M1ON=[];M1OFF=[];
M1ONACC=[];M1OFFACC=[];
M1ONstim=[];M1OFFstim=[];
nrepsON=zeros(numgapdurs, numpulseamps);
nrepsOFF=zeros(numgapdurs, numpulseamps);

xlimits=[-350 350]; %xlimits for storing traces
bl_window=[-100 0];
startle_window=[0 100]; %hard coded integration region for startle response (ms)
fprintf('\nprocessing with xlimits [%d - %d]', xlimits(1), xlimits(2))
fprintf('\nprocessing with startle integration window [%d - %d]', startle_window(1), startle_window(2))

%extract the traces into a big matrix M
j=0;
if flag.plot
    for idur = 1:length(gapdurs)
        Hfig(idur) = figure;
        hold on;
        title(['raw traces for gapdur = ' num2str(gapdurs(idur)) ' ms'])
    end
end

for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') |  strcmp(Events(i).type, 'toneGPIAS')
        %note: gapdelay is the time from the soundcardtrigger (pos) to the
        %gap termination. The time to startle onset should be
        %(gapdelay + soa) after pos, if soaflag=soa
        gapdur=Events(i).gapdur;
        switch soaflag
            case 'isi'
                isi=Events(i).soa;
                soa=isi+gapdur;
            case 'soa'
                soa=Events(i).soa;
                isi=soa-gapdur;
        end
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=LaserTrials(i);
        %since this is all about quantifying startle response, we want a trace
        %locked to startle pulse (not gap)
        startle_onset=pos+gapdelay/1000+isi/1000;
        start=startle_onset + xlimits(1)/1000; %start is in seconds, should be at xlimits relative to startle onset
        stop=startle_onset + xlimits(2)/1000; %stop is in seconds
        if start>0 %(disallow negative or zero start times)
            gdindex= find(gapdur==gapdurs);
            pulseamp=Events(i).pulseamp;
            paindex= find(pulseamp==pulseamps);
            %start=round(pos+xlimits(1)*1e-3*samprate);
            %stop=round(pos+xlimits(2)*1e-3*samprate)-1;
            region=round(start*samprate)+1:round(stop*samprate);
            if isempty(find(region<1))
                if laser
                    nrepsON(gdindex,paindex)=nrepsON(gdindex,paindex)+1;
                    M1ON(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse2scaledtrace(region);
                    M1ONACC(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse2scaledtraceACC(region);
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stimtrace(region);
                else
                    if flag.unwrap
                        temp = mouse2scaledtrace(region);
                        tempD = [0; diff(temp)];
                        indUp = find(tempD>2);
                        indDown = find(tempD<-2)-1;
                        njumps = min1([length(indUp) length(indDown)]);
                        for ijump = 1:njumps
                            temp(indUp(ijump):indDown(ijump)) = temp(indUp(ijump):indDown(ijump)) - (10-.4096);
                        end
                        mouse2scaledtrace(region) = temp;
                    end
                    
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse2scaledtrace(region);
                    % M1OFFACC(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse2scaledtraceACC(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stimtrace(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, mouse2scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
                    end
                    if flag.plot
                        %some sanity check to see if the stimulus played
                        %and there's sensor data
                        figure(8)
                        hold on
                        offset=i*.1;
                        t=region;t=t/samprate; t=t-t(1);
                        plot(t, stimtrace(region)-stimtrace(1)+offset, 'm', t, mouse2scaledtrace(region)-mouse2scaledtrace(1)+offset, 'b')
                        gap_termination=pos+gapdelay/1000;
                        gap_onset=pos+gapdelay/1000-gapdur/1000;
%                         plot(gap_onset, 0, '^', gap_termination,0, 'v')
%                         plot(startle_onset, 0, 'bo')
%                         plot(pos, 0, 'r*')
                    %                     keyboard
                    end
                end
            end
        end
    end
end
if flag.plot
    str = cell(1);
    Hfig(length(gapdurs)+1) = figure; hold on
    for ifig = 1:length(gapdurs)
        figure(Hfig(ifig))
        t=region; t = t-t(1);t=t/samprate;
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,.5]),'linewidth',2)
        figure(Hfig(length(gapdurs)+1))
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,1]),'linewidth',2)
        str{ifig} = num2str(gapdurs(ifig));
    end
    legend(str)
    [~, FN, ~] = fileparts(pwd);
    
    title(FN,'interpreter','none')
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))

PeakON=[];
PeakOFF=[];
SumOFF=[];   %%%Here's what I added
mSumOFF=[];
bl_SumOFF=[];




mM1OFF=mean(M1OFF, 3);
mM1ON=mean(M1ON, 3);
mM1OFFACC=mean(M1OFFACC, 3);
mM1ONACC=mean(M1ONACC, 3);

mM1OFFstim=mean(M1OFFstim, 3);
mM1ONstim=mean(M1ONstim, 3);


% Accumulate startle response across trials using peak (or summed) rectified signal in region
%%%% Added throwing out oddballs 6/27/2018
bl_start=(bl_window(1)-xlimits(1))*samprate/1000;
bl_stop=bl_start+diff(bl_window)*samprate/1000;
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;

SumON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
SumONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

bl_SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

PeakON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
PeakONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

fprintf('\n')
for paindex=1:numpulseamps
    for gdindex=1:numgapdurs; % Hardcoded.
        for k=1:nrepsON(gdindex, paindex);
            %traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
            %PeakON(gdindex, paindex, k) = max(abs(traceON));
            temp = squeeze(M1ON(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakON(gdindex, paindex, k) = nan;
                SumON(gdindex, paindex, k) = nan;
            end
            temp = squeeze(M1ONACC(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                PeakONACC(gdindex, paindex, k) = nan;
                SumONACC(gdindex, paindex, k) = nan;
            end
        end
        for k=1:nrepsOFF(gdindex, paindex);
            temp = squeeze(M1OFF(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFF(gdindex, paindex, k) = nan;
                bl_SumOFF(gdindex, paindex, k) = nan;
                SumOFF(gdindex, paindex, k) = nan;
            end
            % temp = squeeze(M1OFFACC(gdindex,paindex, k, 1:10000));
            % if abs(mean(temp)) <.01 & std(temp) < .1
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % elseif flag.includeALL
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % else
            %     fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
            %     PeakOFFACC (gdindex, paindex, k) = nan;
            %     SumOFFACC (gdindex, paindex, k) = nan;
            % end
        end
        if isempty(PeakON)
            mPeakON=[];
            semPeakON=[];
            mPeakONACC=[];
            semPeakONACC=[];
        else
            mPeakON(gdindex, paindex)=median(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakON(gdindex, paindex)=std(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mPeakONACC(gdindex, paindex)=median(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakONACC(gdindex, paindex)=std(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(SumON)
            mSumON=[];
            semSumON=[];
            mSumONACC=[];
            semSumONACC=[];
        else
            mSumON(gdindex, paindex)=median(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumON(gdindex, paindex)=std(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mSumONACC(gdindex, paindex)=median(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumONACC(gdindex, paindex)=std(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(PeakOFF)
            mPeakOFF=[];
            semPeakOFF=[];
            mPeakOFFACC=[];
            semPeakOFFACC=[];
        else
            mPeakOFF(gdindex, paindex)=median(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFF(gdindex, paindex)=std(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mPeakOFFACC(gdindex, paindex)=median(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFFACC(gdindex, paindex)=std(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
        if isempty(SumOFF)
            mSumOFF=[];
            semSumOFF=[];
            mSumOFFACC=[];
            semSumOFFACC=[];
        else
            mSumOFF(gdindex, paindex)=median(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFF(gdindex, paindex)=std(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mSumOFFACC(gdindex, paindex)=median(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFFACC(gdindex, paindex)=std(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
    end
    
    %sanity check that first gapdur is 0 (i.e. control condition)
    if gapdurs(1)~=0
        error('first gapdur is not 0, what is wrong?')
    end
    
    fprintf('\nusing MEDIAN of peak(abs(trace)) or of sum(abs(trace))  responses\n')
    
    %only makes sense for numgapdurs >= 2
    if isempty(PeakON)
        percentGPIAS_ON=[];
        percentGPIAS_ONACC=[];
        percentGPIAS_ONsum=[];
        percentGPIAS_ONACCsum=[];
        pON=[];
    else
        percentGPIAS_ON(1)=nan;
        pON(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser ON  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakON(1, paindex);
                m2=mPeakON(p, paindex);
                percentGPIAS_ON(p)=((m1-m2)/m1)*100;
                A=PeakON(1,paindex, 1:nrepsON(1, paindex));
                B=PeakON(p,paindex, 1:nrepsON(p, paindex));
                [H,pON(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ON(p),H,pON(p));
                
                m1=mSumON(1, paindex);
                m2=mSumON(p, paindex);
                percentGPIAS_ONsum(p)=((m1-m2)/m1)*100;
                A=SumON(1,paindex, 1:nrepsON(1, paindex));
                B=SumON(p,paindex, 1:nrepsON(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONsum(p),H,temp);
            end
            
            m1=mPeakONACC(1, paindex);
            m2=mPeakONACC(p, paindex);
            percentGPIAS_ONACC(p)=((m1-m2)/m1)*100;
            A=PeakONACC(1,paindex, 1:nrepsON(1, paindex));
            B=PeakONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACC(p),H,temp);
            
            m1=mSumONACC(1, paindex);
            m2=mSumONACC(p, paindex);
            percentGPIAS_ONACCsum(p)=((m1-m2)/m1)*100;
            A=SumONACC(1,paindex, 1:nrepsON(1, paindex));
            B=SumONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACCsum(p),H,temp);
        end
    end
    if isempty(PeakOFF)
        percentGPIAS_OFF=[];
        percentGPIAS_OFFACC=[];
        percentGPIAS_OFFsum=[];
        percentGPIAS_OFFACCsum=[];
        pOFF=[];
    else
        percentGPIAS_OFF(1)=nan;
        pOFF(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser OFF  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakOFF(1, paindex);
                m2=mPeakOFF(p, paindex);
                percentGPIAS_OFF(p)=((m1-m2)/m1)*100;
                A=PeakOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=PeakOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,pOFF(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFF(p),H,pOFF(p));
                
                m1=mSumOFF(1, paindex);
                m2=mSumOFF(p, paindex);
                percentGPIAS_OFFsum(p)=((m1-m2)/m1)*100;
                A=SumOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=SumOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFsum(p),H,temp);
            end
            m1=mPeakOFFACC(1, paindex);
            m2=mPeakOFFACC(p, paindex);
            percentGPIAS_OFFACC(p)=((m1-m2)/m1)*100;
            A=PeakOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=PeakOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACC(p),H,temp);
            
            m1=mSumOFFACC(1, paindex);
            m2=mSumOFFACC(p, paindex);
            percentGPIAS_OFFACCsum(p)=((m1-m2)/m1)*100;
            A=SumOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=SumOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACCsum(p),H,temp);
        end
    end
    
end



%save to outfiles
[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysFolder=EphysFolder;
out.generated_by=mfilename;
out.generated_on=datestr(now);

out.IL=IL;
out.flag=flag;

out.M1ON=M1ON;    % mouse2scaledtrace (depends on flag.accel)
out.M1OFF=M1OFF;    % mouse2scaledtrace
out.mM1ON=mM1ON;    % mouse2scaledtrace
out.mM1OFF=mM1OFF;    % mouse2scaledtrace
% out.M1ONACC=M1ONACC;
% out.M1OFFACC=M1OFFACC;
% out.mM1ONACC=mM1ONACC;
% out.mM1OFFACC=mM1OFFACC;
out.M1ONstim=M1ONstim;
out.M1OFFstim=M1OFFstim;
out.mM1ONstim=mM1ONstim;
out.mM1OFFstim=mM1OFFstim;
out.PeakON=PeakON;
out.PeakOFF=PeakOFF;
out.mPeakON=mPeakON;
out.mPeakOFF=mPeakOFF;
out.semPeakON=semPeakON;
out.semPeakOFF=semPeakOFF;
out.SumON=SumON;
out.SumOFF=SumOFF;
out.bl_SumOFF=bl_SumOFF;
out.mSumON=mSumON;
out.mSumOFF=mSumOFF;
out.semSumON=semSumON;
out.semSumOFF=semSumOFF;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.samprate=samprate;

out.percentGPIAS_OFF=percentGPIAS_OFF;
out.pOFF=pOFF;
out.percentGPIAS_ON=percentGPIAS_ON;
out.pON=pON;


if IL
    out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);
else
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
end


out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;
out.isi=isi;
out.soaflag=soaflag;
out.xlimits=xlimits;
out.startle_window=startle_window;
out.samprate=samprate;
out.mouse_number=2;
try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
    out.mouseID=nb.mouse2ID;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
    out.mouseID='unknown';
end
outfilename=sprintf('outGPIAS_BehaviorMouse2.mat');
out.outfilename=outfilename;
save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)



close all

%%%%Here is the non-looping code for mouse 3.

fprintf('\ncomputing mouse 3 ...');

M1ON=[];M1OFF=[];
M1ONACC=[];M1OFFACC=[];
M1ONstim=[];M1OFFstim=[];
nrepsON=zeros(numgapdurs, numpulseamps);
nrepsOFF=zeros(numgapdurs, numpulseamps);

xlimits=[-350 350]; %xlimits for storing traces
bl_window=[-100 0];
startle_window=[0 100]; %hard coded integration region for startle response (ms)
fprintf('\nprocessing with xlimits [%d - %d]', xlimits(1), xlimits(2))
fprintf('\nprocessing with startle integration window [%d - %d]', startle_window(1), startle_window(2))

%extract the traces into a big matrix M
j=0;
if flag.plot
    for idur = 1:length(gapdurs)
        Hfig(idur) = figure;
        hold on;
        title(['raw traces for gapdur = ' num2str(gapdurs(idur)) ' ms'])
    end
end

for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') |  strcmp(Events(i).type, 'toneGPIAS')
        %note: gapdelay is the time from the soundcardtrigger (pos) to the
        %gap termination. The time to startle onset should be
        %(gapdelay + soa) after pos, if soaflag=soa
        gapdur=Events(i).gapdur;
        switch soaflag
            case 'isi'
                isi=Events(i).soa;
                soa=isi+gapdur;
            case 'soa'
                soa=Events(i).soa;
                isi=soa-gapdur;
        end
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=LaserTrials(i);
        %since this is all about quantifying startle response, we want a trace
        %locked to startle pulse (not gap)
        startle_onset=pos+gapdelay/1000+isi/1000;
        start=startle_onset + xlimits(1)/1000; %start is in seconds, should be at xlimits relative to startle onset
        stop=startle_onset + xlimits(2)/1000; %stop is in seconds
        if start>0 %(disallow negative or zero start times)
            gdindex= find(gapdur==gapdurs);
            pulseamp=Events(i).pulseamp;
            paindex= find(pulseamp==pulseamps);
            %start=round(pos+xlimits(1)*1e-3*samprate);
            %stop=round(pos+xlimits(2)*1e-3*samprate)-1;
            region=round(start*samprate)+1:round(stop*samprate);
            if isempty(find(region<1))
                if laser
                    nrepsON(gdindex,paindex)=nrepsON(gdindex,paindex)+1;
                    M1ON(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse3scaledtrace(region);
                    M1ONACC(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse3scaledtraceACC(region);
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stimtrace(region);
                else
                    if flag.unwrap
                        temp = mouse3scaledtrace(region);
                        tempD = [0; diff(temp)];
                        indUp = find(tempD>2);
                        indDown = find(tempD<-2)-1;
                        njumps = min1([length(indUp) length(indDown)]);
                        for ijump = 1:njumps
                            temp(indUp(ijump):indDown(ijump)) = temp(indUp(ijump):indDown(ijump)) - (10-.4096);
                        end
                        mouse3scaledtrace(region) = temp;
                    end
                    
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse3scaledtrace(region);
                    % M1OFFACC(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse3scaledtraceACC(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stimtrace(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, mouse3scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
                    end
                    if flag.plot
                        %some sanity check to see if the stimulus played
                        %and there's sensor data
                        figure(8)
                        hold on
                        offset=i*.1;
                        t=region;t=t/samprate; t=t-t(1);
                        plot(t, stimtrace(region)-stimtrace(1)+offset, 'm', t, mouse3scaledtrace(region)-mouse3scaledtrace(1)+offset, 'b')
                        gap_termination=pos+gapdelay/1000;
                        gap_onset=pos+gapdelay/1000-gapdur/1000;
%                         plot(gap_onset, 0, '^', gap_termination,0, 'v')
%                         plot(startle_onset, 0, 'bo')
%                         plot(pos, 0, 'r*')
                    %                     keyboard
                    end
                end
            end
        end
    end
end
if flag.plot
    str = cell(1);
    Hfig(length(gapdurs)+1) = figure; hold on
    for ifig = 1:length(gapdurs)
        figure(Hfig(ifig))
        t=region; t = t-t(1);t=t/samprate;
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,.5]),'linewidth',2)
        figure(Hfig(length(gapdurs)+1))
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,1]),'linewidth',2)
        str{ifig} = num2str(gapdurs(ifig));
    end
    legend(str)
    [~, FN, ~] = fileparts(pwd);
    
    title(FN,'interpreter','none')
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))

PeakON=[];
PeakOFF=[];
SumOFF=[];   %%%Here's what I added
mSumOFF=[];
bl_SumOFF=[];




mM1OFF=mean(M1OFF, 3);
mM1ON=mean(M1ON, 3);
mM1OFFACC=mean(M1OFFACC, 3);
mM1ONACC=mean(M1ONACC, 3);

mM1OFFstim=mean(M1OFFstim, 3);
mM1ONstim=mean(M1ONstim, 3);


% Accumulate startle response across trials using peak (or summed) rectified signal in region
%%%% Added throwing out oddballs 6/27/2018
bl_start=(bl_window(1)-xlimits(1))*samprate/1000;
bl_stop=bl_start+diff(bl_window)*samprate/1000;
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;

SumON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
SumONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

bl_SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

PeakON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
PeakONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

fprintf('\n')
for paindex=1:numpulseamps
    for gdindex=1:numgapdurs; % Hardcoded.
        for k=1:nrepsON(gdindex, paindex);
            %traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
            %PeakON(gdindex, paindex, k) = max(abs(traceON));
            temp = squeeze(M1ON(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakON(gdindex, paindex, k) = nan;
                SumON(gdindex, paindex, k) = nan;
            end
            temp = squeeze(M1ONACC(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                PeakONACC(gdindex, paindex, k) = nan;
                SumONACC(gdindex, paindex, k) = nan;
            end
        end
        for k=1:nrepsOFF(gdindex, paindex);
            temp = squeeze(M1OFF(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFF(gdindex, paindex, k) = nan;
                bl_SumOFF(gdindex, paindex, k) = nan;
                SumOFF(gdindex, paindex, k) = nan;
            end
            % temp = squeeze(M1OFFACC(gdindex,paindex, k, 1:10000));
            % if abs(mean(temp)) <.01 & std(temp) < .1
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % elseif flag.includeALL
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % else
            %     fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
            %     PeakOFFACC (gdindex, paindex, k) = nan;
            %     SumOFFACC (gdindex, paindex, k) = nan;
            % end
        end
        if isempty(PeakON)
            mPeakON=[];
            semPeakON=[];
            mPeakONACC=[];
            semPeakONACC=[];
        else
            mPeakON(gdindex, paindex)=median(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakON(gdindex, paindex)=std(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mPeakONACC(gdindex, paindex)=median(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakONACC(gdindex, paindex)=std(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(SumON)
            mSumON=[];
            semSumON=[];
            mSumONACC=[];
            semSumONACC=[];
        else
            mSumON(gdindex, paindex)=median(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumON(gdindex, paindex)=std(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mSumONACC(gdindex, paindex)=median(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumONACC(gdindex, paindex)=std(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(PeakOFF)
            mPeakOFF=[];
            semPeakOFF=[];
            mPeakOFFACC=[];
            semPeakOFFACC=[];
        else
            mPeakOFF(gdindex, paindex)=median(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFF(gdindex, paindex)=std(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mPeakOFFACC(gdindex, paindex)=median(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFFACC(gdindex, paindex)=std(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
        if isempty(SumOFF)
            mSumOFF=[];
            semSumOFF=[];
            mSumOFFACC=[];
            semSumOFFACC=[];
        else
            mSumOFF(gdindex, paindex)=median(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFF(gdindex, paindex)=std(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mSumOFFACC(gdindex, paindex)=median(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFFACC(gdindex, paindex)=std(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
    end
    
    %sanity check that first gapdur is 0 (i.e. control condition)
    if gapdurs(1)~=0
        error('first gapdur is not 0, what is wrong?')
    end
    
    fprintf('\nusing MEDIAN of peak(abs(trace)) or of sum(abs(trace))  responses\n')
    
    %only makes sense for numgapdurs >= 2
    if isempty(PeakON)
        percentGPIAS_ON=[];
        percentGPIAS_ONACC=[];
        percentGPIAS_ONsum=[];
        percentGPIAS_ONACCsum=[];
        pON=[];
    else
        percentGPIAS_ON(1)=nan;
        pON(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser ON  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakON(1, paindex);
                m2=mPeakON(p, paindex);
                percentGPIAS_ON(p)=((m1-m2)/m1)*100;
                A=PeakON(1,paindex, 1:nrepsON(1, paindex));
                B=PeakON(p,paindex, 1:nrepsON(p, paindex));
                [H,pON(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ON(p),H,pON(p));
                
                m1=mSumON(1, paindex);
                m2=mSumON(p, paindex);
                percentGPIAS_ONsum(p)=((m1-m2)/m1)*100;
                A=SumON(1,paindex, 1:nrepsON(1, paindex));
                B=SumON(p,paindex, 1:nrepsON(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONsum(p),H,temp);
            end
            
            m1=mPeakONACC(1, paindex);
            m2=mPeakONACC(p, paindex);
            percentGPIAS_ONACC(p)=((m1-m2)/m1)*100;
            A=PeakONACC(1,paindex, 1:nrepsON(1, paindex));
            B=PeakONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACC(p),H,temp);
            
            m1=mSumONACC(1, paindex);
            m2=mSumONACC(p, paindex);
            percentGPIAS_ONACCsum(p)=((m1-m2)/m1)*100;
            A=SumONACC(1,paindex, 1:nrepsON(1, paindex));
            B=SumONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACCsum(p),H,temp);
        end
    end
    if isempty(PeakOFF)
        percentGPIAS_OFF=[];
        percentGPIAS_OFFACC=[];
        percentGPIAS_OFFsum=[];
        percentGPIAS_OFFACCsum=[];
        pOFF=[];
    else
        percentGPIAS_OFF(1)=nan;
        pOFF(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser OFF  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakOFF(1, paindex);
                m2=mPeakOFF(p, paindex);
                percentGPIAS_OFF(p)=((m1-m2)/m1)*100;
                A=PeakOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=PeakOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,pOFF(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFF(p),H,pOFF(p));
                
                m1=mSumOFF(1, paindex);
                m2=mSumOFF(p, paindex);
                percentGPIAS_OFFsum(p)=((m1-m2)/m1)*100;
                A=SumOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=SumOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFsum(p),H,temp);
            end
            m1=mPeakOFFACC(1, paindex);
            m2=mPeakOFFACC(p, paindex);
            percentGPIAS_OFFACC(p)=((m1-m2)/m1)*100;
            A=PeakOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=PeakOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACC(p),H,temp);
            
            m1=mSumOFFACC(1, paindex);
            m2=mSumOFFACC(p, paindex);
            percentGPIAS_OFFACCsum(p)=((m1-m2)/m1)*100;
            A=SumOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=SumOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACCsum(p),H,temp);
        end
    end
    
end



%save to outfiles
[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysFolder=EphysFolder;
out.generated_by=mfilename;
out.generated_on=datestr(now);

out.IL=IL;
out.flag=flag;

out.M1ON=M1ON;    % mouse3scaledtrace (depends on flag.accel)
out.M1OFF=M1OFF;    % mouse3scaledtrace
out.mM1ON=mM1ON;    % mouse3scaledtrace
out.mM1OFF=mM1OFF;    % mouse3scaledtrace
% out.M1ONACC=M1ONACC;
% out.M1OFFACC=M1OFFACC;
% out.mM1ONACC=mM1ONACC;
% out.mM1OFFACC=mM1OFFACC;
out.M1ONstim=M1ONstim;
out.M1OFFstim=M1OFFstim;
out.mM1ONstim=mM1ONstim;
out.mM1OFFstim=mM1OFFstim;
out.PeakON=PeakON;
out.PeakOFF=PeakOFF;
out.mPeakON=mPeakON;
out.mPeakOFF=mPeakOFF;
out.semPeakON=semPeakON;
out.semPeakOFF=semPeakOFF;
out.SumON=SumON;
out.SumOFF=SumOFF;
out.bl_SumOFF=bl_SumOFF;
out.mSumON=mSumON;
out.mSumOFF=mSumOFF;
out.semSumON=semSumON;
out.semSumOFF=semSumOFF;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.samprate=samprate;

out.percentGPIAS_OFF=percentGPIAS_OFF;
out.pOFF=pOFF;
out.percentGPIAS_ON=percentGPIAS_ON;
out.pON=pON;


if IL
    out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);
else
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
end


out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;
out.isi=isi;
out.soaflag=soaflag;
out.xlimits=xlimits;
out.startle_window=startle_window;
out.samprate=samprate;
out.mouse_number=3;
try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
    out.mouseID=nb.mouse3ID;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
    out.mouseID='unknown';
end
outfilename=sprintf('outGPIAS_BehaviorMouse3.mat');
out.outfilename=outfilename;
save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)



close all

%%%%Here is the non-looping code for mouse 4.

fprintf('\ncomputing mouse 4...');

M1ON=[];M1OFF=[];
M1ONACC=[];M1OFFACC=[];
M1ONstim=[];M1OFFstim=[];
nrepsON=zeros(numgapdurs, numpulseamps);
nrepsOFF=zeros(numgapdurs, numpulseamps);

xlimits=[-350 350]; %xlimits for storing traces
bl_window=[-100 0];
startle_window=[0 100]; %hard coded integration region for startle response (ms)
fprintf('\nprocessing with xlimits [%d - %d]', xlimits(1), xlimits(2))
fprintf('\nprocessing with startle integration window [%d - %d]', startle_window(1), startle_window(2))

%extract the traces into a big matrix M
j=0;
if flag.plot
    for idur = 1:length(gapdurs)
        Hfig(idur) = figure;
        hold on;
        title(['raw traces for gapdur = ' num2str(gapdurs(idur)) ' ms'])
    end
end

for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS') |  strcmp(Events(i).type, 'toneGPIAS')
        %note: gapdelay is the time from the soundcardtrigger (pos) to the
        %gap termination. The time to startle onset should be
        %(gapdelay + soa) after pos, if soaflag=soa
        gapdur=Events(i).gapdur;
        switch soaflag
            case 'isi'
                isi=Events(i).soa;
                soa=isi+gapdur;
            case 'soa'
                soa=Events(i).soa;
                isi=soa-gapdur;
        end
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=LaserTrials(i);
        %since this is all about quantifying startle response, we want a trace
        %locked to startle pulse (not gap)
        startle_onset=pos+gapdelay/1000+isi/1000;
        start=startle_onset + xlimits(1)/1000; %start is in seconds, should be at xlimits relative to startle onset
        stop=startle_onset + xlimits(2)/1000; %stop is in seconds
        if start>0 %(disallow negative or zero start times)
            gdindex= find(gapdur==gapdurs);
            pulseamp=Events(i).pulseamp;
            paindex= find(pulseamp==pulseamps);
            %start=round(pos+xlimits(1)*1e-3*samprate);
            %stop=round(pos+xlimits(2)*1e-3*samprate)-1;
            region=round(start*samprate)+1:round(stop*samprate);
            if isempty(find(region<1))
                if laser
                    nrepsON(gdindex,paindex)=nrepsON(gdindex,paindex)+1;
                    M1ON(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse4scaledtrace(region);
                    M1ONACC(gdindex,paindex, nrepsON(gdindex,paindex),:)=mouse4scaledtraceACC(region);
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stimtrace(region);
                else
                    if flag.unwrap
                        temp = mouse4scaledtrace(region);
                        tempD = [0; diff(temp)];
                        indUp = find(tempD>2);
                        indDown = find(tempD<-2)-1;
                        njumps = min1([length(indUp) length(indDown)]);
                        for ijump = 1:njumps
                            temp(indUp(ijump):indDown(ijump)) = temp(indUp(ijump):indDown(ijump)) - (10-.4096);
                        end
                        mouse4scaledtrace(region) = temp;
                    end
                    
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse4scaledtrace(region);
                    % M1OFFACC(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=mouse4scaledtraceACC(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stimtrace(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, mouse4scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
                    end
                    if flag.plot
                        %some sanity check to see if the stimulus played
                        %and there's sensor data
                        figure(8)
                        hold on
                        offset=i*.1;
                        t=region;t=t/samprate; t=t-t(1);
                        plot(t, stimtrace(region)-stimtrace(1)+offset, 'm', t, mouse4scaledtrace(region)-mouse4scaledtrace(1)+offset, 'b')
                        gap_termination=pos+gapdelay/1000;
                        gap_onset=pos+gapdelay/1000-gapdur/1000;
%                         plot(gap_onset, 0, '^', gap_termination,0, 'v')
%                         plot(startle_onset, 0, 'bo')
%                         plot(pos, 0, 'r*')
                    %                     keyboard
                    end
                end
            end
        end
    end
end
if flag.plot
    str = cell(1);
    Hfig(length(gapdurs)+1) = figure; hold on
    for ifig = 1:length(gapdurs)
        figure(Hfig(ifig))
        t=region; t = t-t(1);t=t/samprate;
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,.5]),'linewidth',2)
        figure(Hfig(length(gapdurs)+1))
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(gapdurs),1,1]),'linewidth',2)
        str{ifig} = num2str(gapdurs(ifig));
    end
    legend(str)
    [~, FN, ~] = fileparts(pwd);
    
    title(FN,'interpreter','none')
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))

PeakON=[];
PeakOFF=[];
SumOFF=[];   %%%Here's what I added
mSumOFF=[];
bl_SumOFF=[];




mM1OFF=mean(M1OFF, 3);
mM1ON=mean(M1ON, 3);
mM1OFFACC=mean(M1OFFACC, 3);
mM1ONACC=mean(M1ONACC, 3);

mM1OFFstim=mean(M1OFFstim, 3);
mM1ONstim=mean(M1ONstim, 3);


% Accumulate startle response across trials using peak (or summed) rectified signal in region
%%%% Added throwing out oddballs 6/27/2018
bl_start=(bl_window(1)-xlimits(1))*samprate/1000;
bl_stop=bl_start+diff(bl_window)*samprate/1000;
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;

SumON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
SumONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

bl_SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

PeakON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
PeakONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

fprintf('\n')
for paindex=1:numpulseamps
    for gdindex=1:numgapdurs; % Hardcoded.
        for k=1:nrepsON(gdindex, paindex);
            %traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
            %PeakON(gdindex, paindex, k) = max(abs(traceON));
            temp = squeeze(M1ON(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                PeakON(gdindex, paindex, k) = max(abs(traceON));
                SumON(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakON(gdindex, paindex, k) = nan;
                SumON(gdindex, paindex, k) = nan;
            end
            temp = squeeze(M1ONACC(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ONACC(gdindex,paindex, k, start:stop));
                PeakONACC(gdindex, paindex, k) = max(abs(traceON));
                SumONACC(gdindex, paindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                PeakONACC(gdindex, paindex, k) = nan;
                SumONACC(gdindex, paindex, k) = nan;
            end
        end
        for k=1:nrepsOFF(gdindex, paindex);
            temp = squeeze(M1OFF(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                bl_traceOFF=squeeze(M1OFF(gdindex,paindex, k, bl_start:bl_stop));
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                bl_SumOFF(gdindex, paindex, k) = sum(abs(bl_traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFF(gdindex, paindex, k) = nan;
                bl_SumOFF(gdindex, paindex, k) = nan;
                SumOFF(gdindex, paindex, k) = nan;
            end
            % temp = squeeze(M1OFFACC(gdindex,paindex, k, 1:10000));
            % if abs(mean(temp)) <.01 & std(temp) < .1
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % elseif flag.includeALL
            %     traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
            %     PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
            %     SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            % else
            %     fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
            %     PeakOFFACC (gdindex, paindex, k) = nan;
            %     SumOFFACC (gdindex, paindex, k) = nan;
            % end
        end
        if isempty(PeakON)
            mPeakON=[];
            semPeakON=[];
            mPeakONACC=[];
            semPeakONACC=[];
        else
            mPeakON(gdindex, paindex)=median(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakON(gdindex, paindex)=std(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mPeakONACC(gdindex, paindex)=median(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semPeakONACC(gdindex, paindex)=std(PeakONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(SumON)
            mSumON=[];
            semSumON=[];
            mSumONACC=[];
            semSumONACC=[];
        else
            mSumON(gdindex, paindex)=median(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumON(gdindex, paindex)=std(SumON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
            mSumONACC(gdindex, paindex)=median(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3,'omitnan');
            semSumONACC(gdindex, paindex)=std(SumONACC(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(PeakOFF)
            mPeakOFF=[];
            semPeakOFF=[];
            mPeakOFFACC=[];
            semPeakOFFACC=[];
        else
            mPeakOFF(gdindex, paindex)=median(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFF(gdindex, paindex)=std(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mPeakOFFACC(gdindex, paindex)=median(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semPeakOFFACC(gdindex, paindex)=std(PeakOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
        if isempty(SumOFF)
            mSumOFF=[];
            semSumOFF=[];
            mSumOFFACC=[];
            semSumOFFACC=[];
        else
            mSumOFF(gdindex, paindex)=median(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFF(gdindex, paindex)=std(SumOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
            mSumOFFACC(gdindex, paindex)=median(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3,'omitnan');
            semSumOFFACC(gdindex, paindex)=std(SumOFFACC(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3,'omitnan')/sqrt(nrepsOFF(gdindex, paindex));
        end
    end
    
    %sanity check that first gapdur is 0 (i.e. control condition)
    if gapdurs(1)~=0
        error('first gapdur is not 0, what is wrong?')
    end
    
    fprintf('\nusing MEDIAN of peak(abs(trace)) or of sum(abs(trace))  responses\n')
    
    %only makes sense for numgapdurs >= 2
    if isempty(PeakON)
        percentGPIAS_ON=[];
        percentGPIAS_ONACC=[];
        percentGPIAS_ONsum=[];
        percentGPIAS_ONACCsum=[];
        pON=[];
    else
        percentGPIAS_ON(1)=nan;
        pON(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser ON  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakON(1, paindex);
                m2=mPeakON(p, paindex);
                percentGPIAS_ON(p)=((m1-m2)/m1)*100;
                A=PeakON(1,paindex, 1:nrepsON(1, paindex));
                B=PeakON(p,paindex, 1:nrepsON(p, paindex));
                [H,pON(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ON(p),H,pON(p));
                
                m1=mSumON(1, paindex);
                m2=mSumON(p, paindex);
                percentGPIAS_ONsum(p)=((m1-m2)/m1)*100;
                A=SumON(1,paindex, 1:nrepsON(1, paindex));
                B=SumON(p,paindex, 1:nrepsON(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONsum(p),H,temp);
            end
            
            m1=mPeakONACC(1, paindex);
            m2=mPeakONACC(p, paindex);
            percentGPIAS_ONACC(p)=((m1-m2)/m1)*100;
            A=PeakONACC(1,paindex, 1:nrepsON(1, paindex));
            B=PeakONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACC(p),H,temp);
            
            m1=mSumONACC(1, paindex);
            m2=mSumONACC(p, paindex);
            percentGPIAS_ONACCsum(p)=((m1-m2)/m1)*100;
            A=SumONACC(1,paindex, 1:nrepsON(1, paindex));
            B=SumONACC(p,paindex, 1:nrepsON(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ONACCsum(p),H,temp);
        end
    end
    if isempty(PeakOFF)
        percentGPIAS_OFF=[];
        percentGPIAS_OFFACC=[];
        percentGPIAS_OFFsum=[];
        percentGPIAS_OFFACCsum=[];
        pOFF=[];
    else
        percentGPIAS_OFF(1)=nan;
        pOFF(1)=nan;
        for p=2:numgapdurs;
            fprintf('Laser OFF  pa:%ddB, \n', pulseamps(paindex));
            if flag.accel ==4
                m1=mPeakOFF(1, paindex);
                m2=mPeakOFF(p, paindex);
                percentGPIAS_OFF(p)=((m1-m2)/m1)*100;
                A=PeakOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=PeakOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,pOFF(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFF(p),H,pOFF(p));
                
                m1=mSumOFF(1, paindex);
                m2=mSumOFF(p, paindex);
                percentGPIAS_OFFsum(p)=((m1-m2)/m1)*100;
                A=SumOFF(1,paindex, 1:nrepsOFF(1, paindex));
                B=SumOFF(p,paindex, 1:nrepsOFF(p, paindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFsum(p),H,temp);
            end
            m1=mPeakOFFACC(1, paindex);
            m2=mPeakOFFACC(p, paindex);
            percentGPIAS_OFFACC(p)=((m1-m2)/m1)*100;
            A=PeakOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=PeakOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACC(p),H,temp);
            
            m1=mSumOFFACC(1, paindex);
            m2=mSumOFFACC(p, paindex);
            percentGPIAS_OFFACCsum(p)=((m1-m2)/m1)*100;
            A=SumOFFACC(1,paindex, 1:nrepsOFF(1, paindex));
            B=SumOFFACC(p,paindex, 1:nrepsOFF(p, paindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_OFFACCsum(p),H,temp);
        end
    end
    
end



%save to outfiles
[~,BonsaiFolder,~]=fileparts(BonsaiPath);
out.BonsaiFolder=BonsaiFolder;
out.BonsaiPath=BonsaiPath;
out.EphysFolder=EphysFolder;
out.generated_by=mfilename;
out.generated_on=datestr(now);

out.IL=IL;
out.flag=flag;

out.M1ON=M1ON;    % mouse4scaledtrace (depends on flag.accel)
out.M1OFF=M1OFF;    % mouse4scaledtrace
out.mM1ON=mM1ON;    % mouse4scaledtrace
out.mM1OFF=mM1OFF;    % mouse4scaledtrace
% out.M1ONACC=M1ONACC;
% out.M1OFFACC=M1OFFACC;
% out.mM1ONACC=mM1ONACC;
% out.mM1OFFACC=mM1OFFACC;
out.M1ONstim=M1ONstim;
out.M1OFFstim=M1OFFstim;
out.mM1ONstim=mM1ONstim;
out.mM1OFFstim=mM1OFFstim;
out.PeakON=PeakON;
out.PeakOFF=PeakOFF;
out.mPeakON=mPeakON;
out.mPeakOFF=mPeakOFF;
out.semPeakON=semPeakON;
out.semPeakOFF=semPeakOFF;
out.SumON=SumON;
out.SumOFF=SumOFF;
out.bl_SumOFF=bl_SumOFF;
out.mSumON=mSumON;
out.mSumOFF=mSumOFF;
out.semSumON=semSumON;
out.semSumOFF=semSumOFF;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.samprate=samprate;

out.percentGPIAS_OFF=percentGPIAS_OFF;
out.pOFF=pOFF;
out.percentGPIAS_ON=percentGPIAS_ON;
out.pON=pON;


if IL
    out.LaserStart=unique(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);
else
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
end


out.numpulseamps = numpulseamps;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.gapdurs = gapdurs;
out.gapdelay = gapdelay;
out.soa=soa;
out.isi=isi;
out.soaflag=soaflag;
out.xlimits=xlimits;
out.startle_window=startle_window;
out.samprate=samprate;
out.mouse_number=4;
try
    out.nb=nb;
    out.stimlog=stimlog;
    out.user=nb.user;
    out.mouseID=nb.mouse4ID;
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
    out.mouseID='unknown';
end
outfilename=sprintf('outGPIAS_BehaviorMouse4.mat');
out.outfilename=outfilename;
save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)



close all




