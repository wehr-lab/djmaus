function ProcessClicktrain_PSTH2(varargin)

%processes clustered spiking tuning curve data from djmaus
%using new OpenEphys and kilosort file formats and hierarchy
%creates only a single outfile per session
%(automatically processes all cells)
%-mike 07.19.2024
%
% usage: ProcessClicktrain_PSTH2(SortedUnits, BonsaiPath, EphysPath, EphysPath_KS, [xlimits],[ylimits])
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
    case 'PlotClicktrain_PSTH'
        %that's fine
    otherwise
        error('This stimulus protcol does not appear to be clicktrain protocol.')
end

%read  Kilosort data
fprintf('\ncomputing tuning curve...');
samprate=Sky.OEsamplerate;

j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'clicktrain')
        j=j+1;
        allicis(i)=Events(i).ici;
        alldurs(i)=Events(i).duration;
        %       if isfield(Events(i), 'nclicks')
        allnclicks(i)=Events(i).nclicks;
        allamps(i)=Events(i).amplitude;
        allclickdurations(i)=Events(i).clickduration;
    end
end
%%%%%%%%%%%%%%
icis=unique(allicis);
nclicks=unique(allnclicks);
nclicks=sort(nclicks, 'descend');
durs=unique(alldurs);
amps=unique(allamps);
numicis=length(icis);
numnclicks=length(nclicks);
numamps=length(amps);
numdurs=length(durs);
nrepsON=zeros(numicis, 1);
nrepsOFF=zeros(numicis, 1);

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
        warning('ProcessClicktrain_PSTH2: Cannot tell if laser button was turned on in djmaus GUI');
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
    warning('ProcessClicktrain_PSTH2: still need to figure out how to find/plot Laser monitor channel ')
else
    LaserRecorded=1;
end
if exist('stimtrace')==1
    %we loaded stimtrace from open ephys Session object
        StimRecorded=1;
elseif isempty(getStimfile('.'))
    StimRecorded=0;
    warning('Stimulus monitor channel not recorded')
    warning('ProcessClicktrain_PSTH2: still need to figure out how to find/plot Stimulus monitor channel ')
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


% Mt: matrix with each complete train, size: Mt(numicis, nreps).spiketimes
% Ms: stimulus matrix in same format as Mt
% Mc: matrix sorted into each click response, size: Mt(numicis, nclicks, nreps).spiketimes

%first sort trains into matrix Mt

MtON=[];MtOFF=[];MtONStim=[];MtOFFStim=[];

numcells=length(SortedUnits);
M1=[];M1ON=[];M1OFF=[];
mM1ON=[];
mM1OFF=[];
nreps=zeros(numcells, numicis);
nrepsON=zeros(numcells, numicis);
nrepsOFF=zeros(numcells, numicis);
nreps_ssON=zeros(numcells);
nreps_ssOFF=zeros(numcells);
%extract the traces into a big matrix M
j=0;
inRange=0;



%extract the traces into a big matrix M
for cellnum=1:numcells;
    spiketimes=[];
    fprintf('\nreading SortedUnits cell %d/%d', cellnum, numcells)
    spiketimes= SortedUnits(cellnum).spiketimes;

    totalnumspikes=length(spiketimes);
    fprintf('\nsuccessfully loaded spike data with %d spikes\n',totalnumspikes)

    for i=1:length(Events)
        if strcmp(Events(i).type, 'clicktrain')

            pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
            laser=LaserTrials(i);
            start=pos + xlimits(1)/1000; %start is in seconds
            stop= pos + xlimits(2)/1000; %stop is in seconds
            region=round(start*samprate)+1:round(stop*samprate);
            if start>0 %(disallow negative or zero start times)
                ici=Events(i).ici;
                iciindex= find(ici==icis);
                st=spiketimes; %are in seconds
                spiketimes1=st(st>start & st<stop); % spiketimes in region
                spikecount=length(spiketimes1); % No. of spikes fired in response to this rep of this stim.
                inRange=inRange+ spikecount; %accumulate total spikecount in region
                spiketimes1=(spiketimes1-pos)*1000;%convert to ms after tone onset
                spont_spikecount=length(find(st<start & st>(start-(stop-start)))); % No. spikes in a region of same length preceding response window


                if laser
                    nrepsON(cellnum, iciindex)=nrepsON(cellnum, iciindex)+1;
                    MtON(cellnum, iciindex, nrepsON(cellnum, iciindex)).spiketimes=spiketimes1; % Spike times
                    MtONspikecount(cellnum, iciindex,nrepsON(cellnum, iciindex))=spikecount; % No. of spikes
                    MtspontON(cellnum, iciindex, nrepsON(cellnum, iciindex))=spont_spikecount; % No. of spikes in spont window, for each presentation.
                    if LaserRecorded
                        MtONLaser(iciindex, nrepsON(cellnum, iciindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        MtONStim(iciindex, nrepsON(cellnum, iciindex),:)=Stimtrace(region);
                    end
                else
                    nrepsOFF(cellnum, iciindex)=nrepsOFF(cellnum, iciindex)+1;

                    MtOFF(cellnum, iciindex, nrepsOFF(cellnum, iciindex)).spiketimes=spiketimes1;
                    MtOFFspikecount(cellnum, iciindex,nrepsOFF(cellnum, iciindex))=spikecount;
                    MtspontOFF(cellnum, iciindex, nrepsOFF(cellnum, iciindex))=spont_spikecount;
                    if LaserRecorded
                        MtOFFLaser(iciindex, nrepsOFF(cellnum, iciindex),:)=Lasertrace(region);
                    end
                    if StimRecorded
                        MtOFFStim(iciindex, nrepsOFF(cellnum, iciindex),:)=Stimtrace(region);
                    end

                    % %%% Individual stimulus stats in chronological order -Nick 8/23/18 %%%
                    % onset=ici*(0:nclicks(iciindex)-1);
                    % phase=[];
                    % for s= MtOFF(cellnum, iciindex, nrepsOFF(iciindex)).spiketimes
                    %     if s>0 & s<onset(end)+ici
                    %         p=s-onset;
                    %         q=p(p>0);
                    %         u=q(end);
                    %         phase=[phase 2*pi*u/ici];
                    %     end
                    % end
                    % n=length(phase);
                    % r=sqrt(sum(cos(phase)).^2 + sum(sin(phase)).^2)/n;
                    % %note: if there is only one spike, Vs=1 artifactually
                    % %therefore I (Mike) am excluding cases where there is only one spike
                    % if n==1 r=nan; end
                    % stimstats(i).CT = iciindex; %ici indices in chronological order
                    % stimstats(i).VS = r;  %Vector strength
                    % stimstats(i).RZ = n*(r^2);  %Rayleigh Z statistic RZ=n*(VS)^2, Zar p.616
                    % stimstats(i).PV = exp(-n*(r^2)); %p-value of Rayleigh Z statistic
                    % stimstats(i).SC = n; %spikecount (between stimulus onset and offset)
                    % %%% End of Nick's addition %%%

                end
            end
        end
    end

    fprintf('\ncell %d, total num spikes: %d', cellnum, length(spiketimes))
    fprintf('\nIn range: %d', inRange)


    %accumulate across trials
    for iciindex=[1:numicis]
        spiketimesON=[];
        for rep=1:nrepsON(cellnum, iciindex)
            spiketimesON=[spiketimesON MtON(cellnum, iciindex, rep).spiketimes];
        end
        mMtON(cellnum, iciindex).spiketimes=spiketimesON;
        spiketimesOFF=[];
        for rep=1:nrepsOFF(cellnum, iciindex)
            spiketimesOFF=[spiketimesOFF MtOFF(cellnum, iciindex, rep).spiketimes];
        end
        mMtOFF(cellnum, iciindex).spiketimes=spiketimesOFF;
    end
    for iciindex=[1:numicis]
        if IL
            %we don't need to save a copy for every cell, there's only one trace
        
            mMtONStim(iciindex,:)=mean(MtONStim(iciindex, 1:nrepsON(cellnum, iciindex),:), 2);
        end
        mMtOFFStim(iciindex,:)=mean(MtOFFStim(iciindex, 1:nrepsOFF(cellnum, iciindex),:), 2);
    end

    % %spikecount matrices
    if IL
        mMtONspikecount=mean(MtONspikecount, 3);
        sMtONspikecount=std(MtONspikecount, 0, 3);
        semMtONspikecount=sMtONspikecount./nrepsON(cellnum, iciindex);
        mMtspontON=mean(MtspontON, 3);
        sMtspontON=std(MtspontON, 0,3);
        semMtspontON=sMtspontON./nrepsON(cellnum, iciindex);
    else
        MtONspikecount=[];
        mMtONspikecount=[];
        sMtONspikecount=[];
        semMtONspikecount=[];
        mMtspontON=[];
        sMtspontON=[];
        semMtspontON=[];
    end
    mMtOFFspikecount=mean(MtOFFspikecount, 3);
    sMtOFFspikecount=std(MtOFFspikecount, 0, 3);
    semMtOFFspikecount=sMtOFFspikecount./nrepsOFF(cellnum, iciindex);
    mMtspontOFF=mean(MtspontOFF, 3);
    sMtspontOFF=std(MtspontOFF, 0,3);
    semMtspontOFF=sMtspontOFF./nrepsOFF(cellnum, iciindex);

    % %WHAT SHOULD TRACELENGTH BE????
     tracelength=50;%25; %ms
     fprintf('\nusing response window of 0-%d ms after tone onset', tracelength);
    %
    fprintf( '\nsorting into click matrix...');
    for iciindex=1:numicis
        for k=1:nclicks(iciindex)
            for rep=1:nrepsON(cellnum, iciindex)
                ici=icis(iciindex);
                start=(0+(k-1)*(ici));
                if start<1 start=1;end
                stop=start+tracelength ;
                region=round(start):round(stop);
                spiketimes=MtON(cellnum, iciindex, rep).spiketimes;
                st=spiketimes(spiketimes>start & spiketimes<stop); % spiketimes in region
                clickstimtrace=squeeze(MtONStim(iciindex, rep, region));
                McON(cellnum, iciindex,  k, rep).spiketimes=st;
                 McONStim(iciindex,  k,rep, :)=clickstimtrace;
            end
            for rep=1:nrepsOFF(cellnum, iciindex)
                ici=icis(iciindex);
                start=(0+(k-1)*(ici));
                if start<1 start=1;end
                stop=start+tracelength ;
                region=round(start):round(stop);
                spiketimes=MtOFF(cellnum, iciindex, rep).spiketimes;
                st=spiketimes(spiketimes>start & spiketimes<stop); % spiketimes in region
                clickstimtrace=squeeze(MtOFFStim(iciindex, rep, region));
                McOFF(cellnum, iciindex,  k, rep).spiketimes=st;
                 McOFFStim(iciindex,  k, rep, :)=clickstimtrace;
            end

            if IL
                mMcON(cellnum, iciindex, k).spiketimes = [McON(cellnum, iciindex, k, 1:nrepsON(cellnum, iciindex)).spiketimes];
            end
            mMcOFF(cellnum, iciindex, k).spiketimes=[McOFF(cellnum, iciindex, k, 1:nrepsOFF(cellnum, iciindex)).spiketimes];
        end
    end

    % Accumulate spiketimes across trials, for psth...
    % Removing unnecessary calls and de-duplicating iteration for speed
    % JLS 080417
    for iciindex=1:numicis
        for k=1:nclicks(iciindex)
            if IL
                mMcON(cellnum, iciindex, k).spiketimes=[McON(cellnum, iciindex, k, 1:nrepsON(cellnum, iciindex)).spiketimes];
                mMcONStim(iciindex, k,:)=mean(McONStim(iciindex, k, 1:nrepsON(cellnum, iciindex), :), 3);
            end
            mMcOFF(cellnum, iciindex, k).spiketimes=[McOFF(cellnum, iciindex, k, 1:nrepsOFF(cellnum, iciindex)).spiketimes];
            mMcOFFStim(iciindex, k,:)=mean(McOFFStim(iciindex, k, 1:nrepsOFF(cellnum, iciindex), :), 3);
        end
    end
    fprintf('done')

    % % compute RRTF as ratio of last5/first click response -- rep-by-rep
    for iciindex=[1:numicis]
        for rep=1:nrepsON(iciindex)
            spiketimes1=(McON(cellnum, iciindex, 1, rep).spiketimes);
            spiketimesn=[];
            for click=nclicks(iciindex)-4:nclicks(iciindex)
                spiketimesn=[spiketimesn (McON(cellnum, iciindex, click, rep).spiketimes)];
            end
            RRTF_ON(cellnum, iciindex, rep)=(length(spiketimesn)/5)/(length(spiketimes1));
            if isinf(RRTF_ON(cellnum, iciindex, rep)), RRTF_ON(cellnum, iciindex, rep)=nan;end
        end
        for rep=1:nrepsOFF(cellnum, iciindex)
            spiketimes1=(McOFF(cellnum, iciindex, 1, rep).spiketimes);
            spiketimesn=[];
            for click=nclicks(iciindex)-4:nclicks(iciindex)
                spiketimesn=[spiketimesn (McOFF(cellnum, iciindex, click, rep).spiketimes)];
            end
            RRTF_OFF(cellnum, iciindex, rep)=(length(spiketimesn)/5)/(length(spiketimes1));
            if isinf(RRTF_OFF(cellnum, iciindex, rep)), RRTF_OFF(cellnum, iciindex, rep)=nan;end
        end
    end

    %compute RRTF as ratio of last5/first click response -- mean across reps
    for iciindex=[1:numicis]
        if IL
            spiketimes1=cat(2,McON(cellnum, iciindex, 1, 1:nrepsON(cellnum, iciindex)).spiketimes);
            spiketimesn=[];
            for click=nclicks(iciindex)-4:nclicks(iciindex)
                spiketimesn=[spiketimesn (McON(cellnum, iciindex, click, 1:nrepsON(cellnum, iciindex)).spiketimes)];
            end
            mRRTF_ON(cellnum, iciindex)=(length(spiketimesn)/5)/(length(spiketimes1));
        else
            mRRTF_ON=[];
        end
        spiketimes1=cat(2,McOFF(cellnum, iciindex, 1, 1:nrepsOFF(cellnum, iciindex)).spiketimes);
        spiketimesn=[];
        for click=nclicks(iciindex)-4:nclicks(iciindex)
            spiketimesn=[spiketimesn (McOFF(cellnum, iciindex, click, 1:nrepsOFF(cellnum, iciindex)).spiketimes)];
        end
        mRRTF_OFF(cellnum, iciindex)=(length(spiketimesn)/5)/(length(spiketimes1));
    end

    %compute P2P1 as ratio of second/first click response -- rep-by-rep
    for iciindex=[1:iciindex]
        for rep=1:nrepsON(cellnum, iciindex)
            spiketimes1=(McON(cellnum, iciindex, 1, rep).spiketimes);
            spiketimes2=(McON(cellnum, iciindex, 2, rep).spiketimes);
            P2P1_ON(cellnum, iciindex,rep)=length(spiketimes2)/(length(spiketimes1));
            if isinf(P2P1_ON(cellnum, iciindex, rep)), P2P1_ON(cellnum, iciindex, rep)=nan;end
        end
        for rep=1:nrepsOFF(cellnum, iciindex)
            spiketimes1=(McOFF(cellnum, iciindex, 1, rep).spiketimes);
            spiketimes2=(McOFF(cellnum, iciindex, 2, rep).spiketimes);
            P2P1_OFF(cellnum, iciindex,rep)=length(spiketimes2)/(length(spiketimes1));
            if isinf(P2P1_OFF(cellnum, iciindex, rep)), P2P1_OFF(cellnum, iciindex, rep)=nan;end
        end
    end

    %compute P2P1 as ratio of second/first click response
    % % mean across reps
    for iciindex=[1:iciindex]
        if IL
            spiketimes1=cat(2, McON(cellnum, iciindex, 1, 1:nrepsON(cellnum, iciindex)).spiketimes);
            spiketimes2=cat(2, McON(cellnum, iciindex, 2, 1:nrepsON(cellnum, iciindex)).spiketimes);
            mP2P1_ON(cellnum, iciindex)=length(spiketimes2)/(length(spiketimes1));
        else
            mP2P1_ON=[];
        end
        spiketimes1=cat(2, McOFF(cellnum, iciindex, 1, 1:nrepsON(cellnum, iciindex)).spiketimes);
        spiketimes2=cat(2, McOFF(cellnum, iciindex, 2, 1:nrepsON(cellnum, iciindex)).spiketimes);
        mP2P1_OFF(cellnum, iciindex)=length(spiketimes2)/(length(spiketimes1));
    end

    % compute vector strength, rayleigh statistic, etc
    for iciindex=[1:numicis]
        ici=icis(iciindex);
        onsets=ici*(0:nclicks(iciindex)-1);
        spiketimesON=mMtON(cellnum, iciindex).spiketimes;
        spiketimesOFF=mMtOFF(cellnum, iciindex).spiketimes;
        phaseON=[];
        phaseOFF=[];
        for s=spiketimesON
            if s>0 & s<onsets(end)+ici
                p=s-onsets;
                q=p(p>0);
                u=q(end);
                phaseON=[phaseON 2*pi*u/ici];
            end
        end
        n=length(phaseON);
        r=sqrt(sum(cos(phaseON)).^2 + sum(sin(phaseON)).^2)/n;
        %note: if there is only one spike, Vs=1 artifactually
        % therefore I am excluding cases where there is only one spike
        if n==1 r=nan; end
        PhaseON(cellnum, iciindex).phase=phaseON;
        VsON(cellnum, iciindex)=r;
        RZON(cellnum, iciindex)=n*(r^2);  %Rayleigh Z statistic RZ=n*(VS)^2, Zar p.616
        p_ON(cellnum, iciindex)=exp(-n*(r^2)); %p-value of Rayleigh Z statistic

        for s=spiketimesOFF
            if s>0 & s<onsets(end)+ici
                p=s-onsets;
                q=p(p>0);
                u=q(end);
                phaseOFF=[phaseOFF 2*pi*u/ici];
            end
        end
        n=length(phaseOFF);
        r=sqrt(sum(cos(phaseOFF)).^2 + sum(sin(phaseOFF)).^2)/n;
        %note: if there is only one spike, Vs=1 artifactually
        % therefore I am excluding cases where there is only one spike
        if n==1 r=nan; end
        PhaseOFF(cellnum, iciindex).phase=phaseOFF;
        VsOFF(cellnum, iciindex)=r;
        RZOFF(cellnum, iciindex)=n*(r^2);  %Rayleigh Z statistic RZ=n*(VS)^2, Zar p.616
        p_OFF(cellnum, iciindex)=exp(-n*(r^2)); %p-value of Rayleigh Z statistic
    end
    %
    % %for some reason the p-value of .001 is not matching up with RZ=13.8
    % %I'm off by a factor of 1000, or p.^2
    %
    % %average laser and stimulus monitor M matrices across trials
    if LaserRecorded
        for iciindex=[1:iciindex]
            if nrepsON(cellnum, iciindex)>0
                mMtONLaser(iciindex, :)=mean(MtONLaser(iciindex, 1:nrepsON(cellnum, iciindex),:), 2);
            else %no reps for this stim, since rep=0
                mMtONLaser(iciindex,:)=zeros(size(region));
            end
            if nrepsOFF(cellnum, iciindex)>0
                mMtOFFLaser(iciindex,:)=mean(MtOFFLaser(iciindex, 1:nrepsOFF(cellnum, iciindex),:), 2);
            else %no reps for this stim, since rep=0
                mMtOFFLaser(iciindex,:)=zeros(size(region));
            end
        end
    end
    if StimRecorded
        for iciindex=[1:iciindex]
            if nrepsON(cellnum, iciindex)>0
                mMtONStim(iciindex,:)=mean(MtONStim(iciindex, 1:nrepsON(cellnum, iciindex),:), 2);
            else %no reps for this stim, since rep=0
                mMtONStim(iciindex,:)=zeros(size(region));
            end
            if nrepsOFF(cellnum, iciindex)>0
                mMtOFFStim(iciindex,:)=mean(MtOFFStim(iciindex, 1:nrepsOFF(cellnum, iciindex),:), 2);
            else %no reps for this stim, since rep=0
                mMtOFFStim(iciindex,:)=zeros(size(region));
            end
        end
    end
end % for cellnum=1:length(SortedUnits);

%save to outfiles

%sizes:
% MtON(numicis, nrepsON).spiketimes
% mMtON(numicis).spiketimes
% mMtONspikecount(numicis)

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
out.icis=icis;
out.numicis=numicis;
out.nclicks=nclicks;
out.numnclicks=numnclicks;
out.durs=durs;
out.numdurs=numdurs;
if IL
    out.MtON=MtON;
    out.mMtON=mMtON;
    out.mMtONspikecount=mMtONspikecount;
    out.sMtONspikecount=sMtONspikecount;
    out.semMtONspikecount=semMtONspikecount;
    out.mMtspontON=mMtspontON;
    out.sMtspontON=sMtspontON;
    out.semMtspontON=semMtspontON;
    out.M_LaserStart=M_LaserStart;
    out.M_LaserWidth=M_LaserWidth;
    out.M_LaserNumPulses=M_LaserNumPulses;
    out.M_LaserISI=M_LaserISI;
    out.LaserStart=uniquelo(LaserStart); %only saving one value for now, assuming it's constant
    out.LaserWidth=unique(LaserWidth);
    out.LaserNumPulses=unique(LaserNumPulses);
    out.P2P1_ON=P2P1_ON;
    out.mP2P1_ON=mP2P1_ON;
    out.RRTF_ON=RRTF_ON;
    out.mRRTF_ON=mRRTF_ON;

else
    out.MtON=[];
    out.mMtONspikecount=[];
    out.sMtONspikecount=[];
    out.semMtONspikecount=[];
    out.mMtON=[];
    out.mMtspontON=[];
    out.sMtspontON=[];
    out.semMtspontON=[];
    out.LaserStart=[];
    out.LaserWidth=[];
    out.Lasernumpulses=[];
    out.M_LaserStart=[];
    out.M_LaserWidth=[];
    out.M_LaserNumPulses=[];
    out.M_LaserISI=[];
end
if isempty(MtOFF)
    out.MtOFF=[];
    out.mMtOFF=[];
    out.mMtOFFspikecount=[];
    out.sMtOFFspikecount=[];
    out.semMtOFFspikecount=[];
    out.mMtspontOFF=[];
    out.sMtspontOFF=[];
    out.semMtspontOFF=[];
else
    out.MtOFF=MtOFF;
    out.mMtOFF=mMtOFF;
     out.mMtOFFspikecount=mMtOFFspikecount;
     out.sMtOFFspikecount=sMtOFFspikecount;
     out.semMtOFFspikecount=semMtOFFspikecount;
     out.mMtspontOFF=mMtspontOFF;
     out.sMtspontOFF=sMtspontOFF;
     out.semMtspontOFF=semMtspontOFF;
    % out.stimstats=stimstats; %Added by Nick 8/23/18
end
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.xlimits=xlimits;
out.samprate=samprate;

out.VsON=VsON; %vector strength
out.RZON=RZON; %rayleigh statistic
out.p_ON=p_ON; %p-value
out.PhaseON=PhaseON; %spike times represented as phases
out.VsOFF=VsOFF;
out.RZOFF=RZOFF;
out.p_OFF=p_OFF;
out.PhaseOFF=PhaseOFF; %spike times represented as phases

out.LaserRecorded=LaserRecorded; %whether the laser signal was hooked up and recorded as a continuous channel
out.StimRecorded=StimRecorded; %%whether the sound stimulus signal was hooked up and recorded as a continuous channel

if LaserRecorded
    if exist('MtONLaser')
        out.MtONLaser=MtONLaser;
        out.mMtONLaser=mMtONLaser;
    else
        out.MtONLaser=[];
        out.mMtONLaser=[];
    end
    if exist('MtOFFLaser')
        out.MtOFFLaser=MtOFFLaser;
        out.mMtOFFLaser=mMtOFFLaser;
    else
        out.MtOFFLaser=[];
        out.mMtOFFLaser=[];
    end
else
    out.MtONLaser=[];
    out.mMtONLaser=[];
    out.MtOFFLaser=[];
    out.mMtOFFLaser=[];
end
if StimRecorded
    if exist('MtONStim')
        out.MtONStim=MtONStim;
        out.mMtONStim=mMtONStim;
    else
        out.MtONStim=[];
        out.mMtONStim=[];
    end
    if exist('MtOFFStim')
        out.MtOFFStim=MtOFFStim;
        out.mMtOFFStim=mMtOFFStim;
    else
        out.MtOFFStim=[];
        out.mMtOFFStim=[];
    end
else
    out.MtONStim=[];
    out.mMtONStim=[];
    out.MtOFFStim=[];
    out.mMtOFFStim=[];
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


