function ProcessPPI_Behavior(datadir,flag_accel)

%processes PPI behavioral data from djmaus using the tilt platform
%(pressure sensor)
%
% usage: ProcessPPI_Behavior(datadir,flag_accel))
% saves to outfile

if nargin < 2 || isempty(flag_accel)
    flag.accel = 4;
else
    flag.accel = flag_accel;
end

flag.includeALL = 1;     % 0) eliminate traces with too much activity before trial
flag.unwrap = 0;
flag.plot = 1;

if nargin==0
    datadir = pwd;
end

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

try
    load notebook.mat
catch
    warning('could not find notebook file')
end

if flag.includeALL
    fprintf('using all traces\n')
else
    fprintf('eliminating traces with too much activity before test\n')
end

if exist('nb', 'var')
    if isfield(nb, 'mouse2ID')
        fprintf('\nnotebook file has a second mouse ID, will try to process that after the first mouse')
    end
end

%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);


%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

if 0
    fprintf('\nrunning SCT_Monitor to examine soundcard triggers, timestamps, and network messages...')
    
    %   I'm running the soundcard trigger (SCT) into ai1 as another sanity check.
    SCTfname=getSCTfile(datadir);
    Stimfname = getStimfile(datadir);
    if isempty(SCTfname)
        warning('could not find soundcard trigger file')
    else
        [SCTtrace, SCTtimestamps, SCTinfo] =load_open_ephys_data(SCTfname);
        ind0 = find(SCTtrace>0);
        ind1 = find(diff(SCTtimestamps(ind0))>1);
        ind1 = [1; ind1+1];
        SCTtimes = SCTtimestamps(ind0(ind1));
        
        [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(Stimfname);
        ind0 = find(Stimtrace>-4.5);
        ind1 = find(diff(Stimtimestamps(ind0))>1);
        ind1 = [1; ind1+1];
        Stimtimes = Stimtimestamps(ind0(ind1));
        
        StimtimesCorr = Stimtimes-1.04905;
        indNearest = nearest_index(SCTtimes,StimtimesCorr);
        indNearest = indNearest(indNearest>1);
        
    end
    end
%get Events and soundcard trigger timestamps
[Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps

%check if this is an appropriate stimulus protocol
if ~strcmp(GetPlottingFunction(datadir), 'PlotPPI_PSTH')
    error('This does not appear to be a PPI stimulus protcol');
end

%accelerometer channels are 33, 34, 35
node='';
NodeIds=getNodes(pwd);
for i=1:length(NodeIds)
    filename=sprintf('%s_AUX2.continuous', NodeIds{i});
    if exist(filename,'file')
        node=NodeIds{i};
    end
end
filename1=sprintf('%s_AUX1.continuous', node);
filename2=sprintf('%s_AUX2.continuous', node);
filename3=sprintf('%s_AUX3.continuous', node);

fprintf('\n')
if exist(filename1, 'file')~=2 %couldn't find it
    fprintf('could not find AUX1 file %s in datadir %s', filename1, datadir)
else
    [scaledtrace1, datatimestamps, datainfo] =load_open_ephys_data(filename1);
    scaledtrace1 = scaledtrace1 - mean(scaledtrace1);
end
if exist(filename2, 'file')~=2 %couldn't find it
    fprintf('could not find AUX2 file %s in datadir %s', filename2, datadir)
else
    [scaledtrace2, datatimestamps, datainfo] =load_open_ephys_data(filename2);
    scaledtrace2 = scaledtrace2 - mean(scaledtrace2);
end
if exist(filename3, 'file')~=2 %couldn't find it
    fprintf('could not find AUX3 file %s in datadir %s', filename3, datadir)
else
    [scaledtrace3, datatimestamps, datainfo] =load_open_ephys_data(filename3);
    scaledtrace3 = scaledtrace3 - mean(scaledtrace3);
end

% get tilt-table recording from ADC4
filename4=sprintf('%s_ADC4.continuous', node);
if exist(filename4, 'file')~=2 %couldn't find it
    error(sprintf('could not find ADC4 file %s in datadir %s', filename4, datadir))
end
[scaledtrace4, datatimestamps, datainfo] =load_open_ephys_data(filename4);
scaledtrace4 = scaledtrace4 - mean(scaledtrace4);

if 0
    fprintf('\nPlotting rectified traces\n')
    switch flag.accel
        case 1
            scaledtrace = sqrt( scaledtrace1.^2);
            fprintf('\nUsing ACCELEROMETER #1 ONLY\n')
        case 2
            scaledtrace = sqrt( scaledtrace2.^2);
            fprintf('\nUsing ACCELEROMETER #2 ONLY\n')
        case 3
            scaledtrace = sqrt( scaledtrace3.^2);
            fprintf('\nUsing ACCELEROMETER #3 ONLY\n')
        case 4
            scaledtrace = sqrt(scaledtrace4.^2);
            fprintf('\nUsing TILT PLATFORM\n')
        otherwise
            scaledtrace = sqrt( scaledtrace1.^2 + scaledtrace2.^2 + scaledtrace3.^2 );
            fprintf('\nUsing ALL ACCELEROMETERS\n')
    end
else
    fprintf('\nPlotting full (unrectified) traces\n')
    scaledtraceACC = scaledtrace1 + scaledtrace2 + scaledtrace3;
    switch flag.accel
        case 1
            scaledtrace = scaledtrace1;
            fprintf('\nUsing ACCELEROMETER #1 ONLY\n')
        case 2
            scaledtrace = scaledtrace2;
            fprintf('\nUsing ACCELEROMETER #2 ONLY\n')
        case 3
            scaledtrace = scaledtrace3;
            fprintf('\nUsing ACCELEROMETER #3 ONLY\n')
        case 4
            scaledtrace = scaledtrace4;
            fprintf('\nUsing TILT PLATFORM\n')
        otherwise
            scaledtrace = scaledtrace1 + scaledtrace2 + scaledtrace3;
            fprintf('\nUsing ALL ACCELEROMETERS\n')
    end
end

SCTfname=getSCTfile(datadir);
stimfile=sprintf('%s_ADC2.continuous', node);
[stim, stimtimestamps, stiminfo] =load_open_ephys_data(stimfile);

%uncomment this to run some sanity checks
%SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)

fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'ASR')
        j=j+1;
        allsoas(j)=Events(i).soa;
        allsoaflags{j}=Events(i).soaflag;
        allprepulsedurs(j)=Events(i).prepulsedur;
        allprepulseamps(j)=Events(i).prepulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allpulseamps(j)=Events(i).pulseamp;
    end
    
end
prepulsedurs=unique(allprepulsedurs);
pulsedurs=unique(allpulsedurs);
soas=unique(allsoas);
soaflags=unique(allsoaflags);
prepulseamps=unique(allprepulseamps);
pulseamps=unique(allpulseamps);
numprepulsedurs=length(prepulsedurs);
numprepulseamps=length(prepulseamps);
nrepsON=zeros( numprepulsedurs, numprepulseamps);
nrepsOFF=zeros( numprepulsedurs, numprepulseamps);

if length(pulseamps)~=1
    error('not able to handle multiple noiseamps')
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
pulseamp=pulseamps;
soa=soas;
pulsedur=pulsedurs;
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
        warning('ProcessPPI_Behavior: Cannot tell if laser button was turned on in djmaus GUI');
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
nrepsON=zeros(numprepulsedurs, numprepulseamps);
nrepsOFF=zeros(numprepulsedurs, numprepulseamps);

xlimits=[-350 350]; %xlimits for storing traces
startle_window=[0 100]; %hard coded integration region for startle response (ms)
fprintf('\nprocessing with xlimits [%d - %d]', xlimits(1), xlimits(2))
fprintf('\nprocessing with startle integration window [%d - %d]', startle_window(1), startle_window(2))

%extract the traces into a big matrix M
j=0;
if flag.plot
    for idur = 1:length(prepulsedurs)
        Hfig(idur) = figure;
        hold on;
        title(['raw traces for prepulsedur = ' num2str(prepulsedurs(idur)) ' ms'])
    end
end

for i=1:length(Events)
    if strcmp(Events(i).type, 'ASR') 
        prepulsedur=Events(i).prepulsedur;
        switch soaflag
            case 'isi'
                isi=Events(i).soa;
                soa=isi+prepulsedur;
            case 'soa'
                soa=Events(i).soa;
                isi=soa-prepulsedur;
        end
        pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        laser=LaserTrials(i);
        %since this is all about quantifying startle response, we want a trace
        %locked to startle pulse (not prepulse)
        startle_onset=pos+isi/1000;
        start=startle_onset + xlimits(1)/1000; %start is in seconds, should be at xlimits relative to startle onset
        stop=startle_onset + xlimits(2)/1000; %stop is in seconds
        if start>0 %(disallow negative or zero start times)
            ppdindex= find(prepulsedur==prepulsedurs);
            prepulseamp=Events(i).prepulseamp;
            ppaindex= find(prepulseamp==prepulseamps);
            %start=round(pos+xlimits(1)*1e-3*samprate);
            %stop=round(pos+xlimits(2)*1e-3*samprate)-1;
            region=round(start*samprate)+1:round(stop*samprate);
            if isempty(find(region<1))
                if laser
                    nrepsON(ppdindex,ppaindex)=nrepsON(ppdindex,ppaindex)+1;
                    M1ON(ppdindex,ppaindex, nrepsON(ppdindex,ppaindex),:)=scaledtrace(region);
                    M1ONACC(ppdindex,ppaindex, nrepsON(ppdindex,ppaindex),:)=scaledtraceACC(region);
                    M1ONstim(ppdindex, ppaindex, nrepsON(ppdindex, ppaindex),:)=stim(region);
                else
                    if flag.unwrap
                        temp = scaledtrace(region);
                        tempD = [0; diff(temp)];
                        indUp = find(tempD>2);
                        indDown = find(tempD<-2)-1;
                        njumps = min1([length(indUp) length(indDown)]);
                        for ijump = 1:njumps
                            temp(indUp(ijump):indDown(ijump)) = temp(indUp(ijump):indDown(ijump)) - (10-.4096);
                        end
                        scaledtrace(region) = temp;
                    end
                    
                    nrepsOFF(ppdindex,ppaindex)=nrepsOFF(ppdindex,ppaindex)+1;
                    M1OFF(ppdindex,ppaindex, nrepsOFF(ppdindex,ppaindex),:)=scaledtrace(region);
                    M1OFFACC(ppdindex,ppaindex, nrepsOFF(ppdindex,ppaindex),:)=scaledtraceACC(region);
                    M1OFFstim(ppdindex, ppaindex, nrepsOFF(ppdindex, ppaindex),:)=stim(region);
                    
                    if flag.plot
                        figure(Hfig(ppdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, scaledtrace(region), 'color',hsv2rgb([ppdindex/length(prepulsedurs),1,1]))
                    end
                  
                end
            end
        end
    end
end
if flag.plot
    str = cell(1);
    Hfig(length(prepulsedurs)+1) = figure; hold on
    for ifig = 1:length(prepulsedurs)
        figure(Hfig(ifig))
        t=region; t = t-t(1);t=t/samprate;
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(prepulsedurs),1,.5]),'linewidth',2)
        figure(Hfig(length(prepulsedurs)+1))
        plot(t,squeeze(mean(M1OFF(ifig,1,:,:),3)),'color',hsv2rgb([ifig/length(prepulsedurs),1,1]),'linewidth',2)
        str{ifig} = num2str(prepulsedurs(ifig));
    end
    legend(str)
    [~, FN, ~] = fileparts(pwd);
    
    title(FN,'interpreter','none')
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))

PeakON=[];
PeakOFF=[];

mM1OFF=mean(M1OFF, 3);
mM1ON=mean(M1ON, 3);
mM1OFFACC=mean(M1OFFACC, 3);
mM1ONACC=mean(M1ONACC, 3);

mM1OFFstim=mean(M1OFFstim, 3);
mM1ONstim=mean(M1ONstim, 3);


% Accumulate startle response across trials using peak (or summed) rectified signal in region
%%%% Added throwing out oddballs 6/27/2018
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;
SumON=nan(numprepulsedurs, numprepulseamps, max(nrepsON(:)));
SumOFF=nan(numprepulsedurs, numprepulseamps, max(nrepsOFF(:)));
SumONACC=nan(numprepulsedurs, numprepulseamps, max(nrepsON(:)));
SumOFFACC=nan(numprepulsedurs, numprepulseamps, max(nrepsOFF(:)));

PeakON=nan(numprepulsedurs, numprepulseamps, max(nrepsON(:)));
PeakOFF=nan(numprepulsedurs, numprepulseamps, max(nrepsOFF(:)));
PeakONACC=nan(numprepulsedurs, numprepulseamps, max(nrepsON(:)));
PeakOFFACC=nan(numprepulsedurs, numprepulseamps, max(nrepsOFF(:)));

fprintf('\n')
for ppaindex=1:numprepulseamps
    for ppdindex=1:numprepulsedurs; % Hardcoded.
        for k=1:nrepsON(ppdindex, ppaindex);
            %traceON=squeeze(M1ON(ppdindex,ppaindex, k, start:stop));
            %PeakON(ppdindex, ppaindex, k) = max(abs(traceON));
            temp = squeeze(M1ON(ppdindex,ppaindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ON(ppdindex,ppaindex, k, start:stop));
                PeakON(ppdindex, ppaindex, k) = max(abs(traceON));
                SumON(ppdindex, ppaindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ON(ppdindex,ppaindex, k, start:stop));
                PeakON(ppdindex, ppaindex, k) = max(abs(traceON));
                SumON(ppdindex, ppaindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out trial#%d of prepulsedur#%d\n',k,ppdindex)
                PeakON(ppdindex, ppaindex, k) = nan;
                SumON(ppdindex, ppaindex, k) = nan;
            end
            temp = squeeze(M1ONACC(ppdindex,ppaindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceON=squeeze(M1ONACC(ppdindex,ppaindex, k, start:stop));
                PeakONACC(ppdindex, ppaindex, k) = max(abs(traceON));
                SumONACC(ppdindex, ppaindex, k) = sum(abs(traceON));
            elseif flag.includeALL
                traceON=squeeze(M1ONACC(ppdindex,ppaindex, k, start:stop));
                PeakONACC(ppdindex, ppaindex, k) = max(abs(traceON));
                SumONACC(ppdindex, ppaindex, k) = sum(abs(traceON));
            else
                fprintf('Throwing out ACC trial#%d of prepulsedur#%d\n',k,ppdindex)
                PeakONACC(ppdindex, ppaindex, k) = nan;
                SumONACC(ppdindex, ppaindex, k) = nan;
            end
        end
        for k=1:nrepsOFF(ppdindex, ppaindex);
            temp = squeeze(M1OFF(ppdindex,ppaindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceOFF=squeeze(M1OFF(ppdindex,ppaindex, k, start:stop));
                PeakOFF(ppdindex, ppaindex, k) = max(abs(traceOFF));
                SumOFF(ppdindex, ppaindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                traceOFF=squeeze(M1OFF(ppdindex,ppaindex, k, start:stop));
                PeakOFF(ppdindex, ppaindex, k) = max(abs(traceOFF));
                SumOFF(ppdindex, ppaindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of prepulsedur#%d\n',k,ppdindex)
                PeakOFF(ppdindex, ppaindex, k) = nan;
                SumOFF(ppdindex, ppaindex, k) = nan;
            end
            temp = squeeze(M1OFFACC(ppdindex,ppaindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceOFF=squeeze(M1OFFACC(ppdindex,ppaindex, k, start:stop));
                PeakOFFACC(ppdindex, ppaindex, k) = max(abs(traceOFF));
                SumOFFACC(ppdindex, ppaindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                traceOFF=squeeze(M1OFFACC(ppdindex,ppaindex, k, start:stop));
                PeakOFFACC(ppdindex, ppaindex, k) = max(abs(traceOFF));
                SumOFFACC(ppdindex, ppaindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out ACC trial#%d of prepulsedur#%d\n',k,ppdindex)
                PeakOFFACC (ppdindex, ppaindex, k) = nan;
                SumOFFACC (ppdindex, ppaindex, k) = nan;
            end
        end
        if isempty(PeakON)
            mPeakON=[];
            semPeakON=[];
            mPeakONACC=[];
            semPeakONACC=[];
        else
            mPeakON(ppdindex, ppaindex)=median(PeakON(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),3,'omitnan');
            semPeakON(ppdindex, ppaindex)=std(PeakON(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsON(ppdindex, ppaindex));
            mPeakONACC(ppdindex, ppaindex)=median(PeakONACC(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),3,'omitnan');
            semPeakONACC(ppdindex, ppaindex)=std(PeakONACC(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsON(ppdindex, ppaindex));
        end
        if isempty(SumON)
            mSumON=[];
            semSumON=[];
            mSumONACC=[];
            semSumONACC=[];
        else
            mSumON(ppdindex, ppaindex)=median(SumON(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),3,'omitnan');
            semSumON(ppdindex, ppaindex)=std(SumON(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsON(ppdindex, ppaindex));
            mSumONACC(ppdindex, ppaindex)=median(SumONACC(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),3,'omitnan');
            semSumONACC(ppdindex, ppaindex)=std(SumONACC(ppdindex,ppaindex, 1:nrepsON(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsON(ppdindex, ppaindex));
        end
        if isempty(PeakOFF)
            mPeakOFF=[];
            semPeakOFF=[];
            mPeakOFFACC=[];
            semPeakOFFACC=[];
        else
            mPeakOFF(ppdindex, ppaindex)=median(PeakOFF(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),3,'omitnan');
            semPeakOFF(ppdindex, ppaindex)=std(PeakOFF(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsOFF(ppdindex, ppaindex));
            mPeakOFFACC(ppdindex, ppaindex)=median(PeakOFFACC(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),3,'omitnan');
            semPeakOFFACC(ppdindex, ppaindex)=std(PeakOFFACC(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsOFF(ppdindex, ppaindex));
        end
        if isempty(SumOFF)
            mSumOFF=[];
            semSumOFF=[];
            mSumOFFACC=[];
            semSumOFFACC=[];
        else
            mSumOFF(ppdindex, ppaindex)=median(SumOFF(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),3,'omitnan');
            semSumOFF(ppdindex, ppaindex)=std(SumOFF(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsOFF(ppdindex, ppaindex));
            mSumOFFACC(ppdindex, ppaindex)=median(SumOFFACC(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),3,'omitnan');
            semSumOFFACC(ppdindex, ppaindex)=std(SumOFFACC(ppdindex,ppaindex, 1:nrepsOFF(ppdindex, ppaindex)),0,3,'omitnan')/sqrt(nrepsOFF(ppdindex, ppaindex));
        end
    end
    
    %sanity check that first prepulsedur is 0 OR prepulseamp is -1000 (i.e. control condition)
    if ~(prepulsedurs(1)==0 | prepulseamps(1) == -1000)
        error('first prepulsedur is not 0, or firstprepulseamp is not -1000, what is wrong?')
    end
    
    
    %we expect that the pure startle control condition has prepulsedur=0
    %and prepulseamp=-1000
    fprintf('\nusing MEDIAN of peak(abs(trace)) or of sum(abs(trace))  responses\n')
    
    %only makes sense for numprepulsedurs >= 2
    if isempty(PeakON)
        percentPPI_ON=[];
        percentPPI_ONACC=[];
        percentPPI_ONsum=[];
        percentPPI_ONACCsum=[];
        pON=[];
    else
        percentPPI_ON(1)=nan;
        pON(1)=nan;
        for p=2:numprepulsedurs;
            fprintf('Laser ON  pa:%ddB, \n', prepulseamps(ppaindex));
            if flag.accel ==4
                m1=mPeakON(1, ppaindex);
                m2=mPeakON(p, ppaindex);
                percentPPI_ON(p)=((m1-m2)/m1)*100;
                A=PeakON(1,ppaindex, 1:nrepsON(1, ppaindex));
                B=PeakON(p,ppaindex, 1:nrepsON(p, ppaindex));
                [H,pON(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_ON(p),H,pON(p));
                
                m1=mSumON(1, ppaindex);
                m2=mSumON(p, ppaindex);
                percentPPI_ONsum(p)=((m1-m2)/m1)*100;
                A=SumON(1,ppaindex, 1:nrepsON(1, ppaindex));
                B=SumON(p,ppaindex, 1:nrepsON(p, ppaindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_ONsum(p),H,temp);
            end
            
            m1=mPeakONACC(1, ppaindex);
            m2=mPeakONACC(p, ppaindex);
            percentPPI_ONACC(p)=((m1-m2)/m1)*100;
            A=PeakONACC(1,ppaindex, 1:nrepsON(1, ppaindex));
            B=PeakONACC(p,ppaindex, 1:nrepsON(p, ppaindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_ONACC(p),H,temp);
            
            m1=mSumONACC(1, ppaindex);
            m2=mSumONACC(p, ppaindex);
            percentPPI_ONACCsum(p)=((m1-m2)/m1)*100;
            A=SumONACC(1,ppaindex, 1:nrepsON(1, ppaindex));
            B=SumONACC(p,ppaindex, 1:nrepsON(p, ppaindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_ONACCsum(p),H,temp);
        end
    end
    if isempty(PeakOFF)
        percentPPI_OFF=[];
        percentPPI_OFFACC=[];
        percentPPI_OFFsum=[];
        percentPPI_OFFACCsum=[];
        pOFF=[];
    else
        percentPPI_OFF(1)=nan;
        pOFF(1)=nan;
        for p=2:numprepulsedurs;
            fprintf('Laser OFF  pa:%ddB, \n', prepulseamps(ppaindex));
            if flag.accel ==4
                m1=mPeakOFF(1, ppaindex);
                m2=mPeakOFF(p, ppaindex);
                percentPPI_OFF(p)=((m1-m2)/m1)*100;
                A=PeakOFF(1,ppaindex, 1:nrepsOFF(1, ppaindex));
                B=PeakOFF(p,ppaindex, 1:nrepsOFF(p, ppaindex));
                [H,pOFF(p)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_OFF(p),H,pOFF(p));
                
                m1=mSumOFF(1, ppaindex);
                m2=mSumOFF(p, ppaindex);
                percentPPI_OFFsum(p)=((m1-m2)/m1)*100;
                A=SumOFF(1,ppaindex, 1:nrepsOFF(1, ppaindex));
                B=SumOFF(p,ppaindex, 1:nrepsOFF(p, ppaindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_OFFsum(p),H,temp);
            end
            m1=mPeakOFFACC(1, ppaindex);
            m2=mPeakOFFACC(p, ppaindex);
            percentPPI_OFFACC(p)=((m1-m2)/m1)*100;
            A=PeakOFFACC(1,ppaindex, 1:nrepsOFF(1, ppaindex));
            B=PeakOFFACC(p,ppaindex, 1:nrepsOFF(p, ppaindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_OFFACC(p),H,temp);
            
            m1=mSumOFFACC(1, ppaindex);
            m2=mSumOFFACC(p, ppaindex);
            percentPPI_OFFACCsum(p)=((m1-m2)/m1)*100;
            A=SumOFFACC(1,ppaindex, 1:nrepsOFF(1, ppaindex));
            B=SumOFFACC(p,ppaindex, 1:nrepsOFF(p, ppaindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%PPI = %.1f%%,  T-test:%d,  p-value:%.3f\n',prepulsedurs(p),percentPPI_OFFACCsum(p),H,temp);
        end
    end
    
end



%save to outfiles
out.IL=IL;
out.flag=flag;

out.M1ON=M1ON;    % scaledtrace (depends on flag.accel)
out.M1OFF=M1OFF;    % scaledtrace
out.mM1ON=mM1ON;    % scaledtrace
out.mM1OFF=mM1OFF;    % scaledtrace
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
out.datadir=datadir;
out.nrepsON=nrepsON;
out.nrepsOFF=nrepsOFF;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.samprate=samprate;

out.percentPPI_OFF=percentPPI_OFF;
out.pOFF=pOFF;
out.percentPPI_ON=percentPPI_ON;
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


out.numprepulseamps = numprepulseamps;
out.numprepulsedurs = numprepulsedurs;
out.prepulseamps = prepulseamps;
out.prepulsedurs = prepulsedurs;
out.soa=soa;
out.isi=isi;
out.soaflag=soaflag;
out.xlimits=xlimits;
out.startle_window=startle_window;
out.samprate=samprate;
out.datadir=datadir;
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
outfilename=sprintf('outPPI_Behavior.mat');
out.outfilename=outfilename;
save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)

if exist('nb', 'var')
    if isfield(nb, 'mouse2ID')
        filename5=sprintf('%s_ADC5.continuous', node);
        
        
        if exist(filename5, 'file')==2
            fprintf('\nnow processing second mouse...')
            
            ProcessPPI_BehaviorMouse2(out);
        end
    end
end
    
    