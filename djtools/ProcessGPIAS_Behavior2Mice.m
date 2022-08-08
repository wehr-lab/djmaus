function ProcessGPIAS_Behavior2Mice(datadir,flag_accel)

%processes accelerometer behavioral data from djmaus
%
% usage: ProcessGPIAS_Behavior(datadir,flag_accel))
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
if ~strcmp(GetPlottingFunction(datadir), 'PlotGPIAS_PSTH')
    error('This does not appear to be a GPIAS stimulus protcol');
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

% get tilt-table recording from ADC4 for mouse 1
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
stimfile=sprintf('%s_ADC1.continuous', node);
[stim, stimtimestamps, stiminfo] =load_open_ephys_data(stimfile);

%uncomment this to run some sanity checks
%SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)

fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

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
                    M1ON(gdindex,paindex, nrepsON(gdindex,paindex),:)=scaledtrace(region);
                    M1ONACC(gdindex,paindex, nrepsON(gdindex,paindex),:)=scaledtraceACC(region);
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stim(region);
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
                    
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=scaledtrace(region);
                    M1OFFACC(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=scaledtraceACC(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stim(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
                    end
                    if flag.plot
                        %some sanity check to see if the stimulus played
                        %and there's sensor data
                        figure(8)
                        hold on
                        offset=i*.1;
                        t=region;t=t/samprate; t=t-t(1);
                        plot(t, stim(region)-stim(1)+offset, 'm', t, scaledtrace(region)-scaledtrace(1)+offset, 'b')
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
% out.SumOFF=[];
% out.mSumOFF=[];



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

SumON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
SumONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

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
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFF(gdindex, paindex, k) = nan;
                SumOFF(gdindex, paindex, k) = nan;
            end
            temp = squeeze(M1OFFACC(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
                PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
                PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFFACC (gdindex, paindex, k) = nan;
                SumOFFACC (gdindex, paindex, k) = nan;
            end
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
out.SumON=SumON;
out.SumOFF=SumOFF;
out.mSumON=mSumON;
out.mSumOFF=mSumOFF;
out.semSumON=semSumON;
out.semSumOFF=semSumOFF;
out.datadir=datadir;
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
outfilename=sprintf('outGPIAS_Behavior.mat');
out.outfilename=outfilename;
save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)


close all


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
if ~strcmp(GetPlottingFunction(datadir), 'PlotGPIAS_PSTH')
    error('This does not appear to be a GPIAS stimulus protcol');
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

% get tilt-table recording from ADC5 for Mouse2
filename4=sprintf('%s_ADC5.continuous', node);
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
stimfile=sprintf('%s_ADC1.continuous', node);
[stim, stimtimestamps, stiminfo] =load_open_ephys_data(stimfile);

%uncomment this to run some sanity checks
%SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)

fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

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
                    M1ON(gdindex,paindex, nrepsON(gdindex,paindex),:)=scaledtrace(region);
                    M1ONACC(gdindex,paindex, nrepsON(gdindex,paindex),:)=scaledtraceACC(region);
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stim(region);
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
                    
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=scaledtrace(region);
                    M1OFFACC(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=scaledtraceACC(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stim(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
                    end
                    if flag.plot
                        %some sanity check to see if the stimulus played
                        %and there's sensor data
                        figure(8)
                        hold on
                        offset=i*.1;
                        t=region;t=t/samprate; t=t-t(1);
                        plot(t, stim(region)-stim(1)+offset, 'm', t, scaledtrace(region)-scaledtrace(1)+offset, 'b')
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
% out.SumOFF=[];
% out.mSumOFF=[];



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

SumON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
SumONACC=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
SumOFFACC=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));

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
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
                PeakOFF(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFF(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFF(gdindex, paindex, k) = nan;
                SumOFF(gdindex, paindex, k) = nan;
            end
            temp = squeeze(M1OFFACC(gdindex,paindex, k, 1:10000));
            if abs(mean(temp)) <.01 & std(temp) < .1
                traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
                PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            elseif flag.includeALL
                traceOFF=squeeze(M1OFFACC(gdindex,paindex, k, start:stop));
                PeakOFFACC(gdindex, paindex, k) = max(abs(traceOFF));
                SumOFFACC(gdindex, paindex, k) = sum(abs(traceOFF));
            else
                fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                PeakOFFACC (gdindex, paindex, k) = nan;
                SumOFFACC (gdindex, paindex, k) = nan;
            end
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
out.SumON=SumON;
out.SumOFF=SumOFF;
out.mSumON=mSumON;
out.mSumOFF=mSumOFF;
out.semSumON=semSumON;
out.semSumOFF=semSumOFF;
out.datadir=datadir;
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
out.datadir=datadir;
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




