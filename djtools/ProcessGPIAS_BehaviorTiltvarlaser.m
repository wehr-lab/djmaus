function ProcessGPIAS_BehaviorTiltvarlaser(datadir,flag_accel)

%processes piezo (tilt platform) and accelerometer behavioral data from djmaus
%
% usage: ProcessGPIAS_BehaviorTiltvarlaser(datadir,flag_accel))
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
if ~(strcmp(GetPlottingFunction(datadir), 'PlotGPIASflashtrain_PSTH') | ...
        strcmp(GetPlottingFunction(datadir), 'PlotGPIASLaserPulse_PSTH') )
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

%check for laser in Events
% LaserScheduled = zeros(1,length(Events));
% LaserOnOffButton = zeros(1,length(Events));
% LaserTrials = zeros(1,length(Events));
%
% allVarLaserstarts=[];
% alltrainnumpulses=[];
% alltrainisis=[];
% alltrainpulsewidths=[];
% allpulsewidths=[];
% for i=1:length(Events)
%     VarLaser(i)=Events(i).VarLaser;
%     if Events(i).VarLaser
%         %        Events(i).VarLaserstart
%         %        Events(i).VarLaserpulsewidth
%         %        Events(i).VarLaserisi
%         %        Events(i).VarLasernumpulses
%         allVarLaserstarts=[allVarLaserstarts Events(i).VarLaserstart];
%         allpulsewidths=[ allpulsewidths Events(i).VarLaserpulsewidth];
%         alltrainisis=[alltrainisis Events(i).VarLaserisi]; %isi in train
%         alltrainnumpulses=[alltrainnumpulses Events(i).VarLasernumpulses]; %isi in train
%         LaserPulsewidth(i)=Events(i).VarLaserpulsewidth;
%     else
%         LaserPulsewidth(i)=0;
%         allpulsewidths=[ allpulsewidths 0];
%
%     end
% end
% laserstarts=unique(allVarLaserstarts);
% trainnumpulses=unique(alltrainnumpulses);
% trainisis=unique(alltrainisis);
% pulsewidths=unique(allpulsewidths);
% numpulsewidths=length(pulsewidths);

samprate=sampleRate;

%get gpias params used
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
nreps=zeros( numgapdurs, numpulseamps);

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

LaserScheduled = zeros(1,length(Events));
LaserOnOffButton = zeros(1,length(Events));
LaserTrials = zeros(1,length(Events));

allVarLaserstarts=[];
alltrainnumpulses=[];
alltrainisis=[];
alltrainpulsewidths=[];
allpulsewidths=[];
for i=1:length(Events)
    if isempty(Events(i).VarLaser)
        VarLaser(i)=0;
    else
        VarLaser(i)=Events(i).VarLaser;
    end
    if VarLaser(i)
        allVarLaserstarts=[allVarLaserstarts Events(i).VarLaserstart];
        allpulsewidths=[ allpulsewidths Events(i).VarLaserpulsewidth];
        alltrainisis=[alltrainisis Events(i).VarLaserisi]; %isi in train
        alltrainnumpulses=[alltrainnumpulses Events(i).VarLasernumpulses]; %isi in train
        LaserPulsewidth(i)=Events(i).VarLaserpulsewidth;
        LaserStart(i)=Events(i).VarLaserstart;
    else
        LaserPulsewidth(i)=0;
        allpulsewidths=[ allpulsewidths 0];
        LaserStart(i)=-666;
    end
end
laserstarts=unique(allVarLaserstarts);
numlaserstarts=length(laserstarts);
trainnumpulses=unique(alltrainnumpulses);
trainisis=unique(alltrainisis);
pulsewidths=unique(allpulsewidths);
numpulsewidths=length(pulsewidths);
fprintf('\n%d laser trials in this Events file', sum(VarLaser))
fprintf('\n%d laser pulsewidths:', numpulsewidths)
fprintf('%g ', pulsewidths)
fprintf('\n%d laser starts:', numlaserstarts)
fprintf('%g ', laserstarts-gapdelay)

%there is only 1 matrix, containing both laser on and off
%pulsewidth 0 is where we put the laser-off trials
M1=[];
nreps=zeros(numgapdurs, numpulsewidths, numlaserstarts); %we don't need to add 1 because pw 0 was inlcuded above for laser-off
M1ACC=[];
M1stim=[];

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
        laser=VarLaser(i);
        %since this is all about quantifying startle response, we want a trace
        %locked to startle pulse (not gap)
        startle_onset=pos+gapdelay/1000+isi/1000;
        start=startle_onset + xlimits(1)/1000; %start is in seconds, should be at xlimits relative to startle onset
        stop=startle_onset + xlimits(2)/1000; %stop is in seconds
        if start>0 %(disallow negative or zero start times)
            gdindex= find(gapdur==gapdurs);
            pulsewidth=LaserPulsewidth(i);
            pwindex= find(pulsewidth==pulsewidths);
            laserstart=LaserStart(i);
            if laserstart==-666
                lsindex=1; %I'll put laser-off trials here for now, consider adding an additional laserstart of -666?
            else
                lsindex= find(laserstart==laserstarts);
            end
            nreps(gdindex,pwindex,lsindex)=nreps(gdindex,pwindex,lsindex)+1;
            %start=round(pos+xlimits(1)*1e-3*samprate);
            %stop=round(pos+xlimits(2)*1e-3*samprate)-1;
            region=round(start*samprate)+1:round(stop*samprate);
            if isempty(find(region<1))
                if 0
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
                    
                    M1(gdindex,pwindex,lsindex, nreps(gdindex,pwindex, lsindex),:)=scaledtrace(region);
                    M1ACC(gdindex,pwindex,lsindex, nreps(gdindex,pwindex, lsindex),:)=scaledtraceACC(region);
                    M1stim(gdindex, pwindex,lsindex, nreps(gdindex,pwindex, lsindex),:)=stim(region);
                    
                    if flag.plot
                        figure(Hfig(gdindex))
                        t=region; t = t-t(1);t=t/samprate;
                        plot(t, scaledtrace(region), 'color',hsv2rgb([gdindex/length(gapdurs),1,1]))
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
        plot(t,squeeze(mean(M1(ifig,1,1,:,:),4)),'color',hsv2rgb([ifig/length(gapdurs),1,.5]),'linewidth',2)
        figure(Hfig(length(gapdurs)+1))
        plot(t,squeeze(mean(M1(ifig,1,1,:,:),4)),'color',hsv2rgb([ifig/length(gapdurs),1,1]),'linewidth',2)
        str{ifig} = num2str(gapdurs(ifig));
    end
    legend(str)
    [~,FN,~]=fileparts(datadir);
    title(FN,'interpreter','none')
end

fprintf('\nmin num reps: %d\nmax num reps: %d', min(nreps(:)), max(nreps(:)))

PeakTilt=[];
mM1=mean(M1, 4); %average across reps
mM1ACC=mean(M1ACC, 4);
mM1stim=mean(M1stim, 4);


% Accumulate startle response across trials using peak (or summed) rectified signal in region
%%%% Added throwing out oddballs 6/27/2018
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;
SumTilt=nan(numgapdurs, numpulsewidths, numlaserstarts,max(nreps(:)));
SumACC=SumTilt;
PeakTilt=SumTilt;
PeakACC=SumTilt;

fprintf('\n')
for lsindex=1:numlaserstarts
    for pwindex=1:numpulsewidths
        for gdindex=1:numgapdurs; % 
            for k=1:nreps(gdindex, pwindex);
                %traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
                %PeakON(gdindex, paindex, k) = max(abs(traceON));
                temp = squeeze(M1(gdindex,pwindex, lsindex, k, 1:10000));
                if abs(mean(temp)) <.01 & std(temp) < .1
                    traceTilt=squeeze(M1(gdindex,pwindex, lsindex, k, start:stop));
                    PeakTilt(gdindex, pwindex, lsindex, k) = max(abs(traceTilt));
                    SumTilt(gdindex, pwindex, lsindex, k) = sum(abs(traceTilt));
                elseif flag.includeALL
                    traceTilt=squeeze(M1(gdindex,pwindex, lsindex, k, start:stop));
                    PeakTilt(gdindex, pwindex, lsindex, k) = max(abs(traceTilt));
                    SumTilt(gdindex, pwindex, lsindex, k) = sum(abs(traceTilt));
                else
                    fprintf('Throwing out trial#%d of gapdur#%d\n',k,gdindex)
                    PeakTilt(gdindex, pwindex, lsindex, k) = nan;
                    SumTilt(gdindex, pwindex, lsindex, k) = nan;
                end
                temp = squeeze(M1ACC(gdindex,pwindex, lsindex, k, 1:10000));
                if abs(mean(temp)) <.01 & std(temp) < .1
                    traceACC=squeeze(M1ACC(gdindex,pwindex, lsindex, k, start:stop));
                    PeakACC(gdindex, pwindex, lsindex, k) = max(abs(traceTilt));
                    SumACC(gdindex, pwindex, lsindex, k) = sum(abs(traceTilt));
                elseif flag.includeALL
                    traceACC=squeeze(M1ACC(gdindex,pwindex, lsindex, k, start:stop));
                    PeakACC(gdindex, pwindex, lsindex, k) = max(abs(traceTilt));
                    SumACC(gdindex, pwindex, lsindex, k) = sum(abs(traceTilt));
                else
                    fprintf('Throwing out ACC trial#%d of gapdur#%d\n',k,gdindex)
                    PeakACC(gdindex, pwindex, lsindex, k) = nan;
                    SumACC(gdindex, pwindex, lsindex, k) = nan;
                end
            end
        end
    end
end

for lsindex=1:numlaserstarts
    for pwindex=1:numpulsewidths
        for gdindex=1:numgapdurs; % Hardcoded.
            
            mPeak(gdindex, pwindex, lsindex)=median(PeakTilt(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),4,'omitnan');
            semPeak(gdindex, pwindex, lsindex)=std(PeakTilt(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),0,4,'omitnan')/sqrt(nreps(gdindex, pwindex, lsindex));
            mPeakACC(gdindex, pwindex, lsindex)=median(PeakACC(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),4,'omitnan');
            semPeakACC(gdindex, pwindex, lsindex)=std(PeakACC(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),0,4,'omitnan')/sqrt(nreps(gdindex, pwindex, lsindex));
            
            mSum(gdindex, pwindex, lsindex)=median(SumTilt(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),4,'omitnan');
            semSum(gdindex, pwindex, lsindex)=std(SumTilt(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),0,4,'omitnan')/sqrt(nreps(gdindex, pwindex, lsindex));
            mSumACC(gdindex, pwindex, lsindex)=median(SumACC(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),4,'omitnan');
            semSumACC(gdindex, pwindex, lsindex)=std(SumACC(gdindex,pwindex, lsindex, 1:nreps(gdindex, pwindex, lsindex)),0,4,'omitnan')/sqrt(nreps(gdindex, pwindex, lsindex));
        end
    end
end


%sanity check that first gapdur is 0 (i.e. control condition)
if gapdurs(1)~=0
    error('first gapdur is not 0, what is wrong?')
end

fprintf('\nusing MEDIAN of peak(abs(trace)) or of sum(abs(trace))  responses\n')

for lsindex=1:numlaserstarts
    for pwindex=1:numpulsewidths
        fprintf('\nLaser pw:%d ms, laserstart:%d ms\n', pulsewidths(pwindex), laserstarts(lsindex));
        
        %only makes sense for numgapdurs >= 2
        percentGPIAS(1, pwindex, lsindex)=nan;
        pTilt(1, pwindex, lsindex)=nan;
        for p=2:numgapdurs;
            if flag.accel ==4
                m1=mPeak(1, pwindex, lsindex);
                m2=mPeak(p, pwindex, lsindex);
                percentGPIAS(p, pwindex, lsindex)=((m1-m2)/m1)*100;
                A=PeakTilt(1,pwindex, lsindex, 1:nreps(1, pwindex, lsindex));
                B=PeakTilt(p,pwindex, lsindex, 1:nreps(p, pwindex, lsindex));
                [H,pTilt(p, pwindex, lsindex)]=ttest2(A,B);
                fprintf('TILT(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS(p),H,pTilt(p));
                
                m1=mSum(1, pwindex, lsindex);
                m2=mSum(p, pwindex, lsindex);
                percentGPIAS_sum(p, pwindex, lsindex)=((m1-m2)/m1)*100;
                A=SumTilt(1,pwindex, lsindex, 1:nreps(1, pwindex, lsindex));
                B=SumTilt(p,pwindex, lsindex, 1:nreps(p, pwindex, lsindex));
                [H,temp]=ttest2(A,B);
                fprintf('TILT(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_sum(p),H,temp);
            end
            
            m1=mPeakACC(1, pwindex, lsindex);
            m2=mPeakACC(p, pwindex, lsindex);
            percentGPIAS_ACC(p, pwindex, lsindex)=((m1-m2)/m1)*100;
            A=PeakACC(1,pwindex, lsindex, 1:nreps(1, pwindex, lsindex));
            B=PeakACC(p,pwindex, lsindex, 1:nreps(p, pwindex, lsindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(peak): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ACC(p),H,temp);
            
            m1=mSumACC(1, pwindex, lsindex);
            m2=mSumACC(p, pwindex, lsindex);
            percentGPIAS_ACCsum(p, pwindex, lsindex)=((m1-m2)/m1)*100;
            A=SumACC(1,pwindex, lsindex, 1:nreps(1, pwindex, lsindex));
            B=SumACC(p,pwindex, lsindex, 1:nreps(p, pwindex, lsindex));
            [H,temp]=ttest2(A,B);
            fprintf('ACCEL(sum): gd: %dms  %%GPIAS = %.1f%%,  T-test:%d,  p-value:%.3f\n',gapdurs(p),percentGPIAS_ACCsum(p),H,temp);
        end
    end
end








%save to outfiles
% (numgapdurs, numpulsewidths, numlaserstarts)
out.IL=1;

out.M1=M1;    % scaledtrace (depends on flag.accel)
out.mM1=mM1;    % scaledtrace
out.M1ACC=M1ACC;
out.mM1ACC=mM1ACC;
out.M1stim=M1stim;
out.mM1stim=mM1stim;

out.Peak=PeakTilt;
out.mPeak=mPeak;
out.semPeak=semPeak;
out.mSumACC=mSumACC;
out.SumACC=SumACC;
out.mPeakACC=mPeakACC;
out.mSum=mSum;
out.Sum=SumTilt;


out.datadir=datadir;
out.nreps=nreps;
out.Events=Events;
out.LaserTrials=LaserTrials;
out.samprate=samprate;

out.percentGPIAS=percentGPIAS;
out.percentGPIAS_sum=percentGPIAS_sum;
out.percentGPIAS_ACCsum=percentGPIAS_ACCsum;



out.p=pTilt;

out.LaserStarts=unique(LaserStart); %this currently includes the -666 laser-off flag, and is not corrected for gapdelay
out.LaserWidths=unique(LaserPulsewidth);
out.LaserNumPulses=trainnumpulses;

out.trainnumpulses=trainnumpulses;
out.numlaserstarts=numlaserstarts;
out.numpulseamps = numpulseamps;
out.numpulsewidths = numpulsewidths;
out.numgapdurs = numgapdurs;
out.pulseamps = pulseamps;
out.pulsewidths = pulsewidths;
out.laserstarts=laserstarts-gapdelay;   %does not include the -666 laser-off flag, and is corrected for gapdelay


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
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
end
outfilename=sprintf('outGPIAS_Behavior.mat');
save (outfilename, 'out')
fprintf('\nsaved outfile %s \nin directory %s\n', outfilename, pwd)