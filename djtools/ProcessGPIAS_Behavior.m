function ProcessGPIAS_Behavior(datadir)

%processes accelerometer behavioral data from djmaus
%
% usage: ProcessGPIAS_Behavior(datadir)
% saves to outfile


if nargin==0
    fprintf('\nno input');
    return;
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

%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);


%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

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
    filename=sprintf('%s_AUX1.continuous', NodeIds{i});
    if exist(filename,'file')
        node=NodeIds{i};
    end
end
filename1=sprintf('%s_AUX1.continuous', node);
filename2=sprintf('%s_AUX2.continuous', node);
filename3=sprintf('%s_AUX3.continuous', node);

if exist(filename1, 'file')~=2 %couldn't find it
    error(sprintf('could not find AUX file %s in datadir %s', filename1, datadir))
end
if exist(filename2, 'file')~=2 %couldn't find it
    error(sprintf('could not find AUX file %s in datadir %s', filename2, datadir))
end
if exist(filename3, 'file')~=2 %couldn't find it
    error(sprintf('could not find AUX file %s in datadir %s', filename3, datadir))
end
fprintf('\n')
[scaledtrace1, datatimestamps, datainfo] =load_open_ephys_data(filename1);
[scaledtrace2, datatimestamps, datainfo] =load_open_ephys_data(filename2);
[scaledtrace3, datatimestamps, datainfo] =load_open_ephys_data(filename3);

%combine X,Y,Z accelerometer channels by RMS
scaledtrace=sqrt(scaledtrace1.^2 + scaledtrace2.^2 + scaledtrace3.^2 );

SCTfname=getSCTfile(datadir);
stimfile=getStimfile; %mw 08.30.2107 old: sprintf('%s_ADC2.continuous', node);
[stim, stimtimestamps, stiminfo] =load_open_ephys_data(stimfile);

%uncomment this to run some sanity checks
SCT_Monitor(datadir, StartAcquisitionSec, Events, all_channels_data, all_channels_timestamps, all_channels_info)

fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS')
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
M1ONstim=[];M1OFFstim=[];
nrepsON=zeros(numgapdurs, numpulseamps);
nrepsOFF=zeros(numgapdurs, numpulseamps);

xlimits=[-200 200]; %xlimits for storing traces
startle_window=[0 150]; %hard coded integration region for startle response

fprintf('\nprocessing with xlimits [%d - %d]', xlimits(1), xlimits(2))
fprintf('\nprocessing with startle integration window [%d - %d]', startle_window(1), startle_window(2))

%extract the traces into a big matrix M
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'GPIAS')
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
                    M1ONstim(gdindex, paindex, nrepsON(gdindex, paindex),:)=stim(region);
                else
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    M1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex),:)=scaledtrace(region);
                    M1OFFstim(gdindex, paindex, nrepsOFF(gdindex, paindex),:)=stim(region);
                    
                    %                     figure(8),clf;hold on
                    %                     t=region;t=t/samprate;
                    %                     plot(t, stim(region), 'm', t, scaledtrace(region), 'b')
                    %                     gap_termination=pos+gapdelay/1000;
                    %                     gap_onset=pos+gapdelay/1000-gapdur/1000;
                    %                     plot(gap_onset, 0, '^', gap_termination,0, 'v')
                    %                     plot(startle_onset, 0, 'bo')
                    %                     plot(pos, 0, 'r*')
                    %                     keyboard
                end
            end
        end
    end
end

fprintf('\nmin num ON reps: %d\nmax num ON reps: %d', min(nrepsON(:)), max(nrepsON(:)))
fprintf('\nmin num OFF reps: %d\nmax num OFF reps: %d',min(nrepsOFF(:)), max(nrepsOFF(:)))

PeakON=[];
PeakOFF=[];

mM1OFF=mean(M1OFF, 3);
mM1ON=mean(M1ON, 3);
mM1OFFstim=mean(M1OFFstim, 3);
mM1ONstim=mean(M1ONstim, 3);


% Accumulate startle response across trials using peak rectified signal in region
start=(startle_window(1)-xlimits(1))*samprate/1000;
stop=start+diff(startle_window)*samprate/1000;
PeakON=nan(numgapdurs, numpulseamps, max(nrepsON(:)));
PeakOFF=nan(numgapdurs, numpulseamps, max(nrepsOFF(:)));
for paindex=1:numpulseamps
    for gdindex=1:numgapdurs; % Hardcoded.
        for k=1:nrepsON(gdindex, paindex);
            traceON=squeeze(M1ON(gdindex,paindex, k, start:stop));
            PeakON(gdindex, paindex, k)=max(abs(traceON));
        end
        for k=1:nrepsOFF(gdindex, paindex);
            traceOFF=squeeze(M1OFF(gdindex,paindex, k, start:stop));
            PeakOFF(gdindex, paindex, k)=max(abs(traceOFF));
        end
        if isempty(PeakON)
            mPeakON=[];
            semPeakON=[];
        else
            mPeakON(gdindex, paindex)=mean(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),3);
            semPeakON(gdindex, paindex)=std(PeakON(gdindex,paindex, 1:nrepsON(gdindex, paindex)),0,3)/sqrt(nrepsON(gdindex, paindex));
        end
        if isempty(PeakOFF)
            mPeakOFF=[];
            semPeakOFF=[];
        else
            mPeakOFF(gdindex, paindex)=mean(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),3);
            semPeakOFF(gdindex, paindex)=std(PeakOFF(gdindex,paindex, 1:nrepsOFF(gdindex, paindex)),0,3)/sqrt(nrepsOFF(gdindex, paindex));
        end
        
        
    end
    
    %sanity check that first gapdur is 0 (i.e. control condition)
    if gapdurs(1)~=0
        error('first gapdur is not 0, what is wrong?')
    end
    
    %only makes sense for numgapdurs==2
    fprintf('\n')
    if isempty(PeakON)
        percentGPIAS_ON=[];
        pON=[];
    else
        percentGPIAS_ON(1)=nan;
        pON(1)=nan;
        for p=2:numgapdurs;
            m1=mPeakON(1, paindex);
            m2=mPeakON(p, paindex);
            percentGPIAS_ON(p)=((m1-m2)/m1)*100;
            A=PeakON(1,paindex, 1:nrepsON(1, paindex));
            B=PeakON(p,paindex, 1:nrepsON(p, paindex));
            [H,pON(p)]=ttest2(A,B);
            fprintf('\nLaser ON  pa:%ddB, ', pulseamps(paindex));
            fprintf('gd: %dms %%GPIAS = %.1f%%, T-test:%d, p-value:%.3f',gapdurs(p),percentGPIAS_ON(p),H,pON(p));
        end
    end
    fprintf('\n')
    if isempty(PeakOFF)
        percentGPIAS_OFF=[];
        pOFF=[];
    else
        percentGPIAS_OFF(1)=nan;
        pOFF(1)=nan;
        for p=2:numgapdurs;
            m1=mPeakOFF(1, paindex);
            m2=mPeakOFF(p, paindex);
            percentGPIAS_OFF(p)=((m1-m2)/m1)*100;
            A=PeakOFF(1,paindex, 1:nrepsOFF(1, paindex));
            B=PeakOFF(p,paindex, 1:nrepsOFF(p, paindex));
            [H,pOFF(p)]=ttest2(A,B);
            fprintf('\nLaser OFF  pa:%ddB, ', pulseamps(paindex));
            fprintf('gd: %dms %%GPIAS = %.1f%%, T-test:%d, p-value:%.3f',gapdurs(p),percentGPIAS_OFF(p),H,pOFF(p));
            
        end
    end
    
end



%save to outfiles
out.IL=IL;

out.M1ON=M1ON;
out.M1OFF=M1OFF;
out.mM1ON=mM1ON;
out.mM1OFF=mM1OFF;
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
catch
    out.nb='notebook file missing';
    out.stimlog='notebook file missing';
    out.user='unknown';
end
outfilename=sprintf('outGPIAS_Behavior.mat');
save (outfilename, 'out')
fprintf('\n\nsaved outfile %s \nin directory %s\n', outfilename, pwd)





