function ProcessTC_2P(varargin)
%processes dF/F mesoscope 2P tuning curve data from djmaus running in meso mode.
%sorts dF/F  data into a big response matrix.
%
%usage: ProcessTC_2P(datapath, [channel], [xlimits], [ylimits])
%datapath is the djmaus-created directory with the notebook file in it
%channel is a number not a string
%saves output in an outfile
%
%notes: whitenoise plotted as freq=-1 kHz, silent sound as -2 kHz

djPrefs;
global pref

if nargin==0
    fprintf('\nno input\n')
    return
else
    datadir=varargin{1};
end
if nargin==1
    xlimits=[0 5000]; %x limits for axis
    ylimits=[-.1 .2];
    channel=1;
elseif nargin==2
    channel=varargin{2};
    xlimits=[0 5000]; %x limits for axis
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
if ischar(channel) channel=str2num(channel);end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

try
    load notebook.mat
catch
    warning('could not find notebook file')
end

%assume this directory structure:
cd ..
cd suite2p_output
cd suite2p
cd combined
iscell=0; %initialize because iscell is a built-in function
load Fall
numcells=sum(iscell(:,1));
fprintf('\n%d cells', numcells)
% sort rows such that cells are at the top
[iscell_sorted, I]=sortrows(iscell, 2, 'descend');
i=1;
f=F(I(i),:);
f0=prctile(f, 10);
dff=(f-f0)/f0;
scaledtrace=dff;

cd(datadir)
cd ..
d=dir('*.smrx')
if length(d)<1 error('no smrx file found')
elseif length(d)>1
    warning('multiple smrx files found, using the first one')
end
sctfilename=replace(d(1).name, '.smrx', '.mat');
load (sctfilename)
if num_ProtocolStarts>1 error('wtf');end %sanity check; there can be only one protocol start

% % %read messages
% % messagesfilename='messages.events';
% % [messages] = GetNetworkEvents(messagesfilename);


% % %read digital Events
% % Eventsfilename='all_channels.events';
% % [all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=ops.fs; %mesoscope framerate in Hz

%get Events and soundcard trigger timestamps
%[Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
%there are some general notes on the format of Events and network messages in help GetEventsAndSCT_Timestamps
numstim=length(stimlog);
if numstim~=num_events error('Number of sound stimuli (from stimlog) does not match Number of hardware triggers (soundcardtrig TTLs from Spike2)');end
   % THERE_IS_A_PROBLEM

 Events=stimlog;
   for i=1:numstim
    Events(i).soundcard_trigger_timestamp_sec=event_times_sec(i)-ProtocolStart_secs;
end

%check if this is an appropriate stimulus protocol
switch GetPlottingFunction(datadir)
    case {'PlotTC_LFP', 'PlotTC_PSTH'}
    otherwise
        error('This does not appear to be a tuning curve stimulus protcol');
end


        LaserRecorded=0;
        StimRecorded=0;

%Here I should load the stimulus trace from the Spike2 file

% try
%     [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(getStimfile('.'));
%     %Stimtimestamps=Stimtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
%     Stimtimestamps=Stimtimestamps-Stimtimestamps(1);
%     %  Stimtrace=Stimtrace./max(abs(Stimtrace));
%     fprintf('\nsuccessfully loaded stim trace')
% catch
%     [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data('105_ADC3.continuous');
%     fprintf('\ncouldnt find a stimtrace, loaded a differnt trace for plotting')
%     
%         LaserRecorded=0;
%     if isempty(getStimfile('.'))
%         StimRecorded=0;
%     else
%         StimRecorded=1;
%     end
%     
% 
%     if StimRecorded
%         try
%             [Stimtrace, Stimtimestamps, Stiminfo] =load_open_ephys_data(getStimfile('.'));
%             Stimtimestamps=Stimtimestamps-StartAcquisitionSec; %zero timestamps to start of acquisition
%             Stimtrace=Stimtrace./max(abs(Stimtrace));
%             fprintf('\nsuccessfully loaded stim trace')
%         catch
%             fprintf('\nfound stim file %s but could not load stim trace', getStimfile('.'))
%         end
%     else
%         fprintf('\nSound stimulus trace not recorded')
%     end
% end






fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

%get freqs/amps
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'silentsound') ...
            |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise')| strcmp(Events(i).type, 'grating')
        j=j+1;
        alldurs(j)=Events(i).param.duration;
        if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, '2tone')
            allamps(j)=Events(i).param.amplitude;
            allfreqs(j)=Events(i).param.frequency;
        elseif strcmp(Events(i).type, 'whitenoise')
            allamps(j)=Events(i).param.amplitude;
            allfreqs(j)=-1000;
        elseif strcmp(Events(i).type, 'silentsound')
            allfreqs(j)=-2000;
            allamps(j)=nan; %flagging silent sound amp as nan
        elseif strcmp(Events(i).type, 'fmtone')
            allamps(j)=Events(i).param.amplitude;
            allfreqs(j)=Events(i).param.carrier_frequency;
        elseif strcmp(Events(i).type, 'grating')
            allfreqs(j)=Events(i).param.angle*1000;
            allamps(j)=Events(i).param.spatialfrequency;
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


M1=[];
M1OFF=[];
 M1OFFLaser=[];
M1OFFStim=[];

nreps=zeros(numfreqs, numamps, numdurs);
nrepsOFF=zeros(numfreqs, numamps, numdurs);



fprintf('\ncomputing tuning curve...');

samprate=sampleRate;



%extract the traces into a big matrix M
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') |  strcmp(Events(i).type, 'silentsound') | ...
            strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating')
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate; %pos is in samples (2P frames)
        end
        start=round(pos+xlimits(1)*1e-3*samprate);
        win=round(diff(xlimits)*1e-3*samprate);
        stop=start+win;
        region=start:stop;


        if isempty(find(region<1)) %(disallow negative or zero start times)
            switch Events(i).type
                case {'tone', '2tone'}
                    freq=Events(i).param.frequency;
                    amp=Events(i).param.amplitude;
                case 'fmtone'
                    freq=Events(i).param.carrier_frequency;%
                    amp=Events(i).param.amplitude;
                case 'whitenoise'
                    freq=-1000;
                    amp=Events(i).param.amplitude;
                case 'silentsound'
                    freq=-2000;
                    amp=min(amps); %put silentsound in it's own column (freq=-2) in the lowest row
                case 'grating'
                    amp=Events(i).param.spatialfrequency;
                    freq=Events(i).param.angle*1000;
            end
            
            
            dur=Events(i).param.duration;
            findex= find(freqs==freq);
            aindex= find(amps==amp);
            dindex= find(durs==dur);
            nreps(findex, aindex, dindex)=nreps(findex, aindex, dindex)+1;
            M1(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=scaledtrace(region);
            %M1stim(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=Stimtrace(region);
                nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
                M1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=scaledtrace(region);
             %   M1OFFStim(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=Stimtrace(region);
            
        end
    end
end

%region=length(M1OFF);
traces_to_keep=[];
if ~isempty(traces_to_keep)
    fprintf('\n using only traces %d, discarding others', traces_to_keep);
    mM1=mean(M1(:,:,:,traces_to_keep,:), 4);
    mM1OFF=mean(M1OFF(:,:,:,traces_to_keep,:), 4);
else
    for aindex=1:numamps
        for findex=1:numfreqs
            for dindex=1:numdurs
                if nreps(findex, aindex, dindex)>0
                    mM1(findex, aindex, dindex,:)=mean(M1(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                   % mM1stim(findex, aindex, dindex,:)=mean(M1stim(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1(findex, aindex, dindex,:)=zeros(size(region));
                    %mM1stim(findex, aindex, dindex,:)=zeros(size(region));
                end
           
                if nrepsOFF(findex, aindex, dindex)>0
                    mM1OFF(findex, aindex, dindex,:)=mean(M1OFF(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    %mM1OFFStim(findex, aindex, dindex,:)=mean(M1OFFStim(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                    %mM1OFFLaser(findex, aindex, dindex,:)=mean(M1OFFLaser(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1OFF(findex, aindex, dindex,:)=zeros(size(region));
                    %mM1OFFStim(findex, aindex, dindex,:)=zeros(size(region));
                    %mM1OFFLaser(findex, aindex, dindex,:)=zeros(size(region));
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




%assign outputs
out.scaledtrace=scaledtrace;
out.M1=M1;
out.M1OFF=M1OFF;
out.M1OFFStim=M1OFFStim;
out.M1OFFLaser=M1OFFLaser;
out.mM1OFF=mM1OFF;
% out.mM1OFFStim=mM1OFFStim;
% out.M1stim=M1stim;
% out.mM1stim=mM1stim;
out.mM1=mM1;
out.datadir=datadir;
out.freqs=freqs;
out.amps=amps;
out.durs=durs;
out.nreps=nreps;
out.nrepsOFF=nrepsOFF;
out.numfreqs=numfreqs;
out.numamps=numamps;
out.numdurs=numdurs;
out.traces_to_keep=traces_to_keep;
out.Events=Events;
out.xlimits=xlimits;
out.ylimits=ylimits;
out.samprate=samprate;
% out.nstd=nstd;
out.channel=channel;
%should probably save header info and stuff like that



outfilename=sprintf('outLFP_ch%d.mat',channel);
save(outfilename, 'out')
outfilename=sprintf('outLFP_ch%d',channel);
save(outfilename, 'out', '-v7.3' )
fprintf('\n saved to %s', outfilename)

