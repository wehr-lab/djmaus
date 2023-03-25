function ProcessTC_2P(varargin)
%processes dF/F mesoscope 2P tuning curve data from djmaus running in meso mode.
%sorts dF/F  data into a response matrix.
%
%usage: ProcessTC_2P([datapath], [xlimits], [ylimits])
%datapath is the djmaus-created directory with the notebook file in it
%datapath defaults to current directory
%saves output in an outfile
%
%I'm switching to including all cells, since it makes no sense to have
%thousands of outfiles.
%
%notes: whitenoise plotted as freq=-1 kHz, silent sound as -2 kHz

djPrefs;
global pref

if nargin==0
    datadir=pwd;
    xlimits=[0 5000]; %x limits for axis
    ylimits=[-.1 .2];
elseif nargin==1
    datadir=varargin{1};
    xlimits=[0 5000]; %x limits for axis
    ylimits=[-.1 .2];
elseif nargin==2
    xlimits=varargin{2};
    ylimits=[-.1 .2];
elseif nargin==3
    xlimits=varargin{2};
    ylimits=varargin{3};
else
    error('wrong number of arguments');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

djPrefs;
global pref
cd (pref.datapath);
cd(datadir)

try
    load notebook.mat
catch
    warning('could not find notebook file in this data directory. For datapath use the djmaus-created directory with the notebook file in it.')
    return
end

%assume this directory structure:
cd ..
cd suite2p_output
cd suite2p
cd combined
iscell=0; %initialize because iscell is a built-in function
load Fall
numcells=sum(iscell(:,1));
numframes=size(F, 2);
fprintf('\n%d cells', numcells)
% sort rows such that best cells (highest cell likelihood) are at the top
%exclude non-cells 
[iscell_sorted, I]=sortrows(iscell, 2, 'descend');
f=F(I(1:numcells),:);
f0=prctile(f', 10);
f0=repmat(f0', 1, numframes);
dff=(f-f0)./f0;

cd(datadir)
cd ..
d=dir('*.smrx');
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
M1Stim=[];

nreps=zeros(numfreqs, numamps, numdurs);




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
            M1(findex,aindex,dindex, nreps(findex, aindex, dindex),:,:)=dff(:,region);
            %M1 is freqs x amps x durs x reps x cells x frames
            %M1stim(findex,aindex,dindex, nreps(findex, aindex, dindex),:)=Stimtrace(region);
            
%             if freq==-1000 & amp == 80
%                 keyboard
%             end

        end
    end
end

traces_to_keep=[];
if ~isempty(traces_to_keep)
    fprintf('\n using only traces %d, discarding others', traces_to_keep);
    mM1=mean(M1(:,:,:,traces_to_keep,:,:), 4);
else
    for aindex=1:numamps
        for findex=1:numfreqs
            for dindex=1:numdurs
                if nreps(findex, aindex, dindex)>0
                    mM1(findex, aindex, dindex,:,:)=mean(M1(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:,:), 4);
                    %mM1 is freqs x amps x durs x cells x frames

                   % mM1stim(findex, aindex, dindex,:)=mean(M1stim(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1(findex, aindex, dindex,:,:)=zeros(numcells, length(region));
                    %mM1stim(findex, aindex, dindex,:)=zeros(size(region));
                end
           
                
            end
        end
    end
end








%assign outputs
out.dff=dff;
out.M1=M1;
% out.M1stim=M1stim;
% out.mM1stim=mM1stim;
out.mM1=mM1;
out.datadir=datadir;
out.freqs=freqs;
out.amps=amps;
out.durs=durs;
out.nreps=nreps;
out.numfreqs=numfreqs;
out.numamps=numamps;
out.numdurs=numdurs;
out.traces_to_keep=traces_to_keep;
out.Events=Events;
out.xlimits=xlimits;
out.ylimits=ylimits;
out.numframes=numframes;
out.samprate=samprate;
out.numcells=numcells;
out.nb=nb;
out.stimlog=stimlog;
out.readme={'M1 is freqs x amps x durs x reps x cells x frames', ...
    'mM1 is freqs x amps x durs x cells x frames'};


outfilename=sprintf('out2P.mat');
save(outfilename, 'out')
fprintf('\n saved to %s', outfilename)

