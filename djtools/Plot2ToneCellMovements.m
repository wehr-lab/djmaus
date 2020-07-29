function Plot2ToneCellMovements(varargin)
%currently just works with stimuli for Yashar, WN and laser pulses
dbstop if error
if nargin==0
    datadir=pwd;
    filename={};
    ylimits=[];
    xlimits=[];
end
if nargin==1
    datadir=varargin{1};
end

xlimits=[];
ylimits=[];
try
    xlimits=varargin{3};
end
if nargin==2
    datadir=varargin{1};
    filename=varargin{2};
    ylimits=[];
    xlimits=[];
end
if nargin==3
    datadir=varargin{1};
    filename=varargin{2};
    xlimits=varargin{3};
    ylimits=[];
end

if isempty(xlimits)
    %     xlimits=[-100 200];
    s=GetStimParams(datadir);
    durs=s.durs;
    dur=max(durs);
    SOAs=s.SOAs;
    SOA=max(SOAs);
    xlimits=[-.5*dur SOA+1.5*dur]; %default x limits for axis
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end

%params

binwidth=50;





load notebook.mat
%read messages
messagesfilename='messages.events';
[messages] = GetNetworkEvents(messagesfilename);
if exist(messagesfilename, 'file')~=2
    error('Could not find messages.events file. Is this the right data directory?')
end
%read digital Events
Eventsfilename='all_channels.events';
[all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
sampleRate=all_channels_info.header.sampleRate; %in Hz

if exist('Events.mat')
    load('Events.mat')
else
    [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
    save('Events.mat','Events')
end

if exist('StartAcquisitionSec.mat')
    load('StartAcquisitionSec.mat')
else
    %get Events and soundcard trigger timestamps
    [~, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
    save('StartAcquisitionSec.mat', 'StartAcquisitionSec')
end


if exist('moves.mat')
    load('moves.mat')
else
    [pulses_ms, puls_ms, motion_on_ms, motion_off_ms] = getBallMotion(pwd);
end
all_cells=dir('*.t');

filename=all_cells(1).name;
[p,f,ext]=fileparts(filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2num(ch{2});
clust=str2num(split{end});

%read MClust .t file
fprintf('\nreading MClust output file %s', filename)
spiketimes=read_MClust_output(filename)'/10000; %spiketimes now in seconds
%correct for OE start time, so that time starts at 0
spiketimes=spiketimes-StartAcquisitionSec;
totalnumspikes=length(spiketimes);
fprintf('\nsuccessfully loaded MClust spike data')

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
if exist(outfilename)
    load(outfilename)
else
    %Process2Tone_PSTH_single2(pwd, filename)
    Process2Tone_PSTH2
end
[a,b,c]=load_open_ephys_data('105_CH32.continuous');
length_of_recording=length(a)/30; %length of recording in ms


spiketimes=spiketimes*1000; %convert spiketimes into ms
%load outfile with spikes and trial info
IL=out.IL; %whether there are any interleaved laser trials
freqs=out.freqs; % only one -1 WN
probefreqs=out.probefreqs; % probe freq -1 WN or 0 silent sound
durs=out.durs; % 25 ms one dur for WN and silent sound

nreps=out.nreps;
numfreqs=out.numfreqs;
numprobefreqs=out.numprobefreqs;
numamps=out.numamps;
numdurs=out.numdurs;
numSOAs=out.numSOAs;
samprate=out.samprate; %in Hz

numLaserNumPulses=out.numLaserNumPulses;
numLaserISIs=out.numLaserISIs;
LaserNumPulses=out.LaserNumPulses;
LaserISIs=out.LaserISIs;

dindex=1;
LaserRecorded=out.LaserRecorded;
StimRecorded=out.StimRecorded;
M1Laser=out.M1Laser;
mM1Laser=out.mM1Laser;
M1Stim=out.M1Stim;
mM1Stim=out.mM1Stim;
st=spiketimes*1000; %are in ms
inRange=[];
load('LaserTrials.mat')
xlimits=out.xlimits;


all_channels_timestamps=all_channels_timestamps-all_channels_timestamps(1); %in s
total_duration=all_channels_timestamps(end)*1000; %in ms

move_trace=nan(round(total_duration),1);
load('puls_s.mat');
ids=puls_s(:,2);
p_on=puls_s(ids==1);
p_off=puls_s(ids==0);

if ids(1)==0
    next=1;
else
    next=0;
end

for i=1:length(p_on)-1
    on= p_on(i)*1000;
    off= p_off(i+next)*1000;
    move=(on-off)/(on-p_on(i+1)*1000);
    move_trace(round(on):round(p_on(i+1)*1000))=move;
    
end
idx=~isnan(move_trace);
move_trace=move_trace-(median(move_trace(idx)));
move_trace1=smoothdata(move_trace,'movmedian');
move_trace1=move_trace1*-1; %make forward positive
save('move_trace.mat','move_trace1')

figure;hold on
xst=1:binwidth:total_duration;
[N,X]=hist(spiketimes,xst); %ms
bar(N);
title('spiking');
figure;
%m=resample(move_trace1,  1 );
bar(move_trace1, 'b');
title('movements');


for i=1:length(Events)
    if strcmp(Events(i).type, '2tone') | ...
            strcmp(Events(i).type, 'silentsound') | strcmp(Events(i).type, 'whitenoise')%%%I should pull out silent sound processing into separate stanza
        if  isfield(Events(i), 'soundcard_trigger_timestamp_sec')
            pos=Events(i).soundcard_trigger_timestamp_sec;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec; %pos is in seconds
        end
        start=(pos+xlimits(1)*1e-3); %in seconds
        stop=(pos+xlimits(2)*1e-3);
        moves=[];
        ind=find(puls_s>start & puls_s<stop); % pulses in region
        pulses=puls_s(ind);
        id=ids(ind);
        pulsecount=length(pulses); % No. of pulses
        pulses=(pulses-pos)*1000;%covert to ms after tone onset
       
        moves=move_trace1(start*1000:stop*1000); %in ms
        
        if strcmp(Events(i).type, '2tone') || strcmp(Events(i).type, 'whitenoise')
            freq=-1;
            probefreq=Events(i).probefreq;
            SOA=Events(i).SOA;
            if strcmp(Events(i).type, 'whitenoise')
                probefreq=0;
                SOA=0;
            end
            findex= find(freqs==freq);
            pfindex= find(probefreqs==probefreq);
            SOAindex= find(SOAs==SOA);
            
            
            P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).pulses=pulses;
            P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).ids=ids;
            P1spikecounts(findex, pfindex, SOAindex,nreps(findex, pfindex, SOAindex))=pulsecount;
            P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).moves=moves;
        elseif strcmp(Events(i).type, 'silentsound')
            LaserNumPulse=Events(i).VarLasernumpulses;
            LaserISI=Events(i).VarLaserisi;
            freq=0;
            findex= find(freqs==freq);
            
            
            %reverse N of laser pulses to match pfindex 1=2tone, 2-single
            %should be lnindex=1 = 2tone if lnindex=2 single
            lnpindex= find(LaserNumPulses== LaserNumPulse);
            if lnpindex==1
                lnpindex=2;
                LaserISI=0;
            elseif lnpindex==2
                lnpindex=1;
            end
            liindex=  find(LaserISIs == LaserISI);
            
            P1(findex,lnpindex, liindex, nreps(findex, lnpindex, liindex)).pulses=pulses; % pulses times
            P1(findex,lnpindex, liindex, nreps(findex, lnpindex, liindex)).ids=ids; % pulses times
            P1spikecounts(findex, lnpindex, liindex,nreps(findex, lnpindex, liindex))=pulsecount; % No. of pulses
            P1(findex, lnpindex, liindex, nreps(findex, lnpindex, liindex)).moves=moves;
        end
        
    end
    
end


out.P1=P1;
out.P1spikecounts=P1spikecounts;


save (outfilename, 'out', '-v7.3')
%%
%plot single WN tials first, evoked
M1=out.M1;
P1=out.P1;
nreps=out.nreps;
figure;
findex=1;
pfindex=2; %single stimulus
SOAindex=1; % 0 SOA
% evoked and spont
evoked_xlimits=[0 200]; %processing xlimit(1) is subtracted from spike times and pulse times during processing, so 0 = onset of stimulus
spont_xlimits=[evoked_xlimits(2) 2000+abs(xlimits(1))];
ev_X=evoked_xlimits(1):binwidth:evoked_xlimits(2); %specify bin centers
sp_X=spont_xlimits(1):binwidth:spont_xlimits(2); %specify bin centers

moves=[];
all_eN=[];
all_emN=[];
all_sN=[];
all_smN=[];
for n=1:nreps(findex, pfindex, SOAindex)
    spiketimes=M1(findex, pfindex, SOAindex, n).spiketimes;
    pulses=P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).pulses;  %times of pulses during movements
    spiketimes1=spiketimes(spiketimes>evoked_xlimits(1) & spiketimes<evoked_xlimits(2)); % spikes in region evoked
    spiketimes2=spiketimes(spiketimes>spont_xlimits(1) & spiketimes<spont_xlimits(2)); % spikes in region spont213
    
    
    
    
    ind1=find(pulses>evoked_xlimits(1) & pulses<evoked_xlimits(2)); % pulses in region evoked
    ind2=find(pulses>spont_xlimits(1) & pulses<spont_xlimits(2)); % pulses in region spont
    pulses1=pulses(ind1); %evoked
    pulses2=pulses(ind2); %sp
    ids1=P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).ids(ind1); %evoked
    ids2=P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).ids(ind2); %sp
    
    pON1=pulses1(ids1==1); %evoked
    pOFF1=pulses1(ids1==0); %evoked
    pON2=pulses2(ids2==1); %sp
    pOFF2=pulses2(ids2==0); %sp
    
    
    for j=1:length(pON1)-1
        if pOFF1(1)>pON1(1)
            moves1(j)=abs((pON1(j)-pOFF1(j)))/(pON1(j+1)-pON1(j)); %on dur divide by total dur
        else
            moves1(j)=abs((pON1(j)-pOFF1(j+1)))/(pON1(j+1)-pON1(j));
        end
    end
    
    for j=1:length(pON2)-1
        if pOFF2(1)>pON2(1)
            moves2(j)=abs((pON2(j)-pOFF2(j)))/(pON2(j+1)-pON2(j)); %on dur divide by total dur
        else
            moves2(j)=abs((pON2(j)-pOFF2(j+1)))/(pON2(j+1)-pON2(j));
        end
    end
    moves1=moves1-.50; %0
    moves2=moves2-.50; %0
    [eN, ex]=hist(spiketimes1, ev_X);
    eN=1000*eN./binwidth; %normalize to spike rate in Hz
    [sN, sx]=hist(spiketimes2, sp_X);
    sN=1000*sN./binwidth; %normalize to spike rate in Hz
    
    move_binwinds=ceil(moves1)/length(eN);
    
    emN=[];
    for j=1:length(eN)
        if j<length(eN)
            emN(j)=mean(moves1(j:j+move_binwinds));
        else
            emN(j)=mean(moves1(j:end)); %whatever is left which might be shorter than mbinwindth
        end
    end
    
    move_binwinds=ceil(moves2)/length(sN);
    
    smN=[];
    for j=1:length(sN)
        if j<length(sN)
            smN(j)=mean(moves2(j:j+move_binwinds));
        else
            smN(j)=mean(moves2(j:end)); %whatever is left which might be shorter than mbinwindth
        end
    end
    
    
    %the end the result should be spiking in Hz bins with movement means
    %for evoked and spont regions to be able to calculate corr between them,
    %matched in time and corrected for sampling differences
    figure(2)
    hold on
    plot(eN,emN, 'r.')
    figure(3)
    plot(sN,smN, 'k.')

   all_eN=[all_eN eN];
   all_emN=[all_emN emN];
   all_sN=[all_sN sN];
   all_smN=[all_smN smN];
    
end
%save('testfiringandmovements.mat','all_eN','all_emN','all_sN','all_smN')
figure;

%plot single laser pulsetials first, evoked
M1=out.M1;
P1=out.P1;
nreps=out.nreps;
figure;
findex=2;
pfindex=2; %single stimulus
SOAindex=1; % 0 SOA
% evoked and spont
evoked_xlimits=[0 200]; %processing xlimit(1) is subtracted from spike times and pulse times during processing, so 0 = onset of stimulus
spont_xlimits=[evoked_xlimits(2) 2000+abs(xlimits(1))];
ev_X=evoked_xlimits(1):binwidth:evoked_xlimits(2); %specify bin centers
sp_X=spont_xlimits(1):binwidth:spont_xlimits(2); %specify bin centers

moves=[];
all_eN=[];
all_emN=[];
all_sN=[];
all_smN=[];
for n=1:nreps(findex, pfindex, SOAindex)
    spiketimes=M1(findex, pfindex, SOAindex, n).spiketimes;
    pulses=P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).pulses;  %times of pulses during movements
    spiketimes1=spiketimes(spiketimes>evoked_xlimits(1) & spiketimes<evoked_xlimits(2)); % spikes in region evoked
    spiketimes2=spiketimes(spiketimes>spont_xlimits(1) & spiketimes<spont_xlimits(2)); % spikes in region spont213
    
    
    
    
    ind1=find(pulses>evoked_xlimits(1) & pulses<evoked_xlimits(2)); % pulses in region evoked
    ind2=find(pulses>spont_xlimits(1) & pulses<spont_xlimits(2)); % pulses in region spont
    pulses1=pulses(ind1); %evoked
    pulses2=pulses(ind2); %sp
    ids1=P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).ids(ind1); %evoked
    ids2=P1(findex, pfindex, SOAindex, nreps(findex, pfindex, SOAindex)).ids(ind2); %sp
    
    pON1=pulses1(ids1==1); %evoked
    pOFF1=pulses1(ids1==0); %evoked
    pON2=pulses2(ids2==1); %sp
    pOFF2=pulses2(ids2==0); %sp
    
    
    for j=1:length(pON1)-1
        if pOFF1(1)>pON1(1)
            moves1(j)=abs((pON1(j)-pOFF1(j)))/(pON1(j+1)-pON1(j)); %on dur divide by total dur
        else
            moves1(j)=abs((pON1(j)-pOFF1(j+1)))/(pON1(j+1)-pON1(j));
        end
    end
    
    for j=1:length(pON2)-1
        if pOFF2(1)>pON2(1)
            moves2(j)=abs((pON2(j)-pOFF2(j)))/(pON2(j+1)-pON2(j)); %on dur divide by total dur
        else
            moves2(j)=abs((pON2(j)-pOFF2(j+1)))/(pON2(j+1)-pON2(j));
        end
    end
    moves1=moves1-.50; %0
    moves2=moves2-.50; %0
    [eN, ex]=hist(spiketimes1, ev_X);
    eN=1000*eN./binwidth; %normalize to spike rate in Hz
    [sN, sx]=hist(spiketimes2, sp_X);
    sN=1000*sN./binwidth; %normalize to spike rate in Hz
    
    move_binwinds=ceil(moves1)/length(eN);
    
    emN=[];
    for j=1:length(eN)
        if j<length(eN)
            emN(j)=mean(moves1(j:j+move_binwinds));
        else
            emN(j)=mean(moves1(j:end)); %whatever is left which might be shorter than mbinwindth
        end
    end
    
    move_binwinds=ceil(moves2)/length(sN);
    
    smN=[];
    for j=1:length(sN)
        if j<length(sN)
            smN(j)=mean(moves2(j:j+move_binwinds));
        else
            smN(j)=mean(moves2(j:end)); %whatever is left which might be shorter than mbinwindth
        end
    end
    
    
    %the end the result should be spiking in Hz bins with movement means
    %for evoked and spont regions to be able to calculate corr between them,
    %matched in time and corrected for sampling differences
    figure(2)
    hold on
    plot(eN,emN, 'r.')
    figure(3)
    plot(sN,smN, 'k.')

   all_eN=[all_eN eN];
   all_emN=[all_emN emN];
   all_sN=[all_sN sN];
   all_smN=[all_smN smN];
    
end
figure;




