function PlotPDR(varargin)
% Plots PDR alone or with data if given
% (datadir, video_file, xlimtis, ylimits)
% ira 08.31.17
djPrefs;
global pref
samprate=30e3;
if nargin==0
    fprintf('\nno input\n')
    return
else
    datadir=varargin{1};
end

if nargin==1
    video_filename=[];
    if isempty(video_filename)
        [video_filename, pathname, ~] = uigetfile('*.mp4', 'Select video file');
    end
elseif nargin==2
    video_filename=varargin{2};
    xlimits=[]; %x limits for axis
    ylimits=[];
elseif nargin==3
    video_filename=varargin{2};
    xlimits=varargin{3};
    ylimits=[];
elseif nargin==4
    video_filename=varargin{2};
    xlimits=varargin{3};
    ylimits=varargin{4};
end

if isempty(video_filename)
    [video_filename, pathname, ~] = uigetfile('*.mp4', 'Select video file');
end
prompt=('Would you like to process data with PDR (Y/N): ');
answer=input(prompt,'s');
if strcmp(answer,'Y')
    try
        [outfilename, pathname, ~] = uigetfile('*.mat', 'please select the outfile');
    catch
        fprintf('\noutfile was not found\n')
    end
    out=load(outfilename);
    xlimits=out.xlimits;
    ylimits=out.ylimits;
end


ts_data_file=sprintf('%s.txt', video_filename(1:19)); %process Pi timestamps
if exist('pupil_data.mat')
    sprintf('\nFound PDR data\n')
    load('pupil_data.mat')
else
    Pupil_Dilation(datadir, video_filename);
    pupil_data = ReadPiTimeStamps(ts_data_file);
end
if isempty(xlimits)
    xlimits=[-100 400];
end


if strcmp(answer,'N')
    load('Events.mat');
    events_unknown=0;
    for i=1:length(Events)
        if strcmp(Events(i).type, 'whitenoise')
            freqs(i)=-1000;
            amps(i)=Events(i).amplitude;
            durs(i)=Events(i).duration;
            lasers(i)=Events(i).laser;
        elseif strcmp(Events(i).type, 'tone')
            freqs(i)=Events(i).frequency;
            amps(i)=Events(i).amplitude;
            durs(i)=Events(i).duration;
            gap_durs(i)=[];
            lasers(i)=Events(i).laser;
        elseif strcmp(Events(i).type, 'silentsound')
            freqs(i)=-2000;
            amps(i)=-1000;
            durs(i)=Events(i).duration;
            lasers(i)=Events(i).laser;
        elseif strcmp(Events(i).type, 'GPIAS') %this might be wrong
            freqs(i)=-1;
            amps(i)=Events(i).amplitude;
            gap_durs(i)=Events(i).gap_duration;
            lasers(i)=Events(i).laser;
        else
            fprintf('\ncould not recognize the stimulus %s\n', Events(i).Type)
            events_unknown=1;
        end
    end
end
numamps=length(unique(amps));
numfreqs=length(unique(freqs));
freqs=unique(freqs);
amps=unique(amps);
try
    numgapdurs=length(unique(gap_durs));
    gap_durs=unique(gap_durs);
catch
    numgapdurs=[];
end
numdurs=length(unique(durs));
durs=unique(durs);


ts=pupil_data.ts;
smooth_data=pupil_data.smooth_data;


if events_unknown==1
    fprintf('\n Could not identify all events, will plot PDR assuming one type of event\n')
    nrepsON=[];
    nrepsOFF=[];
    for i=1:length(Events)
        pos=ts(i+1).ts_zero_samprate;
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        if Events(i).laser==1
            nrepsON=nrepsON+1;
            pM1ON(nrepsON)=smooth_data(region);
        elseif Events(i).laser==0
            nrepsOFF=nrepsOFF+1;
            pM1OFF(nrepsOFF)=smooth_data(region);
        end
    end
end

if events_unknown==0 && strcmp(answer,'N') && ~exist('gap_durs')
    fprintf('\n Found all events, plotting PDR only\n')
    nrepsON=zeros(numfreqs, numamps, numdurs);
    nrepsOFF=zeros(numfreqs, numamps, numdurs);
    for i=1:length(Events)
        pos=ts(i+1).ts_zero_samprate;
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        switch Events(i).type
            case {'tone', '2tone'}
                freq=Events(i).frequency;
                amp=Events(i).amplitude;
            case 'fmtone'
                freq=Events(i).carrier_frequency;%
                amp=Events(i).amplitude;
            case 'whitenoise'
                freq=-1000;
                amp=Events(i).amplitude;
            case 'silentsound'
                freq=-2000;
                amp=min(amps); %put silentsound in it's own column (freq=-2) in the lowest row
            case 'grating'
                amp=Events(i).spatialfrequency;
                freq=Events(i).angle*1000;
        end
        dur=Events(i).duration;
        findex= find(freqs==freq);
        aindex= find(amps==amp);
        dindex= find(durs==dur);
        if Events(i).laser==1
            nrepsON(findex, aindex, dindex)=nrepsON(findex, aindex, dindex)+1;
            pM1ON(findex,aindex,dindex, nrepsON(findex, aindex, dindex),:)=smooth_data(region);
        elseif Events(i).laser==0
            nrepsOFF(findex, aindex, dindex)=nrepsOFF(findex, aindex, dindex)+1;
            pM1OFF(findex,aindex,dindex, nrepsOFF(findex, aindex, dindex),:)=smooth_data(region);
        end
    end
elseif events_unknown==0 && strcmp(answer,'N') && ~isempty(gap_durs)
    j=0;
    if strcmp(Events(i).type, 'GPIAS')
        j=j+1;
        allsoas(j)=Events(i).soa;
        allgapdurs(j)=Events(i).gapdur;
        allgapdelays(j)=Events(i).gapdelay;
        allpulseamps(j)=Events(i).pulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allnoiseamps(j)=Events(i).amplitude;
    elseif strcmp(Events(i).type, 'gapinnoise')
        allsoas(j)=Events(i).soa;
        allgapdurs(j)=Events(i).gapdur;
        allgapdelays(j)=Events(i).gapdelay;
        allpulseamps(j)=Events(i).pulseamp;
        allpulsedurs(j)=Events(i).pulsedur;
        allnoiseamps(j)=Events(i).amplitude;
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
    for i=1:length(Events)
        if strcmp(Events(i).type, 'GPIAS') | strcmp(Events(i).type, 'gapinnoise')
            pos=ts(i+1).ts_zero_samprate; %pos is in seconds
            start=pos + (gapdelay/1000)*samprate +(xlimits(1)/1000)*samprate; %start is in samples
            stop=pos+ (gapdelay/1000)*samprate + (xlimits(2)/1000)*samprate; %stop is in samples
            region=round(start)+1:round(stop);
            if start>0 %(disallow negative or zero start times)
                gapdur=Events(i).gapdur;
                gdindex= find(gapdur==gapdurs);
                pulseamp=Events(i).pulseamp;
                paindex= find(pulseamp==pulseamps);
                if Events(i).laser
                    nrepsON(gdindex,paindex)=nrepsON(gdindex,paindex)+1;
                    pM1ON(gdindex,paindex, nrepsON(gdindex,paindex))=smooth_data(region); % Spike times
                elseif Events(i).laser
                    nrepsOFF(gdindex,paindex)=nrepsOFF(gdindex,paindex)+1;
                    pM1OFF(gdindex,paindex, nrepsOFF(gdindex,paindex))=smooth_data(region);
                end
            end
            
        end
        
    end %events
end

%mean
for aindex=1:numamps
    for findex=1:numfreqs
        for dindex=1:numdurs
            if nrepsON(findex, aindex, dindex)>0
                mpM1ON(findex, aindex, dindex,:)=mean(pM1ON(findex, aindex, dindex, 1:nrepsON(findex, aindex, dindex),:), 4);
            else %no reps for this stim, since rep=0
                mpM1ON(findex, aindex, dindex,:)=zeros(size(region));
            end
            if nrepsOFF(findex, aindex, dindex)>0
                mpM1OFF(findex, aindex, dindex,:)=mean(pM1OFF(findex, aindex, dindex, 1:nrepsOFF(findex, aindex, dindex),:), 4);
            else %no reps for this stim, since rep=0
                mpM1OFF(findex, aindex, dindex,:)=zeros(size(region));
            end
        end
    end
end


%
if isempty(ylimits)
    for aindex=[numamps:-1:1]
        for findex=1:numfreqs
            traceON=mpM1ON(findex, aindex, dindex,:);
            ylimitsON(1)=min(ylimitsON(1), min(traceON));
            ylimitsON(2)=max(ylimitsON(2), max(traceON));

            traceOFF=mpM1ON(findex, aindex, dindex,:);
            ylimitsOFF(1)=min(ylimitsOFF(1), min(traceOFF));
            ylimitsOFF(2)=max(ylimitsOFF(2), max(traceOFF));
        end
    end
    ylimits(1)=min(ylimitsON(1), ylimitsOFF(1));
    ylimits(2)=max(ylimitsOFF(2), ylimitsON(2));
end

%plot pupil data
for dindex=1:numdurs
        figure
        p=0;
        subplot1(numamps,numfreqs)
        for aindex=numamps:-1:1
            for findex=1:numfreqs
                p=p+1;
                subplot1(p)
                axis off
                trace1=squeeze(mpM1OFF(findex, aindex, dindex, :));
                trace1=trace1 -mean(trace1(1:100));
                t=1:length(trace1);
                t=1000*t/out.samprate; %convert to ms
                t=t+xlimits(1); %correct for xlim in original processing call
                line([0 0+durs(dindex)], ylimits(1)+[0 0], 'color', 'm', 'linewidth', 5)
            end
        end
end


if events_unknown==0 && strcmp(answer,'Y')
    fprintf('\n Found all events, plotting PDR and data\n')
    nrepsON=0;
    nrepsOFF=0;
    for i=1:length(Events)
        pos=ts(i+1).ts_zero_samprate;
        start=round(pos+xlimits(1)*1e-3*samprate);
        stop=round(pos+xlimits(2)*1e-3*samprate)-1;
        region=start:stop;
        if Events(i).laser==1
            nrepsON=nrepsON+1;
            pM1ON(nrepsON)=smooth_data(region);
        elseif Events(i).laser==0
            nrepsOFF=nrepsOFF+1;
            pM1OFF(nrepsOFF)=smooth_data(region);
        end
    end
end
