function ProcessTC_2P(varargin)
%processes dF/F 2P tuning curve data from taskcontrol stimulus delivery and scanbox scope data
%sorts dF/F  data into a response matrix.
%
%usage: ProcessTC_2P([datapath], [xlimits], [ylimits])
%datapath is the directory containing the .sbx file, the taskcontrol .h5 file, and a suite2p output folder 
%datapath defaults to current directory
%saves output in an outfile
%
%I'm switching to including all cells, since it makes no sense to have
%thousands of outfiles.
%
%notes: whitenoise plotted as freq=-1 kHz, silent sound as -2 kHz

if nargin==0
    datadir=pwd;
    xlimits=[-1000 5000]; %default x limits for axis
    ylimits=[-.1 .2];
elseif nargin==1
    datadir=varargin{1};
    xlimits=[-1000 5000]; %default x limits for axis
    ylimits=[-.1 .2];
elseif nargin==2
    datadir=varargin{1};
    xlimits=varargin{2};
    ylimits=[-.1 .2];
elseif nargin==3
    datadir=varargin{1};
    xlimits=varargin{2};
    ylimits=varargin{3};
else
    error('wrong number of arguments');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if isempty(xlimits)
        xlimits=[-1000 5000]; %default x limits for axis
end
if isempty(ylimits)
    ylimits=[-.1 .2];
end

%assume this directory structure:
cd suite2p
cd plane0
iscell=0; %initialize because iscell is a built-in function
load Fall
numframes=size(F, 2);
sampleRate=ops.fs; %2P scope framerate in Hz
fprintf('\n%d total frames = %.1f s', numframes, numframes/sampleRate)
% sort rows such that best cells (highest cell likelihood) are at the top
%exclude non-cells 
%[iscell_sorted, I]=sortrows(iscell, 2, 'descend');
iscell_threshold=.7;
f=F(find(iscell(:,2)>iscell_threshold), :);
iscell2=find(iscell(:,2)>iscell_threshold);
%here we could impose a different probability threshold for iscell than whatever was set in suite2p, if we wanted
% f=F(I(1:numcells),:);
f0_prctile=50;
% f0=prctile(f', f0_prctile);
% keep=find(f0>0);
% f=f(keep,:);
f0=prctile(f', f0_prctile); %not sure why the extra iteration, must have been a specific case, should have no effect
f0=repmat(f0', 1, numframes);
dff_unsorted=(f-f0)./f0;
numcells=size(dff_unsorted, 1);
fprintf('\n%d cells (using iscell_threshold=%.1f)', numcells, iscell_threshold)

%another way to compute dff is to compute it for each trial, i.e. use the
%pre-stimulus baseline as f0 for each trial. That's how Xu does it (for SOM cells).
% For PNs, Xu did a rolling baseline removal, by computing the 8th %tile of
% the distribution of F in a 20-s window around each frame, and subtracting
% this from F on that frame. Then for each ROI, F0 is computed as the index of the peak
% of the histogram of F, and dff is computed as usual, (F-F0)/F0

%sort by pca
[pcs, score, latent]=pca(dff_unsorted);
[~, I]=sort(score(:,1));
dff=(dff_unsorted(I,:));
fprintf('\nsorted cells by pca (without regard to stimuli)')
%(if you want to see pca of stim-aligned responses, we have to align to
%stimuli first)

%load stim framestamps from scanbox mat file 
cd(datadir)
d=dir('*.sbx');
if length(d)<1 error('no sbx file found')
elseif length(d)>1
    warning('multiple sbx files found, using the first one')
end
sbxfilename=replace(d(1).name, '.sbx', '.mat');
load(sbxfilename)
% frames = info.frame(1:2:(end-2)); unclear whether we need to discard last 2 triggers or not
frames = info.frame(1:2:end);
numtriggers=length(frames);

%load stimlog from taskcontrol h5 file
d=dir('*.h5');
if length(d)<1 error('no taskcontrol h5 file found')
elseif length(d)>1
    warning('multiple taskcontrol h5 files found, using the first one')
end
taskcontrolfilename=d(1).name;
allfreqs=h5read(taskcontrolfilename, '/resultsData/currentFreq/');
allamps=h5read(taskcontrolfilename, '/resultsData/currentIntensity/');
allisis=h5read(taskcontrolfilename, '/resultsData/isi/');
alldurs=h5read(taskcontrolfilename, '/resultsData/stimDur/');
numstim=length(allfreqs);
fprintf('\n%d sound triggers \n%d logged stimuli', numtriggers, numstim)
if numtriggers==numstim
    fprintf('\t good, they match')
else
    fprintf('\nproblem: numtriggers does not match numstim...')
    fprintf('\n\twill use the first %d triggers, and discard the last %d trigger', numstim, numtriggers-numstim)
end

if numtriggers>=numstim
    for i=1:numstim
        Events(i).type='tone';
        Events(i).freq=allfreqs(i);
        Events(i).amp=allamps(i);
        Events(i).dur=alldurs(i);
        Events(i).isi=allisis(i);
        Events(i).frame=frames(i);
    end
else
    for i=1:numtriggers
        Events(i).type='tone';
        Events(i).freq=allfreqs(i);
        Events(i).amp=allamps(i);
        Events(i).dur=alldurs(i);
        Events(i).isi=allisis(i);
        Events(i).frame=frames(i);
    end
end

        LaserRecorded=0;
        StimRecorded=0;

fprintf('\ncomputing tuning curve...');

samprate=sampleRate;

% %get freqs/amps
% j=0;
% for i=1:length(Events)
%     if strcmp(Events(i).type, '2tone') |strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'silentsound') ...
%             |strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, 'whitenoise')| strcmp(Events(i).type, 'grating')
%         j=j+1;
%         alldurs(j)=Events(i).dur;
%         allisis(j)=Events(i).isi;
%         if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, '2tone')
%             allamps(j)=Events(i).amp;
%             allfreqs(j)=Events(i).freq;
%         elseif strcmp(Events(i).type, 'whitenoise')
%             allamps(j)=Events(i).amp;
%             allfreqs(j)=-1000;
%         elseif strcmp(Events(i).type, 'silentsound')
%             allfreqs(j)=-2000;
%             allamps(j)=nan; %flagging silent sound amp as nan
%         elseif strcmp(Events(i).type, 'fmtone')
%             allamps(j)=Events(i).amp;
%             allfreqs(j)=Events(i).param.carrier_frequency;
% 
%         end
%     end
% end
% allamps=allamps(~isnan(allamps)); %strip out nans from silent sound
freqs=unique(allfreqs);
amps=unique(allamps);
durs=unique(alldurs);
isis=unique(allisis);
numfreqs=length(freqs);
numamps=length(amps);
numdurs=length(durs);
numisis=length(isis);

fprintf('\nfreqs in kHz:');fprintf(' %.1f', freqs/1000)

M1=[];
M1Stim=[];

nreps=zeros(numfreqs, numamps, numdurs);




samprate=sampleRate;



%extract the traces into a big matrix M
j=0;
for i=1:length(Events)
    if strcmp(Events(i).type, 'tone') | strcmp(Events(i).type, 'whitenoise') |  strcmp(Events(i).type, 'silentsound') | ...
            strcmp(Events(i).type, 'fmtone') | strcmp(Events(i).type, '2tone')| strcmp(Events(i).type, 'grating')
        if  isfield(Events(i), 'frame')
            pos=Events(i).frame;
        else
            error('???')
            pos=Events(i).soundcard_trigger_timestamp_sec*samprate; %pos is in samples (2P frames)
        end
        start=round(pos+xlimits(1)*1e-3*samprate);
        win=round(diff(xlimits)*1e-3*samprate);
        stop=start+win;
        region=start:stop;


        if isempty(find(region<1)) & (region(end)< length(dff)) %disallow negative or zero start times & make sure there's enough data after the stimulus
            switch Events(i).type
                case {'tone', '2tone'}
                    freq=Events(i).freq;
                    amp=Events(i).amp;
                case 'whitenoise'
                    freq=-1000;
                    amp=Events(i).amp;
                case 'silentsound'
                    freq=-2000;
                    amp=min(amps); %put silentsound in it's own column (freq=-2) in the lowest row
            end
            
            
            dur=Events(i).dur;
            findex= find(freqs==freq);
            aindex= find(amps==amp);
            dindex= find(durs==dur);
            nreps(findex, aindex, dindex)=nreps(findex, aindex, dindex)+1;
            M1dff(findex,aindex,dindex, nreps(findex, aindex, dindex),:,:)=dff(:,region);
            M1f(findex,aindex,dindex, nreps(findex, aindex, dindex),:,:)=f(:,region);
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
                    mM1f(findex, aindex, dindex,:,:)=mean(M1f(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:,:), 4);
                    mM1dff(findex, aindex, dindex,:,:)=mean(M1dff(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:,:), 4);
                    %mM1 is freqs x amps x durs x cells x frames

                   % mM1stim(findex, aindex, dindex,:)=mean(M1stim(findex, aindex, dindex, 1:nreps(findex, aindex, dindex),:), 4);
                else %no reps for this stim, since rep=0
                    mM1f(findex, aindex, dindex,:,:)=zeros(numcells, length(region));
                    mM1dff(findex, aindex, dindex,:,:)=zeros(numcells, length(region));
                    %mM1stim(findex, aindex, dindex,:)=zeros(size(region));
                end
           
                
            end
        end
    end
end








%assign outputs
out.dff=dff;
out.f0_prctile=f0_prctile;
out.iscell_threshold=iscell_threshold;
out.M1f=M1f;
out.M1dff=M1dff;
% out.M1stim=M1stim;
% out.mM1stim=mM1stim;
out.mM1f=mM1f;
out.mM1dff=mM1dff;
out.datadir=datadir;
out.sbxfilename=sbxfilename;
out.taskcontrolfilename=taskcontrolfilename;
out.freqs=freqs;
out.amps=amps;
out.durs=durs;
out.isis=isis;
out.nreps=nreps;
out.numfreqs=numfreqs;
out.numamps=numamps;
out.numdurs=numdurs;
out.numisis=numisis;
out.traces_to_keep=traces_to_keep;
out.Events=Events;
out.xlimits=xlimits;
out.ylimits=ylimits;
out.numframes=numframes;
out.samprate=samprate;
out.numcells=numcells;
out.iscell=iscell2;
out.stat=stat(iscell2);
out.readme={'M1 is freqs x amps x durs x reps x cells x frames', ...
    'mM1 is freqs x amps x durs x cells x frames'};


outfilename=sprintf('out2P.mat');
save(outfilename, 'out')
fprintf('\n saved to %s', outfilename)

