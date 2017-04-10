function analTemp3(varargin)

% usage: analTemp3(datapath, t_filename, [xlimits],[ylimits], [binwidth], Type)
% (xlimits, ylimits, binwidth are optional),
% argin Type is useful as: 'WN w/o laser', 'GPIAS w/o burst', 'GPIAS w/ burst'
%
% reprocesses data if outfile is not found or xlimits are wrong;
% analTemp1 is outdated  (3/23/17)
% differs from analTemp2 by using global CELL_STATS instead of PV_STATS

global CELL_STATS
global ProcessedCellNum

IDnumber = CELL_STATS(ProcessedCellNum).IDnumber;

% rasters=1;
force_reprocess=0;

if nargin==0
    fprintf('\nno input');
    return;
end
datadir=varargin{1};
t_filename=varargin{2};

try
    xlimits=varargin{3};
catch
    xlimits=[];
end
try
    ylimits=varargin{4};
catch
    ylimits=[];
end
try
    binwidth=varargin{5};
catch
    binwidth=5;
end
try
    Type = varargin{6};
catch
    Type = 'unknown';
end

flag.plot = 0;

flag.tailstr = 'both';
%flag.tailstr = 'right';

[~,f]=fileparts(t_filename);
split=strsplit(f, '_');
ch=strsplit(split{1}, 'ch');
channel=str2double(ch{2});
clust=str2double(split{end});

outfilename=sprintf('outPSTH_ch%dc%d.mat',channel, clust);
cd(datadir)

if exist(outfilename,'file') & ~force_reprocess
    load(outfilename)
else
    switch Type
        case 'WN w/o laser'
            ProcessTC_PSTH_single(datadir,  t_filename, xlimits, ylimits);
        case {'GPIAS w/o burst', 'GPIAS w/ burst'}
            ProcessGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits);
        otherwise
            PlottingFunction=GetPlottingFunction;
            eval([PlottingFunction '(datadir,  t_filename, xlimits, ylimits);']);
    end
    load(outfilename);
end

%if xlimits are requested but don't match those in outfile, force preprocess
if ~isempty(xlimits)
    if out.xlimits(1)>xlimits(1) | out.xlimits(2)<xlimits(2) %xlimits in outfile are too narrow, so reprocess
        switch Type
            case 'WN w/o laser'
                ProcessTC_PSTH_single(datadir,  t_filename, xlimits, ylimits);
            case {'GPIAS w/o burst', 'GPIAS w/ burst'}
                ProcessGPIAS_PSTH_single(datadir,  t_filename, xlimits, ylimits);
            otherwise
                PlottingFunction=GetPlottingFunction;
                eval([PlottingFunction '(datadir,  t_filename, xlimits, ylimits);']);
        end
        load(outfilename);
    end
end

if isempty(xlimits) xlimits=out.xlimits;end

% determine type of test:
% if strfind(out.stimlog(1).protocol_description,'tuning curve, WN only')
%     if strfind(out.stimlog(1).protocol_description,'IL')
%         Type = 'WN with laser';
%     else
%         Type = 'WN w/o laser';
%     end
% elseif strfind(out.stimlog(1).protocol_description,'GPIAS')
%     if strfind(out.stimlog(1).protocol_description,'pulsedur0ms')
%         Type = 'GPIAS w/o burst';
%     elseif strfind(out.stimlog(1).protocol_description,'pulsedur25ms')
%         Type = 'GPIAS w/ burst';
%     else
%         Type = 'unknown';
%     end
% else
%     Type = 'unknown';
% end

IL=out.IL;                  % any interleaved laser trials?
switch Type
    case 'WN w/o laser'
        freqs=out.freqs;
        amps=out.amps;
        durs=out.durs;
        nreps=out.nreps;
        numfreqs=out.numfreqs;
        numamps=out.numamps;
        numdurs=out.numdurs;
        
        M_LaserNumPulses=out.M_LaserNumPulses;
        LaserNumPulses=out.LaserNumPulses;
        LaserISI=out.LaserISI;
    case {'GPIAS w/o burst', 'GPIAS w/ burst'}
        gapdurs = out.gapdurs;
        gapdelay = out.gapdelay;
        numpulseamps = out.numpulseamps;
        numgapdurs = out.numgapdurs;
        pulseamps = out.pulseamps;
end

samprate=out.samprate; %in Hz
mM1ON=out.mM1ON;
mM1OFF=out.mM1OFF;
M1ON=out.M1ON;
M1OFF=out.M1OFF;
nrepsON=out.nrepsON;
nrepsOFF=out.nrepsOFF;

% LaserRecorded=out.LaserRecorded;
% StimRecorded=out.StimRecorded;
% M_LaserStart=out.M_LaserStart;
% M_LaserWidth=out.M_LaserWidth;
% M_LaserISI=out.M_LaserISI;
% LaserStart=out.LaserStart;
% LaserWidth=out.LaserWidth;

% M1ONStim=out.M1ONStim;
% M1ONLaser=out.M1ONLaser; % a crash here means this is an obsolete outfile. Set force_reprocess=1 up at the top of this mfile. (Don't forget to reset it to 0 when you're done)
% mM1ONStim=out.mM1ONStim;
% mM1ONLaser=out.mM1ONLaser;
% M1OFFStim=out.M1OFFStim;
% M1OFFLaser=out.M1OFFLaser;
% mM1OFFStim=out.mM1OFFStim;
% mM1OFFLaser=out.mM1OFFLaser;
% mM1OFFspikecount=out.mM1OFFspikecount;
% sM1OFFspikecount=out.sM1OFFspikecount;
% semM1OFFspikecount=out.semM1OFFspikecount;
% mM1spontOFF=out.mM1spontOFF;
% sM1spontOFF=out.sM1spontOFF;
% semM1spontOFF=out.semM1spontOFF;
% spiketimes=out.spiketimes;

% % % M1ON(findex,aindex,dindex, nrepsON).spiketimes
% % % mM1ON(findex,aindex,dindex).spiketimes
% % % mM1ONspikecount(findex,aindex,dindex)

% Nclusters=out.Nclusters;
% tetrode = out.tetrode;
% channel=out.channel;
% clust=out.cluster; %there are some redundant names here
% cell  = out.cell;
datadir=out.datadir;
% nb = out.nb;
% stimlog = out.stimlog;

try
    nreps_ssON=out.nreps_ssON;
    nreps_ssOFF=out.nreps_ssOFF;
    SilentSoundON=out.SilentSoundON;
    SilentSoundOFF=out.SilentSoundOFF;
    SilentSoundONspikecount=out.SilentSoundONspikecount;
    SilentSoundOFFspikecount=out.SilentSoundOFFspikecount;
    mSilentSoundON=out.mSilentSoundON;
    mSilentSoundOFF=out.mSilentSoundOFF;
    
    SilentSoundONStim=out.SilentSoundONStim;
    SilentSoundOFFStim=out.SilentSoundOFFStim;
    SilentSoundONLaser=out.SilentSoundONLaser;
    SilentSoundOFFLaser=out.SilentSoundOFFLaser;
    SS = 1;         % a flag indicating SilentSounds were played
    
catch
    SS = 0;
end

switch Type
    %%%%%%%%%%%%%%%
    case 'WN w/o laser'
        % WN duration tuning
        analysisWindow = 50;        % ms for analysis
        nrepsOFF = min(nrepsOFF);
        if ~SS  % no silent sound recorded
            nreps_ssOFF = nrepsOFF;
        end
        % set up arrays for #spikes/rep
        ONSETspikesPerRepSS = zeros(1,nreps_ssOFF);
        OFFSETspikesPerRepSS = zeros(numdurs,nreps_ssOFF);
        ONSETspikesPerRep = zeros(numdurs,nrepsOFF);
        OFFSETspikesPerRep = zeros(numdurs,nrepsOFF);
        % ONSET2 tests during entire duration of sound;
        % OFFSET2 tests during 2nd analysis window after offset
        ONSET2spikesPerRepSS = zeros(1,nreps_ssOFF);
        OFFSET2spikesPerRepSS = zeros(numdurs,nreps_ssOFF);
        ONSET2spikesPerRep = zeros(sum(durs>=analysisWindow*2),nrepsOFF);
        OFFSET2spikesPerRep = zeros(numdurs,nrepsOFF);
        
        if SS
            for irep = 1:nreps_ssOFF
                temp = SilentSoundOFF(irep).spiketimes;
                % spiketimes are indexed to sound onset
                ONSETspikesPerRepSS(irep) = length(find(temp > 0 & temp<=analysisWindow));
            end
        end
        % set up arrays for hypothesis-test and P-values
        H_offset = nan(1,numdurs); H_onset = H_offset;
        P_offset = H_offset; P_onset = P_offset;
        
        H_offset2 = nan(1,numdurs); H_onset2 = H_offset2;
        P_offset2 = H_offset2; P_onset2 = P_offset2;
        
        for idur = 1:numdurs
            %%%% OFFSETs (how many spikes in AnalysisWindow)
            % get control data from Silent Sounds or (if no SS) from before sound
            if SS
                for irep = 1:nreps_ssOFF
                    temp = SilentSoundOFF(irep).spiketimes;
                    OFFSETspikesPerRepSS(idur,irep) = length(find(temp>durs(idur)+0*analysisWindow & temp<=durs(idur)+analysisWindow));
                    OFFSET2spikesPerRepSS(idur,irep) = length(find(temp>durs(idur)+1*analysisWindow & temp<=durs(idur)+2*analysisWindow));
                end
            else
                for irep = 1:nrepsOFF
                    temp = M1OFF(1,1,idur,irep).spiketimes;
                    OFFSETspikesPerRepSS(idur,irep) = length(find(temp>-100-0*analysisWindow & temp<=-100-analysisWindow));
                    OFFSET2spikesPerRepSS(idur,irep) = length(find(temp>-100-1*analysisWindow & temp<=-100-2*analysisWindow));
                end
            end
            % get test data
            for irep = 1:nrepsOFF
                temp = M1OFF(1,1,idur,irep).spiketimes;
                OFFSETspikesPerRep(idur,irep) = length(find(temp>durs(idur)+0 & temp<=durs(idur)+analysisWindow));
                OFFSET2spikesPerRep(idur,irep) = length(find(temp>durs(idur)+analysisWindow & temp<=durs(idur)+2*analysisWindow));
            end
            % test for significance
            ntemp = min([size(OFFSETspikesPerRep,2) size(OFFSETspikesPerRepSS,2)]);
            [H_offset(idur),P_offset(idur)] = ttest(OFFSETspikesPerRep(idur,1:ntemp),OFFSETspikesPerRepSS(idur,1:ntemp),'tail',flag.tailstr);
            
            ntemp = min([size(OFFSET2spikesPerRep,2) size(OFFSET2spikesPerRepSS,2)]);
            [H_offset2(idur),P_offset2(idur)] = ttest(OFFSET2spikesPerRep(idur,1:ntemp),OFFSET2spikesPerRepSS(idur,1:ntemp),'tail',flag.tailstr);
            %%%% ONSETs
            % if not calculated above, calculate control values here for each duration
            if SS
                for irep = 1:nreps_ssOFF
                    temp = SilentSoundOFF(irep).spiketimes;
                    ONSET2spikesPerRepSS(idur,irep) = length(find(temp>0 & temp<=durs(idur)));
                end
            else
                for irep = 1:nrepsOFF
                    temp = M1OFF(1,1,idur,irep).spiketimes;
                    ONSETspikesPerRepSS(irep) = length(find(temp>-100-0*analysisWindow & temp<=-100-analysisWindow));
                    ONSET2spikesPerRepSS(irep) = length(find(temp>-durs(idur)-100 & temp<=-100));
                end
            end
            % get test data
            for irep = 1:nrepsOFF
                temp = M1OFF(1,1,idur,irep).spiketimes;
                ONSETspikesPerRep(idur,irep) = length(find(temp>0 & temp<=analysisWindow));
                ONSET2spikesPerRep(idur,irep) = length(find(temp>0 & temp<=durs(idur)));
            end
            % test for significance
            ntemp = min([size(ONSETspikesPerRep,2) size(ONSETspikesPerRepSS,2)]);
            [H_onset(idur),P_onset(idur)] = ttest(ONSETspikesPerRep(idur,1:ntemp),ONSETspikesPerRepSS(1:ntemp),'tail',flag.tailstr);
            
            ntemp = min([size(ONSET2spikesPerRep,2) size(ONSET2spikesPerRepSS,2)]);
            [H_onset2(idur),P_onset2(idur)] = ttest(ONSET2spikesPerRep(idur,1:ntemp),ONSET2spikesPerRepSS(idur,1:ntemp),'tail',flag.tailstr);
        end
        
        if flag.plot
            figure; hold on
            plot(log10(durs),P_onset,'o-','color',[0 .6 0])
            plot(log10(durs),P_offset,'o-','color',[.6 0 .6])
            plot(log10(durs),P_onset2,'o-','color',[0 1 0])
            plot(log10(durs),P_offset2,'o-','color',[1 0 1])
            legend('onset','offset','onset2','offset2')
            xl = xlim;
            plot(xl,[1 1]*.05,'k')
            set(gca,'xtick',log10(durs))
            set(gca,'xticklabel',durs)
            ylabel('P-value')
            xlabel('WN duration')
            
            temp=strsplit(datadir, '\');
            title([ temp{end} '  ' t_filename ],'interpreter','none')
            x = xl(1) + range(xl)*.15;
            text(x,.4, ['WNtc ' num2str(nrepsOFF) ' reps'],'fontsize',14)
            if SS
                str = 'calculated vs SS';
            else
                str = 'calculated vs pre-stim window';
            end
            text(x,.3, str,'fontsize',14)
            text(x,.2, 'offsets for dur<64 contaminated by onsets?','fontsize',14)
        end
        
        CELL_STATS(ProcessedCellNum).WNtc.analysisWindow = analysisWindow;
        CELL_STATS(ProcessedCellNum).WNtc.nrepsOFF = nrepsOFF;   % how many reps?
        CELL_STATS(ProcessedCellNum).WNtc.durs = durs;                % list of durations
        CELL_STATS(ProcessedCellNum).WNtc.OFFSETspikesPerRepSS = OFFSETspikesPerRepSS;
        CELL_STATS(ProcessedCellNum).WNtc.OFFSETspikesPerRep = OFFSETspikesPerRep;
        CELL_STATS(ProcessedCellNum).WNtc.ONSETspikesPerRepSS = ONSETspikesPerRepSS;
        CELL_STATS(ProcessedCellNum).WNtc.ONSETspikesPerRep = ONSETspikesPerRep;
        CELL_STATS(ProcessedCellNum).WNtc.P_onset = P_onset;          % Probability of onset for each duration
        CELL_STATS(ProcessedCellNum).WNtc.P_offset = P_offset;        % Probability of offset for each duration
        
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET2spikesPerRepSS = OFFSET2spikesPerRepSS;
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET2spikesPerRep = OFFSET2spikesPerRep;
        CELL_STATS(ProcessedCellNum).WNtc.ONSET2spikesPerRepSS = ONSET2spikesPerRepSS;
        CELL_STATS(ProcessedCellNum).WNtc.ONSET2spikesPerRep = ONSET2spikesPerRep;
        CELL_STATS(ProcessedCellNum).WNtc.P_onset2 = P_onset2;          % Probability of onset for each duration
        CELL_STATS(ProcessedCellNum).WNtc.P_offset2 = P_offset2;        % Probability of offset for each duration
        
        % summary stats WNtc
        ndurs = length(durs);
        % ONSET (first analysisWindow)
        [~,temp]=find(CELL_STATS(ProcessedCellNum).WNtc.P_onset <=.05);
        SPONT = mean(CELL_STATS(ProcessedCellNum).WNtc.ONSETspikesPerRepSS);
        RESP = mean(CELL_STATS(ProcessedCellNum).WNtc.ONSETspikesPerRep,2);  % response to each duration
        CELL_STATS(ProcessedCellNum).WNtc.ONSET_respSIGN = sign(RESP/SPONT-1);
        CELL_STATS(ProcessedCellNum).WNtc.ONSET_minDur = [];
        CELL_STATS(ProcessedCellNum).WNtc.ONSET_bestDur = [];
        startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
        
        if ~isempty(startIndex)
            CELL_STATS(ProcessedCellNum).WNtc.ONSET_minDur = min(temp);
            % calc RESP to every consecutive pair
            temp2 = nan(size(startIndex));
            for iind = 1:length(startIndex)
                temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
            end
            ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
            % bestDur defined by p<=.05
            ind0 = [ind ind+1];
            if CELL_STATS(ProcessedCellNum).WNtc.ONSET_respSIGN(ind)
                RESPmax = max(RESP(ind0));
                CELL_STATS(ProcessedCellNum).WNtc.ONSET_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
            elseif CELL_STATS(ProcessedCellNum).WNtc.ONSET_respSIGN(ind) == -1
                CELL_STATS(ProcessedCellNum).WNtc.ONSET_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
            end
            
            % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
            [~,temp]=find(CELL_STATS(ProcessedCellNum).WNtc.P_onset <=.05 & RESP'>= RESPmax*.5);
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            temp2 = nan(size(startIndex));
            for iind = 1:length(startIndex)
                temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
            end
            ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
            pmet = fliplr(temp);
            endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
            A = ismember(1:ndurs,startIndex) | ismember(1:ndurs,endIndex);
            [~,endIndex]=find(diff(A(ind:end))<0);
            endIndex=min(endIndex);
            if ~isempty(endIndex)
                endIndex = ind+endIndex;
            else
                endIndex = ndurs;
            end
            [~,startIndex]=find(diff(A(1:ind))>0);
            if ~isempty(startIndex)
                startIndex = 1+startIndex;
            else
                startIndex = 1;
            end
            CELL_STATS(ProcessedCellNum).WNtc.ONSET_tunedDurs = startIndex:endIndex;
        end
        
        % ONSET2 (entire duration)
        [~,ind]=find(CELL_STATS(ProcessedCellNum).WNtc.P_onset2 <=.05);
        SPONT = mean(CELL_STATS(ProcessedCellNum).WNtc.ONSET2spikesPerRepSS,2);
        RESP = mean(CELL_STATS(ProcessedCellNum).WNtc.ONSET2spikesPerRep,2);  % response to each duration
        CELL_STATS(ProcessedCellNum).WNtc.ONSET2_respSIGN = sign(RESP./SPONT-1);
        CELL_STATS(ProcessedCellNum).WNtc.ONSET2_minDur = [];
        CELL_STATS(ProcessedCellNum).WNtc.ONSET2_bestDur = [];
        temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
        if ~isempty(temp)
            CELL_STATS(ProcessedCellNum).WNtc.ONSET2_minDur = min(ind);
            % calc RESP to every consecutive pair
            temp2 = nan(size(temp));
            for iind = 1:length(temp)
                temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
            end
            ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
            if CELL_STATS(ProcessedCellNum).WNtc.ONSET2_respSIGN(ind)
                ind = [ind ind+1];
                CELL_STATS(ProcessedCellNum).WNtc.ONSET2_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
            elseif CELL_STATS(ProcessedCellNum).WNtc.ONSET2_respSIGN(ind) == -1
                ind = [ind ind+1];
                CELL_STATS(ProcessedCellNum).WNtc.ONSET2_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
            end
        end
        
        % OFFSET (1st analysisWindow)
        [~,temp]=find(CELL_STATS(ProcessedCellNum).WNtc.P_offset <=.05);
        SPONT = mean(CELL_STATS(ProcessedCellNum).WNtc.OFFSETspikesPerRepSS,2);
        RESP = mean(CELL_STATS(ProcessedCellNum).WNtc.OFFSETspikesPerRep,2);  % response to each duration
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET_respSIGN = sign(RESP./SPONT-1);
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET_minDur = [];
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET_bestDur = [];
        startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
        if ~isempty(startIndex)
            CELL_STATS(ProcessedCellNum).WNtc.OFFSET_minDur = min(temp);
            % calc RESP to every consecutive pair
            temp2 = nan(size(startIndex));
            for iind = 1:length(startIndex)
                temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
            end
            ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
            
            % bestDur defined by p<=.05
                ind0 = [ind ind+1];
            if CELL_STATS(ProcessedCellNum).WNtc.OFFSET_respSIGN(ind)
                RESPmax = max(RESP(ind0));
                CELL_STATS(ProcessedCellNum).WNtc.OFFSET_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
            elseif CELL_STATS(ProcessedCellNum).WNtc.OFFSET_respSIGN(ind) == -1
                CELL_STATS(ProcessedCellNum).WNtc.OFFSET_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
            end
            
            % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
            [~,temp]=find(CELL_STATS(ProcessedCellNum).WNtc.P_offset <=.05 & RESP'>= RESPmax*.5);
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            temp2 = nan(size(startIndex));
            for iind = 1:length(startIndex)
                temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
            end
            ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
            pmet = fliplr(temp);
            endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
            A = ismember(1:ndurs,startIndex) | ismember(1:ndurs,endIndex);
            [~,endIndex]=find(diff(A(ind:end))<0);
            endIndex=min(endIndex);
            if ~isempty(endIndex)
                endIndex = ind+endIndex;
            else
                endIndex = ndurs;
            end
            [~,startIndex]=find(diff(A(1:ind))>0);
            if ~isempty(startIndex)
                startIndex = 1+startIndex;
            else
                startIndex = 1;
            end
            CELL_STATS(ProcessedCellNum).WNtc.OFFSET_tunedDurs = startIndex:endIndex;
        end
        
        % OFFSET2 (2nd analysisWindow)
        [~,ind]=find(CELL_STATS(ProcessedCellNum).WNtc.P_offset2 <=.05);
        SPONT = mean(CELL_STATS(ProcessedCellNum).WNtc.OFFSET2spikesPerRepSS,2);
        RESP = mean(CELL_STATS(ProcessedCellNum).WNtc.OFFSET2spikesPerRep,2);  % response to each duration
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_respSIGN = sign(RESP./SPONT-1);
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_minDur = [];
        CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_bestDur = [];
        temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
        if ~isempty(temp)
            CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_minDur = min(ind);
            % calc RESP to every consecutive pair
            temp2 = nan(size(temp));
            for iind = 1:length(temp)
                temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
            end
            ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
            if CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_respSIGN(ind)
                ind = [ind ind+1];
                CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
            elseif CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_respSIGN(ind) == -1
                ind = [ind ind+1];
                CELL_STATS(ProcessedCellNum).WNtc.OFFSET2_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case {'GPIAS w/o burst' , 'GPIAS w/ burst'}
        % spiketimes are indexed to gap OFFSET (time==0)
        analysisWindow = 50;
        [minGap, ind_minGap] =min(gapdurs);
        nrepsOFF = min(nrepsOFF);
        %%%%% following GAP ONSET (within gaps)
        % set up arrays
        INGAPspikesPerRep = nan(numgapdurs+2,nrepsOFF);
        H_ingap = nan(1,numgapdurs);
        P_ingap = H_ingap;
        
        % ingap2 uses all spikes within the gap and is only computed if minGap==0
        if minGap==0
            INGAP2spikesPerRep = nan(numgapdurs*2,nrepsOFF);
            H_ingap2 = nan(1,numgapdurs);
            P_ingap2 = H_ingap2;
            for idur = 1:numgapdurs
                for irep = 1:nrepsOFF
                    % controls (minGap==0)
                    temp = M1OFF(ind_minGap,1,irep).spiketimes;
                    INGAP2spikesPerRep(numgapdurs+idur,irep) = length(find(temp>-gapdurs(idur) & temp<=0));
                    % tests
                    temp = M1OFF(idur,1,irep).spiketimes;
                    INGAP2spikesPerRep(idur,irep) = length(find(temp>-gapdurs(idur) & temp<=0));
                end
                % Hypothesis and P-value for each duration
                [~,P_ingap2(idur)] = ttest(INGAP2spikesPerRep(idur,:), INGAP2spikesPerRep(numgapdurs+idur,:),'tail',flag.tailstr);
            end
        end
        
        % ingap control data taken from tests with minimum gap (could be zero)
        for irep = 1:nrepsOFF
            temp = M1OFF(ind_minGap,1,irep).spiketimes;
            % for use with full analysisWindow
            INGAPspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap + analysisWindow) & temp<=-minGap));
            % for use with 32 ms gap
            INGAPspikesPerRep(numgapdurs+2,irep) = length(find(temp>-(minGap + 32) & temp<=-minGap));
        end
        % ingap test data for each duration
        for idur = 1:numgapdurs
            if gapdurs(idur) == 32
                for irep = 1:nrepsOFF
                    temp = M1OFF(idur,1,irep).spiketimes;
                    INGAPspikesPerRep(idur,irep) = length(find(temp>-gapdurs(idur) & temp<=0));
                end
                % Hypothesis and P-value for each duration
                [H_ingap(idur),P_ingap(idur)] = ttest(INGAPspikesPerRep(idur,:), INGAPspikesPerRep(numgapdurs+2,:),'tail',flag.tailstr);
            elseif gapdurs(idur) >=analysisWindow
                for irep = 1:nrepsOFF
                    temp = M1OFF(idur,1,irep).spiketimes;
                    INGAPspikesPerRep(idur,irep) = length(find(temp>-gapdurs(idur) & temp<=-gapdurs(idur)+analysisWindow));
                end
                % Hypothesis and P-value for each duration
                [H_ingap(idur),P_ingap(idur)] = ttest(INGAPspikesPerRep(idur,:), INGAPspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            end
        end
        
        if flag.plot
            figure; hold on
            plot(log10(gapdurs),P_ingap,'o-','color',[.6 0 .6])
            plot(log10(gapdurs),P_ingap2,'*-','color',[1 0 1])
            legend('GAP ONSET', 'GAP ENTIRE')
            ylim([-.05 1]);
            xl = xlim; yl = ylim;
            plot(xl,[1 1]*.05,'k')
            set(gca,'xtick',log10(gapdurs))
            set(gca,'xticklabel',gapdurs)
            ylabel('P-value')
            xlabel('GAP duration')
            temp=strsplit(datadir, '\');
            title([ temp{end} '  ' t_filename ],'interpreter','none')
            text(log10(40),yl(2)*.2,['First 50 ms of Gap ONSET ' num2str(nrepsOFF) ' reps'],'fontsize',14,'color',[.6 0 .6])
            text(log10(40),yl(2)*.3,'ENTIRE Gap','fontsize',14,'color',[1 0 1])
        end
        
        %%%%% following GAP OFFSET (after gap finishes)
        % set up arrays
        PGIspikesPerRep = nan(numgapdurs+1,nrepsOFF);      % entire analysisWindow
        PGIAspikesPerRep = nan(numgapdurs+1,nrepsOFF);     % first half of analysisWindow
        PGIBspikesPerRep = nan(numgapdurs+1,nrepsOFF);     % second half of analysisWindow
        PGICspikesPerRep = nan(numgapdurs+1,nrepsOFF);     % 1:1.5 analysisWindow
        PGIDspikesPerRep = nan(numgapdurs+1,nrepsOFF);     % 1.5:2 analysisWindow
        
        %[minGap, ind] = min(gapdurs);
        % control data taken from tests with minimum gap (could be zero)
        binwidth = 2;
        bins = xlimits(1):binwidth:xlimits(2);
        nbins = length(bins);
        
        SPONT = [];
        if minGap           % if there is no 0 gap, then use data before minGap
            for irep = 1:nrepsOFF
                temp = M1OFF(ind_minGap,1,irep).spiketimes;
                SPONT = [SPONT temp(temp<=-minGap)];
                PGIspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap+analysisWindow) & temp<=-minGap));
                PGIAspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap+analysisWindow) & temp<=-(minGap+analysisWindow/2)));
                PGIBspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap+analysisWindow/2) & temp<=-minGap));
                PGICspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap+analysisWindow*1.5) & temp<=-(minGap+analysisWindow)));
                PGIDspikesPerRep(numgapdurs+1,irep) = length(find(temp>-(minGap+analysisWindow*2) & temp<=-(minGap+analysisWindow*1.5)));
            end
            nbin = find(bins<-minGap,1,'last');
        else        % minGap==0
            for irep = 1:nrepsOFF
                temp = M1OFF(ind_minGap,1,irep).spiketimes;
                SPONT = [SPONT temp];
                PGIspikesPerRep(numgapdurs+1,irep) = length(find(temp>0 & temp<=analysisWindow));
                PGIAspikesPerRep(numgapdurs+1,irep) = length(find(temp>0 & temp<=(analysisWindow/2)));
                PGIBspikesPerRep(numgapdurs+1,irep) = length(find(temp>analysisWindow/2 & temp<=analysisWindow));
                PGICspikesPerRep(numgapdurs+1,irep) = length(find(temp>analysisWindow & temp<=1.5*analysisWindow));
                PGIDspikesPerRep(numgapdurs+1,irep) = length(find(temp>1.5*analysisWindow & temp<=2*analysisWindow));
            end
            nbin=nbins;
        end
        if isempty(SPONT)
            SPONT = 0;
        else
            SPONT = histc(SPONT,bins)/(nrepsOFF*binwidth/1000);     % spikes/second
            SPONT = smooth(SPONT,5);
            SPONT = mean(SPONT(1:nbin));
        end
        
        % test data for all gap durations
        % set up arrays
        H_PGI = nan(1,numgapdurs); H_PGIA = H_PGI; H_PGIB = H_PGI; H_PGIC = H_PGI; H_PGID = H_PGI;
        P_PGI = H_PGI;  P_PGIA = P_PGI; P_PGIB = P_PGI; P_PGIC = P_PGI;  P_PGID = P_PGI;
        
        for idur = 1:numgapdurs
            allSpiketimes = [];
            LATENCY_LOC = [];
            for irep = 1:nrepsOFF
                temp = M1OFF(idur,1,irep).spiketimes;
                PGIspikesPerRep(idur,irep) = length(find(temp>0 & temp<=analysisWindow));
                PGIAspikesPerRep(idur,irep) = length(find(temp>0 & temp<=analysisWindow/2));
                PGIBspikesPerRep(idur,irep) = length(find(temp>0+analysisWindow/2 & temp<=analysisWindow));
                PGICspikesPerRep(idur,irep) = length(find(temp>0+analysisWindow & temp<=1.5*analysisWindow));
                PGIDspikesPerRep(idur,irep) = length(find(temp>0+1.5*analysisWindow & temp<=2*analysisWindow));
                allSpiketimes = [allSpiketimes temp];
            end
            % hypothesis and P-values (two-tailed)
            [H_PGI(idur),P_PGI(idur)] = ttest(PGIspikesPerRep(idur,:), PGIspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            [H_PGIA(idur),P_PGIA(idur)] = ttest(PGIAspikesPerRep(idur,:), PGIAspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            [H_PGIB(idur),P_PGIB(idur)] = ttest(PGIBspikesPerRep(idur,:), PGIBspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            [H_PGIC(idur),P_PGIC(idur)] = ttest(PGICspikesPerRep(idur,:), PGICspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            [H_PGID(idur),P_PGID(idur)] = ttest(PGIDspikesPerRep(idur,:), PGIDspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            
            % trying to get latency to PGI
            if gapdurs(idur)==256 & P_PGI(idur)<= 0.05
                binwidth = 2;
                bins = xlimits(1):binwidth:xlimits(2);
                nbin(1) = find(bins>=-gapdurs(idur), 1 );
                nbin(2) = find(bins>=0, 1 );
                nbin(3) = find(bins>=100,1);
                EVOKED = histc(allSpiketimes,bins)/(nrepsOFF*binwidth/1000);     % spikes/second
                EVOKEDsm = smooth(EVOKED,5);
                SPONT = mean(EVOKED(1:nbin(1)));
                SPONTsm = mean(EVOKEDsm(1:nbin(1)));
                SPONTsmstd = std(EVOKEDsm(1:nbin(1)));
                
                figure; hold on
                bar(bins,EVOKED,'g','edgecolor','g')
                plot(bins,EVOKEDsm,'b','linewidth',2)
                plot([xlimits(1) xlimits(2)],[1 1]*(SPONTsm),'y')
                plot([xlimits(1) xlimits(2)],[1 1]*(SPONTsm+3*SPONTsmstd),'y','linewidth',2)
                LATENCY_LOC = find(EVOKEDsm(nbin(2):end)>=(SPONTsm+3*SPONTsmstd),1)-1;
                if ~isempty(LATENCY_LOC)
                    plot(bins(nbin(2)+LATENCY_LOC),SPONTsm+3*SPONTsmstd,'k*')
                    LATENCY_LOC=LATENCY_LOC*binwidth;
                end
                [PKsm,PKsm_LOC] = max(EVOKEDsm(nbin(2):nbin(3)));
                HALF_PKsm = SPONTsm + (PKsm - SPONTsm)/2;
                HALF_PKsm_LOC = find(EVOKEDsm(nbin(2):nbin(2)+PKsm_LOC-1)>=HALF_PKsm,1);
                if ~isempty(HALF_PKsm_LOC)
                    plot(bins(PKsm_LOC+nbin(2)-1),PKsm,'bo','markersize',12);
                    plot(bins(HALF_PKsm_LOC+nbin(2)-1),HALF_PKsm,'bo','markersize',12);
                end
                if ~isempty(LATENCY_LOC)
                    if LATENCY_LOC<2 | LATENCY_LOC>30
                        fprintf('\nIDnumber: %d   Latency outlier :%d\n',IDnumber,LATENCY_LOC)
                        pause
                    end
                end
            end
            
       end
        
        if flag.plot
            figure; hold on
            plot(log10(gapdurs),P_PGI,'o-','color',[0 0 .6])
            plot(log10(gapdurs),P_PGIA,'o-','color',[0 .6 0])
            plot(log10(gapdurs),P_PGIB,'o-','color',[.6 0 0])
            plot(log10(gapdurs),P_PGIC,'*-','color',[0 0 1])
            plot(log10(gapdurs),P_PGID,'*-','color',[0 1 0])
            legend('PGI', 'PGI-A','PGI-B','PGI-C', 'PGI-D')
            ylim([-.05 1]);
            xl = xlim; yl = ylim;
            plot(xl,[1 1]*.05,'k')
            set(gca,'xtick',log10(gapdurs))
            set(gca,'xticklabel',gapdurs)
            ylabel('P-value')
            xlabel('GAP duration')
            temp=strsplit(datadir, '\');
            title([ temp{end} '  ' t_filename ],'interpreter','none')
            text(log10(gapdurs(3)),yl(2)*.8,['GapTerminationResp ' num2str(nrepsOFF) ' reps'],'fontsize',14)
            text(log10(gapdurs(3)),yl(2)*.7,'For dur<64, contaminated with gap onset','fontsize',14)
            text(log10(gapdurs(3)),yl(2)*.6,['AnalysisWin: ' num2str(analysisWindow)],'fontsize',14)
        end
        
        if strcmp(Type,'GPIAS w/ burst')
            % analyze response to loud burst
            burstDelay = 50;
            if flag.plot
                text(log10(gapdurs(3)),yl(2)*.6,['AnalysisWin: ' num2str(analysisWindow) ' burstDelay: ' num2str(burstDelay)],'fontsize',14)
            end
            
            % set up arrays
            BURSTspikesPerRep = nan(numgapdurs+1,nrepsOFF);
            H_burst = nan(1,numgapdurs);
            P_burst = H_burst;
            BURST2spikesPerRep = nan(numgapdurs+1,nrepsOFF);
            H_burst2 = nan(1,numgapdurs);
            P_burst2 = H_burst2;
            
            % control response (to burst without gap or with minimum gap, using 200% analysisWindow)
            for irep = 1:nrepsOFF
                temp = M1OFF(ind_minGap,1,irep).spiketimes;
                BURSTspikesPerRep(numgapdurs+1,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
                BURST2spikesPerRep(numgapdurs+1,irep) = length(find(temp>burstDelay+analysisWindow*2 & temp<=burstDelay+analysisWindow*4));
            end
            % response to burst with each gap duration (using 200% analysis window)
            for idur = 1:numgapdurs
                for irep = 1:nrepsOFF
                    temp = M1OFF(idur,1,irep).spiketimes;
                    BURSTspikesPerRep(idur,irep) = length(find(temp>burstDelay & temp<=burstDelay+analysisWindow*2));
                    BURST2spikesPerRep(idur,irep) = length(find(temp>burstDelay+analysisWindow*2 & temp<=burstDelay+analysisWindow*4));
                end
                % Hypothesis and P-value
                [H_burst(idur),P_burst(idur)] = ttest(BURSTspikesPerRep(idur,:), BURSTspikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
                [H_burst2(idur),P_burst2(idur)] = ttest(BURST2spikesPerRep(idur,:), BURST2spikesPerRep(numgapdurs+1,:),'tail',flag.tailstr);
            end
            
            if flag.plot
                figure; hold on
                plot(log10(gapdurs),P_burst,'o-','color',[.6 0 .6])
                %                 plot(log10(gapdurs),Prank_burst,'*-','color',[1 0 1])
                %                 legend('burst', 'burst (nonP)')
                plot(log10(gapdurs),Prank_burst2,'*-','color',[1 0 1])
                legend('burst', 'burst2')
                ylim([-.05 1]);
                xl = xlim; %yl = ylim;
                plot(xl,[1 1]*.05,'k')
                set(gca,'xtick',log10(gapdurs))
                set(gca,'xticklabel',gapdurs)
                ylabel('P-value')
                xlabel('GAP duration')
                temp=strsplit(datadir, '\');
                title([ temp{end} '  ' t_filename ],'interpreter','none')
                text(log10(gapdurs(2)),.6,['PostGap Response to BURST. Measured over ' num2str(analysisWindow*2)  ' ms  ' num2str(nrepsOFF) ' reps'],'fontsize',14)
            end
            CELL_STATS(ProcessedCellNum).GPbehav.analysisWindow = analysisWindow;
            CELL_STATS(ProcessedCellNum).GPbehav.nrepsOFF = nrepsOFF;
            CELL_STATS(ProcessedCellNum).GPbehav.gapdurs = gapdurs;
            CELL_STATS(ProcessedCellNum).GPbehav.INGAPspikesPerRep = INGAPspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.PGIspikesPerRep = PGIspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.PGIAspikesPerRep = PGIAspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.PGIBspikesPerRep = PGIBspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.P_ingap = P_ingap;
            CELL_STATS(ProcessedCellNum).GPbehav.P_PGI = P_PGI;
            CELL_STATS(ProcessedCellNum).GPbehav.P_PGIA = P_PGIA;
            CELL_STATS(ProcessedCellNum).GPbehav.P_PGIB = P_PGIB;
            CELL_STATS(ProcessedCellNum).GPbehav.BURSTspikesPerRep = BURSTspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.P_burst = P_burst;
            
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP2spikesPerRep = INGAP2spikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.PGICspikesPerRep = PGICspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.PGIDspikesPerRep = PGIDspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.P_ingap2 = P_ingap2;
            CELL_STATS(ProcessedCellNum).GPbehav.P_PGIC = P_PGIC;
            CELL_STATS(ProcessedCellNum).GPbehav.P_PGID = P_PGID;
            CELL_STATS(ProcessedCellNum).GPbehav.BURST2spikesPerRep = BURST2spikesPerRep;
            CELL_STATS(ProcessedCellNum).GPbehav.P_burst2 = P_burst2;
            
            % summary stats: GPbehav
            ngaps = length(CELL_STATS(ProcessedCellNum).GPbehav.gapdurs);
            % INGAP (first analysisWindow)
            [~,temp]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_ingap <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.INGAPspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.INGAPspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP_bestDur = [];
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            if ~isempty(startIndex)
                CELL_STATS(ProcessedCellNum).GPbehav.INGAP_minDur = min(temp);
                % calc RESP to every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                
                % bestDur defined by p<=.05
                    ind0 = [ind ind+1];
                if CELL_STATS(ProcessedCellNum).GPbehav.INGAP_respSIGN(ind)
                    RESPmax = max(RESP(ind0));
                    CELL_STATS(ProcessedCellNum).GPbehav.INGAP_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.INGAP_respSIGN(ind) == -1
                    CELL_STATS(ProcessedCellNum).GPbehav.INGAP_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
                end
                
                % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
                [~,temp]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_ingap <=.05 & RESP'>= RESPmax*.5);
                startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                pmet = fliplr(temp);
                endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
                A = ismember(1:ngaps,startIndex) | ismember(1:ngaps,endIndex);
                [~,endIndex]=find(diff(A(ind:end))<0);
                endIndex=min(endIndex);
                if ~isempty(endIndex)
                    endIndex = ind+endIndex;
                else
                    endIndex = ngaps;
                end
                [~,startIndex]=find(diff(A(1:ind))>0);
                if ~isempty(startIndex)
                    startIndex = 1+startIndex;
                else
                    startIndex = 1;
                end
                CELL_STATS(ProcessedCellNum).GPbehav.INGAP_tunedDurs = startIndex:endIndex;
            end
            
            % INGAP2 (entire gap)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_ingap2 <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.INGAP2spikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.INGAP2spikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.INGAP2_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP (PGI: one analysis Window)
            [~,temp]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_PGI <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.PGI_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.PGI_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.PGI_bestDur = [];
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            if ~isempty(startIndex)
                CELL_STATS(ProcessedCellNum).GPbehav.PGI_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                
                % bestDur defined by p<=.05
                ind0 = [ind ind+1];
                if CELL_STATS(ProcessedCellNum).GPbehav.PGI_respSIGN(ind)
                    RESPmax = max(RESP(ind0));
                    CELL_STATS(ProcessedCellNum).GPbehav.PGI_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.PGI_respSIGN(ind) == -1
                    CELL_STATS(ProcessedCellNum).GPbehav.PGI_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
                end
                
                % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
                [~,temp]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_PGI <=.05 & RESP'>= RESPmax*.5);
                startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                pmet = fliplr(temp);
                endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
                A = ismember(1:ngaps,startIndex) | ismember(1:ngaps,endIndex);
                [~,endIndex]=find(diff(A(ind:end))<0);
                    endIndex=min(endIndex);
                    if ~isempty(endIndex)
                        endIndex = ind+endIndex;
                    else
                        endIndex = ngaps;
                    end
                    [~,startIndex]=find(diff(A(1:ind))>0);
                    if ~isempty(startIndex)
                        startIndex = 1+startIndex;
                    else
                        startIndex = 1;
                    end
                    CELL_STATS(ProcessedCellNum).WNtc.GPbehav.PGI_tunedDurs = startIndex:endIndex;
            end
            
            % postGAP-A (PGIA: first half analysis window)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_PGIA <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIAspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIAspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.PGIA_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.PGIA_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.PGIA_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPbehav.PGIA_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPbehav.PGIA_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGIA_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.PGIA_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGIA_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP-B (PGIB: second half analysis window)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_PGIB <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIBspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIBspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.PGIB_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.PGIB_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.PGIB_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPbehav.PGIB_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPbehav.PGIB_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGIB_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.PGIB_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGIB_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP-C (PGIC: 2nd full analysisWindow)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_PGIC <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGICspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGICspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.PGIC_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.PGIC_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.PGIC_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPbehav.PGIC_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPbehav.PGIC_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGIC_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.PGIC_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGIC_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP-D ( 3rd full analysisWindow)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_PGID <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIDspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.PGIDspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.PGID_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.PGID_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.PGID_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPbehav.PGID_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPbehav.PGID_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGID_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.PGID_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.PGID_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % burst (first two analysisWindows)
            [~,temp]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_burst <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.BURSTspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.BURSTspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.BURST_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.BURST_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.BURST_bestDur = [];
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            if ~isempty(startIndex)
                CELL_STATS(ProcessedCellNum).GPbehav.BURST_minDur = min(temp);
                % calc RESP to every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                
                % bestDur defined by p<=.05
                ind0 = [ind ind+1];
                if CELL_STATS(ProcessedCellNum).GPbehav.BURST_respSIGN(ind)
                    RESPmax = max(RESP(ind0));
                    CELL_STATS(ProcessedCellNum).GPbehav.BURST_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.BURST_respSIGN(ind) == -1
                    CELL_STATS(ProcessedCellNum).GPbehav.BURST_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
                end
                
                % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
                [~,temp]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_burst <=.05 & RESP'>= RESPmax*.5);
                startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                pmet = fliplr(temp);
                endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
                A = ismember(1:ngaps,startIndex) | ismember(1:ngaps,endIndex);
                [~,endIndex]=find(diff(A(ind:end))<0);
                endIndex=min(endIndex);
                if ~isempty(endIndex)
                    endIndex = ind+endIndex;
                else
                    endIndex = ngaps;
                end
                [~,startIndex]=find(diff(A(1:ind))>0);
                if ~isempty(startIndex)
                    startIndex = 1+startIndex;
                else
                    startIndex = 1;
                end
                CELL_STATS(ProcessedCellNum).WNtc.GPbehav.BURST_tunedDurs = startIndex:endIndex;
            end
            
            % burst2 (3rd and 4th analysisWindows)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPbehav.P_burst2 <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPbehav.BURST2spikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPbehav.BURST2spikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPbehav.BURST2_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPbehav.BURST2_minDur = [];
            CELL_STATS(ProcessedCellNum).GPbehav.BURST2_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPbehav.BURST2_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPbehav.BURST2_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.BURST2_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPbehav.BURST2_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPbehav.BURST2_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
        else
            %%%%%%%% GPtc
            CELL_STATS(ProcessedCellNum).GPtc.gapdurs = gapdurs;
            ngaps = length(CELL_STATS(ProcessedCellNum).GPtc.gapdurs);
            CELL_STATS(ProcessedCellNum).GPtc.nrepsOFF = nrepsOFF;
            CELL_STATS(ProcessedCellNum).GPtc.INGAPspikesPerRep = INGAPspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.PGIspikesPerRep = PGIspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.PGIAspikesPerRep = PGIAspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.PGIBspikesPerRep = PGIBspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.P_ingap = P_ingap;
            CELL_STATS(ProcessedCellNum).GPtc.P_PGI = P_PGI;
            CELL_STATS(ProcessedCellNum).GPtc.PGI_latency = LATENCY_LOC;
            CELL_STATS(ProcessedCellNum).GPtc.P_PGIA = P_PGIA;
            CELL_STATS(ProcessedCellNum).GPtc.P_PGIB = P_PGIB;
            
            CELL_STATS(ProcessedCellNum).GPtc.INGAP2spikesPerRep = INGAP2spikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.PGICspikesPerRep = PGICspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.PGIDspikesPerRep = PGIDspikesPerRep;
            CELL_STATS(ProcessedCellNum).GPtc.P_ingap2 = P_ingap2;
            CELL_STATS(ProcessedCellNum).GPtc.P_PGIC = P_PGIC;
            CELL_STATS(ProcessedCellNum).GPtc.P_PGID = P_PGID;
            
            % summary stats: GPtc
            % INGAP (first analysisWindow)
            [~,temp]=find(CELL_STATS(ProcessedCellNum).GPtc.P_ingap <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.INGAPspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.INGAPspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.INGAP_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.INGAP_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.INGAP_bestDur = [];
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            if ~isempty(startIndex)
                CELL_STATS(ProcessedCellNum).GPtc.INGAP_minDur = min(temp);
                % calc RESP to every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                
                % bestDur defined by p<=.05
                ind0 = [ind ind+1];
                if CELL_STATS(ProcessedCellNum).GPtc.INGAP_respSIGN(ind)
                    RESPmax = max(RESP(ind0));
                    CELL_STATS(ProcessedCellNum).GPtc.INGAP_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.INGAP_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.INGAP_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
                end
                
                % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
                [~,temp]=find(CELL_STATS(ProcessedCellNum).GPtc.P_ingap <=.05 & RESP'>= RESPmax*.5);
                startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                pmet = fliplr(temp);
                endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
                A = ismember(1:ngaps,startIndex) | ismember(1:ngaps,endIndex);
                [~,endIndex]=find(diff(A(ind:end))<0);
                endIndex=min(endIndex);
                if ~isempty(endIndex)
                    endIndex = ind+endIndex;
                else
                    endIndex = ngaps;
                end
                [~,startIndex]=find(diff(A(1:ind))>0);
                if ~isempty(startIndex)
                    startIndex = 1+startIndex;
                else
                    startIndex = 1;
                end
                CELL_STATS(ProcessedCellNum).GPtc.INGAP_tunedDurs = startIndex:endIndex;
            end
            
            % INGAP2 (entire gap)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPtc.P_ingap2 <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.INGAP2spikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.INGAP2spikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.INGAP2_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.INGAP2_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.INGAP2_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPtc.INGAP2_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPtc.INGAP2_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.INGAP2_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.INGAP2_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.INGAP2_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP (PGI: one analysis Window)
            [~,temp]=find(CELL_STATS(ProcessedCellNum).GPtc.P_PGI <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.PGI_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.PGI_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.PGI_bestDur = [];
            startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
            if ~isempty(startIndex)
                CELL_STATS(ProcessedCellNum).GPtc.PGI_minDur = min(temp);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                
                % bestDur defined by p<=.05
                ind0 = [ind ind+1];
                if CELL_STATS(ProcessedCellNum).GPtc.PGI_respSIGN(ind)
                    RESPmax = max(RESP(ind0));
                    CELL_STATS(ProcessedCellNum).GPtc.PGI_bestDur = min(ind0(RESP(ind0)==RESPmax));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.PGI_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGI_bestDur = min(ind0(RESP(ind0)==min(RESP(ind0))));   % index into min firing from max consecutive pair
                end
                
                % tunedDurs defined by p <= .05 && RESP>= max(RESP)*.5
                [~,temp]=find(CELL_STATS(ProcessedCellNum).GPtc.P_PGI <=.05 & RESP'>= RESPmax*.5);
                startIndex = temp(diff(temp)==1);  % indices into first of every consecutive pair
                temp2 = nan(size(startIndex));
                for iind = 1:length(startIndex)
                    temp2(iind) = sum(RESP(startIndex(iind):startIndex(iind)+1));
                end
                ind = max(startIndex(temp2==max(temp2)));    % index into consecutive pair with highest firing
                pmet = fliplr(temp);
                endIndex = pmet(diff(pmet)==-1);  % indices into last of every consecutive pair
                A = ismember(1:ngaps,startIndex) | ismember(1:ngaps,endIndex);
                [~,endIndex]=find(diff(A(ind:end))<0);
                endIndex=min(endIndex);
                if ~isempty(endIndex)
                    endIndex = ind+endIndex;
                else
                    endIndex = ngaps;
                end
                [~,startIndex]=find(diff(A(1:ind))>0);
                if ~isempty(startIndex)
                    startIndex = 1+startIndex;
                else
                    startIndex = 1;
                end
                CELL_STATS(ProcessedCellNum).GPtc.PGI_tunedDurs = startIndex:endIndex;
            end
            
            % postGAP-A (PGIA: first half analysis window)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPtc.P_PGIA <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIAspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIAspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.PGIA_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.PGIA_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.PGIA_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPtc.PGIA_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPtc.PGIA_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGIA_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.PGIA_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGIA_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP-B (PGIB: second half analysis window)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPtc.P_PGIB <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIBspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIBspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.PGIB_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.PGIB_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.PGIB_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPtc.PGIB_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPtc.PGIB_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGIB_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.PGIB_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGIB_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP-C (PGIC: 2nd full analysisWindow)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPtc.P_PGIC <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.PGICspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.PGICspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.PGIC_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.PGIC_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.PGIC_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPtc.PGIC_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPtc.PGIC_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGIC_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.PGIC_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGIC_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
            % postGAP-D ( 3rd full analysisWindow)
            [~,ind]=find(CELL_STATS(ProcessedCellNum).GPtc.P_PGID <=.05);
            SPONT = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIDspikesPerRep(ngaps+1,:));
            RESP = mean(CELL_STATS(ProcessedCellNum).GPtc.PGIDspikesPerRep(1:ngaps,:),2);  % response to each duration
            CELL_STATS(ProcessedCellNum).GPtc.PGID_respSIGN = sign(RESP/SPONT-1);
            CELL_STATS(ProcessedCellNum).GPtc.PGID_minDur = [];
            CELL_STATS(ProcessedCellNum).GPtc.PGID_bestDur = [];
            temp = ind(diff(ind)==1);  % indices into first of every consecutive pair
            if ~isempty(temp)
                CELL_STATS(ProcessedCellNum).GPtc.PGID_minDur = min(ind);
                % calc RESP to every consecutive pair
                temp2 = nan(size(temp));
                for iind = 1:length(temp)
                    temp2(iind) = sum(RESP(temp(iind):temp(iind)+1));
                end
                ind = temp(temp2==max(temp2));    % index into consecutive pair with highest firing
                if CELL_STATS(ProcessedCellNum).GPtc.PGID_respSIGN(ind)
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGID_bestDur = min(ind(RESP(ind)==max(RESP(ind))));   % index into max firing from max consecutive pair
                elseif CELL_STATS(ProcessedCellNum).GPtc.PGID_respSIGN(ind) == -1
                    ind = [ind ind+1];
                    CELL_STATS(ProcessedCellNum).GPtc.PGID_bestDur = min(ind(RESP(ind)==min(RESP(ind))));   % index into min firing from max consecutive pair
                end
            end
            
        end
end


