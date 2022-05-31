function Outfile_Combiner(varargin)

%simple GUI to pick multiple outfiles and combine them.
%developed for ABRs (LFP tuning curves) so that we can split data
%acquisition across multiple sessions
%
%Usage: click "Browse and Add" to select outfiles and add them to the list of
%outfiles that you want to combine. Once you're finished adding files to
%the list, then click "Combine them" which will generate an
%"out_combined.mat" file in the chronologically first directory.
%To plot the combined file, run PlotTC_LFP without any input and use the
%dialog box to select the combined outfile.
% current version only allows outfiles with identical
% frequencies/amplitudes, but in the future we can revise it to combine
% outfiles with different parameters.
%
%You can also run without the GUI by passing a filelist. This should be a
%cell array of absolute filenames (including path) passed as the first
%argument. In this case, you can (optionally) include the destination
%directory for the combined outfile as the second argument (defaults to
%current directory)


global P

if nargin > 0
    if  iscell(varargin{1})
        action='outfilelist'; %run without the GUI by passing a filelist
    else
        action=varargin{1};
    end
else
    action = get(gcbo,'tag');
end
if isempty(action)
    action='Init';
end

switch action
    case 'Init'
        InitializeGUI
    case 'reset'
        InitializeGUI
    case 'Close'
        delete(P.fig)
        clear global P
    case 'CreateNewOutfile'
        CreateNewOutfile
    case 'BrowseAndAdd'
        BrowseAndAdd
    case 'outfilelist'
        %assume varargin{1} is a cell array of filenames
        P.outfilelist=varargin{1};
        fprintf('\nusing passed list of outfiles:')
        fprintf('\n%s', P.outfilelist{:})
        P.numoutfiles=length(P.outfilelist);
        for i=1:P.numoutfiles
            P.dirlist{i}='/'; 
        end
        if nargin==2
            %assume varargin{2} is the target directory
            P.targetdir=varargin{2};
        else
            P.targetdir=pwd;
        end
        fprintf('\nusing %s as targetdir', P.targetdir)
        CreateNewOutfile
    otherwise
        error('unrecognized action')
end


function CreateNewOutfile
global P
fprintf('\nloading and combining outfiles...');
wb=waitbar(0,'loading and combining outfiles...');
sortedoutfilelist=sort (P.outfilelist);
sorteddirlist=sort (P.dirlist);
if isfield(P, 'targetdir')
    targetdir=P.targetdir;
else
    targetdir=sorteddirlist{1};
end
Out.dirlist=sorteddirlist;
Out.targetdir=sorteddirlist{1};
Out.generated_by=mfilename;
Out.generated_on=datestr(now);
Out.outfilelist=sortedoutfilelist;

for i=1:P.numoutfiles
    fprintf('.')
    cd( sorteddirlist{i});
    Out_components(i)=load(sortedoutfilelist{i});
    waitbar(i/(1+P.numoutfiles), wb);
end

%we will handle outfiles differently depending on what kind of data they
%are. We can use some fieldnames as a proxy for data type

if isfield(Out_components(i).out, 'freqs') & ...
        isfield(Out_components(i).out, 'amps') & ...
        isfield(Out_components(i).out, 'durs')
    experiment_type = 'tuningcurve';
elseif isfield(Out_components(i).out, 'gapdurs') & ...
        isfield(Out_components(i).out, 'PeakON')
    experiment_type = 'GPIASbehavior';
elseif isfield(Out_components(i).out, 'amps') & ...
        isfield(Out_components(i).out, 'durs') & ...
        isfield(Out_components(i).out, 'sourcefiles')
    experiment_type = 'SpeechContext';
else
    error('did not find the expected fields in outfile')
end
fprintf('\nthis appears to be a %s experiment', experiment_type)

switch experiment_type
    case 'tuningcurve'
        for i=1:P.numoutfiles
            Freqs(i,:)=Out_components(i).out.freqs;
            Amps(i,:)=Out_components(i).out.amps;
            Durs(i,:)=Out_components(i).out.durs;
        end
        % situation 1: all outfiles have the same params
        if size(unique(Freqs, 'rows'), 1)~=1
            error('frequencies of outfiles don''t match')
        end
        if size(unique(Amps, 'rows'), 1)~=1
            error('amps of outfiles don''t match')
        end
        
        %now that we've verified that all parameters are the same, we can just use one of them
        Out.freqs=Out_components(1).out.freqs;
        Out.numfreqs=Out_components(1).out.numfreqs;
        Out.amps=Out_components(1).out.amps;
        Out.numamps=Out_components(1).out.numamps;
        Out.durs=Out_components(1).out.durs;
        Out.numdurs=Out_components(1).out.numdurs;
        
        Out.samprate=Out_components(1).out.samprate;
        Out.IL=Out_components(1).out.IL;
        Out.xlimits=Out_components(1).out.xlimits;
%         if size(Out_components(1).out.nrepsOFF) ~= [Out.numfreqs Out.numamps]
%             error('nreps is not numfreqs x numamps')
%         end
        
        Out.nreps=Out_components(1).out.nreps;
        Out.nrepsON=Out_components(1).out.nrepsON;
        Out.nrepsOFF=Out_components(1).out.nrepsOFF;
        for i=2:P.numoutfiles
            Out.nreps = Out.nreps + Out_components(i).out.nreps;
            Out.nrepsON = Out.nrepsON + Out_components(i).out.nrepsON;
            Out.nrepsOFF = Out.nrepsOFF + Out_components(i).out.nrepsOFF;
        end
        
        %pre-allocate
        sz=size(Out_components(i).out.M1OFF);
        sz(end-1)=max(Out.nrepsOFF(:));
        Out.M1OFF=nan(sz);
        Out.M1OFFLaser=nan(sz);
        Out.M1OFFStim=nan(sz);
        
        for i=1:P.numoutfiles
            nr=max(Out_components(i).out.nrepsOFF(:));
            start=1+(i-1)*nr;
            stop=nr*i;
            Out.M1OFF(:,:,:,start:stop,:)=Out_components(i).out.M1OFF;
            Out.M1OFFLaser(:,:,:,start:stop,:)=Out_components(i).out.M1OFFLaser;
            Out.M1OFFStim(:,:,:,start:stop,:)=Out_components(i).out.M1OFFStim;
        end
        Out.mM1OFF(:,:,1:Out.numdurs,:)= nanmean(Out.M1OFF, 4);
        Out.mM1OFFLaser(:,:,1:Out.numdurs,:)=nanmean(Out.M1OFFLaser, 4);
        Out.mM1OFFStim(:,:,1:Out.numdurs,:)=nanmean(Out.M1OFFStim, 4);
        
        sz=size(Out_components(i).out.M1ON);
        sz(end-1)=max(Out.nrepsON(:));
        Out.M1ON=nan(sz);
        Out.M1ONLaser=nan(sz);
        Out.M1ONStim=nan(sz);
        
        for i=1:P.numoutfiles
            nr=max(Out_components(i).out.nrepsON(:));
            start=1+(i-1)*nr;
            stop=nr*i;
            fprintf('\n%d-%d', start, stop)
            Out.M1ON(:,:,:,start:stop,:)=Out_components(i).out.M1ON;
            Out.M1ONLaser(:,:,:,start:stop,:)=Out_components(i).out.M1ONLaser;
            Out.M1ONStim(:,:,:,start:stop,:)=Out_components(i).out.M1ONStim;
        end
        Out.mM1ON(:,:,1:Out.numdurs,:)=nanmean(Out.M1ON, 4);
        Out.mM1ONLaser(:,:,1:Out.numdurs,:)=nanmean(Out.M1ONLaser, 4);
        Out.mM1ONStim(:,:,1:Out.numdurs,:)=nanmean(Out.M1ONStim, 4);
        combinedoutfilename='out_combined.mat';
        
    case 'GPIASbehavior'
        for i=1:P.numoutfiles
            Gapdurs(i,:)=Out_components(i).out.gapdurs;
        end
        if size(unique(Gapdurs, 'rows'), 1)~=1
            error('gapdurs of outfiles don''t match')
        end
        
        %now that we've verified that all parameters are the same, we can just use one of them
        Out.gapdurs=Out_components(1).out.gapdurs;
        Out.numgapdurs=Out_components(1).out.numgapdurs;
        Out.pulseamps=Out_components(1).out.pulseamps;
        Out.numpulseamps=Out_components(1).out.numpulseamps;
        Out.gapdelay=Out_components(1).out.gapdelay;
        Out.soa=Out_components(1).out.soa;
        Out.isi=Out_components(1).out.isi;
        Out.soaflag=Out_components(1).out.soaflag;
        Out.mouseID=Out_components(1).out.mouseID;
        Out.outfilename=Out_components(1).out.outfilename;
        Out.samprate=Out_components(1).out.samprate;
        Out.IL=Out_components(1).out.IL;
        Out.xlimits=Out_components(1).out.xlimits;
        if size(Out_components(1).out.nrepsOFF)~=[Out.numgapdurs 1]
            error('nreps is not numgapdurs x 1')
        end
        
        Out.nrepsON=Out_components(1).out.nrepsON;
        Out.nrepsOFF=Out_components(1).out.nrepsOFF;
        for i=2:P.numoutfiles
            Out.nrepsON = Out.nrepsON + Out_components(i).out.nrepsON;
            Out.nrepsOFF = Out.nrepsOFF + Out_components(i).out.nrepsOFF;
        end
        
        %pre-allocate
        sz=size(Out_components(i).out.M1OFF);
        sz(end-1)=max(Out.nrepsOFF(:));
        Out.M1OFF=nan(sz);
        Out.M1OFFLaser=nan(sz);
        Out.M1OFFstim=nan(sz);
        Out.PeakOFF=nan(size(Out_components(i).out.PeakOFF));
        for i=1:P.numoutfiles
            nr=max(Out_components(i).out.nrepsOFF(:));
            start=1+(i-1)*nr;
            stop=nr*i;
            Out.M1OFF(:,:,start:stop,:)=Out_components(i).out.M1OFF;
            Out.M1OFFstim(:,:,start:stop,:)=Out_components(i).out.M1OFFstim;
            Out.PeakOFF(:,:,start:stop)=Out_components(i).out.PeakOFF;
            Out.all_percentGPIAS_OFF(i,:)=Out_components(i).out.percentGPIAS_OFF;
            Out.all_pOFF(i,:)=Out_components(i).out.pOFF;
        end
        Out.mM1OFF(:,1:Out.numpulseamps,:)=nanmean(Out.M1OFF, 3);
        Out.mM1OFFstim(:,1:Out.numpulseamps,:)=nanmean(Out.M1OFFstim, 3);
        Out.mPeakOFF=nanmean(Out.PeakOFF, 3);
        Out.semPeakOFF=nanstd(Out.PeakOFF, 0, 3)/sqrt(length(Out.PeakOFF(:,3)));
        Out.percentGPIAS_OFF=nanmean(Out.all_percentGPIAS_OFF, 1);
        Out.pOFF=nanmean(Out.all_pOFF, 1); %pretty sure it's not kosher to average p-values, we should re-generate them if we will use them
        
        sz=size(Out_components(i).out.M1ON);
        sz(end-1)=max(Out.nrepsON(:));
        Out.M1ON=nan(sz);
        Out.M1ONLaser=nan(sz);
        Out.M1ONstim=nan(sz);
        Out.PeakON=nan(size(Out_components(i).out.PeakON));
        
        for i=1:P.numoutfiles
            nr=max(Out_components(i).out.nrepsON(:));
            start=1+(i-1)*nr;
            stop=nr*i;
            fprintf('\n%d-%d', start, stop)
            Out.M1ON(:,:,:,start:stop,:)=Out_components(i).out.M1ON;
            Out.M1ONstim(:,:,:,start:stop,:)=Out_components(i).out.M1ONstim;
            Out.PeakON(:,:,start:stop)=Out_components(i).out.PeakON;
            if ~isempty(Out_components(i).out.percentGPIAS_ON)
                Out.all_percentGPIAS_ON(i,:)=Out_components(i).out.percentGPIAS_ON;
                Out.all_pON(i,:)=Out_components(i).out.pON;
            else
                Out.all_percentGPIAS_ON=[];
                Out.all_pON=[];
            end
        end
        Out.mM1ON(:,:,:)=nanmean(Out.M1ON, 3);
        Out.mM1ONstim(:,:,:)=nanmean(Out.M1ONstim, 3);
        Out.mPeakON=nanmean(Out.PeakON, 3);
        if ~isempty(Out.PeakON)
            Out.semPeakON=nanstd(Out.PeakON, 0, 3)/sqrt(length(Out.PeakON(:,3)));
        else
            Out.semPeakON=nan;
            Out.mPeakON=[];
        end
        Out.percentGPIAS_ON=nanmean(Out.all_percentGPIAS_ON, 1);
        Out.pON=nanmean(Out.all_pON, 1); %pretty sure it's not kosher to average p-values, we should re-generate them if we will use them
        
        combinedoutfilename='outGPIAS_Behavior_combined.mat';
        
    case 'SpeechContext'
        for i=1:P.numoutfiles
            amps(i,:)=Out_components(i).out.amps;
            durs(i,:)=Out_components(i).out.durs;
            cellnum(i,:)=Out_components(i).out.cell;
        end
        % situation 1: all outfiles have the same params
        if size(unique(amps, 'rows'), 1)~=1
            error('amps of outfiles don''t match')
        end
        if size(unique(cellnum, 'rows'), 1)~=1
            error('cell ids of outfiles don''t match')
        end
        
        %now that we've verified that all parameters are the same, we can just use one of them
        Out.IL=Out_components(1).out.IL;
        Out.amps=Out_components(1).out.amps;
        Out.numamps=Out_components(1).out.numamps;
        Out.durs=Out_components(1).out.durs;
        Out.numdurs=Out_components(1).out.numdurs;
        Out.Nclusters=Out_components(1).out.Nclusters;
        Out.tetrode=Out_components(1).out.tetrode;
        Out.channel=Out_components(1).out.channel;
        Out.cluster=Out_components(1).out.cluster;
        Out.cell=Out_components(1).out.cell;
%         Out.sourcefiles=Out_components(1).out.sourcefiles;
        Out.numsourcefiles=Out_components(1).out.numsourcefiles;
%         Out.SilentSoundOFFStim=Out_components(1).out.SilentSoundOFFStim;
%         Out.SilentSoundOFFLaser=Out_components(1).out.SilentSoundOFFLaser;
%         Out.LaserRecorded=Out_components(1).out.LaserRecorded;
%         Out.StimRecorded=Out_components(1).out.StimRecorded;
        
%         Out.nb.user=Out_components(1).out.nb.user;
%         Out.nb.mouseID=Out_components(1).out.nb.mouseID;
%         Out.nb.Depth=Out_components(1).out.nb.Depth;
%         Out.nb.mouseDOB=Out_components(1).out.nb.mouseDOB;
%         Out.nb.mouseSex=Out_components(1).out.nb.mouseSex;
%         Out.nb.mouseGenotype=Out_components(1).out.nb.mouseGenotype;
        
        
%         Out.samprate=Out_components(1).out.samprate;
%         Out.IL=Out_components(1).out.IL;
%         Out.xlimits=Out_components(1).out.xlimits;
        Out.nreps=Out_components(1).out.nreps;
        Out.nrepsON=Out_components(1).out.nrepsON;
        Out.nrepsOFF=Out_components(1).out.nrepsOFF;
        Out.nreps_ssON=Out_components(1).out.nreps_ssON;
        Out.nreps_ssOFF=Out_components(1).out.nreps_ssOFF;
%         if size(Out_components(1).out.nrepsOFF, 2) ~= Out.numdurs || size(Out_components(1).out.nrepsOFF, 2) ~= Out.numamps
%             error('nreps is not numdurs x numamps')
%         end
        
%         Out.datadir=Out_components(1).out.datadir;
%         Out.datadirs{1}=Out_components(1).out.datadir;
%         Out.nb.notes = Out_components(1).out.nb.notes;
%         Out.nb.datapath = Out_components(1).out.nb.datapath;
%         Out.nb.datapaths{1} = Out_components(1).out.nb.datapath;
%         Out.nb.activedir = Out_components(1).out.nb.activedir;
%         Out.nb.activedirs{1} = Out_components(1).out.nb.activedir;
%         Out.stimlog=Out_components(1).out.stimlog;
%         Out.stimlogs{1}=Out_components(1).out.stimlog;
        for i=2:P.numoutfiles
            Out.nreps = Out.nreps + Out_components(i).out.nreps;
%             Out.nrepsON = Out.nrepsON + Out_components(i).out.nrepsON;
%             Out.nrepsOFF = Out.nrepsOFF + Out_components(i).out.nrepsOFF;
            
%             Out.nreps_ssON = Out.nreps_ssON + Out_components(i).out.nreps_ssON;
%             Out.nreps_ssOFF = Out.nreps_ssOFF + Out_components(i).out.nreps_ssOFF;
%             Out.datadirs{i} =  Out_components(i).out.datadir;
%             Out.nb.notes = Out.nb.notes + Out_components(i).out.nb.notes;
%             Out.nb.datapaths{i} =  Out_components(i).out.nb.datapath;
%             Out.nb.activedirs{i} = Out.nb.activedir + Out_components(i).out.nb.activedir;
%             Out.stimlogs{i} =  Out_components(i).out.stimlog;
            
        end
        
        %pre-allocate
        if exist('Out_components(i).out.M1OFFLaser', 'var') == 1
            sz=size(Out_components(i).out.M1OFFLaser);
            Out.M1OFFLaser=nan(sz);
            Out.M1OFFStim=nan(sz);
        else
        end
        % Added conditional statement for preallocation because some outfiles do not have out.M1OFFLaser (?) - SFM 8/6/21
        
        r_cum=0;
        for i=1:P.numoutfiles
            nr=max(Out_components(i).out.nrepsOFF(:));
         
                r_blockstart=r_cum; %start rep -1 for this block 
                r_blocks(i,:)=[r_blockstart+1, r_blockstart+nr]; %start and stop reps for this block
            for r=1:nr %indexing reps
                r_cum=r_cum+1;
                for aindex=[Out.numamps:-1:1]
                    for sourcefileindex=1:Out.numsourcefiles
                        for dindex=1:Out.numdurs
                            st=Out_components(i).out.M1OFF(sourcefileindex, aindex, dindex, r).spiketimes;
                            Out.M1OFF(sourcefileindex, aindex, dindex, r_cum).spiketimes=st;
                             
%                             fprintf('\ni%d a%d d%d s%d nr%d r%d r_cum%d r_blocks%d-%d', i, aindex, dindex, sourcefileindex, nr, r, r_cum, r_blocks(i,:))
                        end
                    end
                end
            end
             Out.M1OFFLaser(:,:,:,r_blocks(i,1):r_blocks(i, 2),:)=Out_components(i).out.M1OFFLaser;
             Out.M1OFFStim(:,:,:,r_blocks(i,1):r_blocks(i, 2),:)=Out_components(i).out.M1OFFStim;
        end
        Out.mM1OFFLaser(:,:,1:Out.numdurs,:)=nanmean(Out.M1OFFLaser, 4);
        Out.mM1OFFStim(:,:, 1:Out.numdurs,:)=nanmean(Out.M1OFFStim, 4);
        
        
        %repeat for ON
%         sz=size(Out_components(i).out.M1ONLaser);
%         Out.M1ONLaser=nan(sz);
%         Out.M1ONStim=nan(sz);
%         r_cum=0;
%         for i=1:P.numoutfiles
%             nr=max(Out_components(i).out.nrepsON(:));
%         
%                 r_blockstart=r_cum; %start rep -1 for this block 
%                 r_blocks(i,:)=[r_blockstart+1, r_blockstart+nr]; %start and stop reps for this block
%             for r=1:nr %indexing reps
%                 r_cum=r_cum+1;
%                 for aindex=[Out.numamps:-1:1]
%                     for sourcefileindex=1:Out.numsourcefiles
%                         for dindex=1:Out.numdurs
%                             st=Out_components(i).out.M1ON(sourcefileindex, aindex, dindex, r).spiketimes;
%                             Out.M1ON(sourcefileindex, aindex, dindex, r_cum).spiketimes=st;
% 
%                         end
%                     end
%                 end
%             end
%             Out.M1ONLaser(:,:,:,r_blocks(i,1):r_blocks(i, 2),:)=Out_components(i).out.M1ONLaser;
%             Out.M1ONStim(:,:,:,r_blocks(i,1):r_blocks(i, 2),:)=Out_components(i).out.M1ONStim;
%         end
%         if ~isempty(Out.M1ONLaser) %don't bother if there are no laser trials
%             Out.mM1ONLaser(:,:,1:Out.numdurs,:)=nanmean(Out.M1ONLaser, 4);
%             Out.mM1ONStim(:,:,1:Out.numdurs,:)=nanmean(Out.M1ONStim, 4);
%         end
        
        % Accumulate spiketimes across trials, for psth...
        for dindex=1:Out.numdurs; % Hardcoded.
            for aindex=[Out.numamps:-1:1]
                for sourcefileidx=1:Out.numsourcefiles
                    
%                     % on
%                     spiketimesON=[];
%                     for rep=1:Out.nrepsON(sourcefileidx, aindex, dindex)
%                         spiketimesON=[spiketimesON Out.M1ON(sourcefileidx, aindex, dindex, rep).spiketimes];
%                     end
%                     Out.mM1ON(sourcefileidx, aindex, dindex).spiketimes=spiketimesON;
                    
                    % off
                    spiketimesOFF=[];
                    for rep=1:Out.nrepsOFF(sourcefileidx, aindex, dindex)
                        spiketimesOFF=[spiketimesOFF Out.M1OFF(sourcefileidx, aindex, dindex, rep).spiketimes];
                    end
                    Out.mM1OFF(sourcefileidx, aindex, dindex).spiketimes=spiketimesOFF;
                end
            end
        end
        
        combinedoutfilename=sprintf('outPSTH_combined_ch%dc%d.mat', Out.tetrode, Out.cell);
        
        
end %switch experiment type


cd(targetdir)
waitbar(.9, wb, 'saving combined outfile...')
fprintf('\nsaving combined outfile...')
out=Out;
save(combinedoutfilename, 'out', '-v7.3')
close(wb)
fprintf('\nsaved %s in directory %s', combinedoutfilename, targetdir)
fprintf('\n')


%include a field in the outfile saying which outfiles are in it




function BrowseAndAdd
global P
[filename, pathname] = uigetfile('out*.mat', 'select an outfile...')
if filename
    fullfilename= fullfile(pathname, filename);
    P.numoutfiles=P.numoutfiles+1;
    P.outfilelist{P.numoutfiles}=fullfilename;
    P.dirlist{P.numoutfiles}=pathname;
    
    h = P.Messageh;
    str=sprintf('%s\n\n', P.outfilelist{:})
    set(h,'string',str,'backgroundcolor',[1 1 1], 'foregroundcolor', [0 0 0]);
    
end

function WriteToCellList(d, ClustQual, PVcell)
global P
str='';
for s=1:length(d)
    str=sprintf('%s\n\ncell',str);
    str=sprintf('%s\nPath: %s',str, pwd);
    str=sprintf('%s\nFilename: %s',str,  d(s).name);
    str=sprintf('%s\nCluster Quality: %d',str,  ClustQual(s));
    str=sprintf('%s\nPV cell: %d',str,  PVcell(s));
    
    wd=pwd;
    [dd, ff]=fileparts(fullfile(pwd,(d(s).name)));
    cd(dd)
    try
        nb=load('notebook.mat');
        %here we can write out any additional stimulus or notebook info we want
        str=sprintf('%s\n%s',str, nb.stimlog(1).protocol_name);
        str=sprintf('%s\n', str);
    end
    cd(wd)
end
response=questdlg(str, 'Write selected cells to file?', 'Write', 'Cancel', 'Write');
switch response
    case 'Write'
        fid=fopen(P.TargetCellList, 'a'); %absolute path
        fprintf(fid, '%s', str);
        fclose(fid);
end

% Return the name of this function.
function out = me
out = mfilename;

function InitializeGUI

global P
if isfield(P, 'fig')
    try
        close(P.fig)
    end
end
fig = figure;
P.fig=fig;
set(fig,'visible','off');
set(fig,'visible','off','numbertitle','off','name','Outfile Combiner',...
    'doublebuffer','on','menubar','none','closerequestfcn','Outfile_Combiner(''Close'')')
height=500; width=350; e=2; H=e;
w=200; h=25;
set(fig,'pos',[1000 200         width         height],'visible','on');
P.numoutfiles=0;

H=400;
P.Messageh=uicontrol(fig,'tag','message','style','edit','fontweight','bold','units','pixels',...
    'enable','inact','horiz','left','Max', 8, 'pos',[e  10 width H ]);



%Combine them button
H=H+h+e;
uicontrol('parent',fig,'string','Combine them','tag','CreateNewOutfile','units','pixels',...
    'position',[e H w h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%help button
uicontrol('parent',fig,'string','help','tag','help','units','pixels',...
    'position',[e+w+20 H 50 h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback', 'help(mfilename)');

%reset button
uicontrol('parent',fig,'string','reset','tag','reset','units','pixels',...
    'position',[e+w+20 H+h+e 50 h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback', me);

%Browse and Add button
H=H+1*h+e;
P.BrowseAndAddh=uicontrol('parent',fig,'string','Browse and Add','tag','BrowseAndAdd','units','pixels',...
    'position',[e H w h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

