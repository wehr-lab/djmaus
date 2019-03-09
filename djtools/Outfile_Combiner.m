function Outfile_Combiner(varargin)

%simple GUI to pick multiple outfiles and combine them
%developed for ABRs (LFP tuning curves) so that we can split data
%acquisition across multiple sessions

global P

if nargin > 0
    action = varargin{1};
else
    action = get(gcbo,'tag');
end
if isempty(action)
    action='Init';
end


switch action
    case 'Init'
        InitializeGUI
    case 'Close'
        delete(P.fig)
        clear global P
    case 'CreateNewOutfile'
        CreateNewOutfile
    case 'BrowseAndAdd'
        BrowseAndAdd
end

% Return the name of this function.
function out = me
out = mfilename;

function CreateNewOutfile
global P
sortedoutfilelist=sort (P.outfilelist);
sorteddirlist=sort (P.dirlist);
targetdir=sorteddirlist{1};

for i=1:P.numoutfiles
    cd( sorteddirlist{i});
    Out_components(i)=load(sortedoutfilelist{i});
end
for i=1:P.numoutfiles
    Freqs(i,:)=Out_components(i).out.freqs;
    Amps(i,:)=Out_components(i).out.amps;
    Durs(i,:)=Out_components(i).out.durs;
end
if size(unique(Freqs, 'rows'), 1)~=1
    error('frequencies of outfiles don''t match')
end
if size(unique(Amps, 'rows'), 1)~=1
    error('amps of outfiles don''t match')
end

%now that we've verified that all parameters are the same, we can just use one of them
Out.freqs=Out_components(1).out.freqs
Out.amps=Out_components(1).out.amps
Out.durs=Out_components(1).out.durs

%     nreps=out.nreps;
%     numfreqs=out.numfreqs;
%     numamps=out.numamps;
%     numdurs=out.numdurs;
%     samprate=out.samprate; %in Hz
%     M1=out.M1;
%     M1stim=out.M1stim;
%     mM1=out.mM1;
%     mM1ON=out.mM1ON;
%     mM1OFF=out.mM1OFF;
%     M1OFF=out.M1OFF;
%     mM1ONLaser=out.mM1ONLaser;
%     M1ONLaser=out.M1ONLaser;
%     mM1OFFLaser=out.mM1OFFLaser;
%     mM1ONStim=out.mM1ONStim;
%     mM1OFFStim=out.mM1OFFStim;
%     nrepsON=out.nrepsON;
%     nrepsOFF=out.nrepsOFF;
%     IL=out.IL;


cd(targetdir)
combinedoutfilename='out_combined.mat';
save combinedoutfilename Out

keyboard

%include a field in the outfile saying which outfiles are in it


% function SelectExistingCellList
% global P
% [fname, path] = uigetfile('*.txt', 'Select cell list');
% if fname
%     P.TargetCellList=fullfile(path, fname);
%     set(P.TargetCellListDisplay, 'string', {'cell list:',path, fname});
%     set([P.BrowseAndAddh P.AddCurrentDirh], 'enable', 'on')
% end

% function AddCurrentDir
% global P
%     d=dir('out*.mat');
% if isempty(d)
%     h = errordlg('no outfiles found in this directory.', 'no outfiles');
% else
%     SelectOutfiles(d)
% end
%
% function SelectOutfiles(d)
%     [selection, ok, ClustQual, PVcell]=myCellListDlg(d);
%
% if ok
% end


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
height=220; width=350; e=2; H=e;
w=200; h=25;
set(fig,'pos',[1000 800         width         height],'visible','on');
P.numoutfiles=0;

P.Messageh=uicontrol(fig,'tag','message','style','edit','fontweight','bold','units','pixels',...
    'enable','inact','horiz','left','Max', 8, 'pos',[e  H width 5*h ]);


% %TargetCellList display
% P.TargetCellListDisplay= uicontrol('parent',fig,'string','','tag','OutfileListDisplay','units','pixels',...
%     'position',[e H width-e 2*h],'enable','on',...
%     'fontweight','bold','horiz', 'left',...
%     'style','text');

%Combine them button
H=H+5*h+e;
uicontrol('parent',fig,'string','Combine them','tag','CreateNewOutfile','units','pixels',...
    'position',[e H w h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%SelectExistingCellList button
% H=H+1*h+e;
% uicontrol('parent',fig,'string','Select Existing Cell List','tag','SelectExistingCellList','units','pixels',...
%     'position',[e H w h],'enable','on',...
%     'fontweight','bold',...
%     'style','pushbutton','callback',[me ';']);

%Browse and Add button
H=H+2*h+e;
P.BrowseAndAddh=uicontrol('parent',fig,'string','Browse and Add','tag','BrowseAndAdd','units','pixels',...
    'position',[e H w h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

% %AddCurrentDir button
% H=H+h+e;
% P.AddCurrentDirh=uicontrol('parent',fig,'string','Add Current Dir','tag','AddCurrentDir','units','pixels',...
%     'position',[e H w h],'enable','off',...
%     'fontweight','bold',...
%     'style','pushbutton','callback',[me ';']);

% %recursive checkbox
% H=H+h+e;
% P.recursive=uicontrol('parent',fig,'string','Recursive Scan','tag','Recursive','units','pixels',...
%     'position',[e H w h],'enable','on',...
%     'style','checkbox');


% function doOK(src, evt, IncludeCellcheckbox, sl, pvcheckbox)
%
% for i=1:length(sl)
%     IncludeCell(i)=get(IncludeCellcheckbox(i), 'value');
%     ClustQual(i)=get(sl(i), 'value');
%     PVcell(i)=get(pvcheckbox(i), 'value');
% end
% delete(gcbf);
% global P
% P.ok=1;
% P.IncludeCell=IncludeCell;
% P.ClustQual=ClustQual;
% P.PVcell=PVcell;


% function doCancel(src, evt)
% % ad.selection = [];
% % ad.ClustQual=[];
% % ad.PVcell=[];
% % setappdata(0,'ListDialogAppData__',ad)
% delete(gcbf);
% global P
% P.ok=0;



%
% function [selection, ok, ClustQual, PVcell]=myCellListDlg(d)
% global P
%
% fig=figure;
% set(fig, 'pos',[800 150 500 800] )
%
% selection=[];
% ClustQual=[];
% PVcell=[];
% P.ok=0;
%
% btn_wid=50;
% btn_ht=25;
%
% %     top=ffs+uh+4*fus+(smode==2)*(fus+uh)+listsize(2);
% fontsize=12;
% linesize = fontsize*1.4;  % height extent per line of uicontrol text (approx)
% linespacing=linesize*1.5;
% pos=get(fig, 'pos');
% width=pos(3);
% top=pos(4);
% col1width=150; %cell name
% col2width=50; %include cell checkbox
% col3width=120; %cluster quality slider
% col4width=50; %cluster quality numeric indicator
% col5width=50; %pv cell checkbox
% sliderwidth=col3width;
%
% uicontrol('style', 'text', 'pos', [0, top-linesize, col1width, linesize], 'string', 'cell', ...
%     'horizontalalignment', 'center', 'fontsize', fontsize)
% uicontrol('style', 'text', 'pos', [col1width, top-linesize, col2width, linesize], 'string', 'include?', ...
%     'horizontalalignment', 'center', 'fontsize', fontsize)
% uicontrol('style', 'text', 'pos', [col1width+col2width, top-linesize, sliderwidth, linesize],...
%     'string', 'cluster quality', 'fontsize', fontsize)
% uicontrol('style', 'text', 'pos', [col1width+col2width+col3width+col4width, top-linesize, 50, linesize], ...
%     'string', 'PV cell?', 'fontsize', fontsize, 'horizontalalign', 'left')
%
% for i=1:length(d)
%     if ~mod(i,3)
%         u=uipanel('units', 'pixels', 'pos',[ 0,top-(i+2)*linespacing-3, width, 2]);
%     end
%     cellstr(i)=uicontrol('style', 'text', 'pos',[2, top-(i+2)*linespacing, col1width, linesize],...
%         'horizontalalignment', 'right','fontsize', fontsize,'string', d(i).name); %cell name
%     IncludeCellcheckbox(i)=uicontrol('style', 'checkbox', 'pos', ...
%         [col1width, top-(i+2)*linespacing, 30, linesize], 'value', 1); %include cell checkbox
%     slstr(i)=uicontrol('style', 'text', 'pos', ...
%         [col1width+col2width+col3width+2, top-(i+2)*linespacing, 10, linesize], 'string', 0); %cluster quality numeric indicator
%     sl(i)=uicontrol('style', 'slider', 'pos', [col1width+col2width, top-(i+2)*linespacing, sliderwidth, linesize], ...
%         'min', 0, 'max', 5, 'sliderstep', [.2 .2], 'value', 0, 'callback', {@doSlider, slstr(i)} ); %cluster quality slider
%     pvcheckbox(i)=uicontrol('style', 'checkbox', 'pos', ...
%         [col1width+col2width+col3width+col4width, top-(i+2)*linespacing, 30, linesize]); %pvcheckbox
% end
%
% ok_btn = uicontrol('Style','pushbutton',...
%     'String','OK',...
%     'Position',[10 10 btn_wid btn_ht],...
%     'Tag','ok_btn',...
%     'Callback',{@doOK,IncludeCellcheckbox, sl, pvcheckbox });
%
% cancel_btn = uicontrol('Style','pushbutton',...
%     'String','Cancel',...
%     'Position',[20+btn_wid 10 btn_wid btn_ht],...
%     'Tag','cancel_btn',...
%     'Callback',@doCancel);
%
% uiwait(fig)
% ok=P.ok;
% if ok
%     selection=P.IncludeCell;
%     ClustQual=P.ClustQual;
%     PVcell=P.PVcell;
% else
%     selection=[];
%     ClustQual=[];
%     PVcell=[];
% end
