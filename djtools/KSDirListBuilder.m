function KSDirListBuilder(varargin)
% A GUI for building a list of kilosort master data directories
% we don't add invidual cells, just the master data directory for each
% recording session.
% it has a button to add current directory to the list
% another button that opens a dialog box so you can browse to the desired directory
%
% what does it do? It creates/appends a machine readable cell list as a
% .txt file

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
    case 'CreateNewCellList'
        CreateNewCellList
    case 'SelectExistingCellList'
        SelectExistingCellList
    case 'AddCurrentDir'
        AddCurrentDir
    case 'BrowseAndAdd'
        BrowseAndAdd
    case 'TargetCellList'
        TargetCellList
end

% Return the name of this function.
function out = me
out = mfilename;

function CreateNewCellList
global P
[fname, path] = uiputfile('*.txt', 'Select cell list');
if fname
    P.TargetCellList=fullfile(path, fname);
    set(P.TargetCellListDisplay, 'string', {'cell list:',path, fname});
    set([P.BrowseAndAddh P.AddCurrentDirh], 'enable', 'on')
end

function SelectExistingCellList
global P
[fname, path] = uigetfile('*.txt', 'Select cell list');
if fname
    P.TargetCellList=fullfile(path, fname);
    set(P.TargetCellListDisplay, 'string', {'cell list:',path, fname});
    set([P.BrowseAndAddh P.AddCurrentDirh], 'enable', 'on')
end

function AddCurrentDir
global P
d=dir('dirs.mat');

if isempty(d)
    h = errordlg('this does not appear to be a kilosort data directory.', 'not a datadir?');
else
    load dirs
    if ~strcmp(dirs{1}, pwd)
        a = questdlg('This does not appear to be the master directory. Do you want to cd to the master dirrectory?', 'not master?');
        switch a
            case 'Yes'
                try
                    cd(dirs{1});
                    AddCurrentDir
                catch
                    errordlg('could not cd to master directory, sorry')
                end
        end
    else
        %it is a master dir so add it
        WriteToCellList(pwd)
    end
end


function BrowseAndAdd
d = uigetdir(pwd, 'choose directory to scan')
if d
    cd(d)
    AddCurrentDir
end

function WriteToCellList(d)
global P
    
        %here we can write out any additional stimulus or notebook info we want
%     wd=pwd;
%     [dd, ff]=fileparts(fullfile(pwd,(d(s).name)));
%     cd(dd)
%     try
%         nb=load('notebook.mat');
%         %here we can write out any additional stimulus or notebook info we want
%         str=sprintf('%s\n%s',str, nb.stimlog(1).protocol_name);
%         str=sprintf('%s\n', str);
%     end
%     cd(wd)
response=questdlg(str, 'Write selected cells to file?', 'Write', 'Cancel', 'Write');
switch response
    case 'Write'
        fid=fopen(P.TargetCellList, 'a'); %absolute path
        fprintf(fid, '\n\ndir');
        fprintf(fid, '\n%s', d);
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
set(fig,'visible','off','numbertitle','off','name','cell list builder',...
    'doublebuffer','on','menubar','none','closerequestfcn','CellListBuilder(''Close'')')
height=220; width=350; e=2; H=e;
w=200; h=25;
set(fig,'pos',[1000 800         width         height],'visible','on');


%TargetCellList display
P.TargetCellListDisplay= uicontrol('parent',fig,'string','','tag','TargetCellListDisplay','units','pixels',...
    'position',[e H width-e 2*h],'enable','on',...
    'fontweight','bold','horiz', 'left',...
    'style','text');

%CreateNewCellList button
H=H+2*h+e;
uicontrol('parent',fig,'string','Create New Cell List','tag','CreateNewCellList','units','pixels',...
    'position',[e H w h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%SelectExistingCellList button
H=H+1*h+e;
uicontrol('parent',fig,'string','Select Existing Cell List','tag','SelectExistingCellList','units','pixels',...
    'position',[e H w h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%Browse and Add button
H=H+2*h+e;
P.BrowseAndAddh=uicontrol('parent',fig,'string','Browse and Add','tag','BrowseAndAdd','units','pixels',...
    'position',[e H w h],'enable','off',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%AddCurrentDir button
H=H+h+e;
P.AddCurrentDirh=uicontrol('parent',fig,'string','Add Current Dir','tag','AddCurrentDir','units','pixels',...
    'position',[e H w h],'enable','off',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

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


function doCancel(src, evt)
% ad.selection = [];
% ad.ClustQual=[];
% ad.PVcell=[];
% setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);
global P
P.ok=0;

% % slider callback
% function doSlider(sliderh, evd, sliderstrh) %#ok
% set(sliderstrh, 'string', get(sliderh, 'value'))
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
%     u=uipanel('units', 'pixels', 'pos',[ 0,top-(i+2)*linespacing-3, width, 2]); 
%     end
%     cellstr(i)=uicontrol('style', 'text', 'pos',[2, top-(i+2)*linespacing, col1width, linesize],...
%          'horizontalalignment', 'right','fontsize', fontsize,'string', d(i).name); %cell name
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
%  uiwait(fig)
%  ok=P.ok;
%  if ok
%      selection=P.IncludeCell;
%      ClustQual=P.ClustQual;
%      PVcell=P.PVcell;
% else
%     selection=[];
%     ClustQual=[];
%     PVcell=[];
% end
