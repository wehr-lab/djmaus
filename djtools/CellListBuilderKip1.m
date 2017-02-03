function CellListBuilder(varargin)
% button to add current directory to cell list
% will automatically scan for clustered tetrode data
% % gives you a checklist of cells in that directory so you can keep or discard
% another button that opens a dialog box so you can browse to the desired directory
% add a checkbox for recursive scan
%
% what does it do? It creates/appends a machine readable cell list
% this will be a text file with each cell, stimuli presented, other
% stuff.
% It would be nice to have it formatted in such a way that it can be used
% to analyze data.
%
%to append to an existing cell list, we could have a dropdown menu of
%existing cell lists. Not sure where this would live.
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
if get(P.recursive, 'value')
    d=rdir('*/**.t');
else
    d=dir('*.t');
end
if isempty(d)
    h = errordlg('no clustered cells found in this directory.', 'no cells');
else
    SelectCells(d)
end

function SelectCells(d)
listsize=25;
if length(d)<=listsize
    [selection, ok, ClustQual, PVcell, PINP_test, WNtc_test, GPtc_test, GPbehav_test]=myCellListDlg(d);
else %break up into chunks
    selection=[];  ClustQual=[]; PVcell=[];
    uiwait(msgbox(sprintf('there are %d cells in this directory, so breaking into %d chunks of %d', length(d), ceil(length(d)/listsize), listsize), 'modal'));
    for i=1:ceil(length(d)/listsize)
        range=1+(i-1)*listsize:i*listsize;
        range=range(range<length(d));
        [selectionchunk, ok, ClustQualchunk, PVcellchunk]=myCellListDlg(d(range));
        if ~ok break, end
        selection=[selection selectionchunk];  ClustQual=[ClustQual ClustQualchunk]; PVcell=[PVcell PVcellchunk];
    end
end
if ok
    WriteToCellList(d(find(selection)), ClustQual(find(selection)), PVcell(find(selection)), ...
        PINP_test(find(selection)), WNtc_test(find(selection)), GPtc_test(find(selection)), GPbehav_test(find(selection)))
end


function BrowseAndAdd
d = uigetdir(pwd, 'choose directory to scan')
if d
    cd(d)
    AddCurrentDir
end

function WriteToCellList(d, ClustQual, PVcell, PINP_test, WNtc_test, GPtc_test, GPbehav_test)
global P
str='';
for s=1:length(d)
    str=sprintf('%s\n\ncell',str);
    str=sprintf('%s\nPath: %s',str, pwd);
    str=sprintf('%s\nFilename: %s',str,  d(s).name);
    str=sprintf('%s\nCluster Quality: %d',str,  ClustQual(s));
    str=sprintf('%s\nPV cell: %d',str,  PVcell(s));
    
    if PINP_test(s)
        str=sprintf('%s\PINP test',str);
    elseif WNtc_test(s)
        str=sprintf('%s\WN tc test',str);
    elseif GPtc_test(s)
        str=sprintf('%s\GAP tc test',str);
    elseif GPbehav_test(s)
        str=sprintf('%s\Gap behavior test',str);
    end
    
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
set(fig,'visible','off','numbertitle','off','name','cell list builder',...
    'doublebuffer','on','menubar','none','closerequestfcn','CellListBuilder(''Close'')')
ScreenSize = get(0,'screensize');
height=220; width=350; e=2; H=e;
w=200; h=25;
set(fig,'pos',[round(ScreenSize(3)*.25) round(ScreenSize(4)-height*2) width height],'visible','on');


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

%recursive checkbox
H=H+h+e;
P.recursive=uicontrol('parent',fig,'string','Recursive Scan','tag','Recursive','units','pixels',...
    'position',[e H w h],'enable','on',...
    'style','checkbox');


function doOK(src, evt, IncludeCellcheckbox, sl, pv_checkbox, PINP_checkbox, WNtc_checkbox, GPtc_checkbox, GPbehav_checkbox )

for i=1:length(sl)
    IncludeCell(i)=get(IncludeCellcheckbox(i), 'value');
    ClustQual(i)=get(sl(i), 'value');
    PVcell(i)=get(pv_checkbox(i), 'value');
    PINP_test(i)=get(PINP_checkbox(i), 'value');
    WNtc_test(i)=get(WNtc_checkbox(i), 'value');
    GPtc_test(i)=get(GPtc_checkbox(i), 'value');
    GPbehav_test(i)=get(GPbehav_checkbox(i), 'value');
end
delete(gcbf);
global P
P.ok=1;
P.IncludeCell=IncludeCell;
P.ClustQual=ClustQual;
P.PVcell=PVcell;

P.PINP_test=PINP_test;
P.WNtc_test=WNtc_test;
P.GPtc_test=GPtc_test;
P.GPbehav_test=GPbehav_test;


function doCancel(src, evt)
% ad.selection = [];
% ad.ClustQual=[];
% ad.PVcell=[];
% setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);
global P
P.ok=0;

% slider callback
function doSlider(sliderh, evd, sliderstrh) %#ok
set(sliderstrh, 'string', get(sliderh, 'value'))



function [selection, ok, ClustQual, PVcell, PINP_test, WNtc_test, GPtc_test, GPbehav_test]=myCellListDlg(d)
global P

fig=figure;
ScreenSize = get(0,'screensize');
set(fig, 'pos',[ScreenSize(3)*.25 ScreenSize(4)-900 700 800] )

selection=[];
ClustQual=[];
PVcell=[];
PINP_test = [];
WNtc_test = [];
GPtc_test = [];
GPbehav_test = [];
P.ok=0;

btn_wid=50;
btn_ht=25;

%     top=ffs+uh+4*fus+(smode==2)*(fus+uh)+listsize(2);
fontsize=10;
linesize = fontsize*1.4;  % height extent per line of uicontrol text (approx)
linespacing=linesize*1.5;
pos=get(fig, 'pos');
width=pos(3);
top=pos(4);
col1width=150; %cell name
col2width=50; %include cell checkbox
col3width=120; %cluster quality slider
col4width=50; %cluster quality numeric indicator
col5width=50; %pv cell checkbox
colWidth = [150 50 120 50 50 50 50 50 50];
sliderwidth=col3width;

uicontrol('style', 'text', 'pos', [0, top-linesize, col1width, linesize], 'string', 'cell', ...
    'horizontalalignment', 'center', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [colWidth(1), top-linesize, col2width, linesize], 'string', 'include?', ...
    'horizontalalignment', 'center', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [sum(colWidth(1:2)), top-linesize, sliderwidth, linesize],...
    'string', 'cluster quality', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [sum(colWidth(1:4)), top-linesize, 50, linesize], ...
    'string', 'PV cell?', 'fontsize', fontsize, 'horizontalalign', 'left')

uicontrol('style', 'text', 'pos', [sum(colWidth(1:5)), top-linesize, 50, linesize], ...
    'string', 'PINP', 'fontsize', fontsize, 'horizontalalign', 'left')
uicontrol('style', 'text', 'pos', [sum(colWidth(1:6)), top-linesize, 50, linesize], ...
    'string', 'WNtc', 'fontsize', fontsize, 'horizontalalign', 'left')
uicontrol('style', 'text', 'pos', [sum(colWidth(1:7)), top-linesize, 50, linesize], ...
    'string', 'GPtc', 'fontsize', fontsize, 'horizontalalign', 'left')
uicontrol('style', 'text', 'pos', [sum(colWidth(1:8)), top-linesize, 50, linesize], ...
    'string', 'GPbehav', 'fontsize', fontsize, 'horizontalalign', 'left')

for i=1:length(d)
    if ~mod(i,3)
        u=uipanel('units', 'pixels', 'pos',[ 0,top-(i+2)*linespacing-3, width, 2]);
    end
    txt = d(i).name;
    txt = strrep(txt,'_simpleclust_',' Clust# ');
    txt = strrep(txt,'.t',' ');
    cellstr(i)=uicontrol('style', 'text', 'pos',[2, top-(i+2)*linespacing, col1width, linesize],...
        'horizontalalignment', 'right','fontsize', fontsize,'string', txt); %cell name
    IncludeCellcheckbox(i)=uicontrol('style', 'checkbox', 'pos', ...
        [sum(colWidth(1)), top-(i+2)*linespacing, 30, linesize]); %include cell checkbox
    slstr(i)=uicontrol('style', 'text', 'pos', ...
        [sum(colWidth(1:3))+2, top-(i+2)*linespacing, 10, linesize], 'string', 0); %cluster quality numeric indicator
    sl(i)=uicontrol('style', 'slider', 'pos', [sum(colWidth(1:2)), top-(i+2)*linespacing, sliderwidth, linesize], ...
        'min', 0, 'max', 5, 'sliderstep', [.2 .2], 'value', 0, 'callback', {@doSlider, slstr(i)} ); %cluster quality slider
    pv_checkbox(i)=uicontrol('style', 'checkbox', 'pos', ...
        [sum(colWidth(1:4)), top-(i+2)*linespacing, 30, linesize]); %pv_checkbox
    
    PINP_checkbox(i)=uicontrol('style', 'checkbox', 'pos', ...
        [sum(colWidth(1:5)), top-(i+2)*linespacing, 30, linesize]); %PINP_checkbox
    WNtc_checkbox(i)=uicontrol('style', 'checkbox', 'pos', ...
        [sum(colWidth(1:6)), top-(i+2)*linespacing, 30, linesize]); %WNTC_checkbox
    GPtc_checkbox(i)=uicontrol('style', 'checkbox', 'pos', ...
        [sum(colWidth(1:7)), top-(i+2)*linespacing, 30, linesize]); %GPTC_checkbox
    GPbehav_checkbox(i)=uicontrol('style', 'checkbox', 'pos', ...
        [sum(colWidth(1:8)), top-(i+2)*linespacing, 30, linesize]); %GPbehav_checkbox
    
end

ok_btn = uicontrol('Style','pushbutton',...
    'String','OK',...
    'Position',[10 10 btn_wid btn_ht],...
    'Tag','ok_btn',...
    'Callback',{@doOK,IncludeCellcheckbox, sl, pv_checkbox PINP_checkbox WNtc_checkbox GPtc_checkbox GPbehav_checkbox });

cancel_btn = uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Position',[20+btn_wid 10 btn_wid btn_ht],...
    'Tag','cancel_btn',...
    'Callback',@doCancel);

uiwait(fig)
ok=P.ok;
if ok
    selection=P.IncludeCell;
    ClustQual=P.ClustQual;
    PVcell=P.PVcell;
    PINP_test = P.PINP_test;
    WNtc_test = P.WNtc_test;
    GPtc_test = P.GPtc_test;
    GPbehav_test = P.GPbehav_test;
else
    selection=[];
    ClustQual=[];
    PVcell=[];
    PINP_test = [];
    WNtc_test = [];
    GPtc_test = [];
    GPbehav_test = [];
end
