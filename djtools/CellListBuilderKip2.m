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
    [selection, ok, ClustQual, PVcell IDnumberList]=myCellListDlg(d);
else %break up into chunks
    selection=[];  ClustQual=[]; PVcell=[]; IDnumberList=[];
    uiwait(msgbox(sprintf('there are %d cells in this directory, so breaking into %d chunks of %d', length(d), ceil(length(d)/listsize), listsize), 'modal'));
    for i=1:ceil(length(d)/listsize)
        range=1+(i-1)*listsize:i*listsize;
        range=range(range<length(d));
        [selectionchunk, ok, ClustQualchunk, PVcellchunk IDnumberListchunk]=myCellListDlg(d(range));
        if ~ok break, end
        selection=[selection selectionchunk];  ClustQual=[ClustQual ClustQualchunk]; PVcell=[PVcell PVcellchunk]; IDnumberList=[IDnumberList IDnumberListchunk];
    end
end

if ok
    WriteToCellList(d(find(selection)), ClustQual(find(selection)), PVcell(find(selection)), IDnumberList(find(selection)))
end


function BrowseAndAdd
d = uigetdir(pwd, 'choose directory to scan')
if d
    cd(d)
    AddCurrentDir
end

function WriteToCellList(d, ClustQual, PVcell, IDnumberList)
global P
str='';
for s=1:length(d)
    str=sprintf('%s\n\ncell ID#',str);
    txt = ['         ' num2str(IDnumberList(s))];
    str=sprintf(['%s' txt(end-9:end)],str);
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


function doUpdateIDnumber(src, evt, IncludeCellcheckbox, IDnumbox, IDnumber)
for i = 1:length(IncludeCellcheckbox)
    if IncludeCellcheckbox(i).Value
        set(IDnumbox(i),'string',num2str(IDnumber));
        IDnumber = IDnumber+1;
    end
end
set(src,'value',0);

function doOK(src, evt, IncludeCellcheckbox, sl, pvcheckbox, IDnumberbox)

for i=1:length(sl)
    IncludeCell(i)=get(IncludeCellcheckbox(i), 'value');
    ClustQual(i)=get(sl(i), 'value');
    PVcell(i)=get(pvcheckbox(i), 'value');
    try
        IDnumberList(i)=str2num(get(IDnumberbox(i),'string'));
    catch
        IDnumberList(i)=nan;
    end
end
delete(gcbf);
global P
P.ok=1;
P.IncludeCell=IncludeCell;
P.ClustQual=ClustQual;
P.PVcell=PVcell;
P.IDnumberList=IDnumberList;


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

% IDnumbox callback
% function doIDnumbox(checkboxh, evd, IncludeCellcheckboxh, IDnumboxh, IDnumber) %#ok
% %set(sliderstrh, 'string', get(sliderh, 'value'))
% for i = 1:length(IDnumboxh)
%     try
%     temp(i) = str2num(IDnumboxh(i).String);
%     catch
%         temp(i)=nan;
%     end
% end
% IDnumber = max([IDnumber temp])+1;
% set(IDnumboxh(end), 'string', num2str(IDnumber));




function [selection, ok, ClustQual, PVcell, IDnumberList]=myCellListDlg(d)
global P

% find next cell# available
try
    fid=fopen(P.TargetCellList, 'r'); %absolute path
    txt = fscanf(fid, '%c');
    fclose(fid);
    ind = findstr(txt,'ID#');
    IDnumber = str2num(txt(ind(end)+3:ind(end)+12))+1;
catch
    IDnumber=1;
end

fig=figure;
ScreenSize = get(0,'screensize');
set(fig, 'pos',[ScreenSize(3)*.25 ScreenSize(4)-900 400 800] )

selection=[];
ClustQual=[];
PVcell=[];
IDnumberList=[];
P.ok=0;

btn_wid=50;
btn_ht=25;

%     top=ffs+uh+4*fus+(smode==2)*(fus+uh)+listsize(2);
fontsize=12;
linesize = fontsize*1.4;  % height extent per line of uicontrol text (approx)
linespacing=linesize*1.5;
pos=get(fig, 'pos');
width=pos(3);
top=pos(4);
% col1width=150; %cell name
% col2width=50; %include cell checkbox
% col3width=120; %cluster quality slider
% col4width=50; %cluster quality numeric indicator
% col5width=50; %pv cell checkbox
colWidth = [120 40 110 40 40 50 50 50 50];
%sliderwidth=colWidth(3);

Yline = top-linesize;
FN = pwd;
temp = strfind(FN,filesep);
uicontrol('style', 'text', 'pos', [0, Yline, 270, linesize], 'string', FN(temp(end)+1:end), ...
    'horizontalalignment', 'center', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [300, Yline, 100, linesize], 'string', ['nextID# ' num2str(IDnumber)], ...
    'horizontalalignment', 'center', 'fontsize', fontsize)

Yline = Yline-linesize;
uicontrol('style', 'text', 'pos', [0, Yline, colWidth(1), linesize], 'string', 'cell', ...
    'horizontalalignment', 'center', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [colWidth(1), Yline, colWidth(2), linesize], 'string', 'include?', ...
    'horizontalalignment', 'center', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [sum(colWidth(1:2)), Yline, colWidth(3), linesize],...
    'string', 'cluster quality', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [sum(colWidth(1:4)), Yline, colWidth(4), linesize], ...
    'string', 'PV?', 'fontsize', fontsize, 'horizontalalign', 'center') 

uicontrol('style', 'text', 'pos', [sum(colWidth(1:5)), Yline, colWidth(5), linesize], ...
    'string', 'ID#', 'fontsize', fontsize, 'horizontalalign', 'center') 

% check how many clusters for each tetrode
ttNum = nan(length(d)+1,1);
ttNum(end)=100;
for i=1:length(d)
    ttNum(i) = str2num(d(i).name(3));
end

Yline = Yline -linespacing;

for i=1:length(d)
    txt = d(i).name;
    txt = strrep(txt,'_simpleclust_',' Clust# ');
    txt = strrep(txt,'.t',' ');
    Yline = Yline -linespacing;
    cellstr(i)=uicontrol('style', 'text', 'pos',[2, Yline, colWidth(1), linesize],...
        'horizontalalignment', 'right','fontsize', fontsize,'string', txt); %cell name
    IncludeCellcheckbox(i)=uicontrol('style', 'checkbox', 'pos', [sum(colWidth(1))+10, Yline, 30, linesize]); %include cell checkbox
    slstr(i)=uicontrol('style', 'text', 'pos', ...
        [sum(colWidth(1:3))+2, top-(i+2)*linespacing, 10, linesize], 'string', 0); %cluster quality numeric indicator
    sl(i)=uicontrol('style', 'slider', 'pos', [sum(colWidth(1:2)), Yline, colWidth(3), linesize], ...
        'min', 0, 'max', 5, 'sliderstep', [.2 .2], 'value', 0, 'callback', {@doSlider, slstr(i)} ); %cluster quality slider
    pvcheckbox(i)=uicontrol('style', 'checkbox', 'pos', [sum(colWidth(1:4)), Yline, 30, linesize]); %pv_checkbox
        
    IDnumbox(i)=uicontrol('style', 'edit', 'pos', [sum(colWidth(1:5)), Yline, colWidth(5), linesize],'string',''); %IDnumbox
    
    if ttNum(i+1)~=ttNum(i)
         u=uipanel('units', 'pixels', 'pos',[ 0,Yline-3, width, 2]); 
    end
end

ok_btn = uicontrol('Style','pushbutton',...
    'String','OK',...
    'Position',[10 10 btn_wid btn_ht],...
    'Tag','ok_btn',...
    'Callback',{@doOK,IncludeCellcheckbox, sl, pvcheckbox IDnumbox});

cancel_btn = uicontrol('Style','pushbutton',...
    'String','Cancel',...
    'Position',[20+btn_wid 10 btn_wid btn_ht],...
    'Tag','cancel_btn',...
    'Callback',@doCancel);

updateIDnumber_btn = uicontrol('Style','pushbutton',...
    'String','updateID#s',...
    'Position',[30+2*btn_wid 10 btn_wid*2 btn_ht],...
    'Tag','updateIDnumber_btn',...
    'Callback',{@doUpdateIDnumber, IncludeCellcheckbox, IDnumbox, IDnumber} );

 uiwait(fig)
 ok=P.ok;
 if ok
     selection=P.IncludeCell;
     ClustQual=P.ClustQual;
     PVcell=P.PVcell;
     IDnumberList=P.IDnumberList;
else
    selection=[];
    ClustQual=[];
    PVcell=[];
    IDnumberList=[];
end
