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

function TargetCellList
global P
[fname, path] = uiputfile('*.txt', 'Select cell list');
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
    [liststr{1:length(d)}]=deal(d.name);
    [selection, ok, ClustQual, PVcell]=listdlg2('ListString', liststr, 'InitialValue', 1:length(d), 'PromptString', 'select which files to include');
    if ok
        WriteToCellList(d(selection), ClustQual, PVcell)
    end


function BrowseAndAdd
d = uigetdir(pwd, 'choose directory to scan')
if d
    cd(d)
    AddCurrentDir
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
set(fig,'visible','off','numbertitle','off','name','cell list builder',...
    'doublebuffer','on','menubar','none','closerequestfcn','CellListBuilder(''Close'')')
height=200; width=350; e=2; H=e;
w=200; h=25;
set(fig,'pos',[1200 900         width         height],'visible','on');


%TargetCellList display
P.TargetCellListDisplay= uicontrol('parent',fig,'string','','tag','TargetCellListDisplay','units','pixels',...
    'position',[e H width-e 2*h],'enable','on',...
    'fontweight','bold','horiz', 'left',...
    'style','text');

%TargetCellList button
H=H+2*h+e;
uicontrol('parent',fig,'string','Select Cell List','tag','TargetCellList','units','pixels',...
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
