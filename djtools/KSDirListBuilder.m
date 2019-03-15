function KSDirListBuilder(varargin)
% A GUI for building a list of kilosort master data directories
% we don't add invidual cells, just the master data directory for each
% recording session.
% it has a button to add current directory to the list
% another button that opens a dialog box so you can browse to the desired directory
%
% what does it do? It creates/appends a machine readable cell list as a
% .txt file
%
%mw 3.14.2019

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
        delete(P.listviewer)
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
    case 'NextDataDir'
        NextDataDir
        set(P.pwdDisplay, 'string', pwd)
    case 'PreviousDataDir'
        PreviousDataDir
        set(P.pwdDisplay, 'string', pwd)
        
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
    UpdateListviewer
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
global P
d = uigetdir(pwd, 'choose directory to scan')
set(P.pwdDisplay, 'string', pwd)
if d
    cd(d)
    AddCurrentDir
end

function WriteToCellList(d)
global P

%here we can write out any additional stimulus or notebook info we want
%see CellListBuilder.m for examples
response=questdlg(pwd, 'Write this data directory to file?', 'Write', 'Cancel', 'Write');
switch response
    case 'Write'
        fid=fopen(P.TargetCellList, 'a'); %absolute path
        fprintf(fid, '\n\ndir');
        fprintf(fid, '\n%s', d);
        fclose(fid);
        UpdateListviewer
end

function UpdateListviewer
global P
        fid=fopen(P.TargetCellList, 'r'); %absolute path
     i=0;
     while 1
            tline = fgetl(fid);
     if ~ischar(tline), break, end
     i=i+1;
            str{i}=tline;
        end
        fclose(fid);
set(P.listviewerbox,'visible','on', 'string', str);
figure(P.listviewer)

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
    'doublebuffer','on','menubar','none','closerequestfcn', [me '(''Close'');'])
height=220; width=350; e=2; H=e;
w=200; h=25;
set(fig,'pos',[800 800         width         height],'visible','on');

P.listviewer=figure;
set(P.listviewer,'pos',[1000 100 500 600],'visible','on');
P.listviewerbox= uicontrol('parent',P.listviewer,'string','','tag','listviewerbox','units','pixels',...
    'position',[1 1 500 600],'enable','inactive',...
    'horiz', 'left', 'style','edit', 'max', 2 );

%TargetCellList display
P.TargetCellListDisplay= uicontrol('parent',fig,'string','','tag','TargetCellListDisplay','units','pixels',...
    'position',[e H width-e 2*h],'enable','inactive',...
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

% Next and Previous buttons
H=H+h+e;
P.Next=uicontrol('parent',fig,'string','Next Dir','tag','NextDataDir','units','pixels',...
    'position',[e+w/2 H w/2 h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);
P.Previous=uicontrol('parent',fig,'string','Previous Dir','tag','PreviousDataDir','units','pixels',...
    'position',[e H w/2 h],'enable','on',...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%Browse and Add button
H=H+h+e;
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


%pwd display
H=H+h;
P.pwdDisplay= uicontrol('parent',fig,'string',pwd,'tag','pwdDisplay','units','pixels',...
    'position',[e H width-e h],'enable','inactive',...
    'fontweight','bold','horiz', 'left',...
    'style','text');
