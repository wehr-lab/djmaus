function [selection,value, ClustQual, PVcell] = CellListdlg(varargin)
%CellListdlg helper function for CellListBuilder
%modified from LISTDLG  List selection dialog box.
%   [SELECTION,OK] = LISTDLG('ListString',S) creates a modal dialog box
%   which allows you to select a string or multiple strings from a list.
%   SELECTION is a vector of indices of the selected strings (length 1 in
%   the single selection mode).  This will be [] when OK is 0.  OK is 1 if
%   you push the OK button, or 0 if you push the Cancel button or close the
%   figure.
%
%   Double-clicking on an item or pressing <CR> when multiple items are
%   selected has the same effect as clicking the OK button.  Pressing <CR>
%   is the same as clicking the OK button. Pressing <ESC> is the same as
%   clicking the Cancel button.
%
%   Inputs are in parameter,value pairs:
%
%   Parameter       Description
%   'ListString'    cell array of strings for the list box.
%   'SelectionMode' string; can be 'single' or 'multiple'; defaults to
%                   'multiple'.
%   'ListSize'      [width height] of listbox in pixels; defaults
%                   to [160 300].
%   'InitialValue'  vector of indices of which items of the list box
%                   are initially selected; defaults to the first item.
%   'Name'          String for the figure's title; defaults to ''.
%   'PromptString'  string matrix or cell array of strings which appears 
%                   as text above the list box; defaults to {}.
%   'OKString'      string for the OK button; defaults to 'OK'.
%   'CancelString'  string for the Cancel button; defaults to 'Cancel'.
%
%   A 'Select all' button is provided in the multiple selection case.
%
%   Example:
%     d = dir;
%     str = {d.name};
%     [s,v] = listdlg('PromptString','Select a file:',...
%                     'SelectionMode','single',...
%                     'ListString',str)
 %
%  See also DIALOG, ERRORDLG, HELPDLG, INPUTDLG,
%    MSGBOX, QUESTDLG, WARNDLG.

%   Copyright 1984-2014 The MathWorks, Inc.

%   'uh'            uicontrol button height, in pixels; default = 22.
%   'fus'           frame/uicontrol spacing, in pixels; default = 8.
%   'ffs'           frame/figure spacing, in pixels; default = 8.

% simple test:
%
% d = dir; [s,v] = listdlg('PromptString','Select a file:','ListString',{d.name});
% 

% Generate a warning in -nodisplay and -noFigureWindows mode.

narginchk(1,inf)

figname = '';
smode = 2;   % (multiple)
promptstring = {};
liststring = [];
listsize = [160 300];
initialvalue = [];
okstring = getString(message('MATLAB:uistring:popupdialogs:OK'));
cancelstring = getString(message('MATLAB:uistring:popupdialogs:Cancel'));
fus = 8;
ffs = 8;
uh = 22;

if mod(length(varargin),2) ~= 0
    % input args have not com in pairs, woe is me
    error(message('MATLAB:listdlg:InvalidArgument'))
end
for i=1:2:length(varargin)
    switch lower(varargin{i})
     case 'name'
      figname = varargin{i+1};
     case 'promptstring'
      promptstring = varargin{i+1};
     case 'selectionmode'
      switch lower(varargin{i+1})
       case 'single'
        smode = 1;
       case 'multiple'
        smode = 2;
      end
     case 'listsize'
      listsize = varargin{i+1};
     case 'liststring'
      liststring = varargin{i+1};
     case 'initialvalue'
      initialvalue = varargin{i+1};
     case 'uh'
      uh = varargin{i+1};
     case 'fus'
      fus = varargin{i+1};
     case 'ffs'
      ffs = varargin{i+1};
     case 'okstring'
      okstring = varargin{i+1};
     case 'cancelstring'
      cancelstring = varargin{i+1};
     otherwise
      error(message('MATLAB:listdlg:UnknownParameter', varargin{ i }))
    end
end

if ischar(promptstring)
    promptstring = cellstr(promptstring); 
end

if isempty(initialvalue)
    initialvalue = 1;
end

if isempty(liststring)
    error(message('MATLAB:listdlg:NeedParameter'))
end

fontsize=12;
ex = fontsize*1.7;  % height extent per line of uicontrol text (approx)

fp = get(0,'DefaultFigurePosition');
w = 2*(fus+ffs)+listsize(1);
h = 2*ffs+6*fus+ex*length(promptstring)+listsize(2)+uh+(smode==2)*(fus+uh);
fp = [fp(1) fp(2)+fp(4)-h 2*w h];  % keep upper left corner fixed

fig_props = { ...
    'windowstyle'            'normal'...%'modal' ...
    'name'                   figname ...
    'color'                  get(0,'DefaultUicontrolBackgroundColor') ...
    'resize'                 'off' ...
    'numbertitle'            'off' ...
    'menubar'                'none' ...
    'visible'                'on' ...
    'createfcn'              ''    ...
    'position'               fp   ...
    'closerequestfcn'        'delete(gcbf)' ...
            };

liststring=cellstr(liststring);

fig = figure(fig_props{:});

if ~isempty(promptstring)
    prompt_text = uicontrol('Style','text','String',promptstring,...
        'HorizontalAlignment','left',...
        'Position',[ffs+fus fp(4)-(ffs+fus+ex*length(promptstring)) ...
        listsize(1) ex*length(promptstring)]); %#ok
end

btn_wid = (fp(3)-2*(ffs+fus)-fus)/2;

listbox = uicontrol('Style','listbox',...
                    'Position',[ffs+fus ffs+uh+4*fus+(smode==2)*(fus+uh) listsize],...
                    'String',liststring,...
                    'BackgroundColor','w',...
                    'Max',smode,...
                    'Tag','listbox',...
                    'Value',initialvalue, ...
                    'fontsize', fontsize, ...
                    'Callback', {@doListboxClick});

top=ffs+uh+4*fus+(smode==2)*(fus+uh)+listsize(2);
ex2 = fontsize*1.2;  % height extent per line of uicontrol text (approx)
sliderwidth=120;
uicontrol('style', 'text', 'pos', [w, top+2, sliderwidth, ex2], 'string', 'cluster quality', 'fontsize', fontsize)
uicontrol('style', 'text', 'pos', [w+sliderwidth+20, top+2, 50, ex2], ...
    'string', 'PV cell?', 'fontsize', fontsize, 'horizontalalign', 'left')
for slidx=1:length(liststring)
    slstr(slidx)=uicontrol('style', 'text', 'pos', [w+sliderwidth+2, top-slidx*ex2, 10, ex2], 'string', 0);
    sl(slidx)=uicontrol('style', 'slider', 'pos', [w, top-slidx*ex2, sliderwidth, ex2], ...
        'min', 0, 'max', 5, 'sliderstep', [.2 .2], 'value', 0, 'callback', {@doSlider, slstr(slidx)} );
 pvcheckbox(slidx)=uicontrol('style', 'checkbox', 'pos', ...
     [w+sliderwidth+20, top-slidx*ex2, 30, ex2]);
end
                
ok_btn = uicontrol('Style','pushbutton',...
                   'String',okstring,...
                   'Position',[ffs+fus ffs+fus btn_wid uh],...
                   'Tag','ok_btn',...
                   'Callback',{@doOK,listbox, sl, pvcheckbox});

cancel_btn = uicontrol('Style','pushbutton',...
                       'String',cancelstring,...
                       'Position',[ffs+2*fus+btn_wid ffs+fus btn_wid uh],...
                       'Tag','cancel_btn',...
                       'Callback',{@doCancel,listbox});

if smode == 2
    selectall_btn = uicontrol('Style','pushbutton',...
                              'String',getString(message('MATLAB:uistring:popupdialogs:SelectAll')),...
                              'Position',[ffs+fus 4*fus+ffs+uh listsize(1) uh],...
                              'Tag','selectall_btn',...
                              'Callback',{@doSelectAll, listbox});

    if length(initialvalue) == length(liststring)
        set(selectall_btn,'Enable','off')
    end
    set(listbox,'Callback',{@doListboxClick, selectall_btn})
end

set([fig, ok_btn, cancel_btn, listbox], 'KeyPressFcn', {@doKeypress, listbox});

% set(fig,'Position',getnicedialoglocation(fp, get(fig,'Units'))); !!!
% Make ok_btn the default button.
setdefaultbutton(fig, ok_btn);

% make sure we are on screen
movegui(fig)
set(fig, 'Visible','on'); drawnow;

try
    % Give default focus to the listbox *after* the figure is made visible
    uicontrol(listbox);
    c = matlab.ui.internal.dialog.DialogUtils.disableAllWindowsSafely();
    uiwait(fig);
    delete(c);
catch
    if ishghandle(fig)
        delete(fig)
    end
end

if isappdata(0,'ListDialogAppData__')
    ad = getappdata(0,'ListDialogAppData__');
    selection = ad.selection;
    value = ad.value;
    ClustQual=ad.ClustQual;
    PVcell=ad.PVcell;
    rmappdata(0,'ListDialogAppData__')
else
    % figure was deleted
    selection = [];
    value = 0;
end
drawnow; % Update the view to remove the closed figure (g1031998)

%% figure, OK and Cancel KeyPressFcn
function doKeypress(src, evd, listbox) %#ok
switch evd.Key
 case 'escape'
  doCancel([],[],listbox);
end

%% slider callback
function doSlider(sliderh, evd, sliderstrh) %#ok
set(sliderstrh, 'string', get(sliderh, 'value'))

%% OK callback
function doOK(ok_btn, evd, listbox, sl, pvcheckbox) %#ok
if 1%(~isappdata(0, 'ListDialogAppData__'))
    ad.value = 1;
    ad.selection = get(listbox,'Value');
    for i=1:length(sl)
        ad.ClustQual(i)=get(sl(i), 'value');
        ad.PVcell(i)=get(pvcheckbox(i), 'value');
    end
    setappdata(0,'ListDialogAppData__',ad);
    delete(gcbf);
end

%% Cancel callback
function doCancel(cancel_btn, evd, listbox) %#ok
ad.value = 0;
ad.selection = [];
ad.ClustQual=[];
setappdata(0,'ListDialogAppData__',ad)
delete(gcbf);

%% SelectAll callback
function doSelectAll(selectall_btn, evd, listbox) %#ok
set(selectall_btn,'Enable','off')
set(listbox,'Value',1:length(get(listbox,'String')));

%% Listbox callback
function doListboxClick(listbox, evd, selectall_btn) %#ok
% if this is a doubleclick, doOK
if strcmp(get(gcbf,'SelectionType'),'open')
    doOK([],[],listbox);
elseif nargin == 3
    if length(get(listbox,'String'))==length(get(listbox,'Value'))
        set(selectall_btn,'Enable','off')
    else
        set(selectall_btn,'Enable','on')
    end
end

function setdefaultbutton(figHandle, btnHandle)
% WARNING: This feature is not supported in MATLAB and the API and
% functionality may change in a future release.

%SETDEFAULTBUTTON Set default button for a figure.
%  SETDEFAULTBUTTON(BTNHANDLE) sets the button passed in to be the default button
%  (the button and callback used when the user hits "enter" or "return"
%  when in a dialog box.
%
%  This function is used by inputdlg.m, msgbox.m, questdlg.m and
%  uigetpref.m.
%
%  Example:
%
%  f = figure;
%  b1 = uicontrol('style', 'pushbutton', 'string', 'first', ...
%       'position', [100 100 50 20]);
%  b2 = uicontrol('style', 'pushbutton', 'string', 'second', ...
%       'position', [200 100 50 20]);
%  b3 = uicontrol('style', 'pushbutton', 'string', 'third', ...
%       'position', [300 100 50 20]);
%  setdefaultbutton(b2);
%

%  Copyright 2005-2007 The MathWorks, Inc.

%--------------------------------------- NOTE ------------------------------------------
% This file was copied into matlab/toolbox/local/private.
% These two files should be kept in sync - when editing please make sure
% that *both* files are modified.

% Nargin Check
narginchk(1,2)

if (usejava('awt') == 1)
    % We are running with Java Figures
    useJavaDefaultButton(figHandle, btnHandle)
else
    % We are running with Native Figures
    useHGDefaultButton(figHandle, btnHandle);
end

    function useJavaDefaultButton(figH, btnH)
        % Get a UDD handle for the figure.
        fh = handle(figH);
        % Call the setDefaultButton method on the figure handle
        fh.setDefaultButton(btnH);
        
        
        function useHGDefaultButton(figHandle, btnHandle)
            % First get the position of the button.
            btnPos = getpixelposition(btnHandle);
            
            % Next calculate offsets.
            leftOffset   = btnPos(1) - 1;
            bottomOffset = btnPos(2) - 2;
            widthOffset  = btnPos(3) + 3;
            heightOffset = btnPos(4) + 3;
            
            % Create the default button look with a uipanel.
            % Use black border color even on Mac or Windows-XP (XP scheme) since
            % this is in natve figures which uses the Win2K style buttons on Windows
            % and Motif buttons on the Mac.
            h1 = uipanel(get(btnHandle, 'Parent'), 'HighlightColor', 'black', ...
                'BorderType', 'etchedout', 'units', 'pixels', ...
                'Position', [leftOffset bottomOffset widthOffset heightOffset]);
            
            % Make sure it is stacked on the bottom.
            uistack(h1, 'bottom');
            

