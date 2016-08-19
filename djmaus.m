function djmaus(varargin)

% Plays a bunch of different stimuli stored in a .mat file ('protocol')
% Uses PPA sound which requires PsychToolbox. Download from psychtoolbox.org.
% Sends stimulus information to Open Ephys GUI using zeroMQ, Win: download from
% zeromq.org Linux: apt-get install libzmqpp-dev Mac: brew install zeromq


%now has the capability to load files programatically
%call djmaus('load', fullfilename)
%where fullfilename includes the absolute path (e.g. 'D:\lab\exper2.2\protocols\Tuning Curve protocols\tuning-curve-tones-20f_1000-20000Hz-1a_80-80dB-1d_400ms-isi500ms.mat'
%mw 09.24.13

%ISI works by setting the djTimerDelay to
%stimulus.param.duration/1000 + stimulus.param.next/1000

global SP pref djTimer djTimerDelay

if nargin > 0
    action = varargin{1};
else
    action = get(gcbo,'tag');
end
if isempty(action)
    action='Init';
end

%fprintf('\naction: %s', action)

switch action
    case 'Init'
        SP.user='lab';
        djPrefs;
        InitializeGUI;
        PPAdj('init');
        % set the timer
        djTimer=timer('TimerFcn',[me '(''next_stimulus'');'],'StopFcn',[me '(''restart_timer'');'],'ExecutionMode','singleShot');
        
    case 'Close'
        try stop(djTimer); end
        delete(djTimer);
        clear djTimer
        delete(SP.fig)
        PPAdj('close')
        zeroMQwrapper('CloseThread',SP.zhandle);
        clear global SP
        fprintf('\nbye\n')
        
    case 'Run'
        if SP.Run
            %we want to stop
            set(SP.Runh,'backgroundcolor',[0 0.9 0],'String','Play');
            djTimerDelay=-1;
            stop(djTimer);
            SP.Run=0;
        else
            set(SP.Runh, 'backgroundcolor',[0.9 0 0],'String','Playing...');
            start(djTimer);
            SP.Run=1;
        end
        
    case 'Reset'
        djMessage('Reset');
        prot=SP.ProtocolIndex;
        if prot>0
            cstim=SP.CurrentStimulus;
            cstim(prot)=0;
            SP.CurrentStimulus=cstim;
            SP.NRepeats=0;
        end
        
    case 'Reset All'
        djMessage('Reset All');
        cstim=SP.CurrentStimulus;
        cstim(:)=0;
        SP.CurrentStimulus=cstim;
        SP.NRepeats=0;
        
    case 'Repeat'
        if SP.Repeat %it's on, so turn it off
            set(SP.Repeath, 'backgroundcolor',[0 0 1],'foregroundcolor',[0 0 0]);
            set(SP.Repeath, 'string','Repeat Off');
            SP.Repeat=0;
        else
            set(SP.Repeath, 'backgroundcolor',[1 0 0],'foregroundcolor',[0 0 1]);
            set(SP.Repeath, 'string','Repeat On');
            SP.Repeat=1;
        end
        
    case 'next_stimulus'
        NextStimulus;
        
    case 'restart_timer'
        if djTimerDelay>-1
            djTimerDelay= round(1000*djTimerDelay)/1000;
            set(djTimer,'StartDelay',djTimerDelay);
            start(djTimer);
        else
            stop(djTimer);
            set(djTimer,'StartDelay',0);    % next time we push the Play button, it will start immediately
        end
        
    case 'load'
        LoadProtocol
        
        
    case 'ProtocolMenu'
        if isempty(SP.NStimuli) return;end
        protocol_choice=get(SP.ProtocolMenuh, 'value');
        SP.ProtocolIndex=protocol_choice;
        current=SP.CurrentStimulus;
        nstimuli=SP.NStimuli(protocol_choice);
        desc=SP.Description;
        pdh=SP.protocol_descriptionh;
        outstring=textwrap(pdh,{desc{protocol_choice}});
        set(pdh,'String',outstring);
        name=SP.Name;
        pnh=SP.protocol_nameh;
        outstring=textwrap(pnh(1),{name{protocol_choice}});
        set(pnh,'String',outstring);
        SP.NRepeats=0;
        djMessage(['0/' num2str(SP.NStimuli(protocol_choice)) ' stimuli, ' num2str(SP.NRepeats), ' repeats' ]);
        
        
    case 'User'
        users=get(SP.userh, 'string');
        user_index=get(SP.userh, 'value');
        user=users{user_index};
        djMessage(sprintf('Hi %s', user))
        SP.user=user;
        djPrefs;
        SP.datapath=pref.datapath;
        set(SP.pathh, 'string', SP.datapath)
        set(SP.mouseIDMenuh, 'string', pref.allmouseIDs)
        
    case 'mouseIDMenu'
        mouseIDs=get(SP.mouseIDMenuh, 'string');
        if ~iscell(mouseIDs) %skip if only one menu item
        else
            mouseID=mouseIDs{get(SP.mouseIDMenuh, 'value')};
            set(SP.mouseIDh, 'string', mouseID)
        end
        SP.mouseID=mouseID;
        LoadMouse
        
    case 'mouseID'
        SP.mouseID=get(SP.mouseIDh, 'string');
        mouseIDs=get(SP.mouseIDMenuh, 'string');
        if ~iscell(mouseIDs) mouseIDs={mouseIDs};end

        if isempty(mouseIDs)
            mouseIDs=SP.mouseID;
        elseif iscell(mouseIDs) & length(mouseIDs)>1 %more than one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouseID});
        elseif iscell(mouseIDs) & length(mouseIDs)==1 %only one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouseID});
            set(SP.mouseIDMenuh, 'value',1)
        elseif length(mouseIDs)==1 %only one menu item
            mouseIDs={mouseIDs, SP.mouseID};
            set(SP.mouseIDMenuh, 'value',1)
        else
            error('?')
        end
        set(SP.mouseIDMenuh, 'string',mouseIDs)
        SP.allmouseIDs=mouseIDs;
        WriteMouseIDtoPrefs
        LoadMouse
        
    case 'Record'
        Record
        
    case 'mouseGenotype'
        SP.mouseGenotype=get(SP.mouseGenotypeh, 'string');
        cd (SP.datapath)
        load mouseDB
        str=sprintf('%s.mouseGenotype=''%s''', SP.mouseID, SP.mouseGenotype);
        eval(str)
        save('mouseDB.mat', SP.mouseID, '-append')
        
    case 'mouseSex'
        SP.mouseSex=lower(get(SP.mouseSexh, 'string'));
        set(SP.mouseSexh, 'string', SP.mouseSex)
        cd (SP.datapath)
        load mouseDB
        str=sprintf('%s.mouseSex=''%s''', SP.mouseID, SP.mouseSex);
        eval(str)
        save('mouseDB.mat', SP.mouseID, '-append')
        
    case 'mouseDOB'
        SP.mouseDOB=get(SP.mouseDOBh, 'string');
        cd (SP.datapath)
        load mouseDB
        str=sprintf('%s.mouseDOB=''%s''', SP.mouseID, SP.mouseDOB);
        eval(str)
        save('mouseDB.mat', SP.mouseID, '-append')
        SP.Age=(datenum(date)-datenum(SP.mouseDOB))/30; %in months
        pos1=get(SP.mouseSexh, 'pos'); 
        pos=get(SP.mouseDOBh, 'pos'); 
        pos(3)=.6*pos1(3);
        set(SP.mouseDOBh, 'pos',pos);
        pos(1)=pos(1)+pos(3);
        SP.Ageh=uicontrol(gcf,'tag','mouseAgelabel','style','text','units','pixels',...
            'string', sprintf('age\n%.1f mo', SP.Age), 'fontsize', 10,...
            'enable','inact','horiz','left','pos', pos);

        
    case 'Reinforcement'
        SP.Reinforcement=get(SP.Reinforcementh, 'string');
        
    case 'Drugs'
        SP.Drugs=get(SP.Drugsh, 'string');
        
    case 'LaserPower'
        SP.LaserPower=get(SP.LaserPowerh, 'string');
        
    case 'Depth'
        SP.Depth=get(SP.Depthh, 'string');
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMouse
global SP
cd (SP.datapath)
mouseDB=load('mouseDB.mat');
if isfield(mouseDB, SP.mouseID)
    try
        str=sprintf('SP.mouseGenotype=mouseDB.%s.mouseGenotype;', SP.mouseID);
        eval(str);
        set(SP.mouseGenotypeh, 'string', SP.mouseGenotype);
    end
    try
        str=sprintf('SP.mouseDOB=mouseDB.%s.mouseDOB;', SP.mouseID);
        eval(str);
        set(SP.mouseDOBh, 'string', SP.mouseDOB);
    end
    try
        str=sprintf('SP.mouseSex=mouseDB.%s.mouseSex;', SP.mouseID);
        eval(str);
        set(SP.mouseSexh, 'string', SP.mouseSex);
    end
else
    SP.mouseGenotype='genotype unknown';
    set(SP.mouseGenotypeh, 'string', SP.mouseGenotype);
    SP.mouseSex='sex unknown';
    set(SP.mouseSexh, 'string', SP.mouseSex);
    SP.mouseDOB='age unknown';
    set(SP.mouseDOBh, 'string', SP.mouseDOB);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteMouseIDtoPrefs
global SP pref
cd (pref.root)
fid=fopen('djPrefs.m', 'a+');
if ~ iscell(SP.allmouseIDs) SP.allmouseIDs={SP.allmouseIDs};end
allmouseIDstr='{';
for i=1:length(SP.allmouseIDs)
    allmouseIDstr=[allmouseIDstr, '''', SP.allmouseIDs{i}, ''','];
end
allmouseIDstr=allmouseIDstr(1:end-1);
allmouseIDstr=[allmouseIDstr, '}'];
key=sprintf('%%saved mouseIDs for %s', SP.user);

Preftext = regexp( fileread('djPrefs.m'), '\n', 'split');
I=strmatch(key, Preftext);
if ~isempty(I) %found key, overwrite with revised entry
    I=I(1);
    newstr=sprintf('\tpref.allmouseIDs=%s;', allmouseIDstr);
    Preftext{I+3}=newstr; %change entry;
    fid = fopen('djPrefs.m', 'w');
    fprintf(fid, '%s\n', Preftext{:});
    fclose(fid);
else
    fprintf(fid, '\n%s', key);
    fprintf(fid, '\nswitch SP.user\ncase  ''%s''', SP.user);
    fprintf(fid, '\n\tpref.allmouseIDs=%s;', allmouseIDstr);
    fprintf(fid, '\nend\n');
    fclose(fid);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadProtocol
%adding capability to push files programatically
%call dj('load', 'fullfilename')
%where fullfilename includes the absolute path (e.g. 'D:\lab\exper2.2\protocols\Tuning Curve protocols\tuning-curve-tones-20f_1000-20000Hz-1a_80-80dB-1d_400ms-isi500ms.mat'
global pref SP
if nargin==2
    if varargin{2}
        try
            [pathname, filename,ext] = fileparts(varargin{2});
        end
    end
else
    
    currentdir=pwd;
    cd(pref.stimuli);
    [filename, pathname] = uigetfile('*.mat', 'Pick a protocol file');
    cd(currentdir);
end %if nargin==2
if isequal(filename,0) || isequal(pathname,0)
    return;
else
    %try
    stimuli=load(fullfile(pathname,filename));
    f=fieldnames(stimuli);
    stimuli=eval(['stimuli.' f{1}]);        % we take the first field in the loaded structure
    if strcmpi(stimuli(1).type,'exper2 stimulus protocol')
        %convert
        stimuli=convert_stimulus_protocol(stimuli);
    end
    if strcmpi(stimuli(1).version,'djmaus')
        
        AllProtocols=SP.AllProtocols;
        desc=stimuli(1).protocol_description;
        description=SP.Description;
        description={description{:} desc};
        SP.Description=description;
        newname=stimuli(1).protocol_name;
        name=SP.Name;
        name={name{:} newname};
        SP.Name=name;
        
        %stimuli(1)=[];
        stim=SP.StimulusProtocols;
        stim={stim{:} stimuli};
        SP.StimulusProtocols=stim;
        
        cstim=SP.CurrentStimulus;
        cstim=[cstim 0];
        SP.CurrentStimulus=cstim;
        
        nstim=SP.NStimuli;
        nstim=[nstim length(stimuli)];
        SP.NStimuli=nstim;
        
        nprot=SP.NProtocols;
        nprot=nprot+1;
        SP.NProtocols=nprot;
        
        allprot=SP.AllProtocols;
        if isequal(allprot,{''})    % this is the first protocol added
            djMessage('load first stim')
            allprot={newname};
            set(SP.Runh,'enable','on');
            pdh=SP.protocol_descriptionh;
            outstring=textwrap(pdh(1),{desc});
            set(pdh,'String',outstring);
            pnh=SP.protocol_nameh;
            outstring=textwrap(pnh(1),{newname});
            set(pnh,'String',outstring);
        else
            allprot={allprot{:} newname};
        end
        SP.AllProtocols=allprot;
        set(SP.ProtocolMenuh,'String',allprot);
        %                          if numel(allprot)==1    % if this is the first protocol, make it active
        %                              SP.ProtocolIndex=1;
        %                          end
        
        %Make this protocol active
        SP.ProtocolIndex=SP.NProtocols;
        set(SP.ProtocolMenuh, 'value', SP.NProtocols);
        djmaus('ProtocolMenu')
        
        SP.Stimulus=[];
        djMessage('Stimuli loaded', 'append');
        djMessage([ num2str(nstim(nprot)), ' stimuli'], 'append');
        
        SP.NRepeats=0;
        set(SP.Runh, 'enable', 'on')
        
        
    else
        djMessage('Not a valid protocol file');
    end
    %catch
    %    djMessage('Can''t load the protocol file');
    %end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SendStimulus(stimulus)
global SP
%make stimulus
fcn=StimTypes(stimulus.type);
samples=feval(fcn,stimulus.param,SP.SoundFs);

%LoadPPA(type,where,param)
PPAdj('load', 'var', samples, stimulus.param)
str=sprintf('TrialType %s', stimulus.stimulus_description);
if ~isempty(SP.zhandle)
    zeroMQwrapper('Send', SP.zhandle, str)
end
PPAdj('playsound')
UpdateNotebookFile(stimulus)
djMessage(stimulus.stimulus_description, 'append')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function NextStimulus
global SP djTimer djTimerDelay
protocol=SP.ProtocolIndex;
current=SP.CurrentStimulus;
nstimuli=SP.NStimuli;
if current(protocol)==nstimuli(protocol) % we played everything
    current(protocol)=0;                % so let's start from the beginning;
    SP.NRepeats=SP.NRepeats+1;
    if ~SP.Repeat           % if we don't want to repeat,
        djTimerDelay=-1;          % this will cause the timer to stop
        SP.Run=0;
        set(SP.Runh,'backgroundcolor',[0 0.9 0],'String','Play');
        SP.CurrentStimulus=current;
        return;
    end
end
current(protocol)=current(protocol)+1;
SP.CurrentStimulus=current;
stimuli=SP.StimulusProtocols;
stimuli=stimuli{protocol};
stimulus=stimuli(current(protocol));
if isfield(stimulus.param,'next')
    iti=stimulus.param.next/1000;
else
    iti=0.5;    % set default iti to 500ms;
end

if ~isfield(stimulus.param, 'duration')
    error('Improperly designed stimulus: no duration field. ')
end
djTimerDelay=stimulus.param.duration/1000+iti;

djMessage([num2str(current(protocol)) '/' num2str(nstimuli(protocol)) ', ' num2str(SP.NRepeats), ' repeats' ]);
SendStimulus(stimulus);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return the name of this file/module.
function out = me
out = mfilename;
% me
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DisableParamChange
fig=findobj('type','figure','tag',me);
h=findobj(fig,'type','uicontrol','style','edit');
for cnt=1:length(h)
    set(h(cnt),'enable','off')
end
% DisableParamChange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EnableParamChange
fig=findobj('type','figure','tag',me);
h=findobj(fig,'type','uicontrol','style','edit');
for cnt=1:length(h)
    set(h(cnt),'enable','on')
end
%EnableParamChange
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fcn=StimTypes(type)
switch type
    case 'tone'
        fcn='MakeTone';
    case 'whitenoise'
        fcn='MakeWhiteNoise';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Record
%toggle Open Ephys recording state, increment datafile
global SP
fprintf('\n')

if SP.Record
    %we want to stop recording
    zeroMQwrapper('Send',SP.zhandle ,'StopRecord');
    set(SP.Recordh,'backgroundcolor',[0 0.9 0],'String','Record');
    SP.Record=0;
    try
        set(SP.pathh, 'string', {SP.datapath, [SP.activedir, ' finished']})
    end
else
    %we want to start recording;
    startstr=sprintf('StartRecord');
    if ~isfield(SP, 'mouseID')
        SP.mouseID='none';
    end
    zeroMQwrapper('Send',SP.zhandle ,'StartAcquisition'); %shouldn't need to do this unless user stopped acquisition, doesn't hurt anyway
    startstr=sprintf('StartRecord CreateNewDir=1 RecDir=%s AppendText=mouse-%s', SP.datapath, SP.mouseID);
    zeroMQwrapper('Send',SP.zhandle ,startstr);
    set(SP.Recordh, 'backgroundcolor',[0.9 0 0],'String','Recording...');
    SP.Record=1;
    SP.stimcounter=0;
    if isfield(SP, 'stimlog')
        SP=rmfield( SP, 'stimlog');
    end
    InitNotebookFile
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitNotebookFile
global SP 

% find active OE data directory and cd into it
try
    zeroMQwrapper('Send',SP.zhandle ,sprintf('ChangeDirectory %s', SP.datapath))
    pause(.1)
    zeroMQwrapper('Send',SP.zhandle ,'GetRecordingPath');
    pause(.2)
    cd(SP.datapath)
    fid=fopen('RecordingPath.txt', 'r');
    RecordingPath=fgetl(fid);
    fclose(fid);
    fprintf('\nread this Recording Path from file:%s', RecordingPath)
% %     cd(SP.datapath)
% %     d=dir;
% %     [x,y]=sort(datenum(char({d.date})));
% %     i=1;
% %     today=datestr(now, 'yyyy-mm-dd');
% %     while isempty(strfind(d(y(i)).name, today))
% %         i=i+1;
% %     end
    SP.activedir=RecordingPath;
    cd(SP.activedir)
    set(SP.pathh, 'string', {SP.datapath, [SP.activedir, ' recording...']})
    
    %future directions: use hdf5 files for this
    %there aren't any easy utilities to store compound data (i.e. structures)
    %in h5, so I'm just using .mat for now, could write utility in the future
    
    nb.user=SP.user;
    nb.mouseID=SP.mouseID;
    nb.Depth=SP.Depth;
    nb.datapath=SP.datapath ;
    nb.activedir=SP.activedir;
    nb.LaserPower=SP.LaserPower;
    nb.mouseDOB=SP.mouseDOB;
    nb.mouseSex=SP.mouseSex;
    nb.mouseGenotype=SP.mouseGenotype;
    
    save('notebook.mat', 'nb')
    fprintf('\ncreated notebook file in %s', nb.activedir)
catch
    fprintf('\nCould not create notebook file in active data directory')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateNotebookFile(stimulus)
global SP

if SP.Record
    try
        SP.stimcounter=SP.stimcounter+1;
        timestamp=datestr(now,'mmmm-dd-yyyy HH:MM:SS.FFF');
        stimulus.timestamp=timestamp;
        cd(SP.datapath)
        cd(SP.activedir)
        SP.stimlog(SP.stimcounter)=stimulus;
        stimlog=SP.stimlog;
        save('notebook.mat', '-append', 'stimlog')
    catch
        fprintf('\nCould not update notebook file in active data directory')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitializeGUI
global SP pref
if isfield(SP, 'fig')
    try
        close(SP.fig)
    end
end
fig = figure;
SP.fig=fig;
set(fig,'visible','off');
set(fig,'visible','off','numbertitle','off','name','djmaus',...
    'doublebuffer','on','menubar','none','closerequestfcn','djmaus(''Close'')')
height=450; width=420; e=2; H=e;
w=100; h=25;
set(fig,'pos',[1200 600 width height],'visible','on');

%Reset button
uicontrol('parent',fig,'string','Reset','tag','Reset','units','pixels',...
    'position',[e H w h],'enable','on','foregroundcolor',[0.9 0 0],...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);

%Reset All button
uicontrol('parent',fig,'string','Reset All','tag','Reset All','units','pixels',...
    'position',[2*e+w H w h],'enable','on','foregroundcolor',[0.9 0 0],...
    'fontweight','bold',...
    'style','pushbutton','callback',[me ';']);
H=H+h+e;

%description window
SP.protocol_descriptionh=uicontrol(fig,'tag','protocol_description','style','text','units','pixels',...
    'enable','inact','horiz','left','backgroundcolor', [.8 .8 .8],'pos',[e H 2*w 3*h]);
H=H+3*h+e;

%name window
SP.protocol_nameh=uicontrol(fig,'tag','protocol_name','style','text','units','pixels',...
    'enable','inact','horiz','left','backgroundcolor', [.8 .8 .8],'pos', [e H 2*w 2*h]);
H=H+2*h+e;

%Protocol menu
SP.ProtocolMenuh=uicontrol(fig,'tag','ProtocolMenu','style','popupmenu','units','pixels','fontweight','bold',...
    'string', 'empty','enable','on','horiz','left','callback',[me ';'], 'pos',[e H 2*w h]);
SP.ProtocolIndex=[];
SP.AllProtocols={''};
SP.StimulusProtocols={};
H=H+h+e;

%Load button
uicontrol('parent',fig,'string','Load from...','tag','load','units','pixels',...
    'position',[e  H w h],'enable','on',...
    'style','pushbutton','callback',[me ';']);

SP.Repeat=0;

%Repeat button
SP.Repeath=uicontrol(fig,'tag','Repeat','style','togglebutton','units','pixels','fontweight','bold',...
    'string', 'Repeat Off','enable','on','horiz','left', 'callback',[me ';'],'pos',[2*e+w H w h]);
H=H+h+e;

%path display
SP.pathh=uicontrol(fig,'tag','pathdisplay','style','text','units','pixels',...
    'string', pref.datapath, 'enable','inact','horiz','left', 'pos',[e  H 3*w 1.5*h ]);
H=H+1.5*h+e;

%PPAactive message window
SP.PPAactive=uicontrol(fig,'tag','PPAactive','style','text','units','pixels',...
    'string', 'PPA messages', ...
    'enable','inact','horiz','left','backgroundcolor', [.8 .8 .8],'pos', [e H 3*w 1*h]);
H=H+h+e;

%Message window
SP.Messageh=uicontrol(fig,'tag','message','style','edit','fontweight','bold','units','pixels',...
    'enable','inact','horiz','left','Max', 6, 'pos',[e  H 3*w 3*h ]);
H=H+3*h+e;

%Run button
SP.Run=0;
SP.Runh=uicontrol(fig,'tag','Run','style','togglebutton','units','pixels','fontweight','bold',...
    'fontsize',12,'fontname','Arial', ...
    'string', 'Play','callback', [me ';'],'enable','off','horiz','left','pos',[e H w 2*h ]);

%Record button
SP.Record=0;
SP.Recordh=uicontrol(fig,'tag','Record','style','togglebutton','units','pixels','fontweight','bold',...
    'fontsize',12,'fontname','Arial', ...
    'string', 'Record','callback', [me ';'],'enable','off','horiz','left','pos',[e+w H w 2*h ]);

% %New Session button
% SP.NewSessionh=uicontrol(fig,'tag','NewSession','style','pushbutton','units','pixels',...
%     'fontsize',10,'fontname','Arial', ...
%     'string', 'New Session','callback', [me ';'],'enable','on','horiz','left','pos',[e+2*w H w 1*h ]);

H=H+2*h+e;

%User menu
SP.user='lab'; %default user
SP.userh=uicontrol(fig,'tag','User','style','popupmenu','units','pixels','fontweight','bold',...
    'string', pref.users,'enable','on','horiz','left','callback',[me ';'], 'pos',[e H w h ]);
H=H+h;
SP.userlabel=uicontrol(fig,'tag','userlabel','style','text','units','pixels',...
    'string', 'user', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [e H w h/2]);
H=H+h+e;



%%%%%%%%%%%%%%%%%%%%%%%%%%%
%notebook features
H=e;
%Notes edit box
SP.Notesh=uicontrol(fig,'tag','Notes','style','edit','fontweight','bold','units','pixels',...
    'string', '','horiz','left', 'Max', Inf, 'callback',[me ';'],'pos',[2*e+2*w  H 2*w 3*h ]);
H=H+3*h+e;
SP.Noteslabel=uicontrol(fig,'tag','Noteslabel','style','text','units','pixels',...
    'string', 'Notes', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+2*w  H w h/2]);
H=H+h+e;

%mouseID menu
warning('off', 'MATLAB:hg:uicontrol:StringMustBeNonEmpty');
if isfield(pref, 'allmouseIDs') SP.allmouseIDs=pref.allmouseIDs; else SP.allmouseIDs='';end
SP.mouseIDMenuh=uicontrol(fig,'tag','mouseIDMenu','style','popupmenu','units','pixels','fontweight','bold',...
    'string', SP.allmouseIDs,'enable','on','horiz','left','callback',[me ';'], 'pos',[e+2*w H w h]);
H=H+h+e;

%mouseID edit box
SP.mouseIDh=uicontrol(fig,'tag','mouseID','style','edit','fontweight','bold','units','pixels',...
    'string', '','horiz','left', 'callback',[me ';'],'pos',[2*e+2*w  H w h ]);
H=H+h+e;
SP.mouseIDlabel=uicontrol(fig,'tag','mouseIDlabel','style','text','units','pixels',...
    'string', 'mouseID', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+2*w  H w h/2]);
H=H+h/2+e;

H=H-3.5*h;
% Manipulation/conditions details (anesthesia, Any other drugs, Shock, Reward, etc)
%Drugs edit box
SP.Drugsh=uicontrol(fig,'tag','Drugs','style','edit','units','pixels',...
    'string', 'none','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h;
SP.Drugslabel=uicontrol(fig,'tag','Drugslabel','style','text','units','pixels',...
    'string', 'Drugs:', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;

%Reinforcement edit box (shock, reward, etc)
SP.Reinforcementh=uicontrol(fig,'tag','Reinforcement','style','edit','units','pixels',...
    'string', 'none','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h;
SP.Reinforcementlabel=uicontrol(fig,'tag','Reinforcementlabel','style','text','units','pixels',...
    'string', 'Reinforcement:', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;

SP.Reinforcement='none';
SP.Drugs='none';

%mouse details
SP.mouseDOBh=uicontrol(fig,'tag','mouseDOB','style','edit','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h;
SP.mouseDOBlabel=uicontrol(fig,'tag','mouseDOBlabel','style','text','units','pixels',...
    'string', 'mouseDOB:', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;

SP.mouseSexh=uicontrol(fig,'tag','mouseSex','style','edit','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h;
SP.mouseSexlabel=uicontrol(fig,'tag','mouseSexlabel','style','text','units','pixels',...
    'string', 'mouseSex:', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;
SP.mouseGenotypeh=uicontrol(fig,'tag','mouseGenotype','style','edit','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h;
SP.mouseGenotypelabel=uicontrol(fig,'tag','mouseGenotypelabel','style','text','units','pixels',...
    'string', 'mouseGenotype:', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;
SP.mouseSex='unknown';
SP.mouseDOB='unknown';
SP.mouseGenotype='unknown';

%Depth edit box
SP.Depthh=uicontrol(fig,'tag','Depth','style','edit','fontweight','bold','units','pixels',...
    'string', '','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h+e;
SP.Depthlabel=uicontrol(fig,'tag','Depthlabel','style','text','units','pixels',...
    'string', 'Depth', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;

%laser power edit box
SP.LaserPowerh=uicontrol(fig,'tag','LaserPower','style','edit','fontweight','bold','units','pixels',...
    'string', '','horiz','left', 'callback',[me ';'],'pos',[2*e+3*w  H w h ]);
H=H+h+e;
SP.LaserPowerlabel=uicontrol(fig,'tag','LaserPowerlabel','style','text','units','pixels',...
    'string', 'LaserPower', 'fontsize', 10,...
    'enable','inact','horiz','left','pos', [2*e+3*w  H w h/2]);
H=H+h/2+e;
SP.Depth='unknown';
SP.LaserPower='unknown';









%%%%%%%%%%%%%%%%%%%%%%%%%%%

%open tcp/ip connection using zeroMQ  to communicate with open ephys
switch computer
    case 'MACI64'
        url='tcp://localhost:5556'; %seems to work for mac
    case 'GLNXA64'
        url='tcp://127.0.0.1:5556'; %seems to work for linux
    case 'PCWIN64'
        url='tcp://127.0.0.1:5556'; %seems to work for linux
end
try
    SP.zhandle=zeroMQwrapper('StartConnectThread', url);
    zeroMQwrapper('Send',SP.zhandle ,'StartAcquisition');
    djMessage('successful zeroMQ connection')
    set(SP.Recordh, 'enable', 'on');
catch
    djMessage('could not open zeroMQ connection', 'error')
    pause(.5)
    SP.zhandle=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
SP.datapath=pref.datapath;
SP.CurrentStimulus=[];
SP.stimcounter=0;
SP.NStimuli=[];
SP.Name={};
SP.Description={};
SP.NProtocols=0;
SP.NRepeats=0;
SP.PPALaseron=0;
set(fig,'visible','on');
djMessage('djmaus initialized', 'append')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%InitializeGui%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



