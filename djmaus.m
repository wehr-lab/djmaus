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
%test comment

global SP pref djTimer djTimerDelay

if nargin > 0
    action = varargin{1};
else
    action = get(gcbo,'tag');
end
if isempty(action)
    action='Init';
end

% %fprintf('\naction: %s', action)
%  [uv, sv]=memory;
%  jhm = java.lang.Runtime.getRuntime.freeMemory;
%  fid=fopen('C:\lab\djmaus\memlog.txt', 'a');
%  fprintf(fid, '\n%.8f %d %d %d %d %d %d %d', ...
%      datenum(now), ...
%      uv.MaxPossibleArrayBytes, uv.MemAvailableAllArrays, uv.MemUsedMATLAB, ...
%      sv.VirtualAddressSpace.Available, sv.SystemMemory.Available, sv.PhysicalMemory.Available, ...
%      jhm);
%  fclose(fid);

switch action
    case 'Init'
        SP.user='lab';
        SP.Campulse=0;
        djPrefs;
        InitializeGUI;
        if SP.mesoscope_mode
            set(SP.Recordh, 'enable', 'on');
            djMessage('using mesoscope mode, no zmq or OpenEphys', 'append')
            fprintf('using mesoscope mode, no zmq or OpenEphys')
            SP.zhandle=[];
        elseif ~SP.mesoscope_mode
            InitZMQ %initialize zeroMQ connection to open-ephys
        end
        InitParams %initialize some default params in SP structure
        PPAdj('init');
        % set the timer
        djTimer=timer('TimerFcn',[me '(''next_stimulus'');'],'StopFcn',[me '(''restart_timer'');'],'ExecutionMode','singleShot');
        djMessage('djmaus initialized', 'append')
        
    case 'Close'
        try stop(djTimer); end
        delete(djTimer);
        clear djTimer
        PPAdj('close')
        try
            zeroMQwrapper('CloseThread',SP.zhandle);
        end
        try
            clearvars -global SP.arduino
        catch
            warning('could not close arduino')
        end
        pause(0.2)
        delete(SP.fig)
        clear global SP
        fprintf('\nbye\n')
        
    case 'Run' %which is the Play button
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
        
    case 'Stop'
        %stop acquisition, i.e. like clicking the > button on the
        %open-ephys GUI
        % This button is currently disabled
        
        zeroMQwrapper('Send',SP.zhandle ,'StopAcquisition');
        
    case 'launchOE'
        djMessage('launching OE', 'append')
        OEpath=pref.OEpath;
        system(OEpath)
        djMessage('launched', 'append')
        
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
            set(SP.Repeath, 'backgroundcolor',[.5 .5 1],'foregroundcolor',[0 0 0]);
            set(SP.Repeath, 'string','Repeat Off');
            SP.Repeat=0;
        else
            set(SP.Repeath, 'backgroundcolor',[1 .5 .5],'foregroundcolor',[0 0 1]);
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
        djMessage([int2str(current(protocol_choice)), '/' num2str(SP.NStimuli(protocol_choice)) ' stimuli, ' num2str(SP.NRepeats), ' repeats' ]);
        %this is a hack to solve the problem where after running a seamless
        %protocol, then regular protocols don't play any sound
        PPAdj('init');
        
    case 'User'
        users=get(SP.userh, 'string');
        user_index=get(SP.userh, 'value');
        user=users{user_index};
        if strcmp(user, 'add new user')
            AddUser
            user=SP.user;
        end
        djMessage(sprintf('Hi %s', user))
        SP.user=user;
        djPrefs;
        cd (pref.datapath)
        pref.datapath=pwd;
        SP.datapath=pref.datapath;
        set(SP.pathh, 'string', SP.datapath)
        %set(SP.mouseIDMenuh, 'string', pref.allmouseIDs)
        
    case 'mouseIDMenu'
        mouseIDs=get(SP.mouseIDMenuh, 'string');
        if ~iscell(mouseIDs) %skip if only one menu item
        else
            mouseID=mouseIDs{get(SP.mouseIDMenuh, 'value')};
            set(SP.mouseIDh, 'string', mouseID)
        end
        SP.mouseID=mouseID;
        LoadMouse
        
    case 'mouse2IDMenu'
        mouseIDs=get(SP.mouse2IDMenuh, 'string');
        if ~iscell(mouseIDs) %skip if only one menu item
        else
            mouseID=mouseIDs{get(SP.mouse2IDMenuh, 'value')};
            set(SP.mouse2IDh, 'string', mouseID)
        end
        SP.mouse2ID=mouseID;
        LoadMouse2
        
    case 'mouse3IDMenu'
        mouseIDs=get(SP.mouse3IDMenuh, 'string');
        if ~iscell(mouseIDs) %skip if only one menu item
        else
            mouseID=mouseIDs{get(SP.mouse3IDMenuh, 'value')};
            set(SP.mouse3IDh, 'string', mouseID)
        end
        SP.mouse3ID=mouseID;
        LoadMouse3
        
    case 'mouse4IDMenu'
        mouseIDs=get(SP.mouse4IDMenuh, 'string');
        if ~iscell(mouseIDs) %skip if only one menu item
        else
            mouseID=mouseIDs{get(SP.mouse4IDMenuh, 'value')};
            set(SP.mouse4IDh, 'string', mouseID)
        end
        SP.mouse4ID=mouseID;
        LoadMouse4
        
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
        %WriteMouseIDtoPrefs %not using prefs for this anymore 12.05.2018
        %save('mouseDB.mat', sprintf('mouseID_%s',SP.mouseID), '-append')
        LoadMouse
        
    case 'mouse2ID'
        SP.mouse2ID=get(SP.mouse2IDh, 'string');
        mouseIDs=get(SP.mouse2IDMenuh, 'string');
        if ~iscell(mouseIDs) mouseIDs={mouseIDs};end
        
        if isempty(mouseIDs)
            mouseIDs=SP.mouse2ID;
        elseif iscell(mouseIDs) & length(mouseIDs)>1 %more than one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouse2ID});
        elseif iscell(mouseIDs) & length(mouseIDs)==1 %only one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouse2ID});
            set(SP.mouseIDMenu2h, 'value',1)
        elseif length(mouseIDs)==1 %only one menu item
            mouseIDs={mouseIDs, SP.mouse2ID};
            set(SP.mouseIDMenu2h, 'value',1)
        else
            error('?')
        end
        set(SP.mouse2IDMenuh, 'string',mouseIDs)
        SP.allmouseIDs=mouseIDs;
        %WriteMouseIDtoPrefs %not using prefs for this anymore 12.05.2018
        %save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse2ID), '-append')
        LoadMouse2
        
        
    case 'mouse3ID'
        SP.mouse3ID=get(SP.mouse3IDh, 'string');
        mouseIDs=get(SP.mouse3IDMenuh, 'string');
        if ~iscell(mouseIDs) mouseIDs={mouseIDs};end
        
        if isempty(mouseIDs)
            mouseIDs=SP.mouse3ID;
        elseif iscell(mouseIDs) & length(mouseIDs)>1 %more than one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouse3ID});
        elseif iscell(mouseIDs) & length(mouseIDs)==1 %only one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouse3ID});
            set(SP.mouseIDMenu2h, 'value',1)
        elseif length(mouseIDs)==1 %only one menu item
            mouseIDs={mouseIDs, SP.mouse3ID};
            set(SP.mouseIDMenu2h, 'value',1)
        else
            error('?')
        end
        set(SP.mouse3IDMenuh, 'string',mouseIDs)
        SP.allmouseIDs=mouseIDs;
        %WriteMouseIDtoPrefs %not using prefs for this anymore 12.05.2018
        %save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse3ID), '-append')
        LoadMouse3
        
        
    case 'mouse4ID'
        SP.mouse4ID=get(SP.mouse4IDh, 'string');
        mouseIDs=get(SP.mouse4IDMenuh, 'string');
        if ~iscell(mouseIDs) mouseIDs={mouseIDs};end
        
        if isempty(mouseIDs)
            mouseIDs=SP.mouse4ID;
        elseif iscell(mouseIDs) & length(mouseIDs)>1 %more than one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouse4ID});
        elseif iscell(mouseIDs) & length(mouseIDs)==1 %only one menu item
            mouseIDs=unique({mouseIDs{:}, SP.mouse4ID});
            set(SP.mouseIDMenu2h, 'value',1)
        elseif length(mouseIDs)==1 %only one menu item
            mouseIDs={mouseIDs, SP.mouse4ID};
            set(SP.mouseIDMenu2h, 'value',1)
        else
            error('?')
        end
        set(SP.mouse4IDMenuh, 'string',mouseIDs)
        SP.allmouseIDs=mouseIDs;
        %WriteMouseIDtoPrefs %not using prefs for this anymore 12.05.2018
        %save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse4ID), '-append')
        LoadMouse4
        
    case 'SecondmouseIDbutton'
        AddSecondMouse
        
    case 'ThirdmouseIDbutton'
        AddThirdMouse
        
    case 'FourthmouseIDbutton'
        AddFourthMouse
        
    case 'Record'
        Record
        
    case 'mouseGenotype'
        SP.mouseGenotype=get(SP.mouseGenotypeh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseGenotype=''%s'';', SP.mouseID, SP.mouseGenotype);
        eval(str);
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouseID), '-append')
        
    case 'mouseSex'
        SP.mouseSex=lower(get(SP.mouseSexh, 'string'));
        set(SP.mouseSexh, 'string', SP.mouseSex)
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseSex=''%s'';', SP.mouseID, SP.mouseSex);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouseID), '-append')
        
    case 'mouseDOB'
        SP.mouseDOB=get(SP.mouseDOBh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseDOB=''%s'';', SP.mouseID, SP.mouseDOB);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouseID), '-append')
        SP.Age=(datenum(date)-datenum(SP.mouseDOB))/30; %in months
        pos1=get(SP.mouseSexh, 'pos');
        pos=get(SP.mouseDOBh, 'pos');
        pos(3)=.6*pos1(3);
        set(SP.mouseDOBh, 'pos',pos);
        pos(1)=pos(1)+pos(3);
        SP.Ageh=uicontrol(gcf,'tag','mouseAgelabel','style','text','units','pixels',...
            'string', sprintf('age\n%.1f mo', SP.Age), 'fontsize', 10,...
            'enable','inact','horiz','left','pos', pos);
        
        %mouse2
    case 'mouse2Genotype'
        SP.mouse2Genotype=get(SP.mouse2Genotypeh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseGenotype=''%s'';', SP.mouse2ID, SP.mouse2Genotype);
        eval(str);
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse2ID), '-append')
        
    case 'mouse2Sex'
        SP.mouse2Sex=lower(get(SP.mouse2Sexh, 'string'));
        set(SP.mouse2Sexh, 'string', SP.mouse2Sex)
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseSex=''%s'';', SP.mouse2ID, SP.mouse2Sex);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse2ID), '-append')
        
    case 'mouse2DOB'
        SP.mouse2DOB=get(SP.mouse2DOBh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseDOB=''%s'';', SP.mouse2ID, SP.mouse2DOB);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse2ID), '-append')
        SP.Age2=(datenum(date)-datenum(SP.mouse2DOB))/30; %in months
        pos1=get(SP.mouse2Sexh, 'pos');
        pos=get(SP.mouse2DOBh, 'pos');
        pos(3)=.6*pos1(3);
        set(SP.mouse2DOBh, 'pos',pos);
        pos(1)=pos(1)+pos(3);
        SP.Age2h=uicontrol(gcf,'tag','mouse2Agelabel','style','text','units','pixels',...
            'string', sprintf('age\n%.1f mo', SP.Age2), 'fontsize', 10,...
            'enable','inact','horiz','left','pos', pos);
        
        %mouse3
    case 'mouse3Genotype'
        SP.mouse3Genotype=get(SP.mouse3Genotypeh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseGenotype=''%s'';', SP.mouse3ID, SP.mouse3Genotype);
        eval(str);
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse3ID), '-append')
        
    case 'mouse3Sex'
        SP.mouse3Sex=lower(get(SP.mouse3Sexh, 'string'));
        set(SP.mouse3Sexh, 'string', SP.mouse3Sex)
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseSex=''%s'';', SP.mouse3ID, SP.mouse3Sex);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse3ID), '-append')
        
    case 'mouse3DOB'
        SP.mouse3DOB=get(SP.mouse3DOBh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseDOB=''%s'';', SP.mouse3ID, SP.mouse3DOB);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse3ID), '-append')
        SP.Age2=(datenum(date)-datenum(SP.mouse3DOB))/30; %in months
        pos1=get(SP.mouse3Sexh, 'pos');
        pos=get(SP.mouse3DOBh, 'pos');
        pos(3)=.6*pos1(3);
        set(SP.mouse3DOBh, 'pos',pos);
        pos(1)=pos(1)+pos(3);
        SP.Age2h=uicontrol(gcf,'tag','mouse3Agelabel','style','text','units','pixels',...
            'string', sprintf('age\n%.1f mo', SP.Age2), 'fontsize', 10,...
            'enable','inact','horiz','left','pos', pos);
        
        %mouse4
    case 'mouse4Genotype'
        SP.mouse4Genotype=get(SP.mouse4Genotypeh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseGenotype=''%s'';', SP.mouse4ID, SP.mouse4Genotype);
        eval(str);
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse4ID), '-append')
        
    case 'mouse4Sex'
        SP.mouse4Sex=lower(get(SP.mouse4Sexh, 'string'));
        set(SP.mouse4Sexh, 'string', SP.mouse4Sex)
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseSex=''%s'';', SP.mouse4ID, SP.mouse4Sex);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse4ID), '-append')
        
    case 'mouse4DOB'
        SP.mouse4DOB=get(SP.mouse4DOBh, 'string');
        cd (pref.root)
        try load mouseDB
        end
        str=sprintf('mouseID_%s.mouseDOB=''%s'';', SP.mouse4ID, SP.mouse4DOB);
        eval(str)
        save('mouseDB.mat', sprintf('mouseID_%s',SP.mouse4ID), '-append')
        SP.Age2=(datenum(date)-datenum(SP.mouse4DOB))/30; %in months
        pos1=get(SP.mouse4Sexh, 'pos');
        pos=get(SP.mouse4DOBh, 'pos');
        pos(3)=.6*pos1(3);
        set(SP.mouse4DOBh, 'pos',pos);
        pos(1)=pos(1)+pos(3);
        SP.Age2h=uicontrol(gcf,'tag','mouse4Agelabel','style','text','units','pixels',...
            'string', sprintf('age\n%.1f mo', SP.Age2), 'fontsize', 10,...
            'enable','inact','horiz','left','pos', pos);
        
    case 'Reinforcement'
        SP.Reinforcement=get(SP.Reinforcementh, 'string');
        
    case 'Drugs'
        SP.Drugs=get(SP.Drugsh, 'string');
        
    case 'LaserPower'
        SP.LaserPower=get(SP.LaserPowerh, 'string');
        
    case 'Depth'
        SP.Depth=get(SP.Depthh, 'string');
        
    case 'Notes'
        SP.Notes=get(SP.Notesh, 'string');
        
    case 'ResetZMQ'
        InitZMQ;
        
    case 'TestZMQ'
        zeroMQwrapper('Send',SP.zhandle ,'Test');
        djMessage(sprintf('sent "Test" message to %s', pref.zmqurl))
        
    case 'LaserOnOff'
        SP.LaserOnOff=get(SP.LaserOnOffh, 'value');
        if SP.LaserOnOff
            set(SP.LaserNumPulsesh, 'enable', 'on')
            set(SP.LaserWidthh, 'enable', 'on')
            set(SP.LaserStarth, 'enable', 'on')
            set(SP.LaserOnOffh, 'backgroundcolor',[1 0 0],'foregroundcolor',[0 0 0]);
            set(SP.LaserOnOffh, 'string','Laser is ON');
            djmaus('LaserNumPulses')
            
        else
            set(SP.LaserStarth, 'enable', 'off')
            set(SP.LaserNumPulsesh, 'enable', 'off')
            set(SP.LaserISIh, 'enable', 'off')
            set(SP.LaserWidthh, 'enable', 'off')
            set(SP.LaserOnOffh, 'backgroundcolor',[.5 .5 .5],'foregroundcolor',[0 0 1]);
            set(SP.LaserOnOffh, 'string','Laser is OFF');
        end
        
    case 'LaserStart'
        SP.LaserStart=str2num(get(SP.LaserStarth, 'string'));
        
    case 'LaserWidth'
        SP.LaserWidth=str2num(get(SP.LaserWidthh, 'string'));
        
    case 'LaserNumPulses'
        SP.LaserNumPulses=str2num(get(SP.LaserNumPulsesh, 'string'));
        if SP.LaserNumPulses<1
            djMessage('number of laser pulses must be at least 1. Turn Laser off to get 0 pulses');
            SP.LaserNumPulses=1;
            set(SP.LaserNumPulsesh, 'string', '1');
        elseif SP.LaserNumPulses==1
            set(SP.LaserISIh, 'enable', 'off')
        elseif SP.LaserNumPulses>1
            set(SP.LaserISIh, 'enable', 'on')
            if strcmp(get(SP.LaserISIh, 'string'),'0')
                set(SP.LaserISIh, 'string', SP.LaserWidth)
            end
        end
        
    case 'LaserISI'
        SP.LaserISI=str2num(get(SP.LaserISIh, 'string'));
        
    case 'camerapulse'
        if SP.Campulse
            %we want to stop
            set(SP.camerapulse,'backgroundcolor',[0 0.9 0],'String','Camera Stopped');
            SP.Campulse=0;
            PPAdj('camerapulse_off')
        else
            set(SP.camerapulse, 'backgroundcolor',[0.9 0 0],'String','Camera Recording');
            SP.Campulse=1;
            PPAdj('camerapulse_on')
        end
        
    case 'mesoscope_mode' %toggle mesoscope/OpenEphys mode
        if get(SP.MesoModeh, 'value') %1=meso, 0=intan/OE
            set(SP.MesoModeh, 'String', 'Mesoscope');
            SP.mesoscope_mode=1;
            
        elseif ~get(SP.MesoModeh, 'value') %set to OE
            set(SP.MesoModeh, 'String', 'OpenEphys');
            SP.mesoscope_mode=0;
        end
        
end

function AddUser
global SP pref
[username]=inputdlg('please enter new user name. It''s recommended to keep it short.', 'Add new user');
username=username{:};
if any(strcmp(username, pref.users))
    errordlg(sprintf('user name %s already taken', username));
    set(SP.userh, 'value', 1)
else
    SP.user=username;
    cd(pref.datapath)
    cd ..
    mkdir(username)
    cd(username)
    pref.datapath=pwd;
    pref.users{end+1}=username;
    
    [pathstr,name,~] =fileparts(pref.remotedatapath);
    pref.remotedatapath=fullfile(pathstr, username);
    
    %write new pref.users
    cd(pref.root)
    fid=fopen('djPrefs.m', 'a+');
    key=sprintf('pref.users={');
    Preftext = regexp( fileread('djPrefs.m'), '\n', 'split');
    fclose(fid)
    
    I=strmatch(key, Preftext);
    if ~isempty(I) %found key, overwrite with revised entry
        I=I(1);
        newstr=sprintf('''%s'',', pref.users{:});
        newstr=newstr(1:end-1);
        newstr2=sprintf('pref.users={%s};', newstr);
        Preftext{I}=newstr2; %change entry;
        fid = fopen('djPrefs.m', 'w');
        fprintf(fid, '%s\n', Preftext{:});
        fclose(fid);
    else
    end
    
    
    %write new individual user prefs
    fid=fopen('djPrefs.m', 'a+');
    key=sprintf('switch SP.user');
    Preftext = regexp( fileread('djPrefs.m'), '\n', 'split');
    fclose(fid);
    I=strmatch(key, strtrim(Preftext));
    if ~isempty(I) %found key, overwrite with revised entry
        I=I(1);
        str1=sprintf('case ''%s''', username);
        str2=sprintf('pref.datapath=''%s'';', pref.datapath);
        str3=sprintf('pref.remotedatapath=''%s'';', pref.remotedatapath);
        Preftext_copy=Preftext;
        Preftext{I+1}=str1;
        Preftext{I+2}=str2;
        Preftext{I+3}=str3;
        for i=I+4:length(Preftext_copy)+2
            Preftext{i}=Preftext_copy{i-3};
        end
        fid = fopen('djPrefs.m', 'w');
        fprintf(fid, '%s\n', Preftext{:});
        fclose(fid);
        
        
    end
    userstr=pref.users;
    userstr{end+1}='add new user';
    set(SP.userh,'string', userstr);
    
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMouse
global SP pref
cd (pref.root)
set(SP.mouseIDh, 'ForegroundColor', 'k')
try
    mouseDB=load('mouseDB.mat');
    if isfield(mouseDB, ['mouseID_', SP.mouseID])
        try
            str=sprintf('SP.mouseGenotype=mouseDB.mouseID_%s.mouseGenotype;', SP.mouseID);
            eval(str);
            set(SP.mouseGenotypeh, 'string', SP.mouseGenotype);
        end
        try
            str=sprintf('SP.mouseDOB=mouseDB.mouseID_%s.mouseDOB;', SP.mouseID);
            eval(str);
            set(SP.mouseDOBh, 'string', SP.mouseDOB);
        end
        try
            str=sprintf('SP.mouseSex=mouseDB.mouseID_%s.mouseSex;', SP.mouseID);
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
catch
    djMessage('could not load mouse database')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMouse2
global SP pref
cd (pref.root)
try
    mouseDB=load('mouseDB.mat');
    if isfield(mouseDB, ['mouseID_', SP.mouse2ID])
        try
            str=sprintf('SP.mouse2Genotype=mouseDB.mouseID_%s.mouseGenotype;', SP.mouse2ID);
            eval(str);
            set(SP.mouse2Genotypeh, 'string', SP.mouse2Genotype);
        end
        try
            str=sprintf('SP.mouse2DOB=mouseDB.mouseID_%s.mouseDOB;', SP.mouse2ID);
            eval(str);
            set(SP.mouse2DOBh, 'string', SP.mouse2DOB);
        end
        try
            str=sprintf('SP.mouse2Sex=mouseDB.mouseID_%s.mouseSex;', SP.mouse2ID);
            eval(str);
            set(SP.mouse2Sexh, 'string', SP.mouse2Sex);
        end
    else
        SP.mouse2Genotype='genotype unknown';
        set(SP.mouse2Genotypeh, 'string', SP.mouse2Genotype);
        SP.mouse2Sex='sex unknown';
        set(SP.mouse2Sexh, 'string', SP.mouse2Sex);
        SP.mouse2DOB='age unknown';
        set(SP.mouse2DOBh, 'string', SP.mouse2DOB);
    end
catch
    djMessage('could not load mouse database')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMouse3
global SP pref
cd (pref.root)
try
    mouseDB=load('mouseDB.mat');
    if isfield(mouseDB, ['mouseID_', SP.mouse3ID])
        try
            str=sprintf('SP.mouse3Genotype=mouseDB.mouseID_%s.mouseGenotype;', SP.mouse3ID);
            eval(str);
            set(SP.mouse3Genotypeh, 'string', SP.mouse3Genotype);
        end
        try
            str=sprintf('SP.mouse3DOB=mouseDB.mouseID_%s.mouseDOB;', SP.mouse3ID);
            eval(str);
            set(SP.mouse3DOBh, 'string', SP.mouse3DOB);
        end
        try
            str=sprintf('SP.mouse3Sex=mouseDB.mouseID_%s.mouseSex;', SP.mouse3ID);
            eval(str);
            set(SP.mouse3Sexh, 'string', SP.mouse3Sex);
        end
    else
        SP.mouse3Genotype='genotype unknown';
        set(SP.mouse3Genotypeh, 'string', SP.mouse3Genotype);
        SP.mouse3Sex='sex unknown';
        set(SP.mouse3Sexh, 'string', SP.mouse3Sex);
        SP.mouse3DOB='age unknown';
        set(SP.mouse3DOBh, 'string', SP.mouse3DOB);
    end
catch
    djMessage('could not load mouse database')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadMouse4
global SP pref
cd (pref.root)
try
    mouseDB=load('mouseDB.mat');
    if isfield(mouseDB, ['mouseID_', SP.mouse4ID])
        try
            str=sprintf('SP.mouse4Genotype=mouseDB.mouseID_%s.mouseGenotype;', SP.mouse4ID);
            eval(str);
            set(SP.mouse4Genotypeh, 'string', SP.mouse4Genotype);
        end
        try
            str=sprintf('SP.mouse4DOB=mouseDB.mouseID_%s.mouseDOB;', SP.mouse4ID);
            eval(str);
            set(SP.mouse4DOBh, 'string', SP.mouse4DOB);
        end
        try
            str=sprintf('SP.mouse4Sex=mouseDB.mouseID_%s.mouseSex;', SP.mouse4ID);
            eval(str);
            set(SP.mouse4Sexh, 'string', SP.mouse4Sex);
        end
    else
        SP.mouse4Genotype='genotype unknown';
        set(SP.mouse4Genotypeh, 'string', SP.mouse4Genotype);
        SP.mouse4Sex='sex unknown';
        set(SP.mouse4Sexh, 'string', SP.mouse4Sex);
        SP.mouse4DOB='age unknown';
        set(SP.mouse4DOBh, 'string', SP.mouse4DOB);
    end
catch
    djMessage('could not load mouse database')
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
function LoadProtocol(varargin)
% adding capability to push files programatically
% call dj('load', 'fullfilename')
% where fullfilename includes the absolute path (e.g. 'D:\lab\exper2.2\protocols\Tuning Curve protocols\tuning-curve-tones-20f_1000-20000Hz-1a_80-80dB-1d_400ms-isi500ms.mat'
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
        set(SP.Runh, 'enable', 'on','backgroundcolor',[0 0.9 0])
        
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
stimulus.param=CalibrateSound(stimulus.param, stimulus.type);
samples=feval(fcn,stimulus.param,SP.SoundFs);

%LoadPPA(type,where,param)
PPAdj('load', 'var', samples, stimulus.param);
str=sprintf('TrialType %s', stimulus.stimulus_description);
%append a field saying whether the LaserON/OFF button is clicked or not:
str=sprintf('%s %s:%g', str, 'LaserOnOff', SP.LaserOnOff);
%note: the maximum length of this string is 255. Any longer and open ephys
%will fail to write the text to messages.events (although the time stamp is OK)
if length(str)>255 fprintf('\n\n\n\n\n\n');warning('TrialType message string is too long!!!! (try re-making this stimulus protocol)');end

if ~isempty(SP.zhandle)
    zeroMQwrapper('Send', SP.zhandle, str)
end
if SP.mesoscope_mode
    SendArduinoSoundcardtrigger
end
PPAdj('playsound')
UpdateStimlog(stimulus);
djMessage(stimulus.stimulus_description, 'append');
status = PsychPortAudio('GetStatus', SP.PPAhandle);

% if ~status.Active & status.XRuns
% fprintf('status: %d Xruns: %d\n',status.Active,status.XRuns)
% end


% figure(100)
% plot(samples)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  SendArduinoSoundcardtrigger
global SP

%the timer object is created during Init
%here we explicitly set the timer callback to terminate the SCT pulse, and
%then turn the pulse on and start the timer.
SP.arduinoTimer.TimerFcn=@(~,~)writeDigitalPin(SP.arduino, SP.arduinoSCTpin, 0);
writeDigitalPin(SP.arduino, SP.arduinoSCTpin, 1);
start(SP.arduinoTimer)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  SendArduinoProtocolStart
global SP

%the timer object is created during Init
%here we explicitly set the timer callback to terminate the pulse, and
%then turn the pulse on and start the timer.
SP.arduinoTimer.TimerFcn=@(~,~)writeDigitalPin(SP.arduino, SP.arduinoProtocolStartpin, 0);
writeDigitalPin(SP.arduino, SP.arduinoProtocolStartpin, 1);
start(SP.arduinoTimer)



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
        set(SP.Runh,'backgroundcolor',[0 0.9 0],'String','Play', 'value', 0);
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

function fcn=StimTypes(type)
switch type
    case 'square'   % Kip 9/5/24
        fcn = 'MakeSquare';
    case 'tone'
        fcn='MakeTone';
    case 'whitenoise'
        fcn='MakeWhiteNoise';
    case 'AsymGPIAS'
        fcn='MakeAsymGPIAS';
    case 'GPIAS'
        fcn='MakeGPIAS';
    case 'ASR'
        fcn='MakeASR';
    case 'toneGPIAS'
        fcn='MakeToneGPIAS';
    case 'noise'
        fcn='MakeNoise';
    case 'clicktrain'
        fcn='MakeClickTrain';
    case 'silentsound'
        fcn='MakeSilentSound';
    case '2tone'
        fcn='Make2Tone';
    case 'soundfile'
        fcn='LoadSoundfile';
    case 'AMNoise'
        fcn='MakeAMNoise';
    otherwise djMessage(sprintf('%s not recognized', type))
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Record
%toggle Open Ephys recording state, increment datafile
global SP pref
fprintf('\n')

if SP.Record
    %we want to stop recording
    %courtesy message if continuous stimuli haven't finished yet
    try
        PPAhandle=SP.PPAhandle;
        status = PsychPortAudio('GetStatus', PPAhandle);
        
        if  status.Active==1; %device running
            protocol=SP.ProtocolIndex;
            current=SP.CurrentStimulus(protocol);
            inQueue=current-status.SchedulePosition-SP.CurrStimatPPAstart;
            Question=sprintf('stimuli have not finished playing yet (%d remaining).\nDo you really want to stop acquisition?', inQueue);
            ButtonName = questdlg(Question, 'Are you sure?', 'Cancel', 'Really Stop', 'Cancel');
            if strcmp(ButtonName, 'Cancel')
                %unclick the record button
                set(findobj('Tag', 'Record'), 'Value', 1);
                return
            end
            
        end
    end
    
    if ~SP.mesoscope_mode
        zeroMQwrapper('Send',SP.zhandle ,'StopRecord');
    end
    set(SP.Recordh,'backgroundcolor',[0 0.9 0],'String','Record');
    SP.Record=0;
    set(SP.mouseIDh, 'enable', 'on');
    set(SP.mouseIDMenuh, 'enable', 'on');
    if isfield(SP, 'mouse2ID')
        set(SP.mouse2IDh, 'enable', 'on');
        set(SP.mouse2IDMenuh, 'enable', 'on');
    end
    if isfield(SP, 'mouse3ID')
        set(SP.mouse3IDh, 'enable', 'on');
        set(SP.mouse3IDMenuh, 'enable', 'on');
    end
    if isfield(SP, 'mouse4ID')
        set(SP.mouse4IDh, 'enable', 'on');
        set(SP.mouse4IDMenuh, 'enable', 'on');
    end
    if pref.camera_rec==1 %set the camera to be turned off when stopped recording
        SP.Campulse=0;
        set(SP.camerapulse, 'backgroundcolor',[0 0.9 0],'String','Camera Stopped');
        PPAdj('camerapulse_off')
    else
        %fprintf('\n camera not recording\n')
    end
    
    UpdateNotebookFile
    try
        set(SP.pathh, 'string', {SP.datapath, [SP.activedir, ' finished']})
    end
else
    %we want to start recording;
    %try stopping
    %zeroMQwrapper('Send',SP.zhandle ,'StopRecord');
    
    %disable play button here, to avoid delivering stimuli before notebook
    %is initialized (which could result in a skipped stimlog entry)
    set(SP.Runh, 'enable', 'off')
    if SP.mesoscope_mode
        SendArduinoProtocolStart
    elseif ~SP.mesoscope_mode %open-ephys/bonsai specific stuff
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% adding in feature to find bonsai-created folder  - Nick 2/7/2021
        destination = pref.datapath;
        cd(destination)
        listoffolders = dir; thedates = arrayfun(@(x) datetime(x.date), listoffolders);
        [~,sidx] = sort(thedates,'descend');
        for i = 1:6 %take latest 6 (b/c sometimes '.','..',thumbs,DS,camerafiles,etc may be most recent)
            foldername = listoffolders(sidx(i)).name;
            if length(foldername)>15 %this will be the directory
                break
            end
        end
    end
    RecDir = strcat(pref.remotedatapath,'\',foldername);
    if pref.local
        SP.bonsaiFolder = strcat(pref.datahost, '\', RecDir);
    else
        SP.bonsaiFolder = strrep(strcat(pref.datahost, '\', RecDir), ':', '');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%% end of feature - Nick 2/7/2021
    startstr=sprintf('StartRecord');
    if ~isfield(SP, 'mouseID')
        SP.mouseID='none';
    end
    zeroMQwrapper('Send',SP.zhandle ,'StartAcquisition'); %shouldn't need to do this unless user stopped acquisition, doesn't hurt anyway
    % startstr=sprintf('StartRecord CreateNewDir=1 RecDir=%s AppendText=mouse-%s', pref.remotedatapath, SP.mouseID);
    %startstr=sprintf('StartRecord CreateNewDir=1 RecDir=%s AppendText=mouse-%s', [pref.remotedatapath,SP.user], SP.mouseID);
    %trying to fix double user name bug mw 05032018
    %startstr=sprintf('StartRecord CreateNewDir=1 RecDir=%s AppendText=mouse-%s', [pref.remotedatapath], SP.mouseID);
    startstr=sprintf('StartRecord CreateNewDir=1 RecDir=%s AppendText=mouse-%s', RecDir, SP.mouseID); %Nick's change 2/7/2021
    fprintf('\nsending this message to Open Ephys:\n%s\n', startstr)
    zeroMQwrapper('Send',SP.zhandle ,startstr);
    
    set(SP.Recordh, 'backgroundcolor',[0.9 0 0],'String','Recording...');
    set(SP.mouseIDh, 'enable', 'off');
    set(SP.mouseIDMenuh, 'enable', 'off');
    if isfield(SP, 'mouse2ID')
        set(SP.mouse2IDh, 'enable', 'off');
        set(SP.mouse2IDMenuh, 'enable', 'off');
    end
    if isfield(SP, 'mouse3ID')
        set(SP.mouse3IDh, 'enable', 'off');
        set(SP.mouse3IDMenuh, 'enable', 'off');
    end
    if isfield(SP, 'mouse4ID')
        set(SP.mouse4IDh, 'enable', 'off');
        set(SP.mouse4IDMenuh, 'enable', 'off');
    end
    SP.Record=1;
    SP.stimcounter=0;
    if isfield(SP, 'stimlog')
        SP=rmfield( SP, 'stimlog');
    end
    if pref.camera_rec==1 %set the camera to be triggered when you start recording
        SP.Campulse=1;
        set(SP.camerapulse, 'backgroundcolor',[0.9 0 0],'String','Camera Recording');
        PPAdj('camerapulse_on')
    else
        fprintf('\n camera not recording\n')
    end
    InitNotebookFile
    
    %re-enable play button here
    set(SP.Runh, 'enable', 'on')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitNotebookFile
global SP nb pref

if SP.mesoscope_mode %1=mesoscope mode, save notebook locally, no OE or zmq
    dirname=sprintf('%s_mouse-%s', datestr(now, 'yyyy-mm-dd_HH-MM-SS'), SP.mouseID);
    SP.activedir=dirname;
    cd(SP.datapath)
    mkdir(SP.activedir)
elseif ~SP.mesoscope_mode %0 = classic OpenEphys mode, use zmq to get recording directory and save notebook there
    
    
    % find active OE data directory and cd into it
    SP.activedir='unknown';
    try
        %    zeroMQwrapper('Send',SP.zhandle ,sprintf('ChangeDirectory %s', pref.root))
        %    pause(.2)
        zeroMQwrapper('Send',SP.zhandle ,'GetRecordingPath');
        pause(1)
        RecordingPath = zeroMQwrapper('GetReply',SP.zhandle);
        
        wb1=waitbar(0, 'zeroMQwrapper:GetReply did not return a reply, re-trying ...', 'visible', 'off');
        while ~exist('RecordingPath','var')
            %sometimes it doesn't work for some reason, and we get
            %"zmq wrapper GetReply: there is no reply available"
            %so we can re-try and see if that helps
            %mw 04.22.2019
            w=0;
            
            set(wb1, 'units', 'pixels', 'visible', 'on');
            pos=get(wb1, 'pos');
            set(wb1, 'pos', [pref.windowpos(1),pref.windowpos(2)+pref.windowpos(4)-2*pos(4), pos(3), pos(4)]);
            
            fprintf('\nzeroMQwrapper:GetReply did not return a reply, re-trying ...')
            pause(2)
            w=w+.1;
            waitbar(mod(w, 1),wb1)
            
            RecordingPath = zeroMQwrapper('GetReply',SP.zhandle);
        end
        close(wb1)
        
        %     cd(pref.root)
        %     fid=fopen('RecordingPath.txt', 'r');
        %     RecordingPath=fgetl(fid);
        %     RecordingPathSize=str2num(fgetl(fid));
        %     fclose(fid);
        %     %hack: on windows I am still getting extra characters -> trim to size
        %     RecordingPath=RecordingPath(1:RecordingPathSize);
        %     fprintf('\ndjmaus: read this Recording Path from file:%s', RecordingPath)
        
        %     SP.activedir=RecordingPath;
        SP.activedir=fullfile(pref.datahost, RecordingPath); %functional 2/8/2021 Nick
        
        %%%%%%%%%%%%%%%% for the new (6/23, Kip) OEphys GUI, need to add one more directory layer:
        cd(SP.bonsaiFolder)
        listoffolders = dir; thedates = arrayfun(@(x) datetime(x.date), listoffolders);
        [~,sidx] = sort(thedates,'descend');
        for i = 1:6 %take latest 6 (b/c sometimes '.','..',thumbs,DS,camerafiles,etc may be most recent)
            OEFolder = listoffolders(sidx(i)).name;
            if length(OEFolder)>15 %this will be the directory
                break
            end
        end
        SP.activedir = strcat(SP.bonsaiFolder,'\',OEFolder);
        %%%%%%%%%%%%%%%
        
        
        %SP.activedir=strrep(SP.activedir, ':', ''); %commenting out to run on
        %single-machine windows configuration - was this important for
        %2-machine config?? mw 04.11.2017
        
        if strcmp(SP.activedir(1:6), 'o:\d:\') %hack mw 080217
            SP.activedir=SP.activedir([1:3 7:end]);
        end
        if strcmp(SP.activedir(1:6), 'n:\e:\') %hack mw 081518
            SP.activedir=SP.activedir([1:3 7:end]);
        end
        
        if ~pref.local
            SP.activedir=strrep(SP.activedir, ':', '');
        end
        
        d=dir(SP.activedir);
        if isempty(d)
            w=0;
            wb=waitbar(0, 'waiting for data directory to mount...');
            set(wb, 'units', 'pixels');
            pos=get(wb, 'pos');
            set(wb, 'pos', [pref.windowpos(1),pref.windowpos(2)+pref.windowpos(4)-2*pos(4), pos(3), pos(4)]);
            while isempty(d)
                d=dir(SP.activedir);
                pause(.1)
                
                w=w+.1;
                waitbar(mod(w, 1),wb1)
                
                RecordingPath = zeroMQwrapper('GetReply',SP.zhandle )
            end
            close(wb1)
            
            %     cd(pref.root)
            %     fid=fopen('RecordingPath.txt', 'r');
            %     RecordingPath=fgetl(fid);
            %     RecordingPathSize=str2num(fgetl(fid));
            %     fclose(fid);
            %     %hack: on windows I am still getting extra characters -> trim to size
            %     RecordingPath=RecordingPath(1:RecordingPathSize);
            %     fprintf('\ndjmaus: read this Recording Path from file:%s', RecordingPath)
            
            %     SP.activedir=RecordingPath;
            SP.activedir=fullfile(pref.datahost, RecordingPath); %functional 2/8/2021 Nick
            %SP.activedir=strrep(SP.activedir, ':', ''); %commenting out to run on
            %single-machine windows configuration - was this important for
            %2-machine config?? mw 04.11.2017
            
            if strcmp(SP.activedir(1:6), 'o:\d:\') %hack mw 080217
                SP.activedir=SP.activedir([1:3 7:end]);
            end
            if strcmp(SP.activedir(1:6), 'n:\e:\') %hack mw 081518
                SP.activedir=SP.activedir([1:3 7:end]);
            end
            
            if ~pref.local
                SP.activedir=strrep(SP.activedir, ':', '');
            end
            
            d=dir(SP.activedir);
            if isempty(d)
                w=0;
                wb=waitbar(0, 'waiting for data directory to mount...');
                set(wb, 'units', 'pixels');
                pos=get(wb, 'pos');
                set(wb, 'pos', [pref.windowpos(1),pref.windowpos(2)+pref.windowpos(4)-2*pos(4), pos(3), pos(4)]);
                while isempty(d)
                    d=dir(SP.activedir);
                    pause(.1)
                    w=w+.1;
                    waitbar(mod(w, 1),wb)
                end
                close(wb)
            end
        end
    catch
        
        lasterr
        %close(wb1)
        errordlg('Go get Mike', '', 'modal');
        keyboard
        
        fprintf('\nCould not create notebook file in active data directory')
        %ask user if they want to manually save notebook file
        ButtonName = questdlg('Could not create notebook file in active data directory. Do you want to manually save the notebook file?');
        switch ButtonName,
            case 'Yes'
                targetdir = uigetdir(SP.datapath, 'Select directory in which to save notebook file.')
                cd(targetdir)
                SP.activedir=pwd;
                nb.activedir=pwd;
            case 'No'
            case 'Cancel'
                return
        end
    end
end %if mesoscope mode

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
nb.Drugs=SP.Drugs;
nb.notes=SP.Notes;
nb.Reinforcement=SP.Reinforcement;

if isfield(SP, 'mouse2ID')
    nb.mouse2ID=SP.mouse2ID;
    nb.mouse2DOB=SP.mouse2DOB;
    nb.mouse2Sex=SP.mouse2Sex;
    nb.mouse2Genotype=SP.mouse2Genotype;
end
if isfield(SP, 'mouse3ID')
    nb.mouse3ID=SP.mouse3ID;
    nb.mouse3DOB=SP.mouse3DOB;
    nb.mouse3Sex=SP.mouse3Sex;
    nb.mouse3Genotype=SP.mouse3Genotype;
end
if isfield(SP, 'mouse4ID')
    nb.mouse4ID=SP.mouse4ID;
    nb.mouse4DOB=SP.mouse4DOB;
    nb.mouse4Sex=SP.mouse4Sex;
    nb.mouse4Genotype=SP.mouse4Genotype;
end

save('notebook.mat', 'nb')
fprintf('\ncreated notebook file in %s', nb.activedir)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateStimlog(stimulus)
global SP

if SP.Record
    try
        SP.stimcounter=SP.stimcounter+1;
        timestamp=datestr(now,'mmmm-dd-yyyy HH:MM:SS.FFF');
        stimulus.timestamp=timestamp;
        stimulus.LaserOnOff=SP.LaserOnOff;
        if stimulus.LaserOnOff
            stimulus.LaserStart=SP.LaserStart;
            stimulus.LaserWidth=SP.LaserWidth;
            stimulus.LaserNumPulses=SP.LaserNumPulses;
            stimulus.LaserISI=SP.LaserISI;
        else
            stimulus.LaserStart=[];
            stimulus.LaserWidth=[];
            stimulus.LaserNumPulses=[];
            stimulus.LaserISI=[];
        end
        cd(SP.datapath)
        cd(SP.activedir)
        SP.stimlog(SP.stimcounter)=stimulus;
        
        %write json encoded stimulus to stimlog txt file
        fid=fopen('stimlog.txt', 'a');
        str=jsonencode(stimulus);
        fprintf(fid, '\n%s', str);
        fclose(fid);
        
        %old way: save stimlog to notebook.mat (takes progressively longer over session)
        %stimlog=SP.stimlog;
        %save('notebook.mat', '-append', 'stimlog');
        
    catch
        fprintf('\nCould not update stimlog in notebook file in active data directory');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function UpdateNotebookFile
global SP nb
nb.user=SP.user;
nb.mouseID=SP.mouseID;
nb.Depth=SP.Depth;
nb.datapath=SP.datapath ;
nb.activedir=SP.activedir;
nb.LaserPower=SP.LaserPower;
nb.mouseDOB=SP.mouseDOB;
nb.mouseSex=SP.mouseSex;
nb.mouseGenotype=SP.mouseGenotype;
nb.notes=SP.Notes;
nb.Drugs=SP.Drugs;
nb.Reinforcement=SP.Reinforcement;

if isfield(SP, 'mouse2ID')
    nb.mouse2ID=SP.mouse2ID;
    nb.mouse2DOB=SP.mouse2DOB;
    nb.mouse2Sex=SP.mouse2Sex;
    nb.mouse2Genotype=SP.mouse2Genotype;
end
if isfield(SP, 'mouse3ID')
    nb.mouse3ID=SP.mouse3ID;
    nb.mouse3DOB=SP.mouse3DOB;
    nb.mouse3Sex=SP.mouse3Sex;
    nb.mouse3Genotype=SP.mouse3Genotype;
end
if isfield(SP, 'mouse4ID')
    nb.mouse4ID=SP.mouse4ID;
    nb.mouse4DOB=SP.mouse4DOB;
    nb.mouse4Sex=SP.mouse4Sex;
    nb.mouse4Genotype=SP.mouse4Genotype;
end

if ~isfield(SP, 'stimlog') %perhaps no stimuli were even delivered
    stimlog=[];
else
    stimlog=SP.stimlog;
end
try
    cd(nb.datapath)
    cd(nb.activedir)
    
    %read stim log from json file and convert to stim structure
    %not necessary since we can just use SP.stimlog
    %     fid=fopen('stimlog.txt', 'r');
    %     stimidx=0;
    %     while 1
    %         str=fgetl(fid);
    %         if ~ischar(str), break,
    %         else
    %             stimulus=jsondecode(str);
    %             stimidx=stimidx+1;
    %             stimlog(stimidx)=stimulus;
    %         end
    %         fclose(fid)
    %     end
    
    save('notebook.mat', '-append', 'nb', 'stimlog')
    json_nb_text=jsonencode(nb);
    json_stimlog_text=jsonencode(stimlog);
    fid=fopen('notebook.json', 'w');
    fprintf(fid, '%s', json_nb_text);
    fclose(fid);
    fid=fopen('stimlog.json', 'w');
    fprintf(fid, '%s', json_stimlog_text);
    fclose(fid);
    
    
    fprintf('\nupdated notebook file in %s', nb.activedir)
catch
    fprintf('\nCould not update notebook file.')
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitZMQ
global SP pref

%open tcp/ip connection using zeroMQ  to communicate with openephys
cd (pref.root)
%cd(pref.mexpath)
try
    %zeroMQwrapper('CloseThread', url); %crashes matlab
    SP.zhandle=zeroMQwrapper('StartConnectThread', pref.zmqurl);
    pause(.2)
    zeroMQwrapper('Send',SP.zhandle ,'StartAcquisition');
    djMessage(sprintf('successful zeroMQ connection to %s, handle %d', pref.zmqurl, SP.zhandle))
    fprintf('successful zeroMQ connection %d', SP.zhandle)
    set(SP.Recordh, 'enable', 'on');
catch
    djMessage('could not open zeroMQ connection', 'error')
    fprintf('could not open zeroMQ connection')
    pause(.5)
    SP.zhandle=[];
end
fprintf('\n')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function stimparam=CalibrateSound(stimparam, stimtype);
global SP
cal=SP.cal;
if ~isempty(cal) %it will be empty if Init failed to load calibration
    if strcmp(stimtype, '2tone') %special case since 2tone has both a frequency and a probefreq
        try
            findex=find(cal.logspacedfreqs<=stimparam.frequency, 1, 'last');
            atten=cal.atten(findex);
            stimparam.amplitude=stimparam.amplitude-atten;
            
            findex=find(cal.logspacedfreqs<=stimparam.probefreq, 1, 'last');
            atten=cal.atten(findex);
            stimparam.probeamp=stimparam.probeamp-atten;
            
            djMessage( 'calibrated', 'append')
        catch
            djMessage( 'NOT calibrated', 'append')
            pause(.1)
        end
        
        
    elseif isfield(stimparam, 'frequency') %it has a freq and therefore is calibratable by frequency
        try
            findex=find(cal.logspacedfreqs<=stimparam.frequency, 1, 'last');
            atten=cal.atten(findex);
            switch stimtype
                case 'bintone'
                    Ratten=cal.Ratten(findex);
                    Latten=cal.Latten(findex);
                    stimparam.Ramplitude=stimparam.Ramplitude-Ratten;
                    stimparam.Lamplitude=stimparam.Lamplitude-Latten;
                otherwise
                    stimparam.amplitude=stimparam.amplitude-atten;
            end
            djMessage( 'calibrated', 'append')
        catch
            djMessage( 'NOT calibrated', 'append')
            pause(.1)
        end
        
    else
        switch stimtype
            case {'clicktrain', 'whitenoise', 'AMNoise'} %stimuli that consist of white noise
                try
                    findex=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex);
                    switch stimtype
                        case 'binwhitenoise'
                            Ratten=cal.Ratten(findex);
                            Latten=cal.Latten(findex);
                            stimparam.Ramplitude=stimparam.Ramplitude-Ratten;
                            stimparam.Lamplitude=stimparam.Lamplitude-Latten;
                        otherwise
                            stimparam.amplitude=stimparam.amplitude-atten;
                    end
                    djMessage( sprintf('calibrated'), 'append')
                catch
                    djMessage( 'NOT calibrated', 'append');pause(.5)
                end
            case {'fmtone'} %stimuli that have a carrier frequency
                try
                    findex=find(cal.logspacedfreqs<=stimparam.carrier_frequency, 1, 'last');
                    atten=cal.atten(findex);
                    stimparam.amplitude=stimparam.amplitude-atten;
                    djMessage( 'calibrated', 'append')
                catch
                    djMessage( 'NOT calibrated', 'append');pause(.5)
                end
            case {'noise'} %narrow-band noise stimuli (use center frequency calibration)
                try
                    findex=find(cal.logspacedfreqs<=stimparam.center_frequency, 1, 'last');
                    atten=cal.atten(findex);
                    stimparam.amplitude=stimparam.amplitude-atten;
                    djMessage( 'calibrated', 'append')
                catch
                    djMessage( 'NOT calibrated', 'append')
                end
            case {'GPIAS', 'AsymGPIAS'} %startle pulse (use whitenoise calibration)
                %note: GPIAS is now whitenoise-based (not band-limited noise)
                try
                    findex=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex);
                    stimparam.amplitude=stimparam.amplitude-atten;
                    stimparam.pulseamp=stimparam.pulseamp-atten;
                    djMessage( 'calibrated', 'append')
                catch
                    djMessage( 'NOT calibrated', 'append')
                end
            case {'ASR'} %startle pulse (use whitenoise calibration)
                try
                    findex=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex);
                    stimparam.prepulseamp=stimparam.prepulseamp-atten;
                    stimparam.pulseamp=stimparam.pulseamp-atten;
                    
                    djMessage( 'calibrated', 'append')
                catch
                    djMessage( 'NOT calibrated', 'append')
                end
            case {'NBASR'} %startle pulse (use whitenoise calibration)
                %plus narrow-band noise pulse (use center frequency calibration)
                try %
                    findex=find(cal.logspacedfreqs<=stimparam.prepulsefreq, 1, 'last');
                    atten=cal.atten(findex);
                    stimparam.prepulseamp=stimparam.prepulseamp-atten;
                    findex2=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                    atten=cal.atten(findex2);
                    stimparam.pulseamp=stimparam.pulseamp-atten;
                    djMessage( 'calibrated', 'append')
                catch
                    djMessage( 'NOT calibrated', 'append')
                end
            case 'soundfile'
                %best we can do for now is use the whitenoise calibration
                findex=find(cal.logspacedfreqs==-1); %freq of -1 indicates white noise
                atten=cal.atten(findex);
                stimparam.amplitude=stimparam.amplitude-atten;
                djMessage( 'calibrated', 'append')
            case 'silentsound'
                djMessage( 'calibration unnecessary', 'append')
            otherwise
                djMessage( 'NOT calibrated', 'append')
                
        end
    end
else
    djMessage( 'NOT calibrated', 'append')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitParams
global SP pref

SP.datapath=pref.datapath;
SP.mouseID='unknown';
SP.CurrentStimulus=[];
SP.stimcounter=0;
SP.NStimuli=[];
SP.Name={};
SP.Description={};
SP.NProtocols=0;
SP.NRepeats=0;
SP.PPALaseron=0;
SP.cal=[];
try
    cd(pref.root)
    cal=load('calibration.mat');
    SP.cal=cal;
    str=sprintf('successfully loaded calibration: %.0f - %.0f Hz, with %d freqs per octave', cal.minfreq, cal.maxfreq, cal.freqsperoctave);
    djMessage( str, 'append');
catch
    djMessage('failed to load calibration', 'append')
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
e=2; H=e;
w=100; h=25;
set(fig,'pos',pref.windowpos,'visible','on');

labelfs=8; %fontsize for labels

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
    'enable','inact','horiz','left','backgroundcolor', [.8 .8 .8],'pos',[e H 2*w 5*h]);
H=H+5*h+e;

%name window
SP.protocol_nameh=uicontrol(fig,'tag','protocol_name','style','text','units','pixels',...
    'enable','inact','horiz','left','backgroundcolor', [.8 .8 .8],'pos', [e H 2*w 3*h]);
H=H+3*h+e;

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

H=H+1*h+1*e;

%path display
SP.pathh=uicontrol(fig,'tag','pathdisplay','style','text','units','pixels','backgroundcolor', [.8 .8 .8],...
    'string', pref.datapath, 'enable','inact','horiz','left', 'pos',[e  H 3*w 2*h ]);
H=H+2*h+e;

%PPAactive message window
SP.PPAactive=uicontrol(fig,'tag','PPAactive','style','text','units','pixels',...
    'string', 'PPA messages', ...
    'enable','inact','horiz','left','backgroundcolor', [.8 .8 .8],'pos', [e H 3*w 1*h]);
H=H+h+e;

%Message window
SP.Messageh=uicontrol(fig,'tag','message','style','edit','fontweight','bold','units','pixels',...
    'enable','inact','horiz','left','Max', 8, 'pos',[e  H 3*w 5*h ]);
H=H+5*h+e;

%Run button
SP.Run=0;
SP.Runh=uicontrol(fig,'tag','Run','style','togglebutton','units','pixels','fontweight','bold',...
    'fontsize',12,'fontname','Arial','backgroundcolor',[0.75 0.75 0.75], ...
    'string', 'Play','callback', [me ';'],'enable','off','horiz','left','pos',[e H w 2*h ]);

%Record button
SP.Record=0;
SP.Recordh=uicontrol(fig,'tag','Record','style','togglebutton','units','pixels','fontweight','bold',...
    'fontsize',12,'fontname','Arial','backgroundcolor',[0 0.9 0],...
    'string', 'Record','callback', [me ';'],'enable','off','horiz','left','pos',[e+w H w 2*h ]);

%Stop Acquisition button - still buggy, not ready for release
SP.Stoph=uicontrol(fig,'tag','Stop','style','pushbutton','units','pixels',...
    'fontsize',10,'fontname','Arial', 'enable', 'off',...
    'string', 'Stop OE',  'callback', [me ';'],'horiz','left','pos',[2*e+2*w H w h ]);
H=H+h+e;

%launch open-ephys button - - still buggy, not ready for release
SP.LaunchOEh=uicontrol(fig,'tag','launchOE','style','pushbutton','units','pixels',...
    'fontsize',10,'fontname','Arial', 'enable', 'off',...
    'string', 'launch OE',  'callback', [me ';'],'horiz','left','pos',[2*e+2*w H w h ]);
H=H+1*h+e;

%reset zeroMQ button
SP.ResetZMQh=uicontrol(fig,'tag','ResetZMQ','style','pushbutton','units','pixels',...
    'fontname','Arial', ...
    'string', 'ResetZMQ','callback', [me ';'],'enable','on','horiz','left','pos',[3*e+2*w H w h ]);
% H=H+h+e;%reset zeroMQ button
SP.ResetZMQh=uicontrol(fig,'tag','ResetZMQ','style','pushbutton','units','pixels',...
    'fontname','Arial', ...
    'string', 'ResetZMQ','callback', [me ';'],'enable','on','horiz','left','pos',[3*e+2*w H w h ]);
H=H+h+e;

%test zeroMQ button
SP.TestZMQh=uicontrol(fig,'tag','TestZMQ','style','pushbutton','units','pixels',...
    'fontname','Arial', ...
    'string', 'TestZMQ','callback', [me ';'],'enable','on','horiz','left','pos',[3*e+2*w H w h ]);
% H=H+h+e;

%User menu
SP.user='lab'; %default user
userstr=pref.users;
userstr{end+1}='add new user';
SP.userh=uicontrol(fig,'tag','User','style','popupmenu','units','pixels','fontweight','bold',...
    'string', userstr,'enable','on','horiz','left','callback',[me ';'], 'pos',[e H w h ]);
H=H+h+e;
SP.userlabel=uicontrol(fig,'tag','userlabel','style','text','units','pixels',...
    'string', 'user', 'enable','inact','horiz','left','pos', [2*e H w h/2]);
H=H+h+e;

H=H+e;
%%%%%%%%%%%%%%%%%%%%%%%%%
% Laser controls
laserpanel = uipanel( 'Title','Laser','units','pixels', 'Position',[4*e+2*w 3*h+3*e .9*w 7.5*h+7*e]);
ww=.8*w;
%laser edit boxes
SP.LaserStart=0;
SP.LaserWidth=100;
SP.LaserNumPulses=10;
SP.LaserISI=100;
SP.LaserOnOff=0;

H=2*e;
SP.LaserISIh=uicontrol('Parent',laserpanel,'tag','LaserISI','style','edit','units','pixels',...
    'string', num2str(SP.LaserISI),'enable','off','horiz','left', 'callback',[me ';'],'pos',[e  H ww h ]);
H=H+h;
SP.LaserISIlabel=uicontrol('Parent', laserpanel,'tag','ISIlabel','style','text','units','pixels',...
    'string', 'ISI:', 'fontsize', labelfs,...
    'enable','on','horiz','left','pos', [e H ww h/2]);
H=H+h/2+e;

SP.LaserNumPulsesh=uicontrol('Parent',laserpanel,'tag','LaserNumPulses','style','edit','units','pixels',...
    'string', num2str(SP.LaserNumPulses),'horiz','left', 'callback',[me ';'],'pos',[e  H ww h ]);
H=H+h;
SP.LaserNumPulseslabel=uicontrol('Parent', laserpanel,'tag','NumPulseslabel','style','text','units','pixels',...
    'string', 'NumPulses:', 'fontsize', labelfs,...
    'enable','on','horiz','left','pos', [e H ww h/2]);
H=H+h/2+e;

SP.LaserWidthh=uicontrol('Parent',laserpanel,'tag','LaserWidth','style','edit','units','pixels',...
    'string', num2str(SP.LaserWidth),'horiz','left', 'callback',[me ';'],'pos',[e  H ww h ]);
H=H+h;
SP.LaserWidthlabel=uicontrol('Parent', laserpanel,'tag','Widthlabel','style','text','units','pixels',...
    'string', 'Width:', 'fontsize', labelfs,...
    'enable','on','horiz','left','pos', [e H ww h/2]);
H=H+h/2+e;

SP.LaserStarth=uicontrol('Parent',laserpanel,'tag','LaserStart','style','edit','units','pixels',...
    'string', num2str(SP.LaserStart),'horiz','left', 'callback',[me ';'],'pos',[e  H ww h ]);
H=H+h;
SP.LaserStartlabel=uicontrol('Parent', laserpanel,'tag','Startlabel','style','text','units','pixels',...
    'string', 'Start:', 'fontsize', labelfs,...
    'enable','on','horiz','left','pos', [e H ww h/2]);
H=H+h/2+e;

SP.LaserOnOffh=uicontrol('Parent',laserpanel,'tag','LaserOnOff','style','toggle','units','pixels',...
    'string', 'Laser is OFF','horiz','left', 'callback',[me ';'],'pos',[e  H ww h ]);
H=H+h;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%notebook features

H=e;
%Notes edit box
% SP.Notesh=uicontrol(fig,'tag','Notes','style','edit','fontweight','bold','units','pixels',...
%     'string', '','horiz','left', 'Max', Inf, 'callback',[me ';'],'pos',[2*e+2*w  H 2*w 3*h ]);
SP.Notesh=uicontrol(fig,'tag','Notes','style','edit','fontweight','bold','units','pixels',... %mw 07.20.2017
    'string', '','horiz','left', 'Max', 1000, 'callback',[me ';'],'pos',[2*e+2*w  H 2*w 3*h ]);
H=H+3*h+e;
SP.Noteslabel=uicontrol(fig,'tag','Noteslabel','style','text','units','pixels',...
    'string', 'Notes', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [3*e+3*w  H w h/2]);

%Notebook panel
H=H+h/2+e;
hp = uipanel( 'Title','Notebook','units','pixels', 'Position',[3*e+3*w H 2*w+2*e 16*h]);
SP.hp=hp;

% Manipulation/conditions details (anesthesia, Any other drugs, Shock, Reward, etc)
H=e;
%Drugs edit box
SP.Drugsh=uicontrol('Parent',hp,'tag','Drugs','style','edit','units','pixels',...
    'string', 'none','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h;
SP.Drugslabel=uicontrol('Parent',hp,'tag','Drugslabel','style','text','units','pixels',...
    'string', 'Drugs:', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e H w h/2]);
H=H+h/2+e;

%Reinforcement edit box (shock, reward, etc)
SP.Reinforcementh=uicontrol(fig,'Parent',hp,'tag','Reinforcement','style','edit','units','pixels',...
    'string', 'none','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h;
SP.Reinforcementlabel=uicontrol(fig,'Parent',hp,'tag','Reinforcementlabel','style','text','units','pixels',...
    'string', 'Reinforcement:', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);
H=H+h/2+e;

SP.Reinforcement='none';
SP.Drugs='none';

%mouse details
SP.mouseDOBh=uicontrol(fig,'Parent',hp,'tag','mouseDOB','style','edit','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h;
SP.mouseDOBlabel=uicontrol(fig,'Parent',hp,'tag','mouseDOBlabel','style','text','units','pixels',...
    'string', 'mouseDOB:', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);
H=H+h/2+e;

SP.mouseSexh=uicontrol(fig,'Parent',hp,'tag','mouseSex','style','edit','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h;
SP.mouseSexlabel=uicontrol(fig,'Parent',hp,'tag','mouseSexlabel','style','text','units','pixels',...
    'string', 'mouseSex:', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);
H=H+h/2+e;
SP.mouseGenotypeh=uicontrol(fig,'Parent',hp,'tag','mouseGenotype','style','edit','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h;
SP.mouseGenotypelabel=uicontrol(fig,'Parent',hp,'tag','mouseGenotypelabel','style','text','units','pixels',...
    'string', 'mouseGenotype:', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);
H=H+h/2+e;
SP.mouseSex='unknown';
SP.mouseDOB='unknown';
SP.mouseGenotype='unknown';
SP.Notes='';

%Depth edit box
SP.Depthh=uicontrol(fig,'Parent',hp,'tag','Depth','style','edit','fontweight','bold','units','pixels',...
    'string', '','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h+e;
SP.Depthlabel=uicontrol(fig,'Parent',hp,'tag','Depthlabel','style','text','units','pixels',...
    'string', 'Depth', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);
H=H+h/2+e;

%laser power edit box
SP.LaserPowerh=uicontrol(fig,'Parent',hp,'tag','LaserPower','style','edit','fontweight','bold','units','pixels',...
    'string', '','horiz','left', 'callback',[me ';'],'pos',[e  H w h ]);
H=H+h+e;
SP.LaserPowerlabel=uicontrol(fig,'Parent',hp,'tag','LaserPowerlabel','style','text','units','pixels',...
    'string', 'LaserPower', 'fontsize', labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);
H=H+h/2+e;
SP.Depth='unknown';
SP.LaserPower='unknown';

%mouseID menu
warning('off', 'MATLAB:hg:uicontrol:StringMustBeNonEmpty');
%old way was to have a list in djprefs
%new way is to just read the database mw 12.05.2018
%if isfield(pref, 'allmouseIDs') SP.allmouseIDs=pref.allmouseIDs; else SP.allmouseIDs='';end
try
    DB_IDs=whos('-file', 'mouseDB.mat');
catch
    DB_IDs=[];
end
SP.allmouseIDs{1}='';
for i=1:length(DB_IDs)
    idname=DB_IDs(i).name;
    idname=strrep(idname, 'mouseID_', '');
    SP.allmouseIDs{i+1}=idname;
end
SP.mouseIDMenuh=uicontrol(fig,'Parent',hp,'tag','mouseIDMenu','style','popupmenu','units','pixels','fontweight','bold',...
    'string', SP.allmouseIDs,'enable','on','horiz','left','callback',[me ';'], 'pos',[e H w h]);
H=H+h+e;

%mouseID edit box
SP.mouseIDh=uicontrol(fig,'Parent',hp,'tag','mouseID','style','edit','fontweight','bold','units','pixels',...
    'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e  H w h ], 'ForegroundColor', 'r');
H=H+h+e;
SP.mouseIDlabel=uicontrol(fig,'Parent',hp,'tag','mouseIDlabel','style','text','units','pixels',...
    'string', 'mouseID', 'fontsize', 2+labelfs,...
    'enable','inact','horiz','left','pos', [e  H w h/2]);

H=H+h/2;
SP.SecondmouseIDbuttonh=uicontrol(fig,'Parent',hp,'tag','SecondmouseIDbutton',...
    'style','togglebutton','fontweight','bold','units','pixels','value', 0, ...
    'string', '+mouse','horiz','left', 'callback',[me ';'],'pos',[e+w/2  H w/2+e h ]);

%H=H+h/2;
SP.ThirdmouseIDbuttonh=uicontrol(fig,'Parent',hp,'tag','ThirdmouseIDbutton',...
    'style','togglebutton','fontweight','bold','units','pixels','value', 0, ...
    'string', '+mouse','horiz','left', 'callback',[me ';'],'pos',[e+w+w/2  H w/2+e h ], 'visible', 'off');

%H=H+h/2;
SP.FourthmouseIDbuttonh=uicontrol(fig,'Parent',hp,'tag','FourthmouseIDbutton',...
    'style','togglebutton','fontweight','bold','units','pixels','value', 0, ...
    'string', '+mouse','horiz','left', 'callback',[me ';'],'pos',[e+2*w+w/2  H w/2+e h ], 'visible', 'off');

H=H+h;

%check for mouseDB and initialize if not present
cd(pref.root)
if exist('mouseDB.mat')~=2
    temp=[];
    save mouseDB.mat temp
end

%%%%%%%%%%%%%%%%%%
%send pi camera pulse button

H=H+5*h+e;
SP.camerapulse=uicontrol(fig,'tag','camerapulse','style','pushbutton','units','pixels',...
    'fontname','Arial', ...
    'string', 'camera','callback', [me ';'],'enable','on','horiz','left','pos',[4*e+3*w H w h ]);
H=H+h+e;

%  %Run button
% SP.Run=0;
% SP.Runh=uicontrol(fig,'tag','Run','style','togglebutton','units','pixels','fontweight','bold',...
%     'fontsize',12,'fontname','Arial','backgroundcolor',[0.75 0.75 0.75], ...
%     'string', 'Play','callback', [me ';'],'enable','off','horiz','left','pos',[e H w 2*h ]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
% mesoscope mode stuff

% a toggle for mescope mode that reads from prefs but allows you to override
if ~isfield(pref, 'mesoscope_mode') %if not in prefs, default to 0=OpenEphys
    pref.mesoscope_mode=0;
end
SP.mesoscope_mode=pref.mesoscope_mode; %1=meso, 0=intan/OE
SP.MesoModeh=uicontrol(fig,'tag','mesoscope_mode','style','togglebutton','units','pixels','fontweight','bold',...
    'string', '','enable','on','horiz','left', 'callback',[me ';'],'pos',[4*e+3*w H w h]);
set(SP.MesoModeh, 'value', SP.mesoscope_mode);
if SP.mesoscope_mode
    set(SP.MesoModeh, 'String', 'Mesoscope')
else
    set(SP.MesoModeh, 'String', 'OpenEphys');
end
if SP.mesoscope_mode
    
    fprintf('\ninitializing arduino...')
    djMessage('initializing arduino...', 'append')
    try
        clear SP.arduino
        SP.arduino=arduino;
        fprintf(' done\n' )
        djMessage(' done', 'append')
    catch
        warning(lasterr)
        djMessage('failed. Could not connect to arduino.', 'append')
        fprintf(' failed. Could not connect to arduino.\n' )
        SP.arduino=-1; %if there's no arduino connected
    end
    SP.arduinoProtocolStartpin='D10';
    SP.arduinoSCTpin='D11';
    SP.arduinoTimer=timer; %this stops the pulse after a delay without suspending execution
    %it's single shot by default
    SP.arduinoTimer.StartDelay=.001; %1 ms
    SP.arduinoTimer.TimerFcn=@(~,~)writeDigitalPin(SP.arduino, SP.arduinoSCTpin, 0);
end %if SP.mesoscope_mode


set(fig,'visible','on');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%InitializeGui%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function AddSecondMouse
global SP pref
fig=SP.fig;
hp=SP.hp;
pos=get(fig, 'position');
e=2; H=e;
w=100; h=25;
labelfs=8; %fontsize for labels

buttonstate=get( SP.SecondmouseIDbuttonh, 'value');
if buttonstate %add second mouse
    djMessage(' adding a second mouse...')
    %keyboard
    pos(3)=pos(3)+w+e;
    set(fig, 'pos', pos)
    set(SP.SecondmouseIDbuttonh,'string', '-mouse')
    set(SP.ThirdmouseIDbuttonh,'visible', 'on')
    pos=get(SP.Notesh, 'Position');
    set(SP.Notesh, 'Position',[pos(1) pos(2) pos(3)+w pos(4)]);
    
    %mouseID menu
    H=H+10.5*h+9*e;
    warning('off', 'MATLAB:hg:uicontrol:StringMustBeNonEmpty');
    %    if isfield(pref, 'allmouseIDs') SP.allmouseIDs=pref.allmouseIDs; else SP.allmouseIDs='';end
    SP.mouse2IDMenuh=uicontrol(fig,'Parent',hp,'tag','mouse2IDMenu','style','popupmenu','units','pixels','fontweight','bold',...
        'string', SP.allmouseIDs,'enable','on','horiz','left','callback',[me ';'], 'pos',[e+w H w h]);
    H=H+h+e;
    
    %mouseID edit box
    SP.mouse2IDh=uicontrol(fig,'Parent',hp,'tag','mouse2ID','style','edit','fontweight','bold','units','pixels',...
        'string', '','horiz','left', 'callback',[me ';'],'pos',[e+w  H w h ]);
    H=H+h+e;
    SP.mouse2IDlabel=uicontrol(fig,'Parent',hp,'tag','mouse2IDlabel','style','text','units','pixels',...
        'string', 'mouse2ID', 'fontsize', 2+labelfs,...
        'enable','inact','horiz','left','pos', [e+w  H w h/2]);
    
    H=81;
    %mouse details
    SP.mouse2DOBh=uicontrol(fig,'Parent',hp,'tag','mouse2DOB','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+w  H w h ]);
    H=H+h;
    SP.mouse2DOBlabel=uicontrol(fig,'Parent',hp,'tag','mouse2DOBlabel','style','text','units','pixels',...
        'string', 'mouse2DOB:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+w  H w h/2]);
    H=H+h/2+e;
    
    SP.mouse2Sexh=uicontrol(fig,'Parent',hp,'tag','mouse2Sex','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+w  H w h ]);
    H=H+h;
    SP.mouse2Sexlabel=uicontrol(fig,'Parent',hp,'tag','mouse2Sexlabel','style','text','units','pixels',...
        'string', 'mouse2Sex:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+w  H w h/2]);
    H=H+h/2+e;
    SP.mouse2Genotypeh=uicontrol(fig,'Parent',hp,'tag','mouse2Genotype','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+w  H w h ]);
    H=H+h;
    SP.mouse2Genotypelabel=uicontrol(fig,'Parent',hp,'tag','mouse2Genotypelabel','style','text','units','pixels',...
        'string', 'mouse2Genotype:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+w  H w h/2]);
    H=H+h/2+e;
    SP.mouse2Sex='unknown';
    SP.mouse2DOB='unknown';
    SP.mouse2Genotype='unknown';
    
elseif ~buttonstate %remove second mouse
    
    clear SP.mouse2Sex SP.mouse2DOB SP.mouse2Genotype SP.mouse2ID
    delete(SP.mouse2Genotypeh)
    delete(SP.mouse2Genotypelabel)
    delete(SP.mouse2Sexlabel)
    delete(SP.mouse2Sexh)
    delete(SP.mouse2DOBlabel)
    delete(SP.mouse2DOBh)
    delete(SP.mouse2IDMenuh)
    delete(SP.mouse2IDh )
    delete(SP.mouse2IDlabel)
    
    pos(3)=pos(3)-w-e;
    set(SP.fig, 'pos', pos)
    set(SP.SecondmouseIDbuttonh,'string', '+mouse')
    pos=get(SP.Notesh, 'Position');
    set(SP.Notesh, 'Position',[pos(1) pos(2) pos(3)-w pos(4)]);
    
    djMessage('removed second mouse')
end

function AddThirdMouse
global SP pref
fig=SP.fig;
hp=SP.hp;
pos=get(fig, 'position');
e=2; H=e;
w=100; h=25;
labelfs=8; %fontsize for labels

buttonstate=get( SP.ThirdmouseIDbuttonh, 'value');
if buttonstate %add Third mouse
    djMessage(' adding a Third mouse...')
    %keyboard
    pos(3)=pos(3)+w+e;
    set(fig, 'pos', pos)
    set(SP.ThirdmouseIDbuttonh,'string', '-mouse')
    set(SP.FourthmouseIDbuttonh,'visible', 'on')
    
    pos=get(SP.hp, 'Position');
    set(SP.hp, 'Position',[pos(1) pos(2) pos(3)+w pos(4)]);
    pos=get(SP.Notesh, 'Position');
    set(SP.Notesh, 'Position',[pos(1) pos(2) pos(3)+w pos(4)]);
    
    %mouseID menu
    H=H+10.5*h+9*e;
    W=2*w;
    warning('off', 'MATLAB:hg:uicontrol:StringMustBeNonEmpty');
    %if isfield(pref, 'allmouseIDs') SP.allmouseIDs=pref.allmouseIDs; else SP.allmouseIDs='';end
    SP.mouse3IDMenuh=uicontrol(fig,'Parent',hp,'tag','mouse3IDMenu','style','popupmenu','units','pixels','fontweight','bold',...
        'string', SP.allmouseIDs,'enable','on','horiz','left','callback',[me ';'], 'pos',[e+W H w h]);
    H=H+h+e;
    
    %mouseID edit box
    SP.mouse3IDh=uicontrol(fig,'Parent',hp,'tag','mouse3ID','style','edit','fontweight','bold','units','pixels',...
        'string', '','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h+e;
    SP.mouse3IDlabel=uicontrol(fig,'Parent',hp,'tag','mouse3IDlabel','style','text','units','pixels',...
        'string', 'mouse3ID', 'fontsize', 2+labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    
    H=81;
    %mouse details
    SP.mouse3DOBh=uicontrol(fig,'Parent',hp,'tag','mouse3DOB','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h;
    SP.mouse3DOBlabel=uicontrol(fig,'Parent',hp,'tag','mouse3DOBlabel','style','text','units','pixels',...
        'string', 'mouse3DOB:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    H=H+h/2+e;
    
    SP.mouse3Sexh=uicontrol(fig,'Parent',hp,'tag','mouse3Sex','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h;
    SP.mouse3Sexlabel=uicontrol(fig,'Parent',hp,'tag','mouse3Sexlabel','style','text','units','pixels',...
        'string', 'mouse3Sex:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    H=H+h/2+e;
    SP.mouse3Genotypeh=uicontrol(fig,'Parent',hp,'tag','mouse3Genotype','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h;
    SP.mouse3Genotypelabel=uicontrol(fig,'Parent',hp,'tag','mouse3Genotypelabel','style','text','units','pixels',...
        'string', 'mouse3Genotype:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    H=H+h/2+e;
    SP.mouse3Sex='unknown';
    SP.mouse3DOB='unknown';
    SP.mouse3Genotype='unknown';
    
elseif ~buttonstate %remove Third mouse
    
    clear SP.mouse3Sex SP.mouse3DOB SP.mouse3Genotype SP.mouse3ID
    delete(SP.mouse3Genotypeh)
    delete(SP.mouse3Genotypelabel)
    delete(SP.mouse3Sexlabel)
    delete(SP.mouse3Sexh)
    delete(SP.mouse3DOBlabel)
    delete(SP.mouse3DOBh)
    delete(SP.mouse3IDMenuh)
    delete(SP.mouse3IDh )
    delete(SP.mouse3IDlabel)
    
    pos(3)=pos(3)-w-e;
    set(SP.fig, 'pos', pos)
    set(SP.ThirdmouseIDbuttonh,'string', '+mouse')
    set(SP.FourthmouseIDbuttonh, 'visible', 'off')
    pos=get(SP.Notesh, 'Position');
    set(SP.Notesh, 'Position',[pos(1) pos(2) pos(3)-w pos(4)]);
    
    djMessage('removed Third mouse')
end

function AddFourthMouse
global SP pref
fig=SP.fig;
hp=SP.hp;
pos=get(fig, 'position');
e=2; H=e;
w=100; h=25;
labelfs=8; %fontsize for labels

buttonstate=get( SP.FourthmouseIDbuttonh, 'value');
if buttonstate %add fourth mouse
    djMessage(' adding a fourth mouse...')
    %keyboard
    pos(3)=pos(3)+w+e;
    set(fig, 'pos', pos)
    set(SP.FourthmouseIDbuttonh,'string', '-mouse')
    pos=get(SP.hp, 'Position');
    set(SP.hp, 'Position',[pos(1) pos(2) pos(3)+w pos(4)]);
    pos=get(SP.Notesh, 'Position');
    set(SP.Notesh, 'Position',[pos(1) pos(2) pos(3)+w pos(4)]);
    
    %mouseID menu
    H=H+10.5*h+9*e;
    W=3*w;
    warning('off', 'MATLAB:hg:uicontrol:StringMustBeNonEmpty');
    %if isfield(pref, 'allmouseIDs') SP.allmouseIDs=pref.allmouseIDs; else SP.allmouseIDs='';end
    SP.mouse4IDMenuh=uicontrol(fig,'Parent',hp,'tag','mouse4IDMenu','style','popupmenu','units','pixels','fontweight','bold',...
        'string', SP.allmouseIDs,'enable','on','horiz','left','callback',[me ';'], 'pos',[e+W H w h]);
    H=H+h+e;
    
    %mouseID edit box
    SP.mouse4IDh=uicontrol(fig,'Parent',hp,'tag','mouse4ID','style','edit','fontweight','bold','units','pixels',...
        'string', '','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h+e;
    SP.mouse4IDlabel=uicontrol(fig,'Parent',hp,'tag','mouse4IDlabel','style','text','units','pixels',...
        'string', 'mouse4ID', 'fontsize', 2+labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    
    H=81;
    %mouse details
    SP.mouse4DOBh=uicontrol(fig,'Parent',hp,'tag','mouse4DOB','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h;
    SP.mouse4DOBlabel=uicontrol(fig,'Parent',hp,'tag','mouse4DOBlabel','style','text','units','pixels',...
        'string', 'mouse4DOB:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    H=H+h/2+e;
    
    SP.mouse4Sexh=uicontrol(fig,'Parent',hp,'tag','mouse4Sex','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h;
    SP.mouse4Sexlabel=uicontrol(fig,'Parent',hp,'tag','mouse4Sexlabel','style','text','units','pixels',...
        'string', 'mouse4Sex:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    H=H+h/2+e;
    SP.mouse4Genotypeh=uicontrol(fig,'Parent',hp,'tag','mouse4Genotype','style','edit','units','pixels',...
        'string', 'unknown','horiz','left', 'callback',[me ';'],'pos',[e+W  H w h ]);
    H=H+h;
    SP.mouse4Genotypelabel=uicontrol(fig,'Parent',hp,'tag','mouse4Genotypelabel','style','text','units','pixels',...
        'string', 'mouse4Genotype:', 'fontsize', labelfs,...
        'enable','inact','horiz','left','pos', [e+W  H w h/2]);
    H=H+h/2+e;
    SP.mouse4Sex='unknown';
    SP.mouse4DOB='unknown';
    SP.mouse4Genotype='unknown';
    
elseif ~buttonstate %remove fourth mouse
    
    clear SP.mouse4Sex SP.mouse4DOB SP.mouse4Genotype SP.mouse4ID
    delete(SP.mouse4Genotypeh)
    delete(SP.mouse4Genotypelabel)
    delete(SP.mouse4Sexlabel)
    delete(SP.mouse4Sexh)
    delete(SP.mouse4DOBlabel)
    delete(SP.mouse4DOBh)
    delete(SP.mouse4IDMenuh)
    delete(SP.mouse4IDh )
    delete(SP.mouse4IDlabel)
    
    pos(3)=pos(3)-w-e;
    set(SP.fig, 'pos', pos)
    set(SP.FourthmouseIDbuttonh,'string', '+mouse')
    pos=get(SP.Notesh, 'Position');
    set(SP.Notesh, 'Position',[pos(1) pos(2) pos(3)-w pos(4)]);
    
    djMessage('removed fourth mouse')
end



