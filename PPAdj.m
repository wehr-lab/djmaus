function out=PPAdj(varargin)

% djmaus module that initializes, loads, and plays sounds using PsychToolbox PortAudio (PPA) routines
%note: this module requires that PsychToolbox is installed. Freely available from psychtoolbox.org
% you also have to separately install the asio driver, use the link on
% http://psychtoolbox.org/wikka.php?wakka=PsychPortAudio

global SP

action = varargin{1};

restart_timer

% Message(action)
switch action
    case 'init'
        %InitializeGUI;                  % show the gui = Message box:-)
        InitPPA;                        % Initialize soundmachine
        SP.PPATimer=timer('tag', 'PPATimer', 'TimerFcn',[me '(''PPATimer'');'],'ExecutionMode','fixedDelay', 'Period', .5);
        start(SP.PPATimer)
        
        
    case 'reset'
        % Stop playback:
        PPAhandle=SP.PPAhandle;
        PsychPortAudio('Stop', PPAhandle);
        PsychPortAudio('Close', PPAhandle);
        InitPPA;
        
    case 'close'
        try 
            stop(SP.PPATimer);
            delete(SP.PPATimer);
            PPAhandle=SP.PPAhandle;
            PsychPortAudio('Stop', PPAhandle);
            PsychPortAudio('Close', PPAhandle);
        end

    case 'load'
        if nargin<3
            return;
        end
        try
            if nargin==4
                LoadPPA(varargin{2},varargin{3},varargin{4});
            else
                param.channel=1;
                LoadPPA(varargin{2},varargin{3},param); % first channel is the default channel
            end
        catch
            Message('Cannot load sound', 'append');
        end
        
        
        
    case 'ppatimer'
        PPAhandle=SP.PPAhandle;
        try        status = PsychPortAudio('GetStatus', PPAhandle);
            
            %             status.Active
            h=SP.PPAactive;
            if status.Active==0; %device not running
                %                  Message('not running', 'append')
                set(h, 'string', 'PPA not running', 'backgroundcolor', [.5 .5 .5])
            elseif status.Active==1; %device running
                %                  Message(' running', 'append')
                set(h, 'string', 'PPA running', 'backgroundcolor', [1 0 0])
                
            end
        catch
            Message('ppatimer could not check status', 'append')
        end
        
        
    case 'playsound'
        PlaySound;
        
    case 'stop'
        % Stop playback:
        PsychPortAudio('Stop',SP.PPAhandle,0);
        
    case 'close'
        try
            % Stop playback:
            PsychPortAudio('Stop',SP.PPAhandle);
            PsychPortAudio('Close');
        catch
            Message('failed to close device')
            pause(.2)
        end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function InitPPA
global SP pref
if isfield(SP, 'PPAhandle')
    PPAhandle=SP.PPAhandle;
    try
        PsychPortAudio('Stop',PPAhandle );
        % Done, close driver and display:
        Priority(0);
        PsychPortAudio('Close');
    catch
        Message( 'InitPPA: failed to close device')
        pause(.2)
    end
end
% Initialize driver, request low-latency preinit:
%InitializePsychSound(1); %  InitializePsychSound([reallyneedlowlatency=0])
InitializePsychSound(0); %  InitializePsychSound([reallyneedlowlatency=0])

%note: this module requires that PsychToolbox is installed. Freely
%available from psychtoolbox.org
%You might need to reinstall PsychToolbox after a fresh matlab
%install/upgrade
%also, after a fresh PsychToolbox install, you need to download the
%ASIO-enabled PPA driver
%(http://psychtoolbox.org/wikka.php?wakka=PsychPortAudio) and copy the
%enclosed portaudio_x86.dll into C:\toolbox\Psychtoolbox
PsychPortAudio('Verbosity', 5); %nm 09.09.08 turn off all text feedback from PPA

%because it is machine dependent, we now set deviceid in Prefs.m
%use printdevices.m to figure out which device id to use for your soundcard
deviceid = pref.soundcarddeviceID; %32; %11;

numChan = pref.num_soundcard_outputchannels; %set in Prefs.m
reqlatencyclass = pref.reqlatencyclass;
%because it is machine dependent, we now set reqlatencyclass in Prefs.m
%on rig1 use 4; %on rig2, set to 1 (the default) to avoid dropouts mw 051809
%on rig1, 1 seems to cause dropouts but 2/3/4 seem better
% class 2 empirically the best, 3 & 4 == 2
% 'reqlatencyclass' Allows to select how aggressive PsychPortAudio should be about
% minimizing sound latency and getting good deterministic timing, i.e. how to
% trade off latency vs. system load and playing nicely with other sound
% applications on the system. Level 0 means: Don't care about latency, this mode

% works always and with all settings, plays nicely with other sound applications.
% Level 1 (the default) means: Try to get the lowest latency that is possible
% under the constraint of reliable playback, freedom of choice for all parameters
% and interoperability with other applications. Level 2 means: Take full control
% over the audio device, even if this causes other sound applications to fail or
% shutdown. Level 3 means: As level 2, but request the most aggressive settings
% for the given device. Level 4: Same as 3, but fail if device can't meet the
% strictest requirements.

SoundFs = pref.SoundFs;        % Must set this. 96khz, 48khz, 44.1khz.
buffSize = 1024;           % Low latency: 32, 64 or 128. High latency: 512>=
% nm 05.07.09 changed to 32, should fix dropouts.  If not, open LynxMixer.exe
% (in C:\lynx) and Settings->Buffer Size->32
% If Lynx seems not to change buffer size then type "CloseAllSoundDevices" into Matlab.
% You can monitor for dropouts using the LynxMixer as well.
% buffer size of 32 seems to cause dropouts on rig3, and 256 or 1024 do not, and
% the buffersize entered here seems to be ignored (needs to be set in lynx
% mixer and locked there)
buffPos = 0;

% Open audio device for low-latency output:


try PPAhandle = PsychPortAudio('Open', deviceid, [], reqlatencyclass, SoundFs, numChan, buffSize);
catch
    error(sprintf('Could not open soundcard device id %d. Call PrintDevices and confirm that the soundcard DeviceIndex matches pref.soundcarddeviceID (in Prefs)\n', deviceid));
end
%runMode = 0; %default, turns off soundcard after playback
runMode = 1; %leaves soundcard on (hot), uses more resources but may solve dropouts? mw 08.25.09: so far so good.
PsychPortAudio('RunMode', PPAhandle, runMode);

% Unknown system: Assume zero bias. User can override with measured
% values:
%latbias = 0.0013;

% Tell driver about hardwares inherent latency, determined via calibration
% once:
%prelat = PsychPortAudio('LatencyBias', PPAhandle, latbias);
%postlat = PsychPortAudio('LatencyBias', PPAhandle);

if isempty(PPAhandle)
    Message(me,'Can''t create PsychPortAudio object...');
    return;
end
SP.PPAhandle=PPAhandle; % hold the PsychPortAudio object
SP.numChan=numChan; %param to hold number of output channels with which we initialized card (num rows of samples must match this)
SP.SoundFs=SoundFs; %param to hold the sampling rate
SP.Samples=[]; %param to hold the samples, used only for looping
SP.loop_flg=0; %param to store loop flag
SP.seamless=0; %param to store whether transition should be seamless or not
SP.buffers=[]; %param to store pointers to buffers for later deletion
Message( sprintf('Initialized PsychPortAudio with device %d, reqlatencyclass %d, Fs %d, numChan %d, buffersize %d\n\n', deviceid, reqlatencyclass, SoundFs, numChan, buffSize), 'append');

%trying to workaround dropout on first sound after initialization by
%playing dummy tone here
%
param=[]; %do I need to put stuff in param?
LoadPPA('var',zeros(1,200),param)
PlaySound
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadPPA(type,where,param)
global SP pref %mw 081412
% loads data to soundmachine. type can be either 'file' or 'var'
switch type
    case 'file'
        try
            load(where,'samples');
            str=[where ' loaded'];  % string to be displayed in the Message box
        catch
            Message(sprintf('Cannot load %s', where));
            return;
        end
    case 'var'
        samples=where;
        str=sprintf('\nvector loaded %s',datestr(now, 13)); % string to be displayed in the Message box
    otherwise
        return;
end
%Message(str, 'append')

% if isfield(param,'channel')
%     channel=param.channel(1);
% else
%     channel=1;  % default channel
% end

if isfield(param, 'loop_flg')
    loop_flg=param.loop_flg;
else
    loop_flg=0;
end

PPAhandle=SP.PPAhandle; %grab PPAhandle object from param
SoundFs=SP.SoundFs; %sampling rate
numChan=SP.numChan; %number of output channels we initialized the soundcard with
nstimchans=min(size(samples)); %number of channels of requested stimulus (i.e. mono or stereo)
samples=reshape(samples, nstimchans, length(samples)); %ensure samples are a row vector
silence=zeros(numChan, length(samples));
silence(1:nstimchans,:)=samples;
samples=silence;

%last channel serves as trigger
trigsamples=zeros(1, length(samples));
triglength=round(SoundFs/1000); %1 ms trigger

%uncomment this line to make the soundcard trigger last as long as the sound
%for example, if you're using the soundcard trigger to drive an LED pulse
% triglength=length(samples);


trigsamples(1:triglength)=ones(size(1:triglength));
%trigsamples(end-triglength+1:end)=-.5*ones(size(1:triglength));
if loop_flg
    trigsamples=0*trigsamples;
end

if pref.num_soundcard_outputchannels==2  %mw 091812
    %only 2channel soundcard, no way you can do PPALaser
    samples(numChan,:)=trigsamples;
else
    %>2 channels on card, we can leave a channel for PPALaser
    samples(numChan-1,:)=trigsamples; %mw 081412
end
%add in a LaserPPA pulse if requested by that module AND if sound params
%laserON==1
on=SP.PPALaseron;
if on
    if isfield(param, 'AOPulseOn')
        if param.AOPulseOn %then deliver a laser pulse
            str='Laser pulse *';
            start=SP.PPALaserstart; %ms
            pulsewidth=SP.PPALaserwidth; %ms
            numpulses=SP.PPALasernumpulses; %ms
            isi=SP.PPALaserisi; %ms
            if length(pulsewidth)>1
                if length(pulsewidth)~=numpulses;
                    Message('PPALaser: numpulses must match number of widths', 'append')
                end
                if length(isi)~=numpulses-1;
                    Message('PPALaser: numpulses must be 1 more than number of isis', 'append')
                end
            end
            if length(pulsewidth)==1
                pulsewidth=repmat(pulsewidth, 1, numpulses);
            end
            if length(isi)==1
                isi=repmat(isi, 1, numpulses-1);
            end
            % width=pulsewidth*numpulses+isi*(numpulses-1); %width of entire pulse train
            %                 elseif length(pulsewidth)>1
            %                     if length(isi)==1
            %                         isi=repmat(isi, 1, numpulses-1);
            %                     end
            %                     %   width=sum(pulsewidth)+sum(isi); %width of entire pulse train
            %                 end
            width=sum(pulsewidth)+sum(isi); %width of entire pulse train
            
            if start<0 %need to prepend some silence before sound
                temp=zeros(numChan,length(samples)+ -start*SoundFs/1000);
                temp(:,(-start*SoundFs/1000)+1:end)=samples;
                samples=temp;
                
                if width*SoundFs/1000>length(samples)
                    %need to append some silence after sound
                    temp=zeros(numChan, width*SoundFs/1000);
                    temp(:,1:length(samples))=samples;
                    samples=temp;
                end
                laserpulse=zeros(1, length(samples));
                for n=1:numpulses
                    if  n==1
                        pstartsamp=1;
                    else
                        pstartsamp=pstartsamp+(pulsewidth(n-1)+isi(n-1))*SoundFs/1000;
                    end
                    pstopsamp=pstartsamp+pulsewidth(n)*SoundFs/1000-1;
                    laserpulse(pstartsamp:pstopsamp)=1;
                end
                laserpulse(end)=0; %make sure to turn off pulse at end
                samples(numChan,:)=laserpulse;
                
            else %start >=0
                if (start+width)*SoundFs/1000>length(samples)
                    %need to append some silence after sound
                    temp=zeros(numChan, (start+width)*SoundFs/1000);
                    temp(:,1:length(samples))=samples;
                    samples=temp;
                end
                laserpulse=zeros(1, length(samples));
                for n=1:numpulses
                    if  n==1
                        pstartsamp=start*SoundFs/1000+1;
                    else
                        pstartsamp=pstartsamp+(pulsewidth(n-1)+isi(n-1))*SoundFs/1000;
                    end
                    pstopsamp=pstartsamp+pulsewidth(n)*SoundFs/1000-1;
                    laserpulse(pstartsamp:pstopsamp)=1;
                end
                %                      for n=1:numpulses
                %                          pstartsamp=start*SoundFs/1000+1+(n-1)*(pulsewidth+isi)*SoundFs/1000;
                %                          pstopsamp=pstartsamp+pulsewidth*SoundFs/1000-1;
                %                          laserpulse(pstartsamp:pstopsamp)=1;
                %
                %                      end
                %                     laserpulse(start*SoundFs/1000+1:(start+width)*SoundFs/1000)=1;
                laserpulse(end)=0; %make sure to turn off pulse at end
                samples(numChan,:)=laserpulse;
            end
            
            Message(str, 'append')
            
            if isfield(param, 'PulseTrace') % AKH 6/29/14
                samples(numChan,:)=param.PulseTrace;
            end
            
        end
    end
end

str1=str;
% we used to add a silent pad at the end to avoid dropouts, but that didn't really work and the problem is now solved differently anyway

SP.samples= samples; %store samples for re-buffering if we're looping (used only for looping)

if isfield(param, 'seamless')
    if param.seamless==1
        seamless=param.seamless;
        status = PsychPortAudio('GetStatus', PPAhandle);
        str=sprintf('PositionSecs=%g\ndur:%g', status.PositionSecs, param.duration);
        %         Message(str, 'append');
        
        if status.Active==0; %device not running, need to start it
            str=sprintf('%s\nhad to start it !!!', str);
            beep
            currstimuli=SP.CurrentStimulus;
            protocol=SP.ProtocolIndex;
            currstimulus=currstimuli(protocol);
            fid=fopen(fullfile(SP.datapath, 'seamless_starts.txt'), 'a');
            fprintf(fid, '\n\nhad to start device');
            %fprintf(fid, '\nexpstart %s', datestr(getparam('control', 'expstart')));
            %fprintf(fid, '\nexptime %.2f', getparam('control', 'exptime'));
            fprintf(fid, '\ncurrent stimulus %d', currstimulus);
            fprintf(fid, '\nnrepeats %d', SP.NRepeats);
            fclose(fid);
            fid=fopen(fullfile(SP.datapath, 'seamless_restarts.txt'), 'a');
            fprintf(fid, '%d\t', currstimulus);
            fclose(fid);
            PsychPortAudio('UseSchedule', PPAhandle, 1);
            buf = PsychPortAudio('CreateBuffer', [], samples);
            PsychPortAudio('AddToSchedule', PPAhandle, buf, 1, 0.0, [], []);
            nreps=1;
            %           when=GetSecs+.01; %mw 032410
            when=GetSecs+.1;
            %when=GetSecs+2;
            PsychPortAudio('Start', PPAhandle,nreps,when,0);
            
        else %already started, just add to schedule
            %mw 07.03.2014 trying new strategy to try to solve "screwups"
            %instead of gotime/pause, I will check for free slots at
            %AddToSchedule and pause there if necessary
            
            % %             %gotime is in normalized units (fraction of stimulus duration)
            % % %             gotime=.5; %this was pre-06.12.2014
            % %             gotime=.01; %mw 06.13.2014
            % % %            gotime=.1; %mw 06.13.2014 last known goo dvlaue 07.03.2014
            % %             if status.PositionSecs<gotime*param.duration/1000 %less than halfway from start of stimulus (aldis you could try lowering this to, say, .25)
            % %                 str=sprintf('%s\npausing', str);
            % %                 %Message(str, 'append');
            % %                 while status.PositionSecs<(gotime*param.duration/1000) %let's pause until halfway from start of stimulus
            % %                     status = PsychPortAudio('GetStatus', PPAhandle);
            % %                     pause(.01)
            % %                     str=sprintf('%s\nActive=%d,pausing at %g', str,status.Active, status.PositionSecs);
            % %                 Message(str, 'append');
            % %                 end
            % %             end
            
            str=sprintf('%s\nalready started', str);
            %  Message(str, 'append');
            buf = PsychPortAudio('CreateBuffer', [], samples);
            
            [success, freeslots] = PsychPortAudio('AddToSchedule', PPAhandle, buf, 1, 0.0, [], []);
            if success
                str=sprintf('%s\nAddToSchedule success, %d free slots', str, freeslots);
            else
                str=sprintf('%s\nAddToSchedule FAIL', str);
            end
            %freeslots defaults to 128
            if freeslots <20 %arbitrary guess
                %                         pause for to play back slots
                pause(5)
                str=sprintf('%s\nWarning: Free Slots <20 !!!', str);
                str=sprintf('%s\nPaused 5 seconds', str);
                
            end
            %troubleshooting an out-of-memory error with GPIAS 3-2010
            %PsychPortAudio('DeleteBuffer',buf, 1);%mw091710 this waits for
            %buffer to finish playing, which blocks seamless play
            %             result=PsychPortAudio('DeleteBuffer',buf, 0); %this doesn't delete it if it's still playing
            %                         str=sprintf('%s\nresult:%d', str, result);
            %                                     Message(str);
            
        end
        %delete any unused buffers still hanging around
        oldbufs=SP.buffers;
        bufs=[oldbufs buf];
        SP.buffers=bufs; %store buffer pointers
        for b=bufs
            result=PsychPortAudio('DeleteBuffer',b, 0);
            if result
                bufs=setdiff(bufs, b);
            end
        end
        SetParam(me,'buffers', bufs); %store buffer pointers
    end
else %this stimulus is not seamless
    %commenting out on rig2 mw 02-10-2011
    %   PsychPortAudio('UseSchedule', PPAhandle, 0);      %Has no effect unless seamless stimuli were previously delivered, in which case we need to turn scheduling off.
    %above line causes complaint on rig2, works fine commented out.
    %but on rig3 it appears to need the line to transition back to
    %non-seamless stimuli, and causes no complaints. mw 09-11-09
    PsychPortAudio('FillBuffer', PPAhandle, samples); % fill buffer now, start in PlaySound
    
    nreps=0; %0=repeat
    seamless=0;
    %this is a different approach which also works
    %  PsychPortAudio('FillBuffer', PPAhandle, 0*samples); %priming fill and start now, refill in PlaySound
    %  PsychPortAudio('Start', PPAhandle,nreps,0,0);
    
end
SP.loop_flg=loop_flg; %store loop flag
SP.seamless=seamless; %store whether transition should be seamless or not
Message(str, 'append');

% write diagnostic logfile -mw 07.07.2014
%??? 9-3-2015 mw
% try
%     protocol=GetParam('stimulusprotocol', 'protocol');
%     currstims=GetParam('stimulusprotocol', 'currentstimulus');
%     currstim=currstims(protocol);
%     stimuli=GetSharedParam('StimulusProtocols');
%     stimuli=stimuli{protocol};
%     stimulus=stimuli(currstim);
%     paths=Control('GetDataPath');
%     cd(paths{3});
%     expids=Control('GetExpid');
%     fid=fopen(sprintf('%sppasoundlog.txt', expids{3}), 'at');
%     fprintf(fid, '\n\nstim %d', currstim);
%     fprintf(fid, '\n%s', stimulus.type);
%     names=fieldnames(param);
% %     fprintf(fid, '\n');
%     for nameidx=1:length(names)
%         fprintf(fid, '\n%s %g', names{nameidx}, eval(['param.', names{nameidx}]));
%     end
%     fprintf(fid, '%s', str1);
%     fprintf(fid, '\n%s', str);
%     fprintf(fid,'\nGetSecs %.6f', GetSecs);
%     fclose(fid);
% catch
%     Message('failed to write logfile', 'append')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlaySound
global SP
PPAhandle=SP.PPAhandle;
samples=SP.samples; %get samples (simulus vector)
seamless=SP.seamless; %whether transition should be seamless or not
loop_flg=SP.loop_flg; %get loop flag

if seamless
    %    do nothing, since sound was added to schedule in LoadPPA
else %not seamless
    %start device here, it was filled in PPALoad
    if loop_flg
        PsychPortAudio('Start', PPAhandle,0,0,0);
        PsychPortAudio('RescheduleStart', PPAhandle,0,0,0);
    else
        nreps=1;
        %when=GetSecs+.1; %this extra latency prevents dropouts somehow
        when=0; %use this to start immediately
        waitForStart=0;
        
        PsychPortAudio('Start', PPAhandle,nreps,when,waitForStart);
        %PsychPortAudio('RefillBuffer', PPAhandle, 0, samples, 0); %this doesn't work
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function InitializeGUI
% fig = ModuleFigure(me);
% set(fig,'doublebuffer','on','visible','off');
%
% hs = 120;
% h = 5;
% vs = 50;
% n = 1;
% % Message box
% uicontrol('parent',fig,'tag','message','style','text',...
%      'BackgroundColor', [.9 .9 .9], ...
%     'enable','inact','horiz','left','units','normal','pos',[.02 .35 .9 .55]); n=n+1;
% screensize=get(0,'screensize');
%
% uicontrol('parent',fig,'string','Reset','tag','reset','units','normal',...
%     'position',[0.02 0.02 0.40 0.20],'enable','on','foregroundcolor',[0.9 0 0],...
%     'fontweight','bold',...
%     'style','pushbutton','callback',[me ';']);
%
% uicontrol('parent',fig,'string','not playing','tag','active','units','normal',...
%     'position',[0.02 0.24 0.80 0.10],'foregroundcolor',[0.4 .4 .4],...
%     'fontweight','bold',...
%     'style','text');
%
% delete(timerfind('tag', 'PPATimer'))
%
% set(fig,'pos', [screensize(3)-128 screensize(4)-n*vs-100 158 150] ,'visible','on');
%
% Message('Initialized GUI');
% pause(.2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function restart_timer
PPATimer=timerfind('tag', 'PPATimer');
if strcmp(get(PPATimer, 'running'), 'off')
    start(PPATimer)
    Message('restarted timer', 'append')
else
    %    start(PPATimer)
    %    Message('timer already running', 'append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=me
% Simple function for getting the name of this m-file.
out=mfilename;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
