function [out] = PPAdj(varargin)

% djmaus module that initializes, loads, and plays sounds using PsychToolbox PortAudio (PPA) routines
%note: this module requires that PsychToolbox is installed. Freely available from psychtoolbox.org

global SP debugging
debugging=0; %print out helpful info for debugging

action = varargin{1};

restart_timer
if debugging
    fprintf('\naction: %s', action)
end

switch action
    case 'init'
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
            % Stop playback:
            PPAhandle=SP.PPAhandle;
            PsychPortAudio('Stop', PPAhandle);
            PsychPortAudio('Close', PPAhandle);
        catch
            djMessage('PPAdj: failed to close PPA device')
            pause(.5)
        end
        
        try
            %delete timer
            stop(SP.PPATimer);
            delete(SP.PPATimer);
            
            %no really, delete timer. I mean it.
            s=timerfind('TimerFcn', 'PPAdj(''PPATimer'');');
            if ~isempty(s)
                stop(s)
                delete(s)
            end
        catch
            djMessage('PPAdj: failed to stop PPA timer')
            pause(.5)
        end
        
    case 'load'
        if nargin<3
            return;
        end
        if nargin==4
            LoadPPA(varargin{2},varargin{3},varargin{4});
        else
            param.channel=1;
            LoadPPA(varargin{2},varargin{3},param); % first channel is the default channel
        end
        
        
        
        
    case 'PPATimer'
        try
            PPAhandle=SP.PPAhandle;
            status = PsychPortAudio('GetStatus', PPAhandle);
            
            
            h=SP.PPAactive;
            if status.Active==0; %device not running
                set(h, 'string', sprintf('PPA not running, XRuns=%g, CPUload=%.3f', status.XRuns, status.CPULoad), 'backgroundcolor', [.5 .5 .5])
                %                 if ~SP.Run
                %                     set(SP.Runh, 'backgroundcolor', [0 1 0], 'Value', 0)
                %                 end
            elseif status.Active==1; %device running
                protocol=SP.ProtocolIndex;
                current=SP.CurrentStimulus(protocol) + SP.NStimuli(protocol)*SP.NRepeats;
                inQueue=current-status.SchedulePosition-SP.CurrStimatPPAstart;
                set(h, 'string', sprintf('PPA running, inQueue=%d, XRuns=%g, CPUload=%.3f', inQueue, status.XRuns, status.CPULoad), 'backgroundcolor', [1 .5 .5])
                
                if status.XRuns>0
                    fprintf('\nXRun')
                    set(h,'backgroundcolor', [1 0 0])
                end
                

                %                 if ~SP.Run
                %                    set(SP.Runh, 'backgroundcolor', [1 .5 .5])
                %                 end
                
            end
        catch
            djMessage('ppatimer could not check status', 'append')
            h=SP.PPAactive;
            set(h, 'string', '', 'backgroundcolor', [.4 .4 .4])
            
        end
        
        
    case 'playsound'
        PlaySound;
        
    case 'stop'
        % Stop playback:
        PsychPortAudio('Stop',SP.PPAhandle,0);
        PsychPortAudio('Close', SP.PPAhandle);
        
    case 'camerapulse_on'
        SendCameraPulse_on
        
    case 'camerapulse_off'
        SendCameraPulse_off
        
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
        djMessage( 'InitPPA: failed to close device')
        pause(.2)
    end
end
% Verbosity:  0 = Shut up. 1 = Print errors, 2 = Print also warnings,
% 3 = Print also some info, 4 = Print more useful info (default), >5 = Be very
% verbose (mostly for debugging the driver itself).
PsychPortAudio('Verbosity', 0);

% Initialize driver, request low-latency preinit:
%InitializePsychSound(1); %  InitializePsychSound([reallyneedlowlatency=1])
%%% commented out line below Kip 2/21
%InitializePsychSound(0); %  InitializePsychSound([reallyneedlowlatency=0])

%note: this module requires that PsychToolbox is installed. Freely
%available from psychtoolbox.org
%You might need to reinstall PsychToolbox after a fresh matlab
%install/upgrade. In this case you don't need to re-download, just update.


%because it is machine dependent, we now set deviceid in djPrefs.m
%you can use printdevices.m to figure out which device id to use for your soundcard
%The default is to let GetASIOLynxDevice to take care of that for you
deviceid = pref.soundcarddeviceID; %32; %11;

numChan = pref.num_soundcard_outputchannels; %set in djPrefs.m
reqlatencyclass = pref.reqlatencyclass;
suggestedLatency=pref.suggestedLatency; %see post from Mario at https://beta.groups.yahoo.com/neo/groups/PSYCHTOOLBOX/conversations/topics/7131
%because it is machine dependent, we now set reqlatencyclass in djPrefs.m
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
buffSize = 1024;           %this seems to be ignored, you can try setting it in Lynx Mixer
% Low latency: 32, 64 or 128. High latency: 512>=
% nm 05.07.09 changed to 32, should fix dropouts.  If not, open LynxMixer.exe
% (in C:\lynx) and Settings->Buffer Size->32
% If Lynx seems not to change buffer size then type "CloseAllSoundDevices" into Matlab.
% You can monitor for dropouts using the LynxMixer as well.
% buffer size of 32 seems to cause dropouts on rig3, and 256 or 1024 do not, and
% the buffersize entered here seems to be ignored (needs to be set in lynx
% mixer and locked there)
buffPos = 0;

% Open audio device for low-latency output:

playbackonly=1;
try 
    fprintf('open deviceID# %d in PPAdj\n',deviceid)
    PPAhandle = PsychPortAudio('Open', deviceid, playbackonly, reqlatencyclass, SoundFs, numChan, buffSize, suggestedLatency);
catch
    error(sprintf('Could not open soundcard device id %d. Call PrintDevices and confirm that the soundcard DeviceIndex matches pref.soundcarddeviceID (in djPrefs)\n', deviceid));
end
%runMode = 0; %default, turns off soundcard after playback
runMode = 1; %leaves soundcard on (hot), uses more resources but may solve dropouts? mw 08.25.09: so far so good.
PsychPortAudio('RunMode', PPAhandle, runMode);


if isempty(PPAhandle)
    djMessage(me,'Can''t create PsychPortAudio object...');
    return;
end
SP.PPAhandle=PPAhandle; % hold the PsychPortAudio object
% fprintf('set SP.PPAhandle to %d in PPAdj line 214\n',SP.PPAhandle);
SP.numChan=numChan; %param to hold number of output channels with which we initialized card (num rows of samples must match this)
SP.SoundFs=SoundFs; %param to hold the sampling rate
SP.samples=[]; %param to hold the samples, used only for looping
SP.loop_flg=0; %param to store loop flag
SP.seamless=0; %param to store whether transition should be seamless or not
SP.buffers=[]; %param to store pointers to buffers for later deletion
djMessage( sprintf('Initialized PsychPortAudio with device %d, reqlatencyclass %d, Fs %d, numChan %d, buffersize %d suggestedLatency %g\n', deviceid, reqlatencyclass, SoundFs, numChan, buffSize, suggestedLatency), 'append');

%trying to workaround dropout on first sound after initialization by
%playing dummy tone here
% old way: call to LoadPPA followed by PlaySound
% % param=[]; %
% % LoadPPA('var',zeros(1,200),param)
% % PlaySound
%ne way: using low-level call instead, to avoid delivering a stray soundcard trigger
PsychPortAudio('FillBuffer', PPAhandle, zeros(numChan,200)); % fill buffer with brief silence
PlaySound
%commenting out above to see if it's still necessary mw 09.28.17
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SendCameraPulse_on
global SP pref
%sends a TTL pulse on soundcard channel 5, intended to start/stop a pi
%camera
SP.seamless=0;
CampulseOn=pref.CampulseOn;
PPAhandle=SP.PPAhandle; %grab PPAhandle object from param
SoundFs=SP.SoundFs; %sampling rate
numChan=SP.numChan; %number of output channels we initialized the soundcard with
cameratriglength=round(SoundFs/1000)*2; %2 ms trigger
samples=zeros(numChan, 2*cameratriglength);
camera_pulse=zeros(1, 2*cameratriglength);
camera_pulse(1:cameratriglength)=ones(size(1:cameratriglength))*.2;
samples(CampulseOn,:)=camera_pulse;

try
    PsychPortAudio('UseSchedule', PPAhandle, 0);   %mw 09.28.17
    PsychPortAudio('FillBuffer', PPAhandle, samples); % fill buffer now, start in PlaySound
    PlaySound
    djMessage('djPPA: sent camera ON pulse', 'append')
catch
    fprintf('camera button does not work during gap stimuli')
    djMessage('djPPA: camera button disabled during gap stimuli', 'append')
end

function SendCameraPulse_off
global SP pref
%sends a TTL pulse on soundcard channel 6, intended to start/stop a pi camera

SP.seamless=0;
CampulseOff=pref.CampulseOff;
PPAhandle=SP.PPAhandle; %grab PPAhandle object from param
SoundFs=SP.SoundFs; %sampling rate
numChan=SP.numChan; %number of output channels we initialized the soundcard with
cameratriglength=round(SoundFs/1000); %1 ms trigger
samples=zeros(numChan, 2*cameratriglength);
camera_pulse=zeros(1, 2*cameratriglength);
camera_pulse(1:cameratriglength)=ones(size(1:cameratriglength))*.2;
samples(CampulseOff,:)=camera_pulse;

try
    PsychPortAudio('UseSchedule', PPAhandle, 0);   %mw 09.28.17
    PsychPortAudio('FillBuffer', PPAhandle, samples); % fill buffer now, start in PlaySound
    PlaySound
    djMessage('djPPA: sent camera OFF pulse', 'append')
catch
    fprintf('camera button does not work during gap stimuli')
    djMessage('djPPA: camera button disabled during gap stimuli', 'append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function LoadPPA(type,where,param)
global SP pref debugging
% loads data. type can be either 'file' or 'var'
str='';
switch type
    case 'file'
        try
            load(where,'samples');
            str=[where ' loaded'];  % string to be displayed in the Message box
        catch
            djMessage(sprintf('Cannot load %s', where));
            return;
        end
    case 'var'
        samples=where;
        if debugging
            str=sprintf('vector loaded %s',datestr(now, 13)); % string to be displayed in the Message box
        end
    otherwise
        return;
end

Soundchannel1=pref.Soundchannel1;
Soundchannel2=pref.Soundchannel2; %only used for binaural
Soundcardtriggerchannel=pref.Soundcardtriggerchannel;
Laserchannel=pref.Laserchannel;
Shockchannel=pref.Shockchannel;
CampulseOff=pref.CampulseOff;
CampulseOn=pref.CampulseOn;


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
silence(Soundchannel1,:)=samples(1,:); %left channel
if nstimchans==2
    silence(Soundchannel2,:)=samples(2,:); %right channel
end
samples=silence;

%Soundcardtriggerchannel serves as trigger
trigsamples=zeros(1, length(samples));
triglength=round(SoundFs/1000); %1 ms trigger

%uncomment this line to make the soundcard trigger last as long as the sound
%for example, if you're using the soundcard trigger to drive an LED pulse
% triglength=length(samples);

%trigsamples(1:triglength)=ones(size(1:triglength));
%(this produces max output for trig amplitude, which is +10V on the Lynx,
%which is the way we used to do it for TTL triggers)

if GetAsioLynxDevice & isempty(GetXonarDevice)
    trigsamples(1:triglength)=.25*ones(size(1:triglength)); %for lynx soundcard
elseif GetXonarDevice & isempty(GetAsioLynxDevice)
    trigsamples(1:triglength)=ones(size(1:triglength));
elseif ismac
    trigsamples(1:triglength)=ones(size(1:triglength));
elseif GetRealtekDevice &  isempty(GetXonarDevice) &  isempty(GetAsioLynxDevice) &  isempty(GetPreSonusDevice)
    trigsamples(1:triglength)=.1*ones(size(1:triglength)); %for lynx soundcard
elseif GetFocusriteDevice & isempty(GetAsioLynxDevice)
    trigsamples(1:triglength)=ones(size(1:triglength));
elseif GetPreSonusDevice & isempty(GetAsioLynxDevice)
    trigsamples(1:triglength)=ones(size(1:triglength))/2;    
else
    error('cannot determine soundcard type')
end
%Lynx delivers +- 10V, so we cut down,
%Xonar delivers about +-1V, so we can use full range soundcardtrigger
%now, since djmaus is designed for open-ephys which wants 3.3V triggers, we
%use .33 = 3.3/10
%.25 is enough to trigger a digital TTL and is safer, so we will use that
%Reid says the digital lines can take 5V, but the ADCs cannot. Neither
%like negative voltages but can tolerate up to -0.5V before blowing
%anything. The soundcard puts some ringing on the edges but it's less than
%+/- 0.2V so I think we're totally safe for both digital lines and an ADC
%line if we want to use that as a monitor.

%trigsamples(end-triglength+1:end)=-.5*ones(size(1:triglength));
if loop_flg
    trigsamples=0*trigsamples;
end

%new way: specify channels in prefs. mw 08.24.2017
samples(Soundcardtriggerchannel,:)=trigsamples;

%old way
%  if pref.num_soundcard_outputchannels==2  %mw 091812
%     %only 2channel soundcard, no way you can do PPALaser
%     samples(numChan,:)=trigsamples;
% else
%     %>2 channels on card, we can leave a channel for PPALaser
%     samples(numChan-1,:)=trigsamples; %mw 081412
% end

if isfield(param, 'VarLaser') %will use variable Laser params set in protocol, so turn LaserOnOff on and disable djmaus GUI
    SP.LaserOnOff = 1;
    set(SP.LaserStarth, 'enable', 'off')
    set(SP.LaserNumPulsesh, 'enable', 'off')
    set(SP.LaserISIh, 'enable', 'off')
    set(SP.LaserWidthh, 'enable', 'off')
    set(SP.LaserOnOffh, 'enable', 'off', 'backgroundcolor',[.9 .5 .5],'foregroundcolor',[0 0 0]);
    set(SP.LaserOnOffh, 'string','Var Laser ON');
else
    set(SP.LaserStarth, 'enable', 'on')
    set(SP.LaserNumPulsesh, 'enable', 'on')
    set(SP.LaserISIh, 'enable', 'on')
    set(SP.LaserWidthh, 'enable', 'on')
    set(SP.LaserOnOffh, 'enable', 'on')
    djmaus('LaserOnOff')
    
end


%add in a Laser pulse if djmaus laser button is turned on AND if sound param.laser==1
laser_prepend=0;
laser_append=0;
on=SP.LaserOnOff;
if on
    if isfield(param, 'laser')
        if param.laser %then deliver a laser pulse
            djMessage('Laser pulse', 'append');
            
            if isfield(param, 'VarLaser') %use variable Laser params set in protocol, and ignore djmaus GUI
                if param.VarLaser %if Laser params say this stimulus has laser turned on
                    start=param.VarLaserstart; %ms
                    pulsewidth=param.VarLaserpulsewidth; %ms
                    numpulses=param.VarLasernumpulses; %ms
                    isi=param.VarLaserisi; %ms
                    
                    SP.LaserOnOff = 1;
                    set(SP.LaserStarth, 'enable', 'off', 'string', start)
                    set(SP.LaserNumPulsesh, 'enable', 'off', 'string',numpulses)
                    set(SP.LaserISIh, 'enable', 'off', 'string',isi)
                    set(SP.LaserWidthh, 'enable', 'off', 'string',pulsewidth)
                    set(SP.LaserOnOffh, 'enable', 'off', 'backgroundcolor',[.5 .5 .5],'foregroundcolor',[0 0 1]);
                    set(SP.LaserOnOffh, 'string','Laser is ON');
                else
                    %Laser params say this stimulus has laser turned OFF
                end
            else %VarLaser does not exist, so get laser pulse params from djmaus GUI
                set(SP.LaserStarth, 'enable', 'on')
                set(SP.LaserNumPulsesh, 'enable', 'on')
                set(SP.LaserISIh, 'enable', 'on')
                set(SP.LaserWidthh, 'enable', 'on')
                set(SP.LaserOnOffh, 'enable', 'on');
                
                if isfield(param, 'soaflag') %workaround to identify a GPIAS stimulus
                    start=SP.LaserStart+param.gapdelay; %ms relative to gap termination
                else
                    start=SP.LaserStart; %ms relative to sound onset
                end
                pulsewidth=SP.LaserWidth; %ms
                numpulses=SP.LaserNumPulses; %ms
                isi=SP.LaserISI; %ms
            end
            
            %insert laser pulse with specified parameters into samples
            %note that for flashtrains, I am changing isi to actually be soa (mw 1-28-2017)
            if length(pulsewidth)>1
                if length(pulsewidth)~=numpulses;
                    djMessage('djPPA: Laser numpulses must match number of widths', 'append')
                end
                if length(isi)~=numpulses-1;
                    djMessage('djPPA: Laser numpulses must be 1 more than number of isis', 'append')
                end
            end
            if length(pulsewidth)==1
                pulsewidth=repmat(pulsewidth, 1, numpulses);
            end
            if length(isi)==1
                isi=repmat(isi, 1, numpulses-1);
            end
            width=sum(pulsewidth)+sum(isi); %width of entire pulse train
            
            if start<0 %need to prepend some silence before sound
                temp=zeros(numChan,length(samples)+ -start*SoundFs/1000);
                temp(:,(-start*SoundFs/1000)+1:end)=samples;
                samples=temp;
                laser_prepend=(-start*SoundFs/1000)+1;
                
                if width*SoundFs/1000>length(samples)
                    %need to append some silence after sound
                    temp=zeros(numChan, width*SoundFs/1000);
                    temp(:,1:length(samples))=samples;
                    samples=temp;
                    laser_append=width*SoundFs/1000-length(samples);
                end
                laserpulse=zeros(1, length(samples));
                for n=1:numpulses
                    if  n==1
                        pstartsamp=1;
                    else
                        %pstartsamp=pstartsamp+(pulsewidth(n-1)+isi(n-1))*SoundFs/1000;
                        pstartsamp=pstopsamp+isi(n-1)*SoundFs/1000; %mw 1-28-2017
                    end
                    pstopsamp=pstartsamp+pulsewidth(n)*SoundFs/1000-1;
                    laserpulse(pstartsamp:pstopsamp)=1;
                end
                laserpulse(end)=0; %make sure to turn off pulse at end
                samples(Laserchannel,:)=laserpulse; %laser always goes on the last channel?
                
            else %start >=0
                if (start+width)*SoundFs/1000>length(samples)
                    %need to append some silence after sound
                    temp=zeros(numChan, (start+width)*SoundFs/1000);
                    temp(:,1:length(samples))=samples;
                    samples=temp;
                    laser_append=width*SoundFs/1000-length(samples);
                end
                laserpulse=zeros(1, length(samples));
                for n=1:numpulses
                    if  n==1
                        pstartsamp=start*SoundFs/1000+1;
                    else
                        %pstartsamp=pstartsamp+(pulsewidth(n-1)+isi(n-1))*SoundFs/1000;
                        pstartsamp=pstopsamp+isi(n-1)*SoundFs/1000; %mw 1-28-2017
                    end
                    pstartsamp=round(pstartsamp);
                    pstopsamp=pstartsamp+pulsewidth(n)*SoundFs/1000-1;
                    pstopsamp=round(pstopsamp);
                    laserpulse(pstartsamp:pstopsamp)=1;
                end
                laserpulse(end)=0; %make sure to turn off pulse at end
                samples(Laserchannel,:)=laserpulse;
            end
        end
    end
end


% we used to add a silent pad at the end to avoid dropouts, but that didn't really work and the problem is now solved differently anyway

%add in a Shock pulse if the stimulus includes a shock
if isfield(param, 'Shock') %if there is a shock field
    fprintf('\n shock stim')
    if param.Shock %if shock is on for this stim
        Shockstart=param.Shockstart; %ms
        Shockpulsewidth=param.Shockpulsewidth; %ms
        Shocknumpulses=param.Shocknumpulses; %ms
        Shockisi=param.Shockisi; %ms
        %insert laser pulse with specified parameters into samples
        %note that for flashtrains, I am changing isi to actually be soa (mw 1-28-2017)
        if length(Shockpulsewidth)>1
            if length(Shockpulsewidth)~=Shocknumpulses;
                djMessage('djPPA: Shock numpulses must match number of widths', 'append')
            end
            if length(Shockisi)~=Shocknumpulses-1;
                djMessage('djPPA: Shock numpulses must be 1 more than number of isis', 'append')
            end
        end
        if length(Shockpulsewidth)==1
            Shockpulsewidth=repmat(Shockpulsewidth, 1, Shocknumpulses);
        end
        if length(Shockisi)==1
            Shockisi=repmat(Shockisi, 1, Shocknumpulses-1);
        end
        Shockwidth=sum(Shockpulsewidth)+sum(Shockisi); %width of entire pulse train
        
        if Shockstart<0-laser_prepend %need to prepend some silence before sound
            temp=zeros(numChan,length(samples)+ -Shockstart*SoundFs/1000);
            temp(:,(-Shockstart*SoundFs/1000)+1:end)=samples;
            samples=temp;
            
            if Shockwidth*SoundFs/1000>length(samples)
                %need to append some silence after sound
                temp=zeros(numChan, Shockwidth*SoundFs/1000);
                temp(:,1:length(samples))=samples;
                samples=temp;
            end
            
            Shockpulse=zeros(1, length(samples));
            for n=1:Shocknumpulses
                if  n==1
                    pShockstartsamp=1;
                else
                    %pShockstartsamp=pShockstartsamp+(pulsewidth(n-1)+isi(n-1))*SoundFs/1000;
                    pShockstartsamp=pShockstartsamp+Shockisi(n-1)*SoundFs/1000; %mw 1-28-2017
                end
                pShockstopsamp=pShockstartsamp+Shockpulsewidth(n)*SoundFs/1000-1;
                Shockpulse(pShockstartsamp:pShockstopsamp)=1;
            end
            Shockpulse(end)=0; %make sure to turn off pulse at end
            samples(Shockchannel,:)=Shockpulse;
            
        else %Shockstart >=0
            if (Shockstart+Shockwidth)*SoundFs/1000>length(samples)
                %need to append some silence after sound
                temp=zeros(numChan, (Shockstart+Shockwidth)*SoundFs/1000);
                temp(:,1:length(samples))=samples;
                samples=temp;
            end
            
            Shockpulse=zeros(1, length(samples));
            for n=1:Shocknumpulses
                if  n==1
                    pShockstartsamp=Shockstart*SoundFs/1000+1;
                else
                    pShockstartsamp=pShockstartsamp+isi(n-1)*SoundFs/1000; %mw 1-28-2017
                end
                pShockstopsamp=pShockstartsamp+Shockpulsewidth(n)*SoundFs/1000-1;
                Shockpulse(pShockstartsamp:pShockstopsamp)=1;
            end
            Shockpulse(end)=0; %make sure to turn off pulse at end
            samples(Shockchannel,:)=Shockpulse;
        end
    end
end

% t=1:length(samples(1,:));
% t=1000*t/SoundFs;
% figure(100)
% plot(t, samples(1,:), 'o-')
% xlim([-1 1])
% figure(101)
% plot(t, samples(1,:), 'o-')
% xlim([ t(end)-1  t(end)+1])


SP.samples= samples; %store samples for re-buffering if we're looping (used only for looping)

if isfield(param, 'seamless')
    if param.seamless==1;
        seamless=param.seamless;
        status = PsychPortAudio('GetStatus', PPAhandle);
        if debugging
            str=sprintf('PositionSecs=%g\tdur:%g', status.PositionSecs, param.duration);
            fprintf('\nPositionSecs=%g\ndur:%g', status.PositionSecs, param.duration);
            fprintf('\nSchedulePosition=%g', status.PositionSecs);
            djMessage(str, 'append');
        end
        if status.Active==0; %device not running, need to start it
            str=sprintf('%shad to start it !!!');
            beep
            currstimuli=SP.CurrentStimulus;
            protocol=SP.ProtocolIndex;
            currstimulus=currstimuli(protocol);
            SP.CurrStimatPPAstart=currstimulus;
%             fid=fopen(fullfile(SP.datapath, 'seamless_starts.txt'), 'a');
%             fprintf(fid, '\n\nhad to start device');
%             %fprintf(fid, '\nexpstart %s', datestr(getparam('control', 'expstart')));
%             %fprintf(fid, '\nexptime %.2f', getparam('control', 'exptime'));
%             fprintf(fid, '\ncurrent stimulus %d', currstimulus);
%             fprintf(fid, '\nnrepeats %d', SP.NRepeats);
%             fclose(fid);
%             fid=fopen(fullfile(SP.datapath, 'seamless_restarts.txt'), 'a');
%             fprintf(fid, '%d\t', currstimulus);
%             fclose(fid);
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
            % %                 %djMessage(str, 'append');
            % %                 while status.PositionSecs<(gotime*param.duration/1000) %let's pause until halfway from start of stimulus
            % %                     status = PsychPortAudio('GetStatus', PPAhandle);
            % %                     pause(.01)
            % %                     str=sprintf('%s\nActive=%d,pausing at %g', str,status.Active, status.PositionSecs);
            % %                 djMessage(str, 'append');
            % %                 end
            % %             end
            
            %str=sprintf('%s\talready started', str);
            %  djMessage(str, 'append');
            buf = PsychPortAudio('CreateBuffer', [], samples);
            
            [success, freeslots] = PsychPortAudio('AddToSchedule', PPAhandle, buf, 1, 0.0, [], []);
            usedslots=128-freeslots;
            if debugging
                if success
                    str=sprintf('%s\tAddToSchedule success, %d used slots, %d free slots', str, usedslots, freeslots);
                else
                    str=sprintf('%s\tAddToSchedule FAIL', str);
                end
            end
            %freeslots defaults to 128
            if freeslots <20 %arbitrary guess
                %                         pause for to play back slots
                pause(5)
                str=sprintf('%s\tWarning: Free Slots <20 !!!', str);
                str=sprintf('%s\tPaused 5 seconds', str);
                
            end
            %troubleshooting an out-of-memory error with GPIAS 3-2010
            %PsychPortAudio('DeleteBuffer',buf, 1);%mw091710 this waits for
            %buffer to finish playing, which blocks seamless play
            %             result=PsychPortAudio('DeleteBuffer',buf, 0); %this doesn't delete it if it's still playing
            %                         str=sprintf('%s\nresult:%d', str, result);
            %                                     djMessage(str);
            
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
        SP.buffers= bufs; %store buffer pointers
    else %this stimulus is not seamless (there is a seamless flag which is set to 0)
    PsychPortAudio('FillBuffer', PPAhandle, samples); % fill buffer now, start in PlaySound
    currstimuli=SP.CurrentStimulus;
    protocol=SP.ProtocolIndex;
    currstimulus=currstimuli(protocol);
    SP.CurrStimatPPAstart=currstimulus;
    nreps=0; %0=repeat
    seamless=0;

    end
else %this stimulus is not seamless (there is no seamless flag)
    %commenting out on rig2 mw 02-10-2011
    %   PsychPortAudio('UseSchedule', PPAhandle, 0);      %Has no effect unless seamless stimuli were previously delivered, in which case we need to turn scheduling off.
    %above line causes complaint on rig2, works fine commented out.
    %but on rig3 it appears to need the line to transition back to
    %non-seamless stimuli, and causes no complaints. mw 09-11-09
    PsychPortAudio('FillBuffer', PPAhandle, samples); % fill buffer now, start in PlaySound
    currstimuli=SP.CurrentStimulus;
    protocol=SP.ProtocolIndex;
    currstimulus=currstimuli(protocol);
    SP.CurrStimatPPAstart=currstimulus;
    nreps=0; %0=repeat
    seamless=0;
    %this is a different approach which also works
    %  PsychPortAudio('FillBuffer', PPAhandle, 0*samples); %priming fill and start now, refill in PlaySound
    %  PsychPortAudio('Start', PPAhandle,nreps,0,0);
    
end
SP.loop_flg=loop_flg; %store loop flag
SP.seamless=seamless; %store whether transition should be seamless or not
try
    djMessage(str, 'append');
end
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
%     djMessage('failed to write logfile', 'append')
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function PlaySound
global SP
PPAhandle=SP.PPAhandle;
% samples=SP.samples; %get samples (simulus vector)
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
% djMessage('Initialized GUI');
% pause(.2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function restart_timer
PPATimer=timerfind('tag', 'PPATimer');
if strcmp(get(PPATimer, 'running'), 'off')
    start(PPATimer)
    djMessage('restarted timer', 'append')
else
    %    start(PPATimer)
    %    djMessage('timer already running', 'append')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out=me
% Simple function for getting the name of this m-file.
out=mfilename;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
