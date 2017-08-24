function MakeGapShockdjProtocol
% usage: MakeGapShockdjProtocol
%Create Gap Shock Protocol and incorporate a previously
% created GPIAS djmaus protocol. Laser and Shock are now independent
% channels. last updated by Mike 8-24-2017
%
%This is a special-purpose protocol for paired gap-shock conditioning.
%First use MakeGPIASdjProtocol to create a GPIAS protocol. Use IL=0.
% Then call this function and select that protocol to add a laser and shock condition to every GPIAS trial.
% Use the shock soundcard channel output to drive the S88 stimulator to
% deliver shock.
% The soundcard channel used for shock is specified in djPrefs - you might
% need to edit your djPrefs the first time you try this on a new computer,
% and add the following line:
% pref.Shockchannel=5; %or whatever channel you have available
%
%New: in principle you can include whatever type of laser conditions you
%want. For this version we hard-code every trial to have both a shock and a
%laser. Slight modification of this function should allow any combination.
%
%example:
% noiseamp=80; gapdurs=16; gapdelay=1000; poststartle=0;
% pulsedur=0;pulseamps=0;soa=50; soaflag='isi'; ramp=0; isi=30000; isi_var=.33; IL=0; nreps=20;
% MakeGPIASdjProtocol(noiseamp, gapdurs, gapdelay, poststartle, pulsedur, pulseamps, soa, soaflag, ramp, isi, isi_var, IL, nreps)
% MakeGapShockdjProtocol


global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
if exist('GPIAS Protocols', 'dir') ~=7
    error('GPIAS Protocols directory not found. You need to create a GPIAS protocol first')
end
cd ('GPIAS Protocols')

[tcfilename, tcpathname] = uigetfile('*.mat', 'Choose GPIAS protocol to incorporate into gap-shock protocol:');
if isequal(tcfilename,0) || isequal(tcpathname,0)
    disp('User pressed cancel')
    return
else
    disp(['User selected ', fullfile(tcpathname, tcfilename)])
end
tc=load(fullfile(tcpathname, tcfilename));
TotalNumStim=length(tc.stimuli);
protocol_description=tc.stimuli(1).protocol_description;
protocol_name=tc.stimuli(1).protocol_name;



tc_n=0; %TC tone index
st_n=0; %output stimuli index
%note that here we keep tc_n == st_n, so they are totally redundant. But
%for future modifications it may be useful, for example if you want to
%change order or insert additional stimuli
while tc_n+1<=length(tc.stimuli)
    
    
    
    %insert embedded tones
    start_tc_n=tc_n; %store starting tc_n to do pulse-off repeat of tones
    edur=0;
    tc_n=tc_n+1;
    st_n=st_n+1;
    tone=tc.stimuli(tc_n);
    while ~strcmp(tone.type, 'GPIAS') %cycle through  noise stimuli with laser off
        stimuli(st_n)=tone;
        stimuli(st_n).param.laser=0;
        edur=edur+tone.param.duration+tone.param.next;
        tc_n=tc_n+1;
        st_n=st_n+1;
        if tc_n>length(tc.stimuli)
            break
        end
        tone=tc.stimuli(tc_n);
    end
    
    if tc_n<=length(tc.stimuli)
        %then add the GPIAS with laser on
        stimuli(st_n)=tone;
        stimuli(st_n).param.laser=1;
        
        %and the shock
        %Weible 2014 used a 50 ms train of 1 ms pulses at 500 Hz
        stimuli(st_n).param.Shock=1; %whether to deliver a shock on this stimulus
        stimuli(st_n).param.Shockstart=stimuli(st_n).param.gapdelay; %ms relative to sound onset
        stimuli(st_n).param.Shockpulsewidth=1; %ms
        stimuli(st_n).param.Shocknumpulses=1; % for shock trains
        stimuli(st_n).param.Shockisi=nan; %ms, for shock trains
        %note: the soundcard channel to deliver shock TTL is specified in djPrefs as pref.Shockchannel
        stimuli(st_n).stimulus_description=GetParamStr(stimuli(st_n));
    end
end

%update protocol_description and protocol_name in stimuli structure

new_protocol_name=sprintf('GapShockProtocol from %s', protocol_name);
new_file_name=sprintf('GapShockProtocol_from_%s', protocol_name);
new_protocol_description=sprintf('GapShockProtocol from %s', protocol_description);
for n=1:TotalNumStim
    stimuli(n).protocol_description=new_protocol_description;
    stimuli(n).protocol_name=new_protocol_name;
end

cd(pref.stimuli) %where stimulus protocols are saved
cd('GPIAS Protocols')
save(new_file_name, 'stimuli')
fprintf('\nwrote file %s \nin directory %s', new_file_name, pwd)
fprintf('\n')





