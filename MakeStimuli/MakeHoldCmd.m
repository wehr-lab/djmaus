function waveform=MakeHoldCmd(varargin)

global exper pref

% Creates a voltage holdcmd waveform to be sent out via AO
% for external command of holding potential
% first input: (fields of structure sent as first input)
%   start   -  delay to the start of the ramp after the trigger (ms)
%   ramp   -   ramp duration (ms)
%   holdduration-    duration of the holding command after ramp (ms)
%   holdcmd_from  -   the previous commanded holding potential (mV)
%   holdcmd_to  -   the target commanded holding potential (mV)
%  (i.e. the ramp will go from holdcmd_from to holdcmd_to)
% second input:
%   samplerate- sampling rate (Hz)
% Output
%   waveform   -   the holdcmd waveform
%

waveform=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
start=params.start;
ramp=params.ramp;
holdcmd_from=params.holdcmd_from;
holdcmd_to=params.holdcmd_to;
holdduration=params.holdduration;
if isfield(params,'channel')
    channel=params.channel;
else
    channel=1;                  % default channel
end

fprintf('\nMakeHoldCmd: %d to %d', holdcmd_from, holdcmd_to)

waveform_length=(start+ramp+holdduration); % in ms
%include extra 1000 for series pulses in hold_duration when calling MakeHoldCmd
samplerate=samplerate/1000;

%get mode
mode=PatchPreProcess('GetMode');
switch mode{channel}
    case {'Track','V-Clamp'}
        factor=20;
    case {'I=0','I-Clamp Normal','I-Clamp Fast'}
        factor=2000;
    otherwise
        factor=1;
end

%series pulsetrain

PTparams.start= params.pulse_start;
PTparams.width= params.pulse_width;
PTparams.height=params.pulse_height;
PTparams.npulses=params.npulses;
PTparams.isi= params.pulse_isi;
PTparams.duration= params.pulseduration;

pulsetrain=MakePulse(PTparams, 1000*samplerate)*factor; %uncorrect factor
pulsetrain=pulsetrain+holdcmd_to;

%prepare the samples
waveform=zeros(waveform_length*samplerate,1);

sampled_ramp=ramp*samplerate;
sampled_start=start*samplerate;
sampled_holdduration=holdduration*samplerate;
sampled_start=max(1,sampled_start); % if sampled_start==0, we would have problems with indices below

waveform(1:sampled_start)=holdcmd_from*ones(size(1:sampled_start));
waveform(sampled_start+1:(sampled_start+sampled_ramp))= linspace(holdcmd_from, holdcmd_to, sampled_ramp);
waveform((sampled_start+sampled_ramp):end)=holdcmd_to*ones(size(waveform((sampled_start+sampled_ramp):end)));
waveform(end-length(pulsetrain)+1:end)=pulsetrain;





waveform=waveform./factor;    % factor for I-clamp with beta==1

%    waveform(end)=0;   % just quickly make it 0, so that the pulses don't stay high...