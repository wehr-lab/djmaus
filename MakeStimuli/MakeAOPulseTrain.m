function pulsetrain=MakeAOPulseTrain(varargin)

global exper pref

% Creates a pulse train to be delivered to AO

% Input
%   npulses         -   number of pulses
%   isi             -   inter-stimulus interval, i.e. interval between the
%                       start of previous pulse and start of the next pulse
%   pulseduration   -   duration of an individual pulse (ms)
%   amplitude       -   pulse amplitude (volts)
%   start           -   start of the first pulse after the trigger (ms)
%   duration        -   total duration of the pulse train (ms) 
%   next            -   inter-pulse-train-interval, i.e. when the next
%                       pulse train should follow the previous one (ms)
%   samplerate      -   sampling rate (Hz)
% Output
%   pulsetrain   -   the required square waveform
%
%modified from MakeClickTrain, mw 060811
%mw 061711

pulsetrain=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
start           =params.start/1000;              % in s
pulseduration   =params.pulseduration/1000;      % in s
amplitude       =params.amplitude;
npulses         =params.npulses;
isi             =params.isi/1000;                % in s

if isfield(params,'channel')
    channel=params.channel;
else
    channel=1;                  % default channel
end


    train_length=(start+pulseduration+(npulses-1)*isi); % in s
    
    %prepare the samples
    pulsetrain=zeros(ceil(train_length*samplerate),1);
    
    sampled_duration=round(pulseduration*samplerate);
    sampled_start=floor(start*samplerate);
    sampled_isi=round(isi*samplerate);

    pulse=ones(1,sampled_duration);       % corresponds to t=0:1/samplerate:duration;
    
    pulse_starts=[0:(npulses-1)]';
    sampled_start=max(1,sampled_start); % if sampled_start==0, we would have problems with indices below
    pulse_starts=sampled_start+pulse_starts*(sampled_isi);
    
    widths=0:sampled_duration-1;
    
    idx=pulse_starts(:,ones(1,sampled_duration))+widths(ones(1,npulses),:);
    pulse=pulse(ones(1,npulses),:);
    
    pulsetrain(idx)=pulse;
    
    pulsetrain=pulsetrain./(max(abs(pulsetrain)));                 % normalize, so we could fit to +/-10V

    pulsetrain=amplitude.*pulsetrain;
    pulsetrain(end)=0; %make sure we don't end and hold on high %mw 080612
