function tonetrain=MakeCFTrainTC(varargin)
% Creates a CFTrainTC type tone train (a tone train at CF followed by a
% test tone)
% 
% Args:
%   cf: the freq of the tones in train
%   ntones:   number of tones
%   isi: inter-stimulus interval, i.e. interval between the
%        start of previous tone and start of the next tone
%   toneduration:  duration of an individual tone (ms)
%   frequency: frequency of the tone (Hz)
%   amplitude:  tone amplitude (dB)
%   start:  start of the first tone after the trigger (ms)
%   duration: total duration of the tone train (ms) %%???not used? mw 04.13.06
%   next: inter-tone-train-interval, i.e. when the next
%         tone train should follow the previous one (ms)
%   ramp: rising/falling edge of individual tones
%   samplerate: sampling rate (Hz)
%
% Returns:
%     pulse   -   the required  waveform
% 
%mw 111907
global exper pref

tonetrain=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
start           =params.start/1000;              % in s
toneduration   =params.toneduration/1000;      % in s
frequency=params.cf;
%frequency=params.frequency;
amplitude       =params.amplitude;
ntones         =params.ntones;
isi             =params.isi/1000;                % in s
ramp            =params.ramp;

train_length=(start+toneduration+(ntones-1+1)*isi); % in s

%prepare the samples
tonetrain=zeros(ceil(train_length*samplerate),1);

sampled_duration=round(toneduration*samplerate);
sampled_start=floor(start*samplerate);
sampled_isi=round(isi*samplerate);

t=1/samplerate:1/samplerate:toneduration;                  % length of the sampled tone
tone=sin(frequency*2*pi*t);       % the new tone itself
%    tone=randn(1,sampled_duration);       % corresponds to t=0:1/samplerate:duration;
[edge,ledge]=MakeEdge(ramp,samplerate);         % and add the edge
tone(1:ledge)=tone(1:ledge).*fliplr(edge);
tone((end-ledge+1):end)=tone((end-ledge+1):end).*edge;

tone_starts=[0:(ntones-1)]';
sampled_start=max(1,sampled_start); % if sampled_start==0, we would have problems with indices below
tone_starts=sampled_start+tone_starts*(sampled_isi);

widths=0:sampled_duration-1;

idx=tone_starts(:,ones(1,sampled_duration))+widths(ones(1,ntones),:);
tone=tone(ones(1,ntones),:);

tonetrain(idx)=tone;
tonetrain=tonetrain./(max(abs(tonetrain)));                 % normalize, so we could fit to +/-10V
amplitude=10*(10.^((params.ca-pref.maxSPL)/20));
tonetrain=amplitude.*tonetrain;

%next add the test tone at the end
tone=sin(params.frequency*2*pi*t);       % the new tone itself
if params.frequency==-1
    tone=randn(size(tone));
end
tone(1:ledge)=tone(1:ledge).*fliplr(edge);
tone((end-ledge+1):end)=tone((end-ledge+1):end).*edge;
tone=tone./max(abs(tone)); %normalize
amplitude=10*(10.^((params.amplitude-pref.maxSPL)/20));
tone=amplitude.*tone;

tone_starts=sampled_start+(ntones)*(sampled_isi);
idx=tone_starts(:,ones(1,sampled_duration))+widths;
tonetrain(idx)=tone;
shg
%tonetrain=tonetrain./(max(abs(tonetrain)));                 % normalize, so we could fit to +/-10V

%keyboard