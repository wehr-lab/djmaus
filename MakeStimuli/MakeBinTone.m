function tone=MakeBinTone(varargin)


% Note: this was the original header
%function tone=MakeTone(frequency,attenuation,duration,samplerate,ramp)

global pref

% Creates a binaural pure tone of a given frequency, R and L amplitude, duration, at a
% given sample rate, with an ascending/descending ramp of a given length.
% channel 1 is right, channel 2 is left
% All lengths (duration, ramp) are in ms
% frequency is in Hz, amplitude is in absolute dB SPL
% Input:
%  frequency          -   frequency of the tone (Hz)
% % % % % % % % % % % %  attenuation        -   attenuation of the tone (dB) 
% % % % % % % % % % % %                         !!!!!attenuation is relative to the max. sound pressure level 
% % % % % % % % % % % %                         as specified in Prefs.m (pref.maxSPL)
%  amplitude          -   sound pressure level of the tone (dB)            
%  duration           -   duration of the tone (ms)
%  samplerate         -   required sampling rate
%  ramp               -   length of an rising/falling edge (ascending/descending ramp) (ms)
%  
%  Output:
%  tone               -   the specified tone (empty if unsuccessful)
%  
%modified from MakeTone %mw 051409

tone=[];

if nargin<2
    return;
end

params=varargin{1};
samplerate=varargin{2};
frequency=params.frequency;
% attenuation=params.attenuation;
Ramplitude=params.Ramplitude;
Lamplitude=params.Lamplitude;
duration=params.duration;
ramp=params.ramp;



%     amplitude=10*(10.^((-attenuation)/20));
%    amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-10<x<10), i.e. pref.maxSPL=+_10V
    Ramplitude=1*(10.^((Ramplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
    Lamplitude=1*(10.^((Lamplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
    duration=duration/1000;                     % adjust the duration to seconds
    t=0:1/samplerate:duration;                  % length of the sampled trial
    Rtone=Ramplitude*sin(frequency*2*pi*t);       % the new tone itself
    Ltone=Lamplitude*sin(frequency*2*pi*t);       % the new tone itself
    if ramp>0
        [edge,ledge]=MakeEdge(ramp,samplerate);     % prepare the edges
        Rtone(1:ledge)=Rtone(1:ledge).*fliplr(edge);
        Rtone((end-ledge+1):end)=Rtone((end-ledge+1):end).*edge;
        
        Ltone(1:ledge)=Ltone(1:ledge).*fliplr(edge);
        Ltone((end-ledge+1):end)=Ltone((end-ledge+1):end).*edge;
    end

  tone(1,:)=Rtone;
  tone(2,:)=Ltone;
  

