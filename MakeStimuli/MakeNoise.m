function noise = MakeNoise(varargin);

global pref

% NOTE: this was the original header
% function y = MakeNoise(filter_operation, cf, lower_frequency_limit, upper_frequency_limit, duration, amplitude, samplerate, ramp);

% [This function generates noise stimuli]
%
% function y = noise(filter_operation, cf, lower_frequency_limit, upper_frequency_limit, duration, SR, rise_fall_time);
%
% filter_operation =
%    1 = wide-band
%    2 = band-pass
%    3 = band-reject
%    4 = low-pass
%    5 = high-pass
% cf                    - center-frequency (for band-pass, band-reject) or
%                         cutoff-frequency (for low-pass, high-pass) in Hz
% lower_frequency_limit - lower frequency limit in Hz (for bandpass, bandreject)
% upper_frequency_limit - upper frequency limit in Hz (for bandpass, bandreject)
% duration              - duration (ms)
% amplitude             - amplitude (dB)
% samplerate            - sampling rate (Hz)
% ramp                  - rise_fall_time (ms)

% based on orig. from Wang Lab, Johns Hopkins University (Edited on January 13, 2004 by Tom Lu)

noise=[];
if nargin<2
    return;
end

filters={'wideband','bandpass','bandreject','lowpass','highpass'};

params          =varargin{1};
samplerate      =varargin{2};
filter_operation=find(strcmp(filters,params.filter_operation));
duration        =params.duration;
amplitude       =params.amplitude;
ramp            =params.ramp;

% amplitude=10*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-10<x<10), i.e. pref.maxSPL=+_10V
amplitude=1*(10.^((amplitude-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V

switch filter_operation
    case 1;
        cf=0;
        lower_frequency_limit=0;
        upper_frequency_limit=0;
    case 2, 3;
        cf=params.center_frequency;
        lower_frequency_limit=params.lower_frequency;
        upper_frequency_limit=params.upper_frequency;
    case 4, 5;
        cf=params.cutoff_frequency;
        lower_frequency_limit=0;
        upper_frequency_limit=0;
end

samplerate=samplerate/1000;

%onset = ramp; % Wang
noise_length = round(samplerate*duration);
onset_length = round(samplerate*ramp);

noise_data = randn(1,noise_length);
start = lower_frequency_limit/1000; % transform to kHz to comply with the original Wang script
stop  = upper_frequency_limit/1000;
cf=cf/1000;

if onset_length >= noise_length
   error(['Error:  Rise/Fall Time is longer than stimulus Duration']);
end

% our (cos) ramp
    [edge,ledge]=MakeEdge(ramp,samplerate*1000);     % prepare the edges
    noise_data(1:ledge)=noise_data(1:ledge).*fliplr(edge);
    noise_data((end-ledge+1):end)=noise_data((end-ledge+1):end).*edge;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Original rising/falling edge from Wang
% for i=1:onset_length
%    noise_data(i)= noise_data(i) *i/onset_length;
% end
% 
% for i=1:onset_length
%    noise_data(noise_length-i) = noise_data (noise_length-i) * i/onset_length;
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch filter_operation
 case 1                      % wide-band
   noise_filt = noise_data;

 case 2                      % band-pass
   if stop > 49.4 & start >= 0.2
      b = fir2(200,[0 start*.99 start (samplerate/2)]/(samplerate/2),[0 0 1 1]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   elseif stop <= 49.4 & start < 0.2
      b = fir2(200,[0 stop stop*1.01 (samplerate/2)]/(samplerate/2),[1 1 0 0]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   elseif stop >49.4 & start < 0.2
      noise_filt = noise_data;      
   else
      b = fir2(200,[0 start*.99 start stop stop*1.01 samplerate/2]/(samplerate/2),[0 0 1 1 0 0]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   end
if ~range(b)
    warning('MakeNoise: ill conditioned filter operation')
end
 case 3                      % band-reject
   if stop > 49.4 & start >= 0.2
      b = fir2(200,[0 start*.99 start (samplerate/2)]/(samplerate/2),[1 1 0 0]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   elseif stop <= 49.4 & start < 0.2
      b = fir2(200,[0 stop stop*1.01 (samplerate/2)]/(samplerate/2),[0 0 1 1]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   elseif stop > 49.4 & start < 0.2
      
      if cf < (samplerate/10)
         disp(['Rejection Region has overrun cutoff set for 49.5 KHz']);
         stop = 49.5;
         b = fir2(200,[0 stop*.99 stop (samplerate/2)]/(samplerate/2),[0 0 1 1]);
         noise_filt=filtfilt(b,1,noise_data);
      else  % cf > samplerate/4
         disp(['Rejection Region has overrun cutoff set for 0.1 KHz']);
         stop = 0.1;
         b = fir2(200,[0 stop stop*1.01 (samplerate/2)]/(samplerate/2),[1 1 0 0]);
         noise_filt=filtfilt(b,1,noise_data);
      end
   else
      b = fir2(200,[0 start*.99 start stop stop*1.01 samplerate/2]/(samplerate/2),[1 1 0 0 1 1]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   end

 case 4                      % low-pass
   if cf < 49.4
      b = fir2(200,[0 cf*.99 cf (samplerate/2)]/(samplerate/2),[1 1 0 0]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   else
      noise_filt = noise_data;
   end 

 case 5                      % high-pass
   if cf > 0.2
      b = fir2(200,[0 cf cf*1.01 (samplerate/2)]/(samplerate/2),[0 0 1 1]);
      a = 1;
      noise_filt=filtfilt(b,a,noise_data);
   else
      noise_filt = noise_data;
   end   

end % switch


% now 'normalize' and put in some amplitude
y = amplitude.*(noise_filt/max(abs(noise_filt)));
noise=y;

% figure
% subplot(211)
% plot((1:length(y))/samplerate*1000/1000, y)
% xlabel('Time (ms)')
% ylabel('amplitude')
% 
% subplot(212)
% plot((1:length(fft(y)))/length(fft(y))*(samplerate/1000), abs(fft(y)))
% specgram(y,[],samplerate)
% ylabel('Frequency (kHz)');
