function dBPSD=dBPSD(Px, sensitivity)
%usage: dBPSD=dBPSD(Px, sensitivity)
%inputs: 
%   Px: power spectral density of the voltage measured from the output of the B&K microphone
%        this is one element of the Pxx array returned by pwelch
%   sensitivity: sensitivity setting on the B&K amplifier, in volts 
%                   (e.g. .1, .316, 1  V/Pa)
%output: the power of the sound in dB SPL (re: 20 microPa)
%
%algorithm: worked backwards from measurement using the B&K calibrator
%at .316 v/Pa:
%94 dB produces a Pxx of .0076 (.31 Vrms, 1 Vp-p)
%114 dB produces a Pxx of .7620 (3.1 Vrms, 9Vp-p)
%
%at 1 v/Pa:
%94 dB produces a Pxx of .0757 (1 Vrms, 3 Vp-p)
%114 dB produces a Pxx of 5.6720 (9.75 Vrms, 27 Vp-p)
%
%at 3.16 v/Pa:
%94 dB produces a Pxx of .7620 (3.10 Vrms
%114 dB produces a Pxx of 6.96 (clips on scope)

%the Pxx is in units of power(dB)/Hz


switch sensitivity
    case .1
        F= 7.57e-04;
    case .316
        F=.00757;
    case 1
        F=.0757;
    case 3.16
        F=.757;
    otherwise
        error('expected sensitivity of .316, 1, or 3.16')
end
dBPSD=94+10*log10(Px/F);

