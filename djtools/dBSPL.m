function dBSPL=dBSPL(Vrms, sensitivity)
%usage: dBSPL=dBSPL(Vrms, sensitivity)
%inputs: 
%   Vrms: cycle RMS voltage measured from the output of the B&K microphone
%   sensitivity: sensitivity setting on the B&K amplifier, in volts 
%                   (e.g. .1, .316, 1  V/Pa)
%   Vrms = (peak to peak)/sqrt(2)
%output: the power of the sound in dB SPL (re: 20 microPa)

dBSPL=20*log10((Vrms/sensitivity)/20e-6);