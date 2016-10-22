function Vrms=inversedBSPL(dBSPL, sensitivity)
%what Vrms do you need on the B&K mic to achieve a desired dBSPL?
%usage: Vrms=inversedBSPL(dBSPL, sensitivity)
%inputs: 
%   dBSPL: the desired power of the sound in dB SPL (re: 20 microPa)
%   sensitivity: sensitivity setting on the B&K amplifier, in volts 
%                   (e.g. .1, .316, 1  V/Pa)
%output: Vrms: cycle RMS voltage you need to measure from the output of the
%B&Kmicrophone that corresponds to your desired dBSPL

Vrms=sensitivity*20e-6*10^(dBSPL/20);