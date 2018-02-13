function MakeSustainedSuppressionProtocol(duration, next, nrepeats)

% % creates a specific djmaus protocol file for validating sustained
% flashtrain suppression. In behavior we have the arduino flashing the
% laser continuously (5ms on- 5ms off) for 100 seconds, or longer. We just
% want to validate the electrophysiological effect of that.

% uses silent sounds, with interleaved laser pulses. There are also some
% flashtrains included to test for reliability.

%
% inputs:
%   duration   -   how long to have the train on, in seconds
%   next            -   inter-flash-train-interval, i.e. when the next
%                       flash train should follow the previous one (ms)
%   nrepeats: number of repetitions
%
% outputs:
%       - creates a suitably named stimulus protocol in djprefs.stimuli/PINPProtocols
%
%
%example calls:
% MakePINPdjProtocol(100, 10000, 20)

pw=5; %flash pulsewidth in ms
%isi is also set to pw, so it's 50% duty cycle
numpulses= duration*1000/(2*pw);

name= sprintf('Sustained Suppression %d s' ,duration);
description=sprintf('Sustained Suppression, %dms duration, %d repeats', duration, nrepeats);
filename=name;

nn=0;

for rep=1:nrepeats
    nn=nn+1;
    stimuli(nn).type='silentsound';
    stimuli(nn).param.laser=1;
    stimuli(nn).param.VarLaser=1;
    stimuli(nn).param.VarLaserstart=0;
    stimuli(nn).param.VarLaserpulsewidth= pw;
    stimuli(nn).param.VarLasernumpulses=numpulses;
    stimuli(nn).param.VarLaserisi=pw; %
    stimuli(nn).param.next=next;
    stimuli(nn).param.ramp=0;
    stimuli(nn).param.duration=duration*1000;
    stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
    stimuli(nn).protocol_name=name;
    stimuli(nn).protocol_description=description;
    stimuli(nn).PlottingFunction='PlotPINP_PSTH';
    stimuli(nn).version='djmaus';
end

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('PINPProtocols')
cd ('PINPProtocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);
