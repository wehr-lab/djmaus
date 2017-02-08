function MakePINPdjProtocol(flashduration, next, nrepeats)
%usage:  MakePINPdjProtocol(trainduration, icis, flashduration, next, nrepeats)
%
% % creates an djmaus protocol file for PINPing. This consists only of
% silent sounds, with interleaved laser pulses. There are also some
% flashtrains included to test for reliability.

%
% inputs:
%   flashduration   -   duration of an individual flash (ms)
%   next            -   inter-flash-train-interval, i.e. when the next
%                       flash train should follow the previous one (ms)
%   nrepeats: number of repetitions
%
% outputs:
%       - creates a suitably named stimulus protocol in djprefs.stimuli/PINPProtocols
%
%
%example calls:
% MakePINPdjProtocol(30, 500, 50)


name= sprintf('PINP %dms' ,flashduration);
description=sprintf('PINP, %dms flashduration, %d repeats', flashduration, nrepeats);
filename=name;

nn=0;
for rep=1:nrepeats
    nn=nn+1;
    stimuli(nn).type='silentsound';
    stimuli(nn).param.laser=1;
    stimuli(nn).param.VarLaser=1;
    stimuli(nn).param.VarLaserstart=0;
    stimuli(nn).param.VarLaserpulsewidth= flashduration;
    stimuli(nn).param.VarLasernumpulses=1;
    stimuli(nn).param.VarLaserisi=0; %
    stimuli(nn).param.next=next;
    stimuli(nn).param.ramp=0;
    stimuli(nn).param.duration=flashduration;
    stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
    stimuli(nn).protocol_name=name;
    stimuli(nn).protocol_description=description;
    stimuli(nn).PlottingFunction='PlotPINP_PSTH';
    stimuli(nn).version='djmaus';
    
    nn=nn+1;
    stimuli(nn).type='silentsound';
    stimuli(nn).param.laser=0;
    stimuli(nn).param.VarLaser=0;
    stimuli(nn).param.VarLaserstart=0;
    stimuli(nn).param.VarLaserpulsewidth= 0;
    stimuli(nn).param.VarLasernumpulses=0;
    stimuli(nn).param.VarLaserisi=0; %
    stimuli(nn).param.next=next;
    stimuli(nn).param.ramp=0;
    stimuli(nn).param.duration=flashduration;
    stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
    stimuli(nn).protocol_name=name;
    stimuli(nn).protocol_description=description;
    stimuli(nn).PlottingFunction='PlotPINP_PSTH';
    stimuli(nn).version='djmaus';
end

%add some flash trains for testing of direct/indirect activation by
%reliability
for rep=1:nrepeats
    nn=nn+1;
    stimuli(nn).type='silentsound';
    stimuli(nn).param.laser=1;
    stimuli(nn).param.VarLaser=1;
    stimuli(nn).param.VarLaserstart=0;
    stimuli(nn).param.VarLaserpulsewidth= 10;
    stimuli(nn).param.VarLasernumpulses=10;
    stimuli(nn).param.VarLaserisi=100; %
    stimuli(nn).param.next=next;
    stimuli(nn).param.ramp=0;
    stimuli(nn).param.duration=1200;
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
