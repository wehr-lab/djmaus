function MakeFlashtraindjProtocol(trainduration, icis, flashduration, next, nrepeats)
%usage:  MakeFlashtraindjProtocol(trainduration, icis, flashduration, next, nrepeats)
%
% % creates an djmaus protocol file for flash trains of desired rate/duration. No sound.
% modified from MakeClicktraindjProtocol - designed to mimic it optogenetically
% uses fixed train duration
% (with variable number of flashes per train)
%
% inputs:
%   trainduration   -   duration of train, in ms
%   icis            -   inter-flash interval, i.e. interval between the
%                       start of previous flash and start of the next flash
%                       (use an array for multiple ICIs)
%   flashduration   -   duration of an individual flash (ms)
%   next            -   inter-flash-train-interval, i.e. when the next
%                       flash train should follow the previous one (ms)
%   nrepeats: number of repetitions (different pseudorandom orders)
%
% outputs:
%       - creates a suitably named stimulus protocol in djprefs.stimuli/ClicktrainProtocols
%
%
%example calls:
%MakeFlashtraindjProtocol(2e3, [1 2 4 8 16 32 64 128 256], .5, 2000, 20)

numicis=length(icis);
icistring=sprintf('%d-', icis);icistring=icistring(1:end-1);

name= sprintf('Flashtrain, %.1fms/%sms/d=%d' ,flashduration, icistring, trainduration);
description=sprintf('Flashtrain, %.1fms flashduration, %d repeats, %d ms trainduration, %d ICIs (%sms), %ds duration per repeat', flashduration, nrepeats, trainduration, numicis, icistring);
filename=sprintf('Flashtrain-%.1fms-%sms-d%d.mat',flashduration, icistring,trainduration);

nn=0;
for rep=1:nrepeats
    for n=randperm(numicis)
        nn=nn+1;
        stimuli(nn).type='silentsound';
        nflashes=floor(trainduration/icis(n));
        stimuli(nn).param.laser=1;
        stimuli(nn).param.VarLaser=1;
        stimuli(nn).param.VarLaserstart=0;
        stimuli(nn).param.VarLaserpulsewidth= flashduration;
        stimuli(nn).param.VarLasernumpulses=nflashes;
        stimuli(nn).param.VarLaserisi=icis(n); %
        stimuli(nn).param.next=next;
        stimuli(nn).param.ramp=0;
        stimuli(nn).param.duration=trainduration;
        stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
        stimuli(nn).protocol_name=name;
        stimuli(nn).protocol_description=description;
        stimuli(nn).PlottingFunction='PlotFlashtrain_PSTH';
        stimuli(nn).version='djmaus';
    end
end

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('ClicktrainProtocols')
cd ('ClicktrainProtocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);
