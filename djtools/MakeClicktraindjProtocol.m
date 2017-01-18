function MakeClicktraindjProtocol(amplitude, trainduration, icis, clickduration, next, ramp, nrepeats)
%usage:  MakeClicktraindjProtocol(amplitude, trainduration, icis, clickduration, next, ramp, nrepeats)
%
%
% modified from MakeWNTrainProtocol2, and ported to djmaus
% we have a single amplitude and a fixed train duration
% (with variable number of clicks per train)
% creates an djmaus protocol file for a train of white noise bursts (click train)
%
% inputs:
%   amplitude:  amplitude in dB SPL (requires system to be calibrated)
%   trainduration   -   duration of train, in ms
%   icis            -   inter-click interval, i.e. interval between the
%                       start of previous click and start of the next click
%                       (use an array for multiple ICIs)
%   clickduration   -   duration of an individual click (ms)
%   start           -   start of the first click after the trigger (ms)
%   next            -   inter-click-train-interval, i.e. when the next
%                       click train should follow the previous one (ms)
%   ramp            -   rising/falling edge of individual clicks (can be zero)
%   nrepeats: number of repetitions (different pseudorandom orders)
%
% outputs:
%       - creates a suitably named stimulus protocol in
%       djprefs.stimuli/ClicktrainProtocols
%
%
%example calls: 
%MakeClicktraindjProtocol(80, 10e3, [1 2 4 8 16 32 64 128 256], .5, 3000, 0, 20) 
numicis=length(icis);


% % % neworder=randperm( numicis );
% % % interclickintervals=zeros(size(neworder*nrepeats));
% % % 
% % % tdur=0;
% % % total_trainduration=trainduration;
% % % tdur=numicis*(total_trainduration)/1000;%duration per repeat


for nn=1:nrepeats
    neworder=randperm( numicis );
    interclickintervals(  icis( neworder );
end

icistring=sprintf('%d-', icis);icistring=icistring(1:end-1);
%put into stimuli structure
stimuli(1).type='exper2 stimulus protocol';

stimuli(1).param.name= sprintf('WNTrain2, %da/min%d,max%ddB%dms/%sms/d=%d', numamplitudes,minamplitude,maxamplitude ,clickduration, isistring, trainduration);
stimuli(1).param.description=sprintf('White noise train2, %d ampl. (%d-%d dB SPL), %dms clickduration, %.1fms ramp, %d repeats, %d ms trainduration, %d ISIs (%sms), %ds duration per repeat', numamplitudes,minamplitude, maxamplitude, clickduration, ramp, nrepeats, trainduration, numicis, isistring, round(tdur));
filename=sprintf('WNTrain2-%da-%ddB-%ddB-%dms-%sms-d%d.mat',numamplitudes,minamplitude,maxamplitude,clickduration, isistring,trainduration);

nn=0;
for rep=1:nrepeats
    for n=1:numicis
        nn=nn+1;
    stimuli(nn).type='clicktrain';
    stimuli(nn).param.amplitude=amplitude;
   
    nclicks=floor(trainduration/icis(n));
    
    stimuli(nn).param.nclicks=nclicks;
    stimuli(nn).param.isi=interclickintervals(nn-1);
    stimuli(nn).param.clickduration=clickduration;
    stimuli(nn).param.start=start;
    stimuli(nn).param.next=next;
    stimuli(nn).param.ramp=ramp;
    stimuli(nn).param.duration=total_trainduration;
end


global pref
Prefs
cd(pref.stimuli)
cd ('WNTrain2 Protocols')

save(filename, 'stimuli')


% keyboard