%simple low-level script to test soundcard output

deviceid=GetAsioLynxDevice; %note: using the asio device for both input and output
PsychPortAudio('Verbosity', 0); 
InitializePsychSound(0); %  InitializePsychSound([reallyneedlowlatency=0])

reqlatencyclass=0; %Rig1: saw dropouts with 0, fewer 1, fewer with 2, still some with 3
suggestedLatency=.01;
SoundFs=192e3;
buffSize = 1024;
playbackonly=1;
numChan=6;
PPAhandle = PsychPortAudio('Open', deviceid, playbackonly, reqlatencyclass, SoundFs, numChan, buffSize, suggestedLatency);
runMode = 1; %leaves soundcard on (hot), uses more resources but may solve dropouts? mw 08.25.09: so far so good.
nreps=1;
when=0; %use this to start immediately
waitForStart=0;
PsychPortAudio('RunMode', PPAhandle, runMode);
t=1:SoundFs*30;
t=t/SoundFs;
tone=sin(2*pi*1000*t);
for i=1:numChan
    samples(i,:)=tone;
end

%load and play tone
PsychPortAudio('FillBuffer', PPAhandle, samples); % fill buffer now, start in PlaySound
PsychPortAudio('Start', PPAhandle,nreps,when,waitForStart);


%tidy up
% PsychPortAudio('Stop', PPAhandle);
% PsychPortAudio('Close', PPAhandle);
