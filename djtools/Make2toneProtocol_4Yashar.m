% this was written to make stimuli for Yashar. It creates white noise
% bursts and laser pulses, not interleaved. ira 12.05.17
dbstop if error
durs=25;
isis=[50 100 150 200 400 800];
type= {'two tone', 'silentsound'};
amps=80;
iti=3000;
nreps=30;
num_stim=length(isis)*2*nreps+2*nreps; %2 for WN and laser 2 tone and probe alone
N=randperm(num_stim);
description=sprintf('2tone WN only, 1 amplitude(80-80 dB SPL), 1 duration (25 ms), 3000 iti, 6 SOAs(50, 100, 150, 200, 400, 800 ms), laser with silentsound only, 30 reps');
name=sprintf('2tone WN only, 1a(80-80dB)/1d(25ms)/IL silentsound/30reps');
m=0;
%for i=1:length(N)

for i=1:nreps
    for is=1:length(isis)
        m=m+1;
        k=N(m);
        stimuli(k).type='2tone';
        stimuli(k).param.amplitude=80;
        stimuli(k).param.frequency=-1;
        stimuli(k).param.probefreq=-1;
        stimuli(k).param.probeamp=80;
        stimuli(k).param.duration=25;
        stimuli(k).param.ramp=5;
        stimuli(k).param.next=3000;
        stimuli(k).param.SOA=isis(is);
        
        stimuli(k).param.laser=0;
        stimuli(k).param.VarLaser=0;
        stimuli(k).param.VarLaserstart=0;
        stimuli(k).param.VarLaserpulsewidth=0;
        stimuli(k).param.VarLasernumpulses=0;
        stimuli(k).param.VarLaserisi=0;
        
        stimuli(k).stimulus_description=GetParamStr(stimuli(k));
        stimuli(k).protocol_name=name;
        stimuli(k).protocol_description=description;
        stimuli(k).PlottingFunction='Plot2Tone_PSTH_single2';
        stimuli(k).version='djmaus';
    end

    for is=1:length(isis)
        m=m+1;
        k=N(m);
        stimuli(k).type='silentsound';
        
        stimuli(k).param.duration=25;
        stimuli(k).param.amplitude=0;
        stimuli(k).param.frequency=0;
        stimuli(k).param.probefreq=0;
        stimuli(k).param.probeamp=0;
        stimuli(k).param.ramp=0;
        stimuli(k).param.next=3000;
        stimuli(k).param.SOA=isis(is);
      
        stimuli(k).param.laser=1;
        stimuli(k).param.VarLaser=1;
        stimuli(k).param.VarLaserstart=0;
        stimuli(k).param.VarLaserpulsewidth=25;
        stimuli(k).param.VarLasernumpulses=2;
        stimuli(k).param.VarLaserisi=isis(is);
        
        stimuli(k).stimulus_description=GetParamStr(stimuli(k));
        stimuli(k).protocol_name=name;
        stimuli(k).protocol_description=description;
        stimuli(k).PlottingFunction='Plot2Tone_PSTH_single2';
        stimuli(k).version='djmaus';
    end
    
        m=m+1;
        k=N(m);
        stimuli(k).type='whitenoise';
        
        stimuli(k).param.amplitude=80;
        stimuli(k).param.frequency=-1;
        stimuli(k).param.duration=25;
        stimuli(k).param.ramp=5;
        stimuli(k).param.next=3000;
        stimuli(k).param.probefreq=0;
        stimuli(k).param.SOA=0;
        
        stimuli(k).param.laser=0;
        stimuli(k).param.VarLaser=0;
        stimuli(k).param.VarLaserstart=0;
        stimuli(k).param.VarLaserpulsewidth=0;
        stimuli(k).param.VarLasernumpulses=0;
        stimuli(k).param.VarLaserisi=0;
        
        stimuli(k).stimulus_description=GetParamStr(stimuli(k));
        stimuli(k).protocol_name=name;
        stimuli(k).protocol_description=description;
        stimuli(k).PlottingFunction='Plot2Tone_PSTH_single2';
        stimuli(k).version='djmaus';
        
        m=m+1;
        k=N(m);
        stimuli(k).type='silentsound';
        
        stimuli(k).param.duration=25;
        stimuli(k).param.amplitude=0;
        stimuli(k).param.frequency=0;
        stimuli(k).param.probefreq=0;
        stimuli(k).param.probeamp=0;
        stimuli(k).param.ramp=0;
        stimuli(k).param.next=3000;
        stimuli(k).param.SOA=0;
        
        stimuli(k).param.laser=1;
        stimuli(k).param.VarLaser=1;
        stimuli(k).param.VarLaserstart=0;
        stimuli(k).param.VarLaserpulsewidth=25;
        stimuli(k).param.VarLasernumpulses=1;
        stimuli(k).param.VarLaserisi=0;
        
        stimuli(k).stimulus_description=GetParamStr(stimuli(k));
        stimuli(k).protocol_name=name;
        stimuli(k).protocol_description=description;
        stimuli(k).PlottingFunction='Plot2Tone_PSTH_single2';
        stimuli(k).version='djmaus';
end
save('Stimuli_for_Yashar_2toneWN+2tonelaser_3siti.mat', 'stimuli')
%end