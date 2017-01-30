function Make2ToneWNProtocol(freqsperoctave, minfreq, maxfreq, numamplitudes, ...
    minamplitude, maxamplitude, include_whitenoise, duration, ramp, isi, nrepeats, numprobeamplitudes, ...
    minprobeamplitude, maxprobeamplitude, numSOA, minSOA, maxSOA, interleave_laser, include_silent_sound)
%usage: Make2ToneProtocol(freqsperoctave, minfreq, maxfreq, numamplitudes, ...
%     minamplitude, maxamplitude, include_whitenoise, duration, ramp, isi, nrepeats, probefreq, numprobeamplitudes, ...
%     minprobeamplitude, maxprobeamplitude, numSOA, minSOA, maxSOA, interleave_laser, include_silent_sound)
%
% modified from Make2ToneProtocol to use WN only (Make2ToneProtocol does not allow WN)
% mw 10.30.2015
%
% creates an exper2 stimulus protocol file for a 2 tone protocol stimulus
% with probe tones of single frequency but multiple amplitudes
% no white noise used
%
% inputs:
% freqsperoctave: number of frequencies per octave (frequency resolution)
% minfreq: lowest frequency in Hz
% maxfreq: highest frequency in Hz
%           note: if minfreq-to-maxfreq cannot be divided evenly into
%           freqsperoctave, the nearest maxfreq will be used
%           (i.e. the requested freqsperoctave will be exactly enforced)
% numamplitudes: number of masker amplitude steps
% minamplitude: maximum masker amplitude in dB SPL (requires system to be calibrated)
% maxamplitude: maximum masker amplitude in dB SPL (requires system to be calibrated)
% include_whitenoise: 1 for yes and 0 for no
% duration: masker duration in ms
% ramp: on-off ramp duration in ms
% isi: inter stimulus interval (onset-to-onset) in ms
% nrepeats: number of repetitions (different pseudorandom orders)
% numprobeamplitudes: number of probe tone amplitude steps
% minprobeamplitude: minimum probe tone amplitude in dB SPL
% maxprobeamplitude: maximum probe tone amplitude in dB SPL
% numSOA: number of Stimulus Onset Asynchrony in ms = time between masker onset and probe tone onset
% minSOA: the shortest SOA
% maxSOA: the longest SOA

% interleave_laser: 1 for yes and 0 for no, x2 nreps with randomly inserted
% include_silent_sound: 0 or 1 to include an extra stimulus that is silent,
%           useful for explicitly collecting spontaneous and/or laser-only
%           trials as a separate condition
% laser trials
%
% outputs:
% creates a suitably named stimulus protocol in
% home/lab/djmaus/stimuli/2Tone
%
%
%example calls: numSOA
% Make2ToneWNProtocol(0,0,0,1, 70, 70,1, 25, 3,1000, 10, 1, 70,70, 1, 100, 100 1, 0)
% Make2ToneWNProtocol(3, 4e3, 20e3, 1, 70, 70, 0, 400, 3, 1000, 10, 1, 70,70,2, 100, 500, 1, 1)
%
% ira 01.28.17

if nargin==0; fprintf('\nno input');return;end
if nargin>19; error('\nMake2ToneProtocol: wrong number of arguments.'); end

if freqsperoctave==0 & minfreq==0 & maxfreq==0 & include_whitenoise==1
    %white noise only
    WNonly=1;
    logspacedfreqs=[];
    numfreqs=0;
    fprintf('\nwhite noise only')
else
    WNonly=0;
    numoctaves=log2(maxfreq/minfreq);
    logspacedfreqs=minfreq*2.^([0:(1/freqsperoctave):numoctaves]);
    newmaxfreq=logspacedfreqs(end);
    numfreqs=length(logspacedfreqs);
    if maxfreq~=newmaxfreq
        fprintf('\nnote: could not divide %d-%d Hz evenly into exactly %d frequencies per octave', minfreq, maxfreq, freqsperoctave)
        fprintf('\nusing new maxfreq of %d to achieve exactly %d frequencies per octave\n', round(newmaxfreq), freqsperoctave)
        maxfreq=newmaxfreq;
    end
end

linspacedamplitudes = linspace( minamplitude , maxamplitude , numamplitudes );
linspacedprobeamplitudes = linspace( minprobeamplitude , maxprobeamplitude , numprobeamplitudes );
linspacedprobeSOAs=linspace(minSOA, maxSOA, numSOA);

if include_whitenoise==1
    logspacedfreqs=[logspacedfreqs -1]; %add whitenoise as extra freq=-1
    numfreqs=numfreqs+1;
end

if interleave_laser==1
    %[Amplitudes, Freqs, Durations, ProbeAmplitudes, probeSOAs, Lasers]=meshgrid(linspacedamplitudes, logspacedfreqs, duration, linspacedprobeamplitudes, linspacedprobeSOAs,[0 1] );
    numlasers=2;
else
    %[Amplitudes,Freqs, Durations, ProbeAmplitudes, probeSOAs]=meshgrid(linspacedamplitudes, logspacedfreqs, durations, linspacedamplitudes, linspacedprobeSOAs );
    numlasers=1;
    Lasers=zeros(size(Amplitudes));
end

StimPerRepeat=numfreqs *numfreqs * numamplitudes *numprobeamplitudes  * numSOA * numlasers; %these are only 2 tone stimuli per repreat

%add single tones or white noise stimuli

StimPerRepeat=StimPerRepeat+(numfreqs * numamplitudes * numlasers); %add single tones or WN stimuli to 2tonestimuli


if include_silent_sound
    StimPerRepeat=StimPerRepeat+numlasers;
end

TotalNumStim=StimPerRepeat*nrepeats;
DurationPerRepeatSecs=StimPerRepeat*(mean(duration)+isi)/1000;%approx. duration per repeat
TotalDurationSecs=DurationPerRepeatSecs*nrepeats;
%durstring=sprintf('%d-', duration);durstring=durstring(1:end-1);

% %put into stimuli structure
% stimuli(1).param.name= sprintf('2TPWN-%da(%d-%d)-%dms-%dpa(%d-%d)-isi-%d-%d-%d-%d-SOAs',...
%     numamplitudes,minamplitude, maxamplitude, duration,...
%     numprobeamplitudes, minprobeamplitude, maxprobeamplitude, isi, numSOA, minSOA, maxSOA);
% stimuli(1).param.description=sprintf(...
%     '2ToneWN Protocol, %d ampl. (%d-%d dB SPL), %dms duration, %d probe ampl (%d-%d dB SPL), %d-%d-%d-SOAs', ...
%     numamplitudes,minamplitude, maxamplitude, duration, ...
%     numprobeamplitudes, minprobeamplitude, maxprobeamplitude, numSOA, minSOA, maxSOA);

% filename and discreption

if interleave_laser
    interleave_laserstr='IL';
else
    interleave_laserstr='';
end
if include_silent_sound
    include_silent_soundstr='ISS';
else
    include_silent_soundstr='';
end
if include_whitenoise
    include_whitenoisestr='+WN';
else
    include_whitenoisestr='';
end

if WNonly
    name= sprintf('2Tone WN only, %da(%d-%ddB)/%dprobea(%d-%ddB)/%dd(ms)/%s-%s-%dmsisi/%d reps/%d (%d-%dms) SOAs', ...
        numamplitudes,minamplitude,...
        maxamplitude, numprobeamplitudes, minprobeamplitude, maxprobeamplitude, duration, interleave_laserstr,include_silent_soundstr, isi,nrepeats, numSOA, minSOA, maxSOA);
    description=sprintf('2Tone WN only, %d ampl. (%d-%d dB SPL), %d probe amp. (%d-%ddB), %d duration (ms),%s, %s, %d ms isi,%d stim per rep, %d repeats, %d total stimuli, %ds per rep, %d (%d-%dms) SOAs, %d s total dur',...
        numamplitudes,minamplitude, maxamplitude, numprobeamplitudes, minprobeamplitude, maxprobeamplitude, duration, interleave_laserstr,include_silent_soundstr, isi, StimPerRepeat, nrepeats,TotalNumStim, round(DurationPerRepeatSecs), numSOA, minSOA, maxSOA,round(TotalDurationSecs));
    filename=sprintf('2Tone-WNonly-%da_%d-%ddB-%dprobea_%d-%ddB-%dd_ms-%s-%s-isi%dms-%dreps-%ds%d_%d-%dms_SOAs',...
        numamplitudes,minamplitude, maxamplitude, numprobeamplitudes, minprobeamplitude, maxprobeamplitude,duration, ...
        interleave_laserstr,include_silent_soundstr, isi, nrepeats, numSOA, minSOA, maxSOA);
else
    name= sprintf('2Tone %s, %dfpo(%d-%dHz)/%da(%d-%ddB)/%dprobea(%d-%ddB)/%dd(ms)/%s-%s-%dmsisi/%d reps/%d-%d-%d SOAs', ...
        include_whitenoisestr, freqsperoctave,minfreq, round(maxfreq), numamplitudes,minamplitude,...
        maxamplitude, numprobeamplitudes, minprobeamplitude, maxprobeamplitude, duration,interleave_laserstr,include_silent_soundstr, isi,nrepeats, numSOA, minSOA, maxSOA);
    description=sprintf('2Tone, tones%s, %d freqs/oct (%d-%dkHz), %d ampl. (%d-%d dB SPL), %d probe amp. (%d-%ddB), %d duration (ms),%s, %s, %d ms isi,%d stim per rep, %d repeats, %d total stimuli, %ds per rep, %d-(%d-%dms)SOAs, %d s total dur',...
        include_whitenoisestr,freqsperoctave, minfreq, round(maxfreq), numamplitudes,minamplitude, maxamplitude, numprobeamplitudes, minprobeamplitude, maxprobeamplitude, duration, ...
        interleave_laserstr,include_silent_soundstr, isi, StimPerRepeat, nrepeats,TotalNumStim,round(DurationPerRepeatSecs), numSOA, minSOA, maxSOA, round(TotalDurationSecs));
    filename=sprintf('2Tone-tones-%s-%dfpo_%d-%dkHz_%da_%d-%ddB-%dprobea_%d-%ddB-%dd_ms-%s-%s-isi%dms-%dreps-%ds%d_%d-%dms_SOAs',...
        include_whitenoisestr, freqsperoctave,minfreq, round(maxfreq), numamplitudes,minamplitude, maxamplitude, numprobeamplitudes, minprobeamplitude, maxprobeamplitude, duration,...
        interleave_laserstr,include_silent_soundstr, isi, nrepeats, numSOA, minSOA, maxSOA);
end




% Make the stimuli


% ProbeAmps=meshgrid(linspacedprobeamplitudes, 1:nrepeats);
% neworder2=randperm(nrepeats *numprobeamplitudes );
% ProbeAmps=ProbeAmps(neworder2);
%
% Amplitudes=meshgrid(linspacedamplitudes, 1:nrepeats);
% neworder2=randperm(nrepeats * linspacedamplitudes);
% Amplitudes=Amplitudes(neworder2);
%
% neworder2=randperm(nrepeats*numSOA);
% ProbeSOAs=meshgrid(linspacedprobeSOAs, 1:nrepeats);
% ProbeSOAs=ProbeSOAs(neworder2);
%
% Freqs=meshgrid(logspacedfreqs, 1:nrepeats);
% neworder2=randperm(nrepeats*numfreqs);
% Freqs=Freqs(neworder2);



if interleave_laser==1
    [Amplitudes, ProbeAmps, Freqs, Durations, ProbeSOAs, Lasers]=ndgrid( linspacedamplitudes , linspacedprobeamplitudes, logspacedfreqs, duration, linspacedprobeSOAs, [0 1] );
    numlasers=2;
else
    [Amplitudes, ProbeAmps, Freqs, Durations, ProbeSOAs]=ndgrid( linspacedamplitudes , linspacedprobeamplitudes, logspacedfreqs, duration, linspacedprobeSOAs );
    numlasers=1;
    Lasers=zeros(size(Amplitudes));
end
numdurations=1;
numSOAs=length(ProbeSOAs);

kk=0;
reps=0:StimPerRepeat:(StimPerRepeat*nrepeats);
for rep=1:nrepeats
    
%     neworder=randperm( numfreqs * numamplitudes * numprobeamplitudes  * numSOAs * numlasers);
%     amplitudes = Amplitudes( neworder );
%     probeAmps= ProbeAmps( neworder );
%     freqs = Freqs( neworder );
%     lasers = Lasers( neworder );
    stimorder=randperm(StimPerRepeat); %order of single tone and 2Tone, random
    
    nn=0;
    %for nn=1:length(stimorder)
    for a=1:numamplitudes
        for f=1:numfreqs
            for l=1:numlasers
                nn=nn+1;
                kk=stimorder(nn)+reps(rep);
                if Freqs(f)==-1
                    stimuli(kk).type='whitenoise';
                else
                    stimuli(kk).type='tone';
                end
                
                stimuli(kk).param.frequency=Freqs(f);
                stimuli(kk).param.amplitude=Amplitudes(a);
                stimuli(kk).param.duration=duration;
                stimuli(kk).param.ramp=ramp;
                stimuli(kk).param.next=isi;
                stimuli(kk).param.laser=Lasers(l);
                stimuli(kk).stimulus_description=GetParamStr(stimuli(kk));
                stimuli(kk).protocol_name=name;
                stimuli(kk).protocol_description=description;
                stimuli(kk).param.SOA=[];
                stimuli(kk).param.probeamp=[];
                stimuli(kk).param.probefreq=[];
                stimuli(kk).version='djmaus';
                
                
                for aa=1:numprobeamplitudes
                    for s=1:numSOA
                        for ff=1:numfreqs
                            nn=nn+1;
                            kk=stimorder(nn)+reps(rep);
                            stimuli(kk).type='2tone';
                            stimuli(kk).param.frequency=Freqs(f);
                            stimuli(kk).param.amplitude=Amplitudes(a);
                            stimuli(kk).param.probeamp=ProbeAmps(aa);
                            stimuli(kk).param.probefreq=Freqs(ff);
                            stimuli(kk).param.duration=duration;
                            stimuli(kk).param.ramp=ramp;
                            stimuli(kk).param.next=isi;
                            stimuli(kk).param.laser=Lasers(l);
                            stimuli(kk).stimulus_description=GetParamStr(stimuli(kk));
                            stimuli(kk).protocol_name=name;
                            stimuli(kk).protocol_description=description
                            stimuli(kk).param.SOA=ProbeSOAs(s);
                            stimuli(kk).version='djmaus';
                            
                        end
                    end
                end
            end
        end
    end
    
    
    % insert silent sounds
    if include_silent_sound==1
        
        nn=nn+1;
        kk=stimorder(nn)+reps(rep);
        stimuli(kk).param.laser=0;
        stimuli(kk).type='silentsound';
        stimuli(kk).param.duration=duration;
        stimuli(kk).param.ramp=0;
        stimuli(kk).param.next=isi;
        stimuli(kk).stimulus_description=GetParamStr(stimuli(kk));
        stimuli(kk).protocol_name=name;
        stimuli(kk).stimulus_description=description;
        stimuli(kk).version='djmaus';
        
        if interleave_laser==1
            nn=nn+1;
            kk=stimorder(nn)+reps(rep);
            stimuli(kk).param.laser=1;
            stimuli(kk).type='silent sound Laser ON';
            stimuli(kk).param.duration=duration;
            stimuli(kk).param.ramp=0;
            stimuli(kk).param.next=isi;
            stimuli(kk).stimulus_description=GetParamStr(stimuli(kk));
            stimuli(kk).protocol_name=name;
            stimuli(kk).pstimulus_description=description;
            stimuli(kk).version='djmaus';
        end
    end
end


global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('2Tone')
cd ('2Tone')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);


% keyboard