function [filename,path]=MakeTonedjProtocol(freqsperoctave, minfreq, maxfreq, numamplitudes, ...
    minamplitude, maxamplitude, durations, ramp, include_whitenoise, interleave_laser, include_silent_sound, isi, nrepeats)





% Usage: MakeTonedjProtocol(freqsperoctave, minfreq, maxfreq, ...
% numamplitudes, minamplitude, maxamplitude, duration, ramp, ...
% include_whitenoise, interleave_laser, isi, nrepeats)
%
% same as MakeTuningCurveOct from exper, except that you can tell it
% interleave laser trials
%
% INPUTS:
% freqsperoctave: number of frequencies per octave (frequency resolution)
% minfreq: lowest frequency in Hz
% maxfreq: highest frequency in Hz
%           note: if minfreq-to-maxfreq cannot be divided evenly into
%           freqsperoctave, the nearest maxfreq will be used
%           (i.e. the requested freqsperoctave will be exactly enforced)
% numamplitudes: number of amplitude steps
% minamplitude: maximum amplitude in dB SPL (requires system to be calibrated)
% maxamplitude: maximum amplitude in dB SPL (requires system to be calibrated)
% durations: vector of different tone durations (in ms) (can be a single duration)
% ramp: on-off ramp duration in ms
% include_whitenoise: 0 or 1 to include white noise bursts at each amplitude
% interleave_laser: 0 or 1 to duplicate all stimuli and interleave laser
%            and non-laser trials in random order
% include_silent_sound: 0 or 1 to include an extra stimulus that is silent,
%           useful for explicitly collecting spontaneous and/or laser-only
%           trials as a separate condition
% isi: inter stimulus interval (onset-to-onset) in ms
% nrepeats: number of repetitions (different pseudorandom orders)
% OUTPUTS:
%       - creates a suitably named stimulus protocol in djprefs.stimuli\Tone Protocols
%       - returns name & path to protocol (AKH 6/19/13)
% ------------------------------------------------------------------------
%
% example call: MakeTonedjProtocol(4, 1000, 32000, 3, 50, 80, 200, 10, 1, 1, 1, 500, 10)
%
% example call with multiple durations:
% MakeTonedjProtocol(4, 1000, 32000, 3, 50, 80, [200 400],10,1,1, 1, 500, 10)
%
% MakeTonedjProtocol(0, 0, 0, 1, 70, 70, [1 2 4 8 16 32 64 128 256],0,1,0, 1, 1000, 10)

if nargin==0; fprintf('\nno input');return;end

numdurations=length(durations);


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
if numfreqs==0; logspacedfreqs=[]; end

if include_whitenoise==1
    logspacedfreqs=[logspacedfreqs -1]; %add whitenoise as extra freq=-1
    numfreqs=numfreqs+1;
end
if interleave_laser==1
    [Amplitudes,Freqs, Durations, Lasers]=ndgrid( linspacedamplitudes , logspacedfreqs, durations, [0 1] );
    numlasers=2;
else
    [Amplitudes,Freqs, Durations]=meshgrid( linspacedamplitudes , logspacedfreqs, durations );
    numlasers=1;
    Lasers=zeros(size(Amplitudes));
end

neworder=randperm( numfreqs * numamplitudes * numdurations * numlasers);

amplitudes=zeros(size(neworder));
freqs=zeros(size(neworder));
durs=zeros(size(neworder));
lasers=zeros(size(neworder));

StimPerRepeat=length(neworder);
if include_silent_sound
    StimPerRepeat=StimPerRepeat+numlasers;
end
TotalNumStim=StimPerRepeat*nrepeats;
DurationPerRepeatSecs=StimPerRepeat*(mean(durations)+isi)/1000;%approx. duration per repeat
TotalDurationSecs=DurationPerRepeatSecs*nrepeats;




durstring=sprintf('%d-', durations);durstring=durstring(1:end-1);
%put into stimuli structure
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
        name= sprintf('Tuning curve WN only, %da(%d-%ddB)/%dd(%sms)/%s-%s-%dmsisi/%d reps', ...
            numamplitudes,minamplitude,...
            maxamplitude, numdurations, durstring,interleave_laserstr,include_silent_soundstr, isi,nrepeats);
        description=sprintf('tuning curve, WN only, %d ampl. (%d-%d dB SPL), %d durations (%sms),%s, %s, %d ms isi,%d stim per rep, %d repeats, %d total stimuli, %ds per rep, %d s total dur',...
            numamplitudes,minamplitude, maxamplitude, numdurations, durstring, ...
            interleave_laserstr,include_silent_soundstr, isi, StimPerRepeat, nrepeats,TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs));
        filename=sprintf('tuning-curve-WNonly-%da_%d-%ddB-%dd_%sms-%s-%s-isi%dms-%dreps',...
             numamplitudes,minamplitude, maxamplitude, numdurations, durstring,...
            interleave_laserstr,include_silent_soundstr, isi, nrepeats);
    else
        name= sprintf('Tuning curve %s, %dfpo(%d-%dHz)/%da(%d-%ddB)/%dd(%sms)/%s-%s-%dmsisi/%d reps', ...
            include_whitenoisestr, freqsperoctave,minfreq, round(maxfreq), numamplitudes,minamplitude,...
            maxamplitude, numdurations, durstring,interleave_laserstr,include_silent_soundstr, isi,nrepeats);
        description=sprintf('tuning curve, tones%s, %d freqs/oct (%d-%dkHz), %d ampl. (%d-%d dB SPL), %d durations (%sms),%s, %s, %d ms isi,%d stim per rep, %d repeats, %d total stimuli, %ds per rep, %d s total dur',...
            include_whitenoisestr,freqsperoctave, minfreq, round(maxfreq), numamplitudes,minamplitude, maxamplitude, numdurations, durstring, ...
            interleave_laserstr,include_silent_soundstr, isi, StimPerRepeat, nrepeats,TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs));
        filename=sprintf('tuning-curve-tones%s-%dfpo_%d-%dHz-%da_%d-%ddB-%dd_%sms-%s-%s-isi%dms-%dreps',...
            include_whitenoisestr, freqsperoctave,minfreq, round(maxfreq), numamplitudes,minamplitude, maxamplitude, numdurations, durstring,...
            interleave_laserstr,include_silent_soundstr, isi, nrepeats);
    end

nn=0;
for rep=1:nrepeats
    
    neworder=randperm( numfreqs * numamplitudes * numdurations * numlasers);
    amplitudes = Amplitudes( neworder );
    freqs = Freqs( neworder );
    durs = Durations( neworder );
    lasers = Lasers( neworder );
    
    for n=1:length(amplitudes)
        nn=nn+1;
        if freqs(n)==-1
            stimuli(nn).type='whitenoise'; %use nn because stimuli(1) is name/description
            stimuli(nn).param.amplitude=amplitudes(n);
            stimuli(nn).param.duration=durs(n);
            stimuli(nn).param.laser=lasers(n);
            stimuli(nn).param.ramp=ramp;
            stimuli(nn).param.next=isi;
            stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
            stimuli(nn).protocol_name=name;
            stimuli(nn).protocol_description=description;
            stimuli(nn).version='djmaus';
            
        else
            stimuli(nn).type='tone';
            stimuli(nn).param.frequency=freqs(n);
            stimuli(nn).param.amplitude=amplitudes(n);
            stimuli(nn).param.duration=durs(n);
            stimuli(nn).param.laser=lasers(n);
            stimuli(nn).param.ramp=ramp;
            stimuli(nn).param.next=isi;
            stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
            stimuli(nn).protocol_name=name;
            stimuli(nn).protocol_description=description;
            stimuli(nn).version='djmaus';
            
        end
        
        
    end
    % insert silent sounds
    if include_silent_sound==1
        nn=nn+1;
        stimuli(nn).param.laser=0;
        stimuli(nn).type='silentsound';
        stimuli(nn).param.duration=durs(1);
        stimuli(nn).param.ramp=0;
        stimuli(nn).param.next=isi;
        stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
        stimuli(nn).protocol_name=name;
        stimuli(nn).protocol_description=description;
        stimuli(nn).version='djmaus';
        
        if interleave_laser==1
            nn=nn+1;
            stimuli(nn).param.laser=1;
            paramstring='silent sound Laser ON';
            stimuli(nn).type='silentsound';
            stimuli(nn).param.duration=durs(1);
            stimuli(nn).param.ramp=0;
            stimuli(nn).param.next=isi;
            stimuli(nn).stimulus_description=GetParamStr(stimuli(nn));
            stimuli(nn).protocol_name=name;
            stimuli(nn).protocol_description=description;
            stimuli(nn).version='djmaus';
        end
    end
end



global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('Tuning Curve protocols')
cd ('Tuning Curve protocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);

