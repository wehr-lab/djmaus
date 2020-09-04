function [filename,path]=MakeSoundfiledjProtocol(amplitude, dur, include_whitenoise, interleave_laser, include_silent_sound, isi, nrepeats)

% Usage: MakeSoundfiledjProtocol(amplitude, dur, include_whitenoise, ...
%        interleave_laser, include_silent_sound, isi, nrepeats)
%
%
%
% loads sounds from sound files (e.g., WAV files). Opens a dialog box for
% you to select the sound files and provide a descriptive name for the protocol. 
%
%
% see help audioread for list of supported file formats. WAV works for
% sure. if stereo only the right channel will be used
%
% INPUTS:
% amplitude: in dB SPL. All sounds will be normalized to peak level across all sounds,
%   then scaled to this amplitude
% duration (in seconds): how much of the soundfile to use
%   (use [] to default to total duration of the soundfile)
% include_whitenoise: 0 or 1 to include white noise bursts (at amplitude)
% interleave_laser: 0 or 1 to duplicate all stimuli and interleave laser
%            and non-laser trials in random order
% include_silent_sound: 0 or 1 to include an extra stimulus that is silent,
%           useful for explicitly collecting spontaneous and/or laser-only
%           trials as a separate condition
% isi: inter stimulus interval (onset-to-onset) in ms
% nrepeats: number of repetitions (different pseudorandom orders)
% OUTPUTS:
%       - creates a suitably named stimulus protocol in djprefs.stimuli\Soundfile Protocols
%       - copies resampled/truncated/amplitude-adjusted sourcefiles to subfolder
%       - returns name & path to protocol (AKH 6/19/13)
% ------------------------------------------------------------------------
%
% example call: 
% amp= 80; dur=[]; include_whitenoise= 1; interleave_laser= 0; include_silent_sound= 1; isi= 800; nrepeats= 20;
%MakeSoundfiledjProtocol(amp, dur, include_whitenoise, interleave_laser, include_silent_sound, isi, nrepeats)



if nargin==0; fprintf('\nno input');return;end
[filename_ext, wavpath] = uigetfile( '*', 'please choose source files','MultiSelect','on');
if isequal(filename_ext,0) || isequal(wavpath,0)
    disp('User pressed cancel')
    return
end
[descriptname] = inputdlg('Please name the protocol');
descriptname = char(descriptname);


cd(wavpath)

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


if ischar(filename_ext) %only one filename
    StimPerRepeat=1;
    fprintf('\n%d sound', StimPerRepeat)
    numsoundfiles=1;
    filename_ext={filename_ext};
elseif iscell(filename_ext)
    StimPerRepeat=length(filename_ext);
    fprintf('\n%d sounds', StimPerRepeat)
    numsoundfiles=StimPerRepeat;
else
    error
end

if isempty(dur)
    durstring='full_duration';
    for i = 1:length(filename_ext)
        fn=filename_ext{i};
        info = audioinfo(fn);
        durations(i)=info.Duration;
    end
else durstring=sprintf('%dms', dur);
    durations=dur*ones(1, numsoundfiles)/1000;
end

if include_silent_sound
    StimPerRepeat=StimPerRepeat+1;
end
if include_whitenoise
    StimPerRepeat=StimPerRepeat+1;
end
if interleave_laser
    StimPerRepeat=StimPerRepeat*2;
end

TotalNumStim=StimPerRepeat*nrepeats;
if isempty(dur) %full duration
    DurationPerRepeatSecs=sum(durations)+StimPerRepeat*(isi)/1000+(StimPerRepeat-numsoundfiles)*durations(1);
else %truncated
    DurationPerRepeatSecs=StimPerRepeat*(dur+isi)/1000;%approx. duration per repeat
end
TotalDurationSecs=DurationPerRepeatSecs*nrepeats;


name= sprintf('soundfile %s %s, %d dB, %s, %s-%s-%dmsisi/%d reps', ...
    descriptname, include_whitenoisestr, amplitude, durstring,interleave_laserstr,include_silent_soundstr, isi,nrepeats);
description= sprintf('soundfile name: %s, %s, %d dB, dur: %s, %s, %s, %dms isi, %d stim per repeat, %d repeats, %d total stimuli, %ds per repeat, %d s total dur',...
    descriptname, include_whitenoisestr, amplitude, durstring,interleave_laserstr,include_silent_soundstr, isi, StimPerRepeat, nrepeats,TotalNumStim, round(DurationPerRepeatSecs), round(TotalDurationSecs));
filename=sprintf('soundfile-%s%s%ddB-%s-%s-%s-isi%dms-%dreps',...
    descriptname, include_whitenoisestr, amplitude, durstring,...
    interleave_laserstr,include_silent_soundstr, isi, nrepeats);

djPrefs
global pref
if isempty(pref) djPrefs; end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('Soundfile Protocols')
cd ('Soundfile Protocols')
sourcepath=[filename, '_sourcefiles'];
mkdir(sourcepath)
cd(sourcepath)
sourcepath_wd=pwd; %full path

amp=amplitude;
n=0;
for nreps=1:nrepeats
    for i = 1:numsoundfiles
        cd(wavpath)
        fn=filename_ext{i};
        n=n+1;
        
        [s, Fs]=audioread(fn);
        s=s(1:durations(i)*Fs); %truncate to requested durations or default to full length
        s=resample(s, pref.SoundFs , Fs); %resample to soundcard samprate
        
        %normalize and set to requested SPL;
        s=s./max(abs(s));
        %amplitude=1*(10.^((amp-pref.maxSPL)/20)); %in volts (-1<x<1), i.e. pref.maxSPL=+_1V
        %s=amplitude.*s;
        
        sample.param.description='soundfile stimulus';
        sample.param.duration=durations(i);
        sample.param.sourcepath=sourcepath;
        sample.param.sourcefilename=fn;
        sample.sample=s;
        
        sourcefilename=sprintf('soundfile_%s_sourcefile_%s_%d_%ddB_%.1fs.mat', descriptname, fn, i, amp, durations(i));
        cd(sourcepath_wd)
        save(sourcefilename, 'sample');
        
        %Make stim structure
        stimuli(n).type='soundfile';
        stimuli(n).param.sourcefile=sourcefilename;
        stimuli(n).param.sourcepath=sourcepath;
        stimuli(n).param.duration=durations(i)*1e3; %in ms
        stimuli(n).param.amplitude=amp; %
        stimuli(n).param.next=isi;
        stimuli(n).param.ramp=0;
        stimuli(n).param.laser=0;
        stimuli(n).stimulus_description=GetParamStr(stimuli(n));
        stimuli(n).protocol_name=name;
        stimuli(n).protocol_description=description;
        stimuli(n).PlottingFunction='PlotSoundfile';
        stimuli(n).version='djmaus';
        
        if interleave_laser
            n=n+1;
            stimuli(n).type='soundfile';
            stimuli(n).param.sourcefile=sourcefilename;
            stimuli(n).param.sourcepath=sourcepath;
            stimuli(n).param.duration=durations(i)*1e3; %in ms
            stimuli(n).param.amplitude=amp; %
            stimuli(n).param.next=isi;
            stimuli(n).param.ramp=0;
            stimuli(n).param.laser=1;
            stimuli(n).stimulus_description=GetParamStr(stimuli(n));
            stimuli(n).protocol_name=name;
            stimuli(n).protocol_description=description;
            stimuli(n).PlottingFunction='PlotSoundfile';
            stimuli(n).version='djmaus';
        end
        
    end
    % insert silent sounds
    if include_silent_sound==1
        n=n+1;
        stimuli(n).param.laser=0;
        stimuli(n).type='silentsound';
        stimuli(n).param.duration=durations(1)*1e3; %in ms;
        stimuli(n).param.ramp=0;
        stimuli(n).param.next=isi;
        stimuli(n).stimulus_description=GetParamStr(stimuli(n));
        stimuli(n).protocol_name=name;
        stimuli(n).protocol_description=description;
        stimuli(n).PlottingFunction='PlotSoundfile';
        stimuli(n).version='djmaus';
        
        if interleave_laser==1
            n=n+1;
            stimuli(n).param.laser=1;
            stimuli(n).type='silentsound';
            stimuli(n).param.duration=durations(1)*1e3; %in ms;
            stimuli(n).param.ramp=0;
            stimuli(n).param.next=isi;
            stimuli(n).stimulus_description=GetParamStr(stimuli(n));
            stimuli(n).protocol_name=name;
            stimuli(n).protocol_description=description;
            stimuli(n).PlottingFunction='PlotSoundfile';
            stimuli(n).version='djmaus';
        end
    end
    % insert white noise
    if include_whitenoise==1
        n=n+1;
        stimuli(n).param.laser=0;
        stimuli(n).type='whitenoise';
        stimuli(n).param.duration=durations(1)*1e3; %in ms;
        stimuli(n).param.amplitude=amp;
        stimuli(n).param.ramp=0;
        stimuli(n).param.next=isi;
        stimuli(n).stimulus_description=GetParamStr(stimuli(n));
        stimuli(n).protocol_name=name;
        stimuli(n).protocol_description=description;
        stimuli(n).PlottingFunction='PlotSoundfile';
        stimuli(n).version='djmaus';
        
        if interleave_laser==1
            n=n+1;
            stimuli(n).param.laser=1;
            stimuli(n).type='whitenoise';
            stimuli(n).param.duration=durations(1)*1e3; %in ms;
            stimuli(n).param.amplitude=amp;
            stimuli(n).param.ramp=0;
            stimuli(n).param.next=isi;
            stimuli(n).stimulus_description=GetParamStr(stimuli(n));
            stimuli(n).protocol_name=name;
            stimuli(n).protocol_description=description;
            stimuli(n).PlottingFunction='PlotSoundfile';
            stimuli(n).version='djmaus';
        end
    end
end

%shuffle stimuli
if n~=TotalNumStim
    error
end

stimorder=randperm(TotalNumStim); %order of single tone and 2Tone, random
for n=1:TotalNumStim
    shuffledstimuli(n)=stimuli(stimorder(n));
end
stimuli=shuffledstimuli;



cd('E:\Stimuli')
cd ('Soundfile Protocols')

path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);

