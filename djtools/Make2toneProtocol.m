function [filename,path]=Make2toneProtocol(masker_freqsperoctave, masker_minfreq, masker_maxfreq, masker_numamps, masker_minamp, masker_maxamp, masker_durs, ...
    probe_freqsperoctave, probe_minfreq, probe_maxfreq, probe_numamps, probe_minamp, probe_maxamp, probe_durs, SOAs, ramp, iti, laser, include_WN, include_silentsound, nrepeats)
%for djmaus
%use Make2toneProtocol(masker_freqsperoctave, masker_minfreq, masker_maxfreq, masker_numamps, masker_minamp, masker_maxamp, masker_durs ...
%probe_freqsperoctave, probe_minfreq, probe_maxfreq, probe_numamps, probe_minamp, probe_maxamp, probe_durs, SOAs, ramp, itis, laser, incluse_WN, include_silent_sound, nrepeats)
%SOAs (ms) from onset of masker to onset of probe
%itis (ms) from onset of masker to onset of the second masker (i think),
%freqs are in Hz, amps are in dB
%if laser=1 will interleave all stimuli including silentsound
%silentsound will be an average of masker + probe + SOAs duration
%creates probe alone stimuli as well
%recommend to only vary one parameter (either dur or freq or amp) varying
%more than one will make the protocol very long
% Make2toneProtocol(1, 200, 12000, 1, 80, 80, [10 25 50], 1, 6000, 6000, 1,
% 70, 70, 30, [50 80 100], 5, 1000, 1, 1, 1, 30) = 142 min
%ira 12.06.17
dbstop if error
if nargin==0; fprintf('\nno input');return;end

%durations for masker and probe
masker_numdurs=length(masker_durs);
probe_numdurs=length(probe_durs);


%calculating masker freqs
if masker_freqsperoctave==0 && makser_minfreq==0 && masker_maxfreq==0 && include_WN==1
    %white noise only
    masker_WNonly=1;
    masker_logspacedfreqs=[];
    masker_numfreqs=0;
    fprintf('\nwhite noise only for masker')
else
    masker_WNonly=0;
    masker_numoctaves=log2(masker_maxfreq/masker_minfreq);
    masker_logspacedfreqs=masker_minfreq*2.^([0:(1/masker_freqsperoctave):masker_numoctaves]);
    masker_newmaxfreq=masker_logspacedfreqs(end);
    masker_numfreqs=length(masker_logspacedfreqs);
    if masker_maxfreq~=masker_newmaxfreq
        fprintf('\nnote: could not divide %d-%d Hz evenly into exactly %d frequencies per octave', masker_minfreq, masker_maxfreq, masker_freqsperoctave)
        fprintf('\nusing new maxfreq of %d to achieve exactly %d frequencies per octave\n', round(masker_newmaxfreq), masker_freqsperoctave)
        masker_maxfreq=masker_newmaxfreq;
    end
end

%calculating probe freqs
if probe_freqsperoctave==0 && probe_minfreq==0 && probe_maxfreq==0 && include_WN==1
    %white noise only
    probe_WNonly=1;
    probe_logspacedfreqs=[];
    probe_numfreqs=0;
    fprintf('\nwhite noise only for masker')
else
    probe_WNonly=0;
    probe_numoctaves=log2(probe_maxfreq/probe_minfreq);
    probe_logspacedfreqs=probe_minfreq*2.^([0:(1/probe_freqsperoctave):probe_numoctaves]);
    probe_newmaxfreq=probe_logspacedfreqs(end);
    probe_numfreqs=length(probe_logspacedfreqs);
    if probe_maxfreq~=probe_newmaxfreq
        fprintf('\nnote: could not divide %d-%d Hz evenly into exactly %d frequencies per octave',probe_minfreq, probe_maxfreq, probe_freqsperoctave)
        fprintf('\nusing new maxfreq of %d to achieve exactly %d frequencies per octave\n', round(probe_newmaxfreq), probe_freqsperoctave)
        probe_maxfreq=probe_newmaxfreq;
    end
end

%calculate amplitudes for masker and for probe
masker_linspacedamplitudes = linspace( masker_minamp , masker_maxamp , masker_numamps );
if masker_numfreqs==0; masker_logspacedfreqs=[]; end

probe_linspacedamplitudes = linspace( probe_minamp , probe_maxamp , probe_numamps );
if probe_numfreqs==0; probe_logspacedfreqs=[]; end

if include_WN==1
    masker_logspacedfreqs=[masker_logspacedfreqs -1]; %add whitenoise as extra freq=-1
    masker_numfreqs=masker_numfreqs+1;
    probe_logspacedfreqs=[probe_logspacedfreqs -1]; %add whitenoise as extra freq=-1
    probe_numfreqs=probe_numfreqs+1;
end

num_SOAs=length(SOAs);
masker_numamps=length(masker_linspacedamplitudes);
probe_numamps=length(probe_linspacedamplitudes);

if laser==1
    numlaser=2;
    lasers=[0 1];
elseif laser==0
    numlaser=1;
    lasers=0;
end

if include_silentsound==1
    num_silentsound=1;
elseif include_silentsound==0
    num_silentsound=0;
end
num_probealone=probe_numamps*probe_numdurs*probe_numfreqs;
%randomize
StimPerRepeat=masker_numfreqs*probe_numfreqs*masker_numamps*probe_numamps*masker_numdurs*probe_numdurs*num_SOAs*numlaser+num_silentsound*numlaser+num_probealone*numlaser; % number of stimuli in each repetition

DurationPerRepeatSecs=StimPerRepeat*(mean(masker_durs)+mean(SOAs)+mean(probe_durs)+iti)/1000;%approx. duration per repeat
TotalDurationMins=(DurationPerRepeatSecs*nrepeats)/60;

%put into stimuli structure
if laser
    interleave_laserstr='IL';
    
else
    interleave_laserstr='';
    
end
if include_silentsound
    include_silent_soundstr='ISS';
else
    include_silent_soundstr='';
end
if include_WN
    include_whitenoisestr='+WN';
else
    include_whitenoisestr='';
end
masker_durstring=sprintf('%d-', masker_durs); masker_durstring=masker_durstring(1:end-1);
probe_durstring=sprintf('%d-', probe_durs); probe_durstring=probe_durstring(1:end-1);
SOA_durstring=sprintf('%d-', SOAs); SOA_durstring=SOA_durstring(1:end-1);

%create filenames and descriptions
name= sprintf('2tone %s, %dmfpo(%d-%dHz) %dpfpo(%d-%dHz)/%dma(%d-%ddB) %dpa(%d-%ddB)/%dmd(%sms) %dpd(%sms)/ %dSOAs(%sms)/ %dmsiti/ %s- %s, %d reps', ...
    include_whitenoisestr, masker_freqsperoctave, masker_minfreq, round(masker_maxfreq), probe_freqsperoctave, probe_minfreq, round(probe_maxfreq), masker_numamps, masker_minamp,...
    masker_maxamp, probe_numamps, probe_minamp, probe_maxamp, masker_numdurs,  masker_durstring, probe_numdurs, probe_durstring,num_SOAs, SOA_durstring, iti, interleave_laserstr, include_silent_soundstr,nrepeats);
description=sprintf('2tone, %s, %d freqs/oct masker (%d-%dkHz), %d freqs/oct probe (%d-%dkHz), %d ampl. masker (%d-%d dB SPL), %d ampl. probe (%d-%d dB SPL), %d durations masker (%sms), %d durations probe (%sms), %s, %s, %d ms soa (%sms), intertrial interval %d, %d stim per rep, %d repeats, %d total stimuli, %ds per rep, %d min total dur',...
    include_whitenoisestr, masker_freqsperoctave, masker_minfreq, round(masker_maxfreq), probe_freqsperoctave, probe_minfreq, round(probe_maxfreq), masker_numamps, masker_minamp,...
    masker_maxamp, probe_numamps, probe_minamp, probe_maxamp, masker_numdurs,  masker_durstring, probe_numdurs, probe_durstring, ...
    interleave_laserstr,include_silent_soundstr, num_SOAs, SOA_durstring, iti, StimPerRepeat, nrepeats, StimPerRepeat*nrepeats, round(DurationPerRepeatSecs), round(TotalDurationMins));
filename=sprintf('2tone-%s-%dmfpo_%d-%dHz_%dpfpo_%d-%dHz-%dma_%d-%ddB_%dpa_%d-%ddB_%dmd_%sms_%dpd_%sms_%dSOAs_%sms_%dmsiti_%s-%s_%d reps',...
    include_whitenoisestr, masker_freqsperoctave, masker_minfreq, round(masker_maxfreq), probe_freqsperoctave, probe_minfreq, round(probe_maxfreq), masker_numamps, masker_minamp,...
    masker_maxamp, probe_numamps, probe_minamp, probe_maxamp, masker_numdurs,  masker_durstring, probe_numdurs, probe_durstring,num_SOAs, SOA_durstring, iti, interleave_laserstr, include_silent_soundstr,nrepeats);

neworder=randperm(StimPerRepeat*nrepeats); %randomize sequence
j=0;
for rep=1:nrepeats
    for l=1:numlaser
        for pf=1:probe_numfreqs
            for pd=1:probe_numdurs
                for pa=1:probe_numamps
                    for i=1:num_SOAs
                        for mf=1:masker_numfreqs
                            for ma=1:masker_numamps
                                for md=1:masker_numdurs
                                    j=j+1;
                                    k=neworder(j);
                                    stimuli(k).type='2tone';
                                    stimuli(k).param.amplitude=masker_linspacedamplitudes(ma);
                                    stimuli(k).param.frequency=masker_logspacedfreqs(mf);
                                    stimuli(k).param.probefreq=probe_logspacedfreqs(pf);
                                    stimuli(k).param.probeamp=probe_linspacedamplitudes(pa);
                                    stimuli(k).param.duration=masker_durs(md);
                                    stimuli(k).param.probe_duration= probe_durs(pd);
                                    stimuli(k).param.laser=lasers(l);
                                    stimuli(k).param.ramp=ramp;
                                    stimuli(k).param.next=iti;
                                    stimuli(k).param.SOA=SOAs(i);
                                    stimuli(k).stimulus_description=sprintf('2tone masker frequency: %d probe frequency:%d masker amplitude: %d probe amplitude: %d masker duration: %d probe duration: %d isi %d laser %d ramp %d next %d', masker_logspacedfreqs(mf), probe_logspacedfreqs(pf), ...
                                        masker_linspacedamplitudes(ma), probe_linspacedamplitudes(pa), masker_durs(md), probe_durs(pd), SOAs(i), lasers(l), ramp, iti);
                                    stimuli(k).protocol_name=name;
                                    stimuli(k).protocol_description=description;
                                    stimuli(k).PlottingFunction='Plot2tone_PSTH';
                                    stimuli(k).version='djmaus';
                                end
                            end
                        end
                    end
                    if probe_logspacedfreqs(pf)==-1
                        j=j+1;
                        k=neworder(j);
                        stimuli(k).type='whitenoise';
                        stimuli(k).param.frequency=probe_logspacedfreqs(pf);
                        stimuli(k).param.amplitude=probe_linspacedamplitudes(pa);
                        stimuli(k).param.duration= probe_durs(pd);
                        stimuli(k).param.laser=lasers(l);
                        stimuli(k).param.ramp=ramp;
                        stimuli(k).param.next=iti;
                        stimuli(k).stimulus_description=sprintf('WN frequency:%d amplitude: %d duration: %d laser %d ramp %d next %d',  probe_logspacedfreqs(pf), ...
                            probe_linspacedamplitudes(pa),  probe_durs(pd), lasers(l), ramp, iti);
                        stimuli(k).protocol_name=name;
                        stimuli(k).protocol_description=description;
                        stimuli(k).PlottingFunction='Plot2tone_PSTH';
                        stimuli(k).version='djmaus';
                        
                    else
                        j=j+1;
                        k=neworder(j);
                        stimuli(k).type='tone';
                        stimuli(k).param.frequency=probe_logspacedfreqs(pf);
                        stimuli(k).param.amplitude=probe_linspacedamplitudes(pa);
                        stimuli(k).param.duration= probe_durs(pd);
                        stimuli(k).param.laser=lasers(l);
                        stimuli(k).param.ramp=ramp;
                        stimuli(k).param.next=iti;
                        stimuli(k).stimulus_description=sprintf('tone frequency:%d amplitude: %d duration: %d laser %d ramp %d next %d',  probe_logspacedfreqs(pf), ...
                            probe_linspacedamplitudes(pa),  probe_durs(pd), lasers(l), ramp, iti);
                        stimuli(k).protocol_name=name;
                        stimuli(k).protocol_description=description;
                        stimuli(k).PlottingFunction='Plot2tone_PSTH';
                        stimuli(k).version='djmaus';
                    end
                    
                end
            end
        end
        j=j+1;
        k=neworder(j);
        stimuli(k).type='silentsound';
        stimuli(k).param.laser=lasers(l);
        stimuli(k).param.ramp=0;
        stimuli(k).param.duration=mean(masker_durs)+mean(probe_durs)+mean(SOAs); %make an average duration silent sound
        stimuli(k).param.next=iti;
        stimuli(k).param.laser=lasers(l);
        %        stimuli(k).param.VarLaser=0;
        %         stimuli(k).param.VarLaserstart=0;
        %         stimuli(k).param.VarLaserpulsewidth=[];
        %         stimuli(k).param.VarLasernumpulses=[];
        %         stimuli(k).param.VarLaserisi=[];
        stimuli(k).stimulus_description=sprintf('2tone silentsound %d duration',stimuli(k).param.duration);
        stimuli(k).protocol_name=name;
        stimuli(k).protocol_description=description;
        stimuli(k).PlottingFunction='Plot2tone_PSTH';
        stimuli(k).version='djmaus';
    end
    
end

global pref
if isempty(pref) djPrefs;end
cd(pref.stimuli)
warning off MATLAB:MKDIR:DirectoryExists
mkdir('2tone protocols')
cd ('2tone protocols')
path=pwd;
save(filename, 'stimuli')
fprintf('\ncreated file %s \nin directory %s\n\n', filename, path);



