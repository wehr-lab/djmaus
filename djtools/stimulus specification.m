% this is a specification for what djmaus stimulus protocols should be like
% mw 2.24.2017
%
% a stimulus protocol is nothing more than a structure array, with each
% element corresponding to a stimulus. djmaus simply plays them in order.
%
% the name of the structure should be "stimuli" 
% the fields are:
% type 
% param
% stimulus_description 
% protocol_name 
% protocol_description 
% PlottingFunction
% version

% type is used to determine which MakeStimuli function is used to generate
% the stimulus 
% for example, MakeWhitenoise is used for type='whitenoise'

% param is a substructure that contains fields for all the params. 
% param contains some general parameters that all stimuli should have, and some
% specific params for certain stimuli 
% all stimuli should have these params:
%   duration (duration of the stimulus in ms) 
%   next (inter-stimulus interval % in ms, time to wait before loading next stimulus) 
%   laser: 0 or 1, whether the laser is on or off on that trial 
% specific params could be things like 
%   frequency
%   amplitude 
%   ramp 
%   seamless (used when there should be no interval between stimuli, e.g. for continuous background noise) 
%   var_laser params can go here to override the GUI settings in djmaus
 
% stimulus_description should be automatically generated as follows:
% stimuli(n).stimulus_description=GetParamStr(stimuli(n)); 
% this is because this field contains a list of ordered pairs of fields and values,
% which are sent directly out by djmaus as a network message. This is read
% by open ephys and stored with the data. We use this to sort data into
% matrices for tuning curves, psths, etc. so it needs to be in a fixed format. 
% the stimulus_description is also displayed in djmaus as the stimulus is presented 
% Important: as you are creating each stimulus, call 
% stimuli(n).stimulus_description=GetParamStr(stimuli(n)); 
% after creating all the params but before you create the name, description, plottingfunction, and version,
% so that the latter fields don't get inserted into the stimulus_description
% (see below example)

% protocol_name - the name that will be displayed in
% the name field in djmaus. It's helpful (but not necessary) to have it
% match the filename
 
% protocol_description - the description displayed in the description
% field in djmaus. It's helpful (but not necessary) to include information
% like the total duration of the stimulus, the duration per repeat.

% PlottingFunction - which plotting function might reasonably be used to
% plot the data. This is used for example by JustPlotIt and
% GetPlottingFunction version - should be 'djmaus'
 
% Here is a very simple example of code to create a stimulus protocol
%   for more complex examples, such as pseudorandom interleaving, or multiple
%   stimulus types, see Make ... djProtocol files in djtools
% 
% amplitude=80;
% dur=25;
% nrepeats=10;
% isi=500;
% name= sprintf('my protocol, %ddB%.1fms/%sms/d=%d/n=%d', amplitude, dur, nrepeats);
% description=sprintf('just white noise, %d dB, %d repeats', amplitude, nrepeats);
% filename=sprintf('my_protocol-%ddB-d%d', amplitude, dur);
% 
% n=0;
% for rep=1:nrepeats
%         n=n+1;
%             stimuli(n).type='whitenoise'; 
%             stimuli(n).param.amplitude=amplitude;
%             stimuli(n).param.duration=dur;
%             stimuli(n).param.laser=0;
%             stimuli(n).param.ramp=1;
%             stimuli(n).param.next=isi;
%             stimuli(n).stimulus_description=GetParamStr(stimuli(n));
%             stimuli(n).protocol_name=name;
%             stimuli(n).protocol_description=description;
%             stimuli(n).PlottingFunction='PlotTC_PSTH';
%             stimuli(n).version='djmaus';
%     end
% 
% 
% global pref
% if isempty(pref) djPrefs;end
% cd(pref.stimuli)
% warning off MATLAB:MKDIR:DirectoryExists
% mkdir('Tuning Curve protocols')
% cd ('Tuning Curve protocols')
% save(filename, 'stimuli')
% fprintf('\ncreated file %s \nin directory %s\n\n', filename, pwd);
% 
