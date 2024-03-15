function [event_times_sec, num_events, ProtocolStart_secs, Stimtrace, StimFs] = GetSCTriggersFromSpike2file(filename)
%reads soundcard triggers and protocol start events and stimulus monitor from a 
% Spike2 .mat file exported from Spike 2
%
%
%usage: 
%   [event_times_sec,num_events, ProtocolStart_secs, Stimtrace, StimFs] = GetSCTriggersFromSpike2file(filename)

SCTchan=18;
ProtocolStartchan=19;
Soundchan=11;

load(filename ); 
%each channel is a separate struct variable in the workspace loaded from filename

if contains(filename, '1414_1') %special case from old session without native spike2-mat-conversion
return
end

%discover the name of the spike2 variable Data n name which changes every time
% SoundCardTrigger
SCT_structname=whos(sprintf('Data*__Stopped__Ch%d', SCTchan))
SCT_struct=eval(SCT_structname.name);
if ~contains(SCT_struct.title, 'SoundTri')
    error ('SCT trigger channel not found')
end

event_times_sec=SCT_struct.times;
num_events=SCT_struct.length;

% Protocol start
ProtocolStart_structname=whos(sprintf('Data*__Stopped__Ch%d', ProtocolStartchan));
ProtocolStart_struct=eval(ProtocolStart_structname.name);
if ~contains(ProtocolStart_struct.title, 'ProtStar')
    error ('Protocol Start  channel not found')
end

ProtocolStart_secs=ProtocolStart_struct.times;
num_ProtocolStarts=ProtocolStart_struct.length;
if num_ProtocolStarts > 1
    error ('more than one protocol start? that''s not supposed to happen')
end

Sound_structname=whos(sprintf('Data*__Stopped__Ch%d', Soundchan));
Sound_struct=eval(Sound_structname.name);
if ~contains(Sound_struct.title, 'sound')
    error ('Sound  channel not found')
end
Stimtrace=Sound_struct.values;
Stimtrace=Stimtrace./max(abs(Stimtrace));
StimFs=1/Sound_struct.interval;

%     sound monitor?
% 
% Data5__Stopped__Ch11 = 
% 
%   struct with fields:
% 
%        title: 'sound'
%      comment: ''
%     interval: 2.0000e-05
%        scale: 1.5259e-04
%       offset: 0
%        units: ''
%        start: 3.2614
%       length: 12135653
%       values: [12135653×1 double]
%       but this channel appears to have no signal