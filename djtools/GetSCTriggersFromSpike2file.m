function GetSCTriggersFromSpike2file(filename)
%reads soundcard triggers and protocol start events from Spike2 smrx file
%and writes to a .mat file with the same name
%WINDOWS ONLY (thanks Spike2)
%usage: 
%    GetSCTriggersFromSpike2file(filename)
if ~ispc error('sorry, Windows only'); end

SCTchan=18;
ProtocolStartchan=19;
cedpath ='C:\CEDMATLAB\CEDS64ML';
addpath( cedpath );
CEDS64LoadLib( cedpath );
fid = CEDS64Open(filename );
[ num_events, event_times_ticks ] = CEDS64ReadEvents( fid, SCTchan, 1e6, 0 );
event_times_sec=CEDS64TicksToSecs(fid, event_times_ticks);

[ num_ProtocolStarts, ProtocolStart_ticks ] = CEDS64ReadEvents( fid, ProtocolStartchan, 1e6, 0 );
event_times_sec=CEDS64TicksToSecs(fid, event_times_ticks);
ProtocolStart_secs=CEDS64TicksToSecs(fid, ProtocolStart_ticks);

outfilename=replace(filename, '.smrx', '.mat');
save (outfilename)
fprintf('\nwrote events file %s\n', outfilename)

    
