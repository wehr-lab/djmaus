function [ events ] = GetNetworkEvents(filename )
%reads open ephys network events messages from messages.events file

if ~strcmp(filename, 'messages.events')
    error('this file should only used for messages.events');
end

fid = fopen(filename);
%filesize = getfilesize(fid);

event_index=0;
while 1
    tline=fgetl(fid);
    if ~ischar(tline), break, end
    %fprintf('\n%s', tline);
    event_index=event_index+1;
    events{event_index}=tline;
end

fclose(fid);

