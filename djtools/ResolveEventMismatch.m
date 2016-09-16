function [ Events, all_SCTs, stimlog ] = ResolveEventMismatch(Events, all_SCTs, stimlog  )
%try to resolve the situation when the number of sound events (from network
%messages) does not match Number of hardware triggers (soundcardtrig TTLs

nstimlog=length(stimlog);
nEvents=length(Events);
nSCTs=length(all_SCTs);
fprintf('\nNumber of stimuli logged to notebook (stimlog): %d', nstimlog);
fprintf('\nNumber of sound events (from network messages): %d', nEvents);
fprintf('\nNumber of hardware triggers (soundcardtrig TTLs): %d', nSCTs);

if nstimlog~=nEvents
    error('ResolveEventMismatch: number of logged stim does not match number of Events')
end

%first double-check that the stimlog and Events match
if nstimlog==nEvents
    for i=1:nstimlog
        if ~strcmp(stimlog(i).type, Events(i).type)
            error('ResolveEventMismatch: stimlog does not match Events')
        end
        if ~stimlog(i).param.duration==Events(i).duration
            error('ResolveEventMismatch: stimlog does not match Events')
        end
        if strcmp(stimlog(i).type,'tone')
            if ~stimlog(i).param.frequency==Events(i).frequency
                error('ResolveEventMismatch: stimlog does not match Events')
            end
            if ~stimlog(i).param.amplitude==Events(i).amplitude
                error('ResolveEventMismatch: stimlog does not match Events')
            end
        end
    end
end
fprintf('\nGood. Logged stimuli match network Events.')
    
if nstimlog==nEvents & nSCTs==1+ nEvents
    %there's one extra soundcard trigger, find it and remove it
    
    %if the extra one comes before any events, don't need to do anything
    if all_SCTs(1) < Events(1).message_timestamp_sec
        %we're fine, because extra SCTs before any message is automatically ignored
        fprintf('\nMismatched resolved. Extra soundcard trigger discarded') 
    else
        error('ResolveEventMismatch: this case is not handled yet')
    end
else
    error('ResolveEventMismatch: this case is not handled yet')
end

