function [Events, StartAcquisitionSec, StartAcquisitionSamples] = ProcessEvents(varargin)
    
    try
        datadir = varargin{1};
    catch
        error('Could not parse directory, check to see if the proper dir name is an input')
    end
    cd(datadir)

    messagesfilename='messages.events';
    [messages] = GetNetworkEvents(messagesfilename);
    Eventsfilename='all_channels.events';
    [all_channels_data, all_channels_timestamps, all_channels_info] = load_open_ephys_data(Eventsfilename);
    load('notebook.mat');
    sampleRate=all_channels_info.header.sampleRate;
    [Events, StartAcquisitionSec] = GetEventsAndSCT_Timestamps(messages, sampleRate, all_channels_timestamps, all_channels_data, all_channels_info, stimlog);
    
   
end