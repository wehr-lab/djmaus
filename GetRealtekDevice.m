function deviceID = GetRealtekDevice
%looks in soundcard device list and 
%returns the deviceID of the 'Speakers (Realtek High Definition Audio)'
% (because sometimes they can change deviceIDs)
%modified from PrintDevices
%mw 10-11-2016 // KK 10-23-2020

%AssertOpenGL;
% InitializePsychSound(1);
%this is throwing a message every stimulus but I
% don't think we need it. mw 08.23.2017
%(we use GetSoundMapperDevice to determine which soundcard is used and
%therefore how to set amplitude of soundcard trigger)

devs = PsychPortAudio('GetDevices');
deviceID=[];
for n = 1:length(devs) 
    if strcmp(devs(n).DeviceName, 'Speakers (Realtek High Definition Audio)')
        deviceID=devs(n).DeviceIndex;
    end
end