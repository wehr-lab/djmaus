function deviceID=GetXonarDevice
%looks in soundcard device list and 
%returns the deviceID of the Xonar STX II+H6: Multichannel adaptor
% (because sometimes they can change deviceIDs)
%modified from PrintDevices
%mw 10-11-2016

%AssertOpenGL;
% InitializePsychSound(1);
%this is throwing a message every stimulus but I
% don't think we need it. mw 08.23.2017
%(we use GetXonarDevice to determine which soundcard is used and
%therefore how to set amplitude of soundcard trigger)

devs = PsychPortAudio('GetDevices');
deviceID=[];
for n = 1:length(devs) 
    if strncmp(devs(n).DeviceName, 'Xonar STX II+H6: Multichannel', 29)
        deviceID=devs(n).DeviceIndex;
    end
end