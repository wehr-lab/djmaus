function deviceID=GetLynxTwoBDevice
%looks in Windows soundcard device list and 
%returns the deviceID of the ASIO Lynx adaptor
%
%modified from PrintDevices
%mw 01-24-2012

%AssertOpenGL;
% InitializePsychSound(1); %this is throwing a message every stimulus but I
% don't think we need it. mw 08.23.2017
%(we use GetAsioLynxDevice to determine which soundcard is used and
%therefore how to set amplitude of soundcard trigger)

try
    devs = PsychPortAudio('GetDevices');
catch
    InitializePsychSound(1);
    devs = PsychPortAudio('GetDevices');
end
deviceID=[];
for n = 1:length(devs) 
    if strcmp(devs(n).DeviceName, 'Speakers (LynxTWO-B 01+02)') & strcmp(devs(n).HostAudioAPIName, 'Windows WDM-KS')
        deviceID=devs(n).DeviceIndex;
        break;
    end
end