function deviceID=GetPreSonusDevice
%looks in Windows soundcard device list and 
%returns the deviceID of the PreSonus adaptor
%
%modified from PrintDevices
%mw 01-24-2012

%AssertOpenGL;
%%InitializePsychSound(1); %this is throwing a message every stimulus but I
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
% for n = 1:length(devs) 
%     if strcmp(devs(n).DeviceName, 'Line Out 1/2 (Studio 1824c)') & ...
%             strcmp(devs(n).HostAudioAPIName, 'Windows DirectSound')
%         deviceID=devs(n).DeviceIndex;
%     end
% end

for n = 1:length(devs) 
    if strcmp(devs(n).DeviceName, 'Studio USB ASIO Driver') & ...
            strcmp(devs(n).HostAudioAPIName, 'ASIO')
        deviceID=devs(n).DeviceIndex;
    end
end