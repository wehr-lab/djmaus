function deviceID=GetAsioLynxDevice
%looks in Windows soundcard device list and 
%returns the deviceID of the ASIO Lynx adaptor
%
%modified from PrintDevices
%mw 01-24-2012

%AssertOpenGL;
InitializePsychSound(1);

devs = PsychPortAudio('GetDevices');
deviceID=[];
for n = 1:length(devs) 
    if strcmp(devs(n).DeviceName, 'ASIO Lynx')
        deviceID=devs(n).DeviceIndex;
    end
end