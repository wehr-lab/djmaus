function PrintDevices(deviceid, apiid, apiname)
AssertOpenGL;
InitializePsychSound(1);

if nargin < 1 
    deviceid = 0;
    apiid = 0;
    apiname = 0;
end

devs = PsychPortAudio('GetDevices');

for idx = 1:length(devs) disp(devs(idx));
    if devs(idx).DeviceIndex == deviceid-1
        break;
    end
end