% compile 0MQ matlab wrapper
%
% NOTE: On OS X, you may need to implement this patch
%       for compilation to work:
%       http://www.mathworks.com/matlabcentral/answers/94092
%


% 1. make sure mex is setup properly and a compiler is available
mex -setup

% it helps to execute this from within the directory
% cd ('/Users/mikewehr/Documents/Analysis/plugin GUI/plugin-GUI/Resources/Matlab')
%cd('D:\lab\djmaus\')
djhome

%GUIfolder = 'C:\Users\Shay\Documents\GitHub\GUI';
%GUIfolder = '/home/jsiegle/Programming/GUI/';
GUIfolder = '/home/lab/plugin-GUI';



%headerFolder = [GUIfolder, '/Resources/windows-libs/ZeroMQ/include'];
%headerFolder = [GUIfolder, '\Resources\windows-libs\ZeroMQ\include'];
%if you're trying to figure out what your headerFolder is try which zmq
headerFolder = '/usr/include';


if strcmp(computer,'PCWIN')
    libFolder = [GUIfolder, '/Resources/ZeroMQ/lib_x86'];
    libraryName = 'libzmq-v110-mt-3_2_2';
    cppFile = 'windows/zeroMQwrapper.cpp';
elseif strcmp(computer,'PCWIN64')
    libFolder = [GUIfolder, '\Resources\windows-libs\ZeroMQ\lib_x64'];
    %libraryName = 'libzmq-v110-mt-3_2_2';
    libraryName = 'libzmq-v120-mt-4_0_4.lib';
    cppFile = 'windows/zeroMQwrapper.cpp';
    
elseif strcmp(computer,'GLNX86') || strcmp(computer,'GLNXA64')
    libFolder = '/usr/local/lib';
    libraryName = 'zmq';
    cppFile = 'unix/zeroMQwrapper.cpp';
elseif strcmp(computer,'MACI64')
    libFolder = '/opt/local/lib';
    libraryName = 'zmq';
    %cppFile = 'unix/zeroMQwrapper.cpp';
    cppFile = 'mac/zeroMQwrapper.cpp';
    
end

% 2. compile
eval(['mex ' cppFile ' -I',headerFolder,' -L',libFolder,' -l',libraryName] )
