% compile 0MQ matlab wrapper 
%
% NOTE: On OS X, you may need to implement this patch
%       for compilation to work:
%       http://www.mathworks.com/matlabcentral/answers/94092
%

% note: you will have to change the file locations below to match your
% local machine configuration

% note: on windows, check the output of mex -setup. YOu might need to install the mingw c++ compiler.
% install the
% MinGW g++ compiler by clicking the Add-Ons button on the matlab toolbar,
% checking the Features box, and installing MinGW support. Then run mex
% -setup and click the MinGW compiler link.
% make sure it's set to use the C++ compiler, not the C compiler (which might be the default)
%
%note: it's possible that restarting windows after re-compiling is helpful.

% 1. make sure mex is setup properly and a compiler is available
mex -setup C++

% it helps to execute this from within the directory
cd('C:\Users\WehrLab\Documents\GitHub\djmaus')

%Where is the plugin-GUI folder?
GUIfolder = 'C:\Users\WehrLab\Documents\GitHub\plugin-GUI'; %wehrrig3


headerFolder = [GUIfolder, '\Resources\windows-libs\ZeroMQ\include'];% for compiled version of OE

if strcmp(computer,'PCWIN')
    libFolder = [GUIfolder, '\Resources\windows-libs\ZeroMQ\lib_x86'];
    libraryName = 'libzmq-v110-mt-4_0_4';
    cppFile = 'windows/zeroMQwrapper.cpp';
elseif strcmp(computer,'PCWIN64')
    libFolder = [GUIfolder, '\Resources\windows-libs\ZeroMQ\lib_x64']; % for compiled version of OE
    %libFolder =   'C:\lab\ZeroMQ4.0.4\lib'; %rig2, binary download
    %libraryName = 'libzmq-v110-mt-3_2_2';
    libraryName = 'libzmq-v120-mt-4_0_4.lib';
    cppFile=  'C:\Users\WehrLab\Documents\GitHub\djmaus\windows\zeroMQwrapper.cpp';
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