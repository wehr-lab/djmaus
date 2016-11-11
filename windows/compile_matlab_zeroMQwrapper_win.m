% compile 0MQ matlab wrapper 
%
% NOTE: On OS X, you may need to implement this patch
%       for compilation to work:
%       http://www.mathworks.com/matlabcentral/answers/94092
%

% note: you will have to change the file locations below to match your
% local machine configuration

% note: on windows, check the output of mex -setup. If it's set up to use
% the Visual Studio compiler -- "Building with 'Microsoft Visual C++ 2013
% Professional'", you will probably not be able to compile because of VS
% peculiarities with variable array size allocation. Instead, install the
% MinGW g++ compiler by clicking the Add-Ons button on the matlab toolbar,
% checking the Features box, and installing MinGW support. Then run mex
% -setup and click the MinGW compiler link.
%note: it's possible that restarting windows after re-compiling is helpful.

% 1. make sure mex is setup properly and a compiler is available
mex -setup

% it helps to execute this from within the directory
% cd ('/Users/mikewehr/Documents/Analysis/plugin GUI/plugin-GUI/Resources/Matlab')
cd('D:\lab\djmaus\')

%Where is the plugin-GUI folder?
%GUIfolder = 'C:\Users\Shay\Documents\GitHub\GUI';
%GUIfolder = '/home/jsiegle/Programming/GUI/';
%GUIfolder = '/Users/mikewehr/Documents/Analysis/plugin GUI';
GUIfolder = 'D:\lab\plugin-GUI'; %wehrrig1
%GUIfolder = 'C:\lab\plugin-GUI';


headerFolder = [GUIfolder, '\Resources\windows-libs\ZeroMQ\include'];
%headerFolder = '/usr/local/include';

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
