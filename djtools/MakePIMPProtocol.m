function MakePIMPProtocol(varargin)
%takes an existing protocol and inserts N laser pulses randomly, specify N

djPrefs
global pref
dbstop if error
if nargin==1
    stimuli=varargin{1};
elseif nargin==0 %select file to convert
    currentdir=pwd;
    cd(pref.stimuli);
    [filename, pathname] = uigetfile('*.mat', 'Pick a protocol file to convert');
    cd(currentdir);
    
    if isequal(filename,0) || isequal(pathname,0)
        return;
    else
        load(fullfile(pathname,filename));
    end
else error('convert_stimulus_protocol: wrong number of arguments')
end

L=length(stimuli);
if nargin>0
    nn=floor(L/N);
else
    N=20;
    nn=floor(L/N);
end

left=stimuli;

NN=nn;
for i=1:N
    left=stimuli(NN:end);
    stimuli(NN).type='whitenoise';
    stimuli(NN).param;
    if strcmp(stimuli(1).type,'tone') || strcmp(stimuli(1).type,'WN')
        stimuli(NN).param.frequency=-1;
        stimuli(NN).param.amplitude=-1000;
        stimuli(NN).param.duration=25;
        stimuli(NN).param.laser=1;
        stimuli(NN).param.ramp=5;
        stimuli(NN).param.next=stimuli(1).param.next;
    else
    end
    stimuli(NN).stimulus_description='silent laser pulse';
    stimuli(NN).protocol_name=stimuli(1).protocol_name;
    stimuli(NN).protocol_description=stimuli(1).protocol_description;
    stimuli(NN).version=stimuli(1).version;
    stimuli(NN+1:L+i)=left;
    NN=NN+nn+1;
end



cd('/home/lab/djmaus/stimuli/PIMP')
filename=[stimuli(1).protocol_name '_PIMP']
save(filename, 'stimuli')
end

