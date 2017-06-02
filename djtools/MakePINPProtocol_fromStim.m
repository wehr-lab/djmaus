function MakePINPProtocol_fromStim(varargin)
%takes an existing protocol and inserts N laser pulses randomly, specify N
%only have to specify how many silent laser pulses to insert, will split
%them and insert them throughout the stimulus protocol
%works with tuning curves and wn stimuli only
%ira 04.18.2017

djPrefs
global pref
dbstop if error

    currentdir=pwd;
    cd(pref.stimuli);
    [filename, pathname] = uigetfile('*.mat', 'Pick a protocol file to convert');
    cd(currentdir);
    
    if isequal(filename,0) || isequal(pathname,0)
        return;
    else
        load(fullfile(pathname,filename));
    end


L=length(stimuli);
if nargin>0
    N=varargin{1};
    nn=floor(L/N);
else
    N=20;
    nn=floor(L/N);
end

left=stimuli;

NN=nn;
for i=1:N
    left=stimuli(NN:end);
    stimuli(NN).type='silentsound';
    if strcmp(stimuli(1).type,'tone') || strcmp(stimuli(1).type,'whitenoise')
        %stimuli(NN).param.frequency=-1;
        %stimuli(NN).param.amplitude=-1000;
        stimuli(NN).param.duration=stimuli(1).param.duration;
        stimuli(NN).param.laser=1;
        stimuli(NN).param.ramp=0;
        stimuli(NN).param.next=stimuli(1).param.next;
    else
    end
    stimuli(NN).stimulus_description='silentsound laser:1 duration:200 ramp:0 next:1200';
    stimuli(NN).protocol_name=stimuli(1).protocol_name;
    stimuli(NN).protocol_description=stimuli(1).protocol_description;
    stimuli(NN).version=stimuli(1).version;
    stimuli(NN+1:L+i)= left;
    NN=NN+nn+1;
end



cd('C:\Users\lab\Documents\GitHub\djmaus\Stimuli\PINP_fromStim')
filename=strtok(filename,'.');
filename1=[filename '_PINP.mat'];
save(filename1, 'stimuli')
end

