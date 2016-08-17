function varargout=convert_stimulus_protocol(varargin)
%converts old exper2 stimulus protocols to new djmaus protocols
%the difference in formats is that exper2 protocols have a stimulus(1) that has name
%and description, with the stimuli starting at 2. The new format is that
%all stimuli have a copy of name and description, and stimuli start at 1.
%usage:
% convert_stimulus_protocol
%    select exper2 protocol with dialog box, writes converted file
% stimuli=convert_stimulus_protocol(stimuli)
%    converts the input stimuli and returns it as output variable

djPrefs
global pref

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
if strcmpi(stimuli(1).type,'exper2 stimulus protocol')
    %needs to be converted
    in=stimuli;
    clear stimuli;
    name=in(1).param.name;
    description=in(1).param.description;
    for i=1:length(in)-1
        a=in(i+1);
        paramstring=sprintf('%s',a.type);
        fns=fieldnames(a.param);
        for fn=1:length(fns)
            paramstring=[paramstring sprintf(' %s:%d',fns{fn},getfield(a.param, fns{fn}))];
        end
        stimuli(i).type=a.type;
        stimuli(i).param=a.param;
        stimuli(i).protocol_name=name;
        stimuli(i).stimulus_description=paramstring;
        stimuli(i).protocol_description=description;
        stimuli(i).version='djmaus';
    end
    
else
    fprintf('this does not appear to be an exper2 stimulus protocol')
end
if nargout==0
        fname=[fullfile(pathname,filename)];
    fname=strrep(fname, '.mat', 'dj.mat');
    save(fname, 'stimuli')
    fprintf('\n wrote file %s\n', fname)
elseif nargout==1
varargout{1}=stimuli;
else
    error('convert_stimulus_protocol: wrong number of output arguments')
end