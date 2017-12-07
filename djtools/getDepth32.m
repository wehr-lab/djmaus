function depth= gecDepch32(varargin)

% Figure ouc che depch of cells in recordings wich silicon probes
% dir, ch, c, depch (micromanipulacor)

if nargin==0
    datadir=pwd;
    ch=[];
    c=[];
    d=[];
elseif nargin==1
    datadir=varargin{1};
    ch=[];
    c=[];
    d=[];
elseif nargin==2
    datadir=varargin{1};
    ch=varargin{2};
    c=[];
    d=[];
elseif nargin==3
    datadir=varargin{1};
    ch=varargin{2};
    c=varargin{3};
    d=[];
elseif nargin==4
    datadir=varargin{1};
    ch=varargin{2};
    c=varargin{3};
    d=varargin{4};
end
cd(datadir)
depth=[];

if isempty(d)
    load('notebook.mat')
    d=nb.Depth;
    fprintf(' \nFound recorded depth - %s \n', d)
    try
        if isempty(d)
            prompt= 'Could not find depth of the recording in the notebook. Please enter depth (um) here:';
            d=input(prompt);
        end
    catch
        fprintf('\nCannot calculate cell depth\n')
    end
end
D=sqrt(d^2+d^2);
if ~isempty(d) && isempty(ch) && isempty(c)
    fnames=dir('*-wv.mat');
    for i=1:length(fnames)
        [p,f,ext]=fileparts(fnames(i).name);
        split=strsplit(f, '_');
        channel=strsplit(split{1}, 'ch');
        ch(i)=str2num(channel{2});
        split=strsplit(split{end}, '-');
        c(i)=str2num(split{1});
        load(f)
        plot(mWV);
        max1_peak_ch=find(min(min(mWV))==min(mWV)); %channel
        pos=32-((ch(i)-1)*4+max1_peak_ch)+1;
        d=D-25*pos; %all site on 32 channel are 25um apart
        depth(i)=d;
        cells(i).ch=ch(i);
        cells(i).c=c(i);
        cells(i).depth=depth(i);
    end
    
end
filename='depth.txt';
for i=1:length(depth)
    Depth=sprintf('%.4f', depth(i));

 save(filename, 'Depth', '-ASCII', '-append');
end

save('depth.mat','cells')
