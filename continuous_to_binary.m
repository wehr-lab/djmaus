function continuous_to_binary(varargin)

if nargin==0
    datadir=pwd;
    cd(datadir)
    file=dir('*_CH*.continuous');
elseif nargin==1
    datadir=varargin{1};
    cd(datadir)
    file=dir('*_CH*.continuous');
end
[filename_ext, datadir] = uigetfile('*.continuous', 'please choose source files','MultiSelect','on');

% for i=1:length(file)
% filename_ext{i}=file(i).name;
% end

% Convert .continuous from open ephys to a binary file that klusta can use

if ischar(filename_ext) % In the case of a single file...
    filename_ext = {filename_ext};
end

%Prepend path to files

%files = cellfun(@(x) fullfile(path,x),filename_ext,'UniformOutput',false);
files = filename_ext;
 
%%%%%%%%%%%%%%
% Make binary file from raw waveform data

data = [];
for i = 1:length(files)
    [data(:,i), timestamps, info_continuous] = load_open_ephys_data(files{i});
end

data = int32(data.*1000);

data_reshape = int32(zeros(length(data)*length(files),1));
for i = 1:length(files)
    data_reshape(i:length(files):length(data)*length(files)) = data(:,i);
end


savefname = sprintf('%s%s_%sCH_test_binary.dat',path,filename_ext{1}(1:3),num2str(length(files)));
bin_file = fopen(savefname,'w');
fwrite(bin_file,data_reshape,'int32');
fclose(bin_file);





