%% Resort Tool - SFM 8/11/21

% The following files are raw data and never to be moved/deleted from any
% djmaus directory:
%
% settings.xml
% stimlog.txt
% 101*.continuous
% TT*.spikes
% notebook.mat
% messages.events
% Continuous_Data.openephys
% all_channels.events
% *fig

masterdirs = {'D:\lab\djmaus\Data\sfm\2021-01-18_11-55-16_mouse-0295'};    % Enter list of MASTER directories (can be a list in the format {'dir1'; 'dir2'; 'etc';}),  a master directory is a data directory (generated by djmaus)
                                                                           % where Kilosort has been run, either for a single session or for an entire mouse
for i = 1:length(masterdirs)
    cd(masterdirs{i});
    numpys = dir('*.npy');
    dirlist = dir('dirs.mat');      % If you need to resort the same group of directories (like an entire mouse's sessions) that won't change, this doesn't have to be deleted
    channelmap = dir('chanMap.mat');
    datfile = dir('OEtetrodes.dat');
    params = dir('params.py');
    recordinglengths = dir('RecLengths.mat');
    clusterlog = dir('cluster_groups.csv');
    phylog = dir('phy.log');
    rez = dir('rez.mat');
    if ~isempty(dir('out*.mat'))
        outlist = dir('out*.mat');
        KilosortGeneratedData = vertcat(numpys, dirlist, channelmap, datfile, params, recordinglengths, clusterlog, phylog, rez, outlist);
    else
        KilosortGeneratedData = vertcat(numpys, dirlist, channelmap, datfile, params, recordinglengths, clusterlog, phylog, rez);
        fprintf('No outfiles detected in %s', masterdirs{i})
    end
    
    phyautodir = strcat(masterdirs{1}, '\', '.phy');
    
%     if isempty(KilosortGeneratedData)
%         error('There is no data in this directory that would be generated by Kilosort/Phy, it has either not been processed yet or has already been moved')
%     else
%     end
    
    if ~exist('Previously Sorted Data', 'dir')
        dirname = 'Previously Sorted Data';
        mkdir(dirname);
        for j = 1:length(KilosortGeneratedData)
            movefile(KilosortGeneratedData(j).name, dirname);
        end
        movefile(phyautodir, dirname);
        fprintf('Previously processed data has been moved to %s in %s', dirname, masterdirs{i});
    else
        fprintf('There is already another folder in this directory containing previously processed data, making a new one!')
        tempdirlist = dir('dir');
        dirnameindex = num2str(length(tempdirlist) + 1);
        dirname = strcat('Previously Sorted Data', dirnameindex);
        mkdir(dirname);
        for j = 1:length(KilosortGeneratedData)
            movefile(KilosortGeneratedData(j).name, dirname);
        end
        movefile(phyautodir, dirname);
    end
end

