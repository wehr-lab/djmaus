function clean_kilosort
%this is a simple utility that deletes all kilosort output from a data
%directory. It runs on the current working directory (pwd). This allows you to start from scratch and re-kilosort the
%directory. If there is a dirs.mat, it will also delete all kilosort output
%from the other directories listed in dirs.mat. With some "are you sure"
%windows.

wd=pwd;
dirs=[];
if exist('dirs.mat')==2
    load dirs.mat
end

f=uifigure;
str=sprintf('delete all kilosort output files in %s ?', wd);
title='Confirm delete kilosort files?';
answer=uiconfirm(f,str, title) ;
close(f)
switch answer
    case 'OK'
        clean_it
    otherwise
        fprintf('\ndid nothing')
end

for i=1:length(dirs)
    str=sprintf('following dirs.mat to companion directory\ndelete all kilosort output files in %s ?', dirs{i});
    f=uifigure;
    answer=uiconfirm(f,str, title) ;
    close(f);
    switch answer
        case 'OK'
            cd(dirs{i})
            clean_it
            cd(wd)
    end
end


function clean_it
filelist={...
    'amplitudes.npy', ...
    'whitening_mat_inv.npy', ...
    'whitening_mat.npy', ...
    'templates_unw.npy', ...
    'templates_ind.npy', ...
    'templates.npy', ...
    'template_features.npy', ...
    'template_feature_ind.npy', ...
    'spike_times.npy', ...
    'spike_templates.npy', ...
    'spike_clusters.npy', ...
    'similar_templates.npy', ...
    'rez.mat', ...
    'phy.log', ...
    'pc_features.npy', ...
    'pc_feature_ind.npy', ...
    'params.py', ...
    'OEtetrodes.dat', ...
    'dirs.mat', ...
    'cluster_groups.csv', ...
    'channel_positions.npy', ...
    'channel_map.npy', ...
    'chanMap.mat'};
count=0;
for i=1:length(filelist)
    d=dir(filelist{i});
    if ~isempty(d)
        delete(filelist{i})
        fprintf('\ndeleted %s', filelist{i})
        count=count+1;
    end
end

if exist(fullfile(pwd, '.phy'),'dir')
    rmdir('.phy', 's')
    fprintf('\ndeleted .phy directory')
    count=count+1;
end
if exist(fullfile(pwd, 'Figs'),'dir')
    cd Figs
    delete *.fig
    cd ..
    rmdir Figs
    fprintf('\ndeleted Figs directory')
    count=count+1;
end
fprintf('\nfound and deleted %d kilosort files/folders',count);
fprintf('\ndeleted all kilosort output files\n in %s ', pwd);
