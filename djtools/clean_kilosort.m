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
end

for i=1:length(dirs)
    str=sprintf('delete all kilosort output files in %s ?', dirs{i});
    f=uifigure;
    answer=uiconfirm(f,str, title) ;
    close(f);
    switch answer
        case 'OK'
            cd(dirs(i})
            clean_it
            cd(wd)
    end
end


function clean_it
delete 'whitening_mat_inv.npy'
delete 'whitening_mat.npy'
delete 'templates_unw.npy'
delete 'templates_ind.npy'
delete 'templates.npy'
delete 'template_features.npy'
delete 'template_feature_ind.npy'
delete 'spike_times.npy'
delete 'spike_templates.npy'
delete 'spike_clusters.npy'
delete 'similar_templates.npy'
delete 'rez.mat'
delete 'phy.log'
delete 'pc_features.npy'
delete 'pc_feature_ind.npy'
delete 'params.py'
delete 'OEtetrodes.dat'
delete 'dirs.mat'
delete 'cluster_groups.csv'
delete 'channel_positions.npy'
delete 'channel_map.npy'
delete 'chanMap.mat'

if exist('.phy','dir')
    cd .phy
    delete 'joblib'
    delete 'memcache'
    delete 'new_cluster_id.json'
    cd ..
    rmdir .phy
end
if exist('Figs','dir')
    cd Figs
    delete *.fig
    cd ..
    rmdir Figs
end
fprintf('deleted all kilosort output files in %s ', pwd);
