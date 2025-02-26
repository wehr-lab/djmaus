function [spiketimes, cell_ID]=readKiloSortOutput2(clust, sp, currentdir_indx, MasterDir) %output spiketimes and Kilosort ID number

%reads Kilosort output, finds cell's spiking time
% plot's statistics about this cell
%
%usage: spiketimes=read_MClust_output(cell)
%input: cell - cell number (different from cell ID assigned by kilosort
%output: spiketimes in seconds
%
%ira 7.25.2018

% load all cells clustered by Kilosort

% load('dirs.mat') %find all directories that were clustered in one session of kilosort. all clusters are saved in the first directory
% try
% masterdir=dirs{1};
% catch
%  masterdir = pwd;
% end
% currentdir=pwd; %remember which directory you are in now
% currentdir_indx=find(strcmp(currentdir, dirs)==1); %which dir are we trying to plot?
if currentdir_indx==0
    error('This directory cannot be found on the list of clustered directories. \n Either this data has not been clustered or something bad happened')
end

cell_ID=clust; %get Kilosort id
cellnum=find(sp.cids==clust);
cg=sp.cgs(cellnum); %whats the group
if cg==0
    qual='noise';
elseif cg==1
    qual='MUA';
elseif cg==2
    qual='good';
elseif cg==3
    qual='unassigned';
end

% fprintf('\nthis cell was saved as %s cluster', qual);

%find spiketimes
spiketimes=sp.st(sp.clu==cell_ID); % in seconds, start at 0

try
    load(fullfile(MasterDir,'RecLengths.mat'))
catch

    try    % load rez, which contains number of samples of each recording 1=1, 2=1+2,
        % 3=1+2+3, etc
        load(fullfile(MasterDir,'rez.mat'))
        %L=(rez.ops.recLength)/sp.sample_rate;
        %8.31.2023 rez.ops.recLength does not exist, I assume it's a relic of
        %KS1, using rez.ops.sampsToRead instead, as a wild guess -mike
        L=(rez.ops.sampsToRead)/sp.sample_rate;

        save(fullfile(MasterDir,'RecLengths.mat'),'L')



        stop=L(currentdir_indx);
        if currentdir_indx==1
            start=0;
        else
            start=L(currentdir_indx-1);
        end
        spiketimes=spiketimes(spiketimes>start & spiketimes<stop); %find spiketimes for this recording;
        spiketimes=spiketimes-start; %all spiketimes will start at 0
        spiketimes=spiketimes';
        % cd(currentdir); %go to the original directory

    catch
        %with kilosort4, there is no longer a rez.mat and I'm not sure how to
        %get the recording durations. So punt for now and just include all
        %spiketimes
        warning('could not find rez.mat, maybe because this is kilosort4 data. Using all spiketimes which will probably work if this is a single session')
        spiketimes=spiketimes';

    end
end

