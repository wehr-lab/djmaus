function spiketimes=readKiloSortOutput(cellnum)
%reads Kilosort output, finds cell's spiking time
% plot's statistics about this cell
%
%usage: spiketimes=read_MClust_output(cell)
%input: cell - cell number (different from cell ID assigned by kilosort
%output: spiketimes in seconds
%
%ira 7.25.2018

% load all cells clustered by Kilosort
sp = loadKSdir(pwd);

[cluster_ID,group_ID]=readClusterGroupsCSV('cluster_groups.csv');
cell_ID=cluster_ID(cellnum); %get KS ID for this cell to find its spiketimes
%find spiketimes
spiketimes=sp.st(sp.clu==cell_ID); % in seconds, start at 0












