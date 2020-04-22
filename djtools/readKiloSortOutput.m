function [spiketimes, cell_ID]=readKiloSortOutput(cellnum, sampleRate) %output spiketimes and Kilosort ID number
%reads Kilosort output, finds cell's spiking time
% plot's statistics about this cell
%
%usage: spiketimes=read_MClust_output(cell)
%input: cell - cell number (different from cell ID assigned by kilosort
%output: spiketimes in seconds
%
%ira 7.25.2018

% load all cells clustered by Kilosort

load('dirs.mat') %find all directories that were clustered in one session of kilosort. all clusters are saved in the first directory
masterdir=dirs{1};
currentdir=pwd; %remember which directory you are in now
currentdir_indx=find(strcmp(currentdir, dirs)==1); %which dir are we trying to plot?
if currentdir_indx==0
    error('This directory cannot be found on the list of clustered directories. \n Either this data has not been clustered or something bad happened')
end
cd(masterdir) %go to the first directory to load clustered data, we can call this master directory
sp = loadKSdir(pwd); %load all cells, all spikes

cell_ID=sp.cids(cellnum); %get Kilosort id
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

fprintf('\nthis cell was saved as %s cluster', qual);

%find spiketimes
spiketimes=sp.st(sp.clu==cell_ID); % in seconds, start at 0
try
    load('RecLengths.mat')
    L=recLengths/sampleRate;
catch
    % load rez, which contains number of samples of each recording 1=1, 2=1+2,
    % 3=1+2+3, etc
    load('rez.mat')
    
    L=(rez.ops.recLength)/sampleRate;
    save('RecLengths.mat','L')
end

stop=L(currentdir_indx);
if currentdir_indx==1
    start=0;
else
    start=L(currentdir_indx-1);
end
spiketimes=spiketimes(spiketimes>start & spiketimes<stop); %find spiketimes for this recording;
spiketimes=spiketimes-start; %all spiketimes will start at 0
spiketimes=spiketimes';
cd(currentdir); %go to the original directory



