useGPU = 1; % do you have a GPU? Kilosorting 1000sec of 32chan simulated data takes 55 seconds on gtx 1080 + M2 SSD.

% fpath    = pwd; % where on disk do you want the simulation? ideally and SSD...
% if ~exist(fpath, 'dir'); mkdir(fpath); end
[fpath] = SettingYourStage(); %function to choose, save, and/or load the directories you wish to process together. fpath is your 'master directory', which is used when you run your configuration file a few lines down.

% This part adds paths
addpath(genpath('C:\Users\wehrlab\Documents\GitHub\djmaus\djtools\Kilosort8TT')) % path to kilosort folder
addpath(genpath('C:\Users\lab\Documents\npy-matlab')) % path to npy-matlab scripts
pathToYourConfigFile = 'C:\Users\wehrlab\Documents\GitHub\djmaus\djtools\Kilosort8TT'; % for this example it's ok to leave this path inside the repo, but for your own config file you *must* put it somewhere else!  

% Run the configuration file, it builds the structure of options (ops)
run(fullfile(pathToYourConfigFile, 'config8TT.m'))

% This part makes the channel map for this simulation
make_8TTChannelMap(fpath);

% This part simulates and saves data. There are many options you can change inside this 
% function, if you want to vary the SNR or firing rates, or number of cells etc. 
% You can vary these to make the simulated data look more like your data.
% Currently it is set to relatively low SNR for illustration purposes in Phy. 
%
% This part runs the normal Kilosort processing on the simulated data
[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

% This runs the benchmark script. It will report both 1) results for the
% clusters as provided by Kilosort (pre-merge), and 2) results after doing the best
% possible merges (post-merge). This last step is supposed to
% mimic what a user would do in Phy, and is the best achievable score
% without doing splits. 
% benchmark_simulation(rez, fullfile(fpath, 'eMouseGroundTruth.mat'));

% save python results file for Phy
rezToPhy(rez, fpath);

fprintf('Kilosort took %2.2f seconds vs 72.77 seconds on GTX 1080 + M2 SSD \n', toc)

% now fire up Phy and check these results. There should still be manual
% work to be done (mostly merges, some refinements of contaminated clusters). 
%% AUTO MERGES 
% after spending quite some time with Phy checking on the results and understanding the merge and split functions, 
% come back here and run Kilosort's automated merging strategy. This block
% will overwrite the previous results and python files. Load the results in
% Phy again: there should be no merges left to do (with the default simulation), but perhaps a few splits
% / cleanup. On realistic data (i.e. not this simulation) there will be drift also, which will usually
% mean there are merges left to do even after this step. 
% Kilosort's AUTO merges should not be confused with the "best" merges done inside the
% benchmark (those are using the real ground truth!!!)

rez = merge_posthoc2(rez);
%benchmark_simulation(rez, fullfile(fpath, 'eMouseGroundTruth.mat'));

% save python results file for Phy
rezToPhy(rez, fpath);

%% save and clean up
% save matlab results file for future use (although you should really only be using the manually validated spike_clusters.npy file)
save(fullfile(fpath,  'rez.mat'), 'rez', '-v7.3');

% remove temporary file
delete(ops.fproc);
%%

% %%%%%%%% functions %%%%%%%%%% -Nick 2/9/2021
function [fpath] = SettingYourStage() %choose, save, and/or load the directories you wish to process together
try
    load('dirs.mat')
catch                   
    try
        load('Bdirs.mat') %For compatible bonsai experiments - Nick 2/9/2021
    catch
        dirs = uigetfile_n_dir(pwd);
        dirs = sort(dirs); %orders dirs chronologically
    end
end
%% Determine if the OE folder is inside of each trial's directory (i.e. if bonsai was used for this experiment) - added by Nick 2/9/2021
chosendirs = dirs; %the dirs you just selected manually
for i =1:length(dirs)
    cd(dirs{i})
    test = dir('*.continuous');
    if isempty(test) %the ephys files are not in the main experiment directory, so let's find the OE folders and reassign the dirs paths so kilosort can successfully navigate to them
        test = dir();
        narrowdown = find([test.isdir]); %identifies folders present in our main experiment folder
        for k = narrowdown
            testing = strsplit(test(k).name,'_mouse-');
            if length(testing) > 1 %then it is a folder with '_mouse-' in its name, so it's almost certainly the OE folder...
                ephysfolder = strcat(test(k).folder,'\',test(k).name);
            end
        end
        
        if exist('ephysfolder','var') %we found the OE folder
            dirs{i} = ephysfolder; %let's reassign the dirs path to point to the OE folder so kilosort can successfully navigate to it
        else %we didn't find an OE folder
            print('no OE files or folder found... what-what-what?!?')
        end
        clear test; clear testing; clear narrowdown;
    end
end

if exist('ephysfolder','var') %Then OE folders were found inside the trial folder directories (therefore: bonsai was used for these experiments)
    %That means we also just properly redefined the dirs variable to point to the OE folders so kilosort can successfully navigate between them.
    %In a few lines down, we'll save this dirs variable in each of the OE folders as usual, but in the meantime let's also save the variable in the bonsai directories 
    %in a file called: 'Bdirs.mat' so kilosort can easily load it in the future. We'll also include a Bdirs variable that contains the bonsai folders
    %we just selected for batch processing:
    Bdirs = chosendirs;
    for i =1:length(chosendirs)
        cd(chosendirs{i})
        save('Bdirs.mat','dirs','Bdirs');
    end
    clear ephysfolder; clear chosendirs; clear i; clear k;
end

for i =1:length(dirs)
    cd(dirs{i})
    save('dirs.mat','dirs');
end
fpath    = dirs{1}; %lastly, we need to set the fpath (the 'master directory') which is used by your configuration file to generate a few critical fields in the 'ops' structure
end
