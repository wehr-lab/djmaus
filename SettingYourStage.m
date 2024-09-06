function [fpath] = SettingYourStage() 
%choose, save, and/or load the directories you wish to run kilosort on together
%tries to identify the bonsai and openephys directories
%each session has the bonsai directory as the parent, and the ephys directory within it
%
%outputs: creates bdirs.mat and dirs.mat files.
% Bdirs is a cell array that contains a list of the bonsai directories that you selected
% dirs is a cell array that contains a list of the openephys directories corresponding to your selection
% these are absolute paths as they appear on the machine you run this function from. For example, if you're running it over the network, they will contain the network path as it appears from your machine 
% directories you selected
% also saves the DataRoot (parent folder to the bonsai directories) as it
% appears on your machine
%   original - Nick 2/9/2021
%   modified by mike 4.18.2024 to work with OpenEphys 0.6 and added to djmaus repository 
%
% can we get it to run from either ephys or bonsai dir?
%
% 
                
try
    load('Bdirs.mat') 
catch
    Bdirs = uigetfile_n_dir(pwd, 'Select one or more Bonsai directories');
    for idx = 1:length(Bdirs) %this just extracts a timestamp from the folder name to sort dirs chronologically
        test = strsplit(Bdirs{idx}, filesep);
        test = test{end}; test = strsplit(test,'_mouse'); test = test{1};
        test = strsplit(test,'mouse'); test = test{1}; %sometimes there's no _ before mouse
        if ~isequal(length(test),19)
            test = strsplit(test,'_');
            test = strcat(test{1},'_0',test{2});
        end
        dirDT(idx) = datetime(test, 'Format','yyyy-MM-dd_HH-mm-ss');
    end
    [~,I] = sort(dirDT);
    for idx = 1:length(Bdirs)
        SortedDirs{idx} = Bdirs{I(idx)};
    end
    Bdirs = SortedDirs; %orders dirs chronologically
end


%% Determine if the OE folder is inside of each trial's directory (i.e. if bonsai was used for this experiment) - added by Nick 2/9/2021
chosendirs = Bdirs; %the dirs you just selected manually
for i =1:length(Bdirs)
    cd(Bdirs{i})
    %test = dir('*.continuous'); this doesn't work anymore since new OE doesn't save .continuos files
    test = dir('notebook.mat'); %notebook.mat indicates we're in an ephys dir 
    if isempty(test) %the ephys files are not in the main experiment directory, so let's find the OE folders and reassign the dirs paths so kilosort can successfully navigate to them
        test = dir();
        narrowdown = find([test.isdir]); %identifies folders present in our main experiment folder
        for k = 1:length(narrowdown)
            testing = strsplit(test(k).name,'mouse');
            if length(testing) > 1 %then it is a folder with 'mouse' in its name, so it's almost certainly the OE folder...
                ephysfolder = strcat(test(k).folder,filesep,test(k).name);
            end
        end
        
        if exist('ephysfolder','var') %we found the OE folder
            dirs{i} = ephysfolder; %let's reassign the dirs path to point to the OE folder so kilosort can successfully navigate to it
        else %we didn't find an OE folder
            error('no OE files or folder found... what-what-what?!?')
        end
        clear test; clear testing; clear narrowdown;
        
        %Lastly, we've confirmed the chosendirs were bonsai folders, so
        %lets declare the bonsaifolder
        bonsaifolder = chosendirs{1};
    end
end

if exist('ephysfolder','var') %Then OE folders were found inside the trial folder directories (therefore: bonsai was used for these experiments)
    %That means we also just properly redefined the dirs variable to point to the OE folders so kilosort can successfully navigate between them.
    %In a few lines down, we'll save this dirs variable in each of the OE folders as usual, but in the meantime let's also save the variable in the bonsai directories 
    %in a file called: 'Bdirs.mat' so kilosort can easily load it in the future. We'll also include a Bdirs variable that contains the bonsai folders
    %we just selected for batch processing:
    Bdirs = chosendirs;
    
    % sanity check to make sure mouseID matches for Bdirs and dirs:
    idx=strfind(ephysfolder, 'mouse');
    mouseidephys=ephysfolder(idx(2):end);
    idx2=strfind(bonsaifolder, 'mouse');
    mouseidbonsai=bonsaifolder(idx2:end);
    if ~isequal(mouseidbonsai, mouseidephys)
       error('mouse ids dont match between ephys and bonsai folders') 
    end
    [DataRoot,~, ~]=fileparts(chosendirs{1});
    if ~strcmp(DataRoot(end), filesep)
        DataRoot=[DataRoot, filesep];
    end
    for i =1:length(chosendirs)
        cd(chosendirs{i})
        save('Bdirs.mat','dirs','Bdirs','DataRoot');
    end
    clear ephysfolder; clear chosendirs; clear i; clear k;
else error('what happened?')
end

for i =1:length(dirs)
    cd(dirs{i})
    save('dirs.mat', 'Bdirs', 'dirs','DataRoot');
end
fpath    = dirs{1}; %lastly, we need to set the fpath (the 'master directory') which is used by your configuration file to generate a few critical fields in the 'ops' structure
end