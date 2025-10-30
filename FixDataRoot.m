function FixedDataRoot = FixDataRoot(DataRoot, datadir)
%  Converts  NAS DataRoot path to what it looks like from local machine based on datadir
% should convert to/from Mac/Linux/Windows.
%
% This function finds the last common folder name between the two input
% paths and uses that anchor point to calculate the Windows root path
% based on the provided Windows file path example.
%
% Inputs:
%   DataRoot      - The  NAS root path as it looks like on another computer
%                      (e.g., '/Volumes/Projects/5XFAD/Rig3Phys').
%   datadir - A full path to a directory *inside* the NAS volume
%                        (e.g.,
%                        'G:\5XFAD\Rig3Phys\2025-05-06_9-14-44_mouse-3630').
%
% Output:
%   FixedDataRoot      - The converted root path as it should look on the local machine
%                      (e.g., 'G:\5XFAD\Rig3Phys').
%
% Example Usage:
%  FixDataRoot('/Volumes/Projects/5XFAD/Rig3Phys', 'G:\5XFAD\Rig3Phys\2025-05-06_9-14-44_mouse-3630')
%    returns 'G:\5XFAD\Rig3Phys'
%

    % --- 1. Identify the Anchor Folder Name ---
    % The last folder name in DataRoot is the anchor that must exist
    % in datadir. 
    
    % Get the name of the final folder (e.g., 'Rig3Phys')
    if ispc
        [~, anchorFolder, ~] = fileparts(DataRoot);
    else
        [~, anchorFolder, ~] = fileparts(macifypath(DataRoot));
    end
    % --- 2. Find the Anchor Folder in the Windows Example Path ---
    % Find the starting index of the anchor folder in the Windows example string.
    % The comparison is case-insensitive ('ignorecase' flag) for robustness,
    % as file systems can vary in case sensitivity.
    idxStart = strfind(lower(datadir), lower(anchorFolder));

    if isempty(idxStart)
        error('ConversionError: Could not find anchor folder "%s" in datadir.', anchorFolder);
    end

    % Use the *last* occurrence of the anchor folder, in case the folder name
    % appears multiple times in the path (e.g., 'C:\Data\Data_Project\Data').
    idxStart = idxStart(end);

    % --- 3. Calculate the End Index and Truncate the Path ---
    % The end index of the desired Windows root path is the start index of the
    % anchor folder plus its length, minus one.
    idxEnd = idxStart + length(anchorFolder) - 1;

    % Truncate the Windows example path to include up to and including the anchor folder
    FixedDataRoot = datadir(1:idxEnd);
    
end
