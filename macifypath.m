function newpath=macifypath(path, varargin)
%convert a  windows path to a remote mac path
%usage: newpath=macifypath(path, [rig])
% you might need to specify rig as rig3, rig4, etc

%to-do: this could be split into a macify function separate from a
%local-to-remote function. For now mac is the only remote use case.

% catalina mounts each drive (D:, E:) as Volumes/wehrrig3.uoregon.edu,
% Volumes/wehrrig3.uoregon.edu-1 based on the order they are mounted.

if ~ismac 
    newpath=path;
    return 
end
if nargin==1
    rig='wehr-nas'; %default 
else
    rig=varargin{1};
end

newpath= strrep(path, '\', '/');
switch rig
    case 'wehr-nas'
        %you probably need to edit this for each use case
        % newpath= strrep(newpath, 'F:/Data/sfm', '/Volumes/Projects/temp');
        newpath= strrep(newpath, '//wehr-nas.uoregon.edu/projects/', '/Volumes/Projects/');

    case 'rig1'
        newpath= strrep(newpath, 'E:', '/Volumes/wehrrig1b.uoregon.edu');
    case 'rig3'
        newpath= strrep(newpath, 'D:', '/Volumes/wehrrig3.uoregon.edu-1');
        newpath= strrep(newpath, 'E:', '/Volumes/wehrrig3.uoregon.edu');
    case 'rig4'
        newpath= strrep(newpath, 'D:', '/Volumes/wehrrig4.uoregon.edu');
    
    case 'ion-nas'
        
        newpath= strrep(newpath, 'ion-nas.uoregon.edu', 'Volumes');
    otherwise
        newpath= strrep(newpath, 'D:', '/Volumes/wehrrig4.uoregon.edu');
end
