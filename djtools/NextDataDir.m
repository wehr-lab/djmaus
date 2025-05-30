function NextDataDir
%simple utility to cd to the alphanumerically next data directory
%I wrote this because I was tired of scrolling through a long list of data
%directories every time you cd ..

if exist('notebook.mat')==2
    %we're in an ephys directory, so first cd .. up to the bonsai directory
    %comment out if your data doesn't have nested bonsai/OE folders
   
    %cd ..
end
wd=pwd;
[path,name,ext]=fileparts(wd);

cd ..
d=dir;

for i=1:length(d)
    if (strcmp(d(i).name, name))
        wdi=i;
    end
end
%wrap around if at end (comment out if you don't want to do this)
if wdi==length(d), fprintf('at last directory in root folder, wrapping around...'); wdi=3; end
cd(d(wdi+1).name)