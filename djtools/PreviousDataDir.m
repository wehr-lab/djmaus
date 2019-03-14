function PreviousDataDir
%simple utility to cd to the alphanumerically previous data directory
%I wrote this because I was tired of scrolling through a long list of data
%directories every time you cd ..

wd=pwd;
[path,name,ext]=fileparts(wd);

cd ..
d=dir;

for i=1:length(d)
    if (strcmp(d(i).name, name))
        wdi=i;
    end
end
cd(d(wdi-1).name)