%enter desired file name in NAME
%need to figure out a way to make it sense desired tetrode 

% usage:    g = number of graphs per cluster in given plotting program
%               (i.e. n = 3 for PINP)
%           c = number of clusters across tetrodes 
%           

function add2ps(g,c)
n = g* c; 
for i = 1:n
    figure(i)
    cd d:\lab\emily
    print -dpsc2 pwd.ps -append
    close
end