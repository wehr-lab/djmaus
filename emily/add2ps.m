%enter desired file name in NAME

% usage:    g = number of graphs per cluster in given plotting program
%               (i.e. g = 3 for PINP)
%           c = number of clusters across tetrodes 
% 
% there must be an easier way to get it to sense how many figures
% there are up and/or how many total clusters there are going to be

% there needs to be a way for it to sense the file that it is dealing with
% and use that as the naem 

% there should be a way so that it can be specific for each kind of
% plotting (like it can sense GPIAS versus PINP so that it already knows
% how many graphs there will be (automatic g) 

%need to figure out how I can get the input to be for just one cell (T and
%cluser) 


function add2ps
for i = 1:57
    figure(i)
    cd D:\lab\djmaus\Data\apw\2017-06-08_10-01-54_mouse-7682
    print -dpsc2 10-01-54.ps -append
    close
end