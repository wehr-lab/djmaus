function spiketimes=read_MClust_output(filename)
%reads MClust timestamp output files
%each "neuron" (or cluster) has its own output file
%
%usage: spiketimes=read_MClust_output(filename)
%input: filename of .t file in curent directory
%output: spiketimes in ms
%
%mw 12-18-2011

fid=fopen(filename, 'rb', 'b');
ReadHeader(fid);
spiketimes = fread(fid,inf,'uint32');	%read as 32 bit ints
fclose(fid);
% spiketimes=spiketimes; %already in ms



