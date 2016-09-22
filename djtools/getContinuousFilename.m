function fname = getContinuousFilename( datadir, channel )
%figures out the filename for the continuous channel
cd(datadir)

NodeIds=getNodes(datadir);
for i=1:length(NodeIds)
    filename=sprintf('%s_CH%d.continuous', NodeIds{i}, channel);
    absfilename=fullfile(datadir, filename);
    if exist(absfilename,'file')
        fname=filename;
    end
end

