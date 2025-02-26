function NodeIds = getNodes( datadir )
%reads settings.xml in datadir to figure out which nodes are being used
% cd(datadir)
settings=djxml2struct(fullfile(datadir,'settings.xml'));

k=0;
NodeIds={};
%this helpfully prints out all the nodes and their names
signalchain=settings.SETTINGS.SIGNALCHAIN;
for i=1:length(signalchain)
    processors=signalchain{i}.PROCESSOR;
    for j=1:length(processors)
        if iscell(processors)
            k=k+1;
            NodeIds{k}=    processors{j}.Attributes.NodeId;
        else
            if length(processors)==1
                k=k+1;
                NodeIds{k}=  processors.Attributes.NodeId;
            else
                error('wtf')
            end
        end
    end
end