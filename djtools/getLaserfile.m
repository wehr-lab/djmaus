function [ fname ] = getLaserfile( datadir )
%reads settings.xml in datadir to figure out which node the Laser
%was recorded on, and hence what the filename is
%usage: [ fname ] = getLaserfile( datadir )

try
    cd(datadir)
    settings=xml2struct('settings.xml');
    Laserchannel=3; %filenames and ADClines are both 1-indexed
catch
    fname=[];
    return
end

% ADC1: sound monitor
% ADC2: soundcard trigger monitor
% ADC3: laser monitor

% As far as we can tell, channel numbers are 0-indexed in the settings.xml
% list. (8.30.17)
% So:
% 0-31 = tetrode channels
% 32 = aux1 accelerometer
% 33 = aux2 accelerometer
% 34 = aux3 accelerometer
% 35 = adc1 sound monitor
% 36 = adc2 soundcard trigger
% 37 = adc3 laser monitor

if (0)
    %this helpfully prints out all the nodes and their names
    signalchain=settings.SETTINGS.SIGNALCHAIN;
    for i=1:length(signalchain)
        processors=signalchain{i}.PROCESSOR;
        for j=1:length(processors)
            if iscell(processors)
                fprintf('\n%s: %s',    processors{j}.Attributes.NodeId, processors{j}.Attributes.name)
            else
                if length(processors)==1
                    fprintf('\n%s: %s',    processors.Attributes.NodeId, processors.Attributes.name)
                else
                    error('wtf')
                end
            end
        end
    end
end


%this searches all the nodes to see which has "record" turned on
%we are looking for the node that has ch37 turned on, which is ADC1 which
%should be recording the laser
NodeId={};
node_idx=0;
signalchain=settings.SETTINGS.SIGNALCHAIN;
for i=1:length(signalchain)
    processors=signalchain{i}.PROCESSOR;
    for j=1:length(processors)
        if iscell(processors)
            if isfield(processors{j}, 'CHANNEL')
                channels=processors{j}.CHANNEL;
                for ch=1:length(channels)
                    if    str2num(channels{ch}.SELECTIONSTATE.Attributes.record)
%                         fprintf('\n%s: %s ch %s is being recorded',    processors{j}.Attributes.NodeId, processors{j}.Attributes.name, channels{ch}.Attributes.number)
                        if  strcmp(channels{ch}.Attributes.number, '37')
                            node_idx=node_idx+1;
                            NodeId{node_idx}=processors{j}.Attributes.NodeId;
                        end
                    end
                end
            end
        else
            if length(processors)==1
                if isfield(processors, 'CHANNEL')
                    
                    channels=processors.CHANNEL;
                    for ch=1:length(channels)
                        if    str2num(channels{ch}.SELECTIONSTATE.Attributes.record)
                            if   strcmp(channels{ch}.Attributes.number, '37')
                                node_idx=node_idx+1;
                                NodeId{node_idx}=processors.Attributes.NodeId;
                            end
                            
                           % fprintf('\n%s: %s ch %s is being recorded',    processors.Attributes.NodeId, processors.Attributes.name, channels{ch}.Attributes.number)
                            %  don't bother, this is the network events
                            %except that on Rig 2 as of Aug 30, 2017, it's FPGA, so we need
                            %to check for NodeID here too

                        end
                    end
                end
            else
                error('wtf')
            end
        end
    end
end

% If there are multiple NodeIds, we just take the first one (hopefully the rhythm FPGA)

filename=sprintf('%s_ADC%d.continuous', NodeId{1}, Laserchannel);
absfilename=fullfile(datadir, filename);
if exist(absfilename,'file')
fname=filename;
else
    fname=[];
end
