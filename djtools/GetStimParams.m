function stimparams=GetStimParams(datapath)
% usage: stimparams=GetStimParams(datapath)
%
%returns structure with stimulus params for stimuli that were presented for
%that data session
%
%
cd(datapath)
try
    load notebook.mat
catch
    fprintf('no notebook file found in %s', datapath)
    stimparams=[];
    return
end



allstimlogs=[];
alldurs=[];
allfreqs=[];
allcarrierfreqs=[];
allamps=[];
allpotentials=[];
allSOAs=[];
allISIs=[];
allprepulseamps=[];
allprepulsefreqs=[];
allprepulsebandwidths=[];
allpulsedurs=[];
allpulseamps=[];
allfilter_operations={};
allcenter_frequencies=[];
alllower_frequencies=[];
allupper_frequencies=[];
allgapdelays=[];
allgapdurs=[];
alldescriptions=[];
allfiles=[];
allfields={};

for i=1:length(stimlog)
    allstimlogs{i}=stimlog(i).type;
    fields=fieldnames(stimlog(i).param);
    allfields={allfields{:}, fields{:}};
    
    
    if isfield(stimlog(i).param, 'duration')
        alldurs=[alldurs stimlog(i).param.duration];
    end
    if isfield(stimlog(i).param, 'frequency')
        allfreqs=[allfreqs stimlog(i).param.frequency];
    end
    if isfield(stimlog(i).param, 'carrier_frequency')
        allcarrierfreqs=[allcarrierfreqs stimlog(i).param.carrier_frequency];
    end
    if strcmp(stimlog(i).type, 'whitenoise')
        allfreqs=[allfreqs -1];
    end
    if isfield(stimlog(i).param, 'amplitude')
        allamps=[allamps stimlog(i).param.amplitude];
    end
    if isfield(stimlog(i).param, 'prepulseamp')
        allprepulseamps=[allprepulseamps stimlog(i).param.prepulseamp];
    end
    if isfield(stimlog(i).param, 'prepulsefreq')
        allprepulsefreqs=[allprepulsefreqs stimlog(i).param.prepulsefreq];
    end
    if isfield(stimlog(i).param, 'prepulsebandwidth')
        allprepulsebandwidths=[allprepulsebandwidths stimlog(i).param.prepulsebandwidth];
    end
    if isfield(stimlog(i).param, 'pulsedur')
        allpulsedurs=[allpulsedurs stimlog(i).param.pulsedur];
    end
    if isfield(stimlog(i).param, 'pulseamp')
        allpulseamps=[allpulseamps stimlog(i).param.pulseamp];
    end
    if isfield(stimlog(i).param, 'filter_operation')
        allfilter_operations={allfilter_operations{:} stimlog(i).param.filter_operation};
    end
    if isfield(stimlog(i).param, 'upper_frequency')
        allupper_frequencies=[allupper_frequencies stimlog(i).param.upper_frequency];
    end
    if isfield(stimlog(i).param, 'lower_frequency')
        alllower_frequencies=[alllower_frequencies stimlog(i).param.lower_frequency];
    end
    if isfield(stimlog(i).param, 'center_frequency')
        allcenter_frequencies=[allcenter_frequencies stimlog(i).param.center_frequency];
    end
    if isfield(stimlog(i).param, 'SOA')
        allSOAs=[allSOAs stimlog(i).param.SOA];
    end
    if isfield(stimlog(i).param, 'gapdelay')
        allgapdelays=[allgapdelays stimlog(i).param.gapdelay];
    end
    if isfield(stimlog(i).param, 'gapdur')
        allgapdurs=[allgapdurs stimlog(i).param.gapdur];
    end
    if isfield(stimlog(i).param, 'next')
        allISIs=[allISIs stimlog(i).param.next];
    end
    if isfield(stimlog(i).param, 'description')
        
        alldescriptions{i}=stimlog(i).param.description;
    end
    if isfield(stimlog(i).param, 'file')
        allfiles{i}= stimlog(i).param.file;
    end
end

stimparams.stimtypes=unique(allstimlogs);
stimparams.durs=unique(alldurs);
stimparams.amps=unique(allamps);
stimparams.freqs=unique(allfreqs);
stimparams.carrier_freqs=unique(allcarrierfreqs);
stimparams.potentials=unique(allpotentials);
stimparams.SOAs=unique(allSOAs);
stimparams.gapdelays=unique(allgapdelays);
stimparams.gapdurs=unique(allgapdurs);
stimparams.ISIs=unique(allISIs);

stimparams.prepulseamps=unique(allprepulseamps);
stimparams.prepulsefreqs=unique(allprepulsefreqs);
stimparams.prepulsebandwidths=unique(allprepulsebandwidths);
stimparams.pulsedurs=unique(allpulsedurs);
stimparams.pulseamps=unique(allpulseamps);
stimparams.filter_operations=unique(allfilter_operations);
stimparams.center_frequencies=unique(allcenter_frequencies);
stimparams.lower_frequencies=unique(alllower_frequencies);
stimparams.upper_frequencies=unique(allupper_frequencies);
stimparams.allfields=unique(allfields);
stimparams.descriptions=unique(alldescriptions);


