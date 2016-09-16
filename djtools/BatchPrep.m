% function BatchPrep
% user options
figure
s_opt = [];

s_opt.mex_intersect=1; % use mex fn for poly intersect (much faster)
s_opt.dont_save_noise=1; % exclude noise from final output, not from the simpleclust file though

s_opt.auto_overlap = 1; % automatically loads other channels from same recording and computes spike overlap
s_opt.auto_overlap_dontask = 1; % dont ask if others should be loaded
s_opt.auto_overlap_max = 1; %if >0, limits how many other channels are loaded

s_opt.auto_noise = 1; % automatically assign channels with high overlap into noise cluster
s_opt.auto_noise_trs = 2; %proportion of channels a spike must co-occur in within .2ms in order to be classified noise

s_opt.auto_number = 1; % if set to 1, simpleclust will assume that there is ONLY ONE number in the MUA filenames and use is to designate the source channel for the resulting data

s_opt.invert= 1; % invert waveforms?

s_opt.skipevery_wf_display = 16; % skip every Nth waveform in the waveform display - this only sets an upper bound! smaller clusters will not be affected by this.

% specify what features to compute:
s_opt.features.pca=0;
s_opt.features.wavelet=0;
s_opt.features.nonlinear_energy=1;
s_opt.features.max_derivative=1;

%% init
run = 1;
dataloaded = 0;
s_opt.batch = 0;

global pref SP
try
    cd (pref.datapath)
    cd(SP.activedir)
end
[FileName,PathName,FilterIndex] = uigetfile({'*.wf;*.nse;*.nst;*.ntt;*.spikes;','all base electrode file types';...
    '*_simpleclust.mat', 'simpleclust file';'*.spikes', 'open ephys file';'*.mat', 'matlab file';...
    '*_extracted.mat', 'extracted matlab file';'*.wf','Waveform file';'*.nse' ,'neuralynx single electrode file';...
    '*.nst',  'neuralynx stereotrode file'; '*.ntt',  'neuralynx tetrode file'},['choose files to process'],'MultiSelect','on');


if FilterIndex(1)~=0
    for b=1:numel(FileName)
        
        features.muafile =[PathName,FileName{b}];
        
        fprintf('processing file %d of %d \n',b,numel(FileName));
        
        sc_load_mua_dialog;
        sc_save_dialog;
        
    end;
end;
dataloaded=0;