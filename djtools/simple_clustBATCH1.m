%{
   -----------------------------------------------------------------------
  
    Simple Clust v0.5
   
    (c) jan. 2012, Jakob Voigts (jvoigts@mit.edu)
 
    Beta version, use at your own risk

    -----------------------------------------------------------------------

    This is a program for manual clustering of spikes in matlab.
    I dont recommend this software for scientific use by
    anyone without full understanding of the methods.

    Please post any issues you encounter, or any improvements or additions
    to github at:

    http://github.com/moorelab/simpleclust

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%}
%%%%
% Kip's batch processing version. Other BS commented out or eliminated

global debugstate;

PathStart = pwd;
%cd 'C:\Kip\Wehr_data';


% user options

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

features.updatezoom=1; % so that new features are checked for zoom state and visibility, set 1 after each zoom operation or new feature

global debugstate
debugstate = 0; % 0: do nothing, 1: go trough following states
debuginput = [0 0 0];

%addpath(pwd); % warn instead
%addpath(fullfile(pwd,'read_cheetah'));

if numel(strfind(path,'read_cheetah')) ==0
    error('make sure the read_cheetah dir is in your matlab path');
end

if s_opt.mex_intersect
    if numel(strfind(path,'InPolygon-MEX')) ==0
        error('make sure the InPolygon-MEX dir is in your matlab path, or disable s_opt.mex_intersect');
    end
end

%% main loop

%while run

F.Xind_ = 20;   F.Xfig_ht = 500; F.Xfig_wd = 400;
H.fig = figure('Units','pixels',...
    'Position',[50 100 F.Xfig_wd F.Xfig_ht],...
    'Name','simple clust v0.5',...
    'NumberTitle','off');
F.Xtot_ht = 0;
F.Xtot_ht = F.Xtot_ht + F.Xind_;
uicontrol('Parent',H.fig,...
    'Style','text',...
    'Units','pixels',...
    'Position',[0 F.Xfig_ht-F.Xtot_ht F.Xfig_wd F.Xind_],...
    'ForegroundColor','green',...
    'BackgroundColor','black',...
    'FontWeight','bold',...
    'FontSize',12,...
    'String', 'batch simple clust v0.5');
%Test types: popupmenu
F.Xtot_ht = F.Xtot_ht + F.Xind_ + 30;
% uicontrol('Parent',H.fig,...
%     'Style','text',...
%     'Units','pixels',...
%     'Position',[0 F.Xfig_ht-F.Xtot_ht 175 F.Xind_],...
%     'ForegroundColor','blue',...
%     'FontWeight','bold',...
%     'String', 'Choose action');
% F.Xtot_ht = F.Xtot_ht + F.Xind_;
% H.action = uicontrol('Parent',H.fig,...
%     'Style','popup',...
%     'Units','pixels',...
%     'Position',[0 F.Xfig_ht-F.Xtot_ht 175 F.Xind_],...
%     'BackgroundColor','white',...
%     'value', 1,...
%     'String', ('choose|Open|save / exit|batch prep|batch run'));
%run button
H.go_pb = uicontrol('Style','togglebutton',...
    'Units','pixels',...
    'Position',[201 F.Xfig_ht-F.Xtot_ht 105 F.Xind_],...
    'BackgroundColor','green',...
    'ForegroundColor','yellow',...
    'FontWeight','bold',...
    'String', 'go');

F.Xtot_ht = F.Xtot_ht + F.Xind_;
%run button
H.exit_pb = uicontrol('Style','togglebutton',...
    'Units','pixels',...
    'Position',[301 F.Xfig_ht-F.Xtot_ht 105 F.Xind_],...
    'BackgroundColor','red',...
    'ForegroundColor','yellow',...
    'FontWeight','bold',...
    'String', 'exit');



%while ~get(H.go_pb,'value')
%drawnow
%end
while ~get(H.exit_pb,'value')
    drawnow
    if dataloaded
        
%         features = sc_plotclusters(features);
%         
%         sc_plotfeatureselection(features);
%         
%         sc_plotclusterselection(features);
%         
%         sc_plotallclusters(features,mua,s_opt);
%         
%         sc_extramenu(features);
%         
%         sc_timeline(features,mua,x,y,b);
%         
%         features.highlight = 0; % remove highlight on each click
    else
%         switch get(H.action,'value')
%             case 2
%                 button='Open';
%                 % load MUa data
%                 if debugstate > 0
%                     PathName = '/home/jvoigts/Dropbox/em003/good/';
%                     FileName =  'ST11.nse';
%                 else
%                     [FileName,PathName,FilterIndex] = uigetfile({'*.wf;*.nse;*.nst;*.ntt;*.spikes;','all base electrode file types';'*_simpleclust.mat', 'simpleclust file';'*.mat', 'matlab file';'*.wf','Waveform file';'*.nse' ,'neuralynx single electrode file'; '*.nst',  'neuralynx stereotrode file'; '*.ntt',  'neuralynx tetrode file'},'choose input file');
%                 end
%                 features.muafile =[PathName,FileName];
%                 sc_load_mua_dialog;
%                 
%             case  3
%                 clf; drawnow;
%                 run=0;
%                 
%             case 4
                [FileName,PathName,FilterIndex] = uigetfile({'*.wf;*.nse;*.nst;*.ntt;*.spikes;','all base electrode file types';...
                    '*_simpleclust.mat', 'simpleclust file';'*.spikes', 'open ephys file';'*.mat', 'matlab file';...
                    '*_extracted.mat', 'extracted matlab file';'*.wf','Waveform file';'*.nse' ,'neuralynx single electrode file';...
                    '*.nst',  'neuralynx stereotrode file'; '*.ntt',  'neuralynx tetrode file'},'choose files to process','MultiSelect','on');
                
                if FilterIndex(1)~=0
                    for b=1:numel(FileName)
                        figure
                        features.muafile =[PathName,FileName{b}];
                        
                        fprintf('processing file %d of %d \n',b,numel(FileName));
                        
                        sc_load_mua_dialog;
                        sc_save_dialog;
                        close
                    end
                end
                %dataloaded=0;
%             case 5
%                 [multifiles,PathName,FilterIndex] = uigetfile({'*.wf;*.nse;*.nst;*.ntt;*.spikes;','all base electrode file types';...
%                     '*_simpleclust.mat', 'simpleclust file';'*.spikes', 'open ephys file';'*.mat', 'matlab file';...
%                     '*_extracted.mat', 'extracted matlab file';'*.wf','Waveform file';'*.nse' ,'neuralynx single electrode file';...
%                     '*.nst',  'neuralynx stereotrode file'; '*.ntt',  'neuralynx tetrode file'},['choose files to cluster'],'MultiSelect','on');
%                 
%                 if FilterIndex(1)~=0
%                     s_opt.batch=1; % indicate we're doing a batch
%                     multi_N=1; % cycle through many files
%                     
%                     features.muafile =[PathName,multifiles{multi_N}];
%                     sc_load_mua_dialog;
%                 else
%                     dataloaded=0;
%                 end
%         end
    end
    %    close
end                     % while ~exit

close
eval(['cd ' PathStart])