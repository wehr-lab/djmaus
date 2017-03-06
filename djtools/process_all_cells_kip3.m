% process_all_cells_kip3
%this is an example script that you can use as a template to process all
%cells in a cell list that you have built using CellListBuilder.
%mike 12-2016
% kip3 opens and displays data from cell_list1 and cell_list2


cd('Z:\lab\djmaus\Data\Kip')
path_home = pwd;
cell_list1='listKK01_pinp.txt';
cell_list2='listKK01_WNtc.txt';
cell_list3='listKK01_GPtc.txt';
cell_list4='listKK01_GPbehav.txt';

temp = [8  9    11    12    14    15    16    20    22    23    24   27    28    29    30    31    32 ...
    36    37    38    43    44    47    48    49    50    51    53    56    58 65    68    70    72 ...
    74    82    84    86    87    90    94    96   105   107  124   129   136   138   141   143   144  ...
    151   153   159   163   164   166   168   173   179   182   189   191   195   198   204   206   210 ...
    216   219   227   233   236   239   242   247   250   255   263   267   272   276   282   287   290 ...
    291   294   295   300   302   303   304   308   311   328   330   331   333   334   335  336 ];        %105  107 335 336 stimlog error (fixed)
PVlist = [ 16   20    22    23    24   27    28    29    30    31    32 ...
    36    37    38    43    44    47    48    49    50    51    53    56    58 65    68    70    72 ...
    74    82    84    86    87    90    94    96   105   107  124   129   136   138   141   143   144  ...
    151   153   159   163   164   166   168   173   179   182   189   191   195   198   204   206   210 ...
    216   219   227   233   236   239   242   247   250   255   263   267   272   276   282   287   290 ...
    291   294   295   300   302   303   304   308   311   328   330   331   333   334   335  336   337 ...
    338];
flag.IDnumber = PVlist;
flag.findSpecificCells = 1;
xlimits=[-350 350];

fid1=fopen(cell_list1,'r+');
fseek(fid1, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid1);
fseek(fid1, 0, -1); %rewind to start

fid2=fopen(cell_list2,'r+');
fseek(fid2, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid2);
fseek(fid2, 0, -1); %rewind to start

fid3=fopen(cell_list3,'r+');
fseek(fid3, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid3);
fseek(fid3, 0, -1); %rewind to start

fid4=fopen(cell_list4,'r+');
fseek(fid4, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid4);
fseek(fid4, 0, -1); %rewind to start


%here you could fast forward to specific spot in your cell list.

iPV=1;
iNonPV = 1;
IDnumber2 = 0;
IDnumber3 = 0;
IDnumber4 = 0;

while 1         % processes until end of file is reached, then breaks
    line=fgetl(fid1);
    if  ~ischar(line), break, end %break at end of file
    while isempty(line)
        line=fgetl(fid1);
    end
    if strfind(line,'cell ID#')
        clc
        %fprintf('\n%s\n',line);
        ind = findstr(line,'ID#');
        IDnumber = str2num(line(ind(end)+3:ind(end)+12));
        posPath = ftell(fid1);           pathstr=fgetl(fid1);
        posFN = ftell(fid1);             filenamestr=fgetl(fid1);
        posClustQual = ftell(fid1);      clusterqualstr=fgetl(fid1);
        posPV = ftell(fid1);             pvstr=fgetl(fid1);
        posProtocol = ftell(fid1);       protocolstr=fgetl(fid1);
        posNow = ftell(fid1);
        
        if flag.findSpecificCells
            
            if ismember(IDnumber,flag.IDnumber)        % process a particular cell
                fprintf('\n%s\n',pathstr);
                fprintf('\n%s\n',filenamestr);
                fprintf('\n%s\n',clusterqualstr);
                
                % plot the pinp response
                datadir=strsplit(pathstr, ': ');
                newdatadir=datadir{2};
                filename=strsplit(filenamestr, ': ');
                filename=filename{2};
                [h,p] = PlotPINP_PSTH_single(newdatadir, filename, xlimits);
                Handles.pinpLaserOFF_fig = figure(1);
                set(Handles.pinpLaserOFF_fig,'position',[10 660 560 420])
                Handles.pinpLaserON_fig = figure(2);
                set(Handles.pinpLaserON_fig,'position',[10 140 560 420])
                
                % characterize the pinp response histogram
                Handles.pinpLaserON_bar = findobj(Handles.pinpLaserON_fig,'BarWidth',1);
                XdataON = get(Handles.pinpLaserON_bar,'Xdata');
                YdataON = get(Handles.pinpLaserON_bar,'Ydata');
                Handles.pinpLaserOFF_bar = findobj(Handles.pinpLaserOFF_fig,'BarWidth',1);
                XdataOFF = get(Handles.pinpLaserOFF_bar,'Xdata');
                YdataOFF = get(Handles.pinpLaserOFF_bar,'Ydata');
                LaserStart = loadOutfileValue(datadir{2}, filename, xlimits,'stimlog(1).LaserStart');
                LaserWidth = loadOutfileValue(datadir{2}, filename, xlimits,'stimlog(1).LaserWidth');
                %spiketimes = loadOutfileValue(datadir{2}, filename, xlimits,'spiketimes');
                nbin = find(XdataON<=LaserStart,1,'last');
                nbin(2) = find(XdataON<=LaserStart+LaserWidth,1,'last');
                binwidth = LaserWidth/range(nbin);
                spont = round(mean(YdataOFF(nbin(1):nbin(2))),2);
                postLaser = round(mean(YdataON(nbin(2):nbin(2)+50/binwidth)),2);
                Xfast = interp(XdataON,10);
                Yfast = interp(YdataON,10); 
                [M,Mind] = max(Yfast - Yfast(nbin(2)*10));
                hMind = [Mind Mind+find((Yfast(Mind+1:nbin(2)*10) - Yfast(nbin(2)*10))>=M/exp(1))];
                tau = round(range(Xfast(hMind)),1);
                figure(Handles.pinpLaserON_fig)
                plot(Xfast(hMind),Yfast(hMind),'r','linewidth',3)
                text(XdataON(10),mean(ylim),['tau (adapt): ' num2str(tau) ' ms'],'fontsize',12)
                text(XdataON(10),mean(ylim)*.9,['spont: ' num2str(spont) ' postL: ' num2str(postLaser)],'fontsize',12)
                
                
                % load the waveforms for the PINP test
                str = [newdatadir '\' filename(1:end-2) '-wv.mat'];
                load(str)
                Y(1,:,:) = mWV';
                S(1,:,:) = sWV'; 
                
                % go to the same position in cell_list2 (if it exists)
                while IDnumber2 ~= IDnumber & IDnumber2 < IDnumber+1
                    line=fgetl(fid2);
                    if  ~ischar(line), break, end %break at end of file
                    while isempty(line)
                        line=fgetl(fid2);
                    end
                    if strfind(line,'cell ID#')
                        ind = findstr(line,'ID#');
                        IDnumber2 = str2num(line(ind(end)+3:ind(end)+12));
                    end
                end
                if IDnumber2 == IDnumber
                    posPath2 = ftell(fid2);           pathstr2=fgetl(fid2);
                    posFN2 = ftell(fid2);             filenamestr2=fgetl(fid2);
                    posClustQual2 = ftell(fid2);      clusterqualstr2=fgetl(fid2);
                    posPV2 = ftell(fid2);             pvstr2=fgetl(fid2);
                    posProtocol2 = ftell(fid2);       protocolstr2=fgetl(fid2);
                    posNow2 = ftell(fid2);
                    
                    fprintf('\n%s\n',line);
                    fprintf('\n%s\n',pathstr2);
                    fprintf('\n%s\n',filenamestr2);
                    fprintf('\n%s\n',clusterqualstr2);
                    % load the waveforms for WNtc
                    datadir=strsplit(pathstr2, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr2, ': ');
                    filename=filename{2}(1:end-2);
                    str = [newdatadir '\' filename '-wv.mat'];
                    load(str)
                    Y(2,:,:) = mWV';
                    S(2,:,:) = sWV';
                    
                    flag.WNtc_exists = 1;
                else
                    fprintf('\n%s\n','WNtc file does not exist');
                    flag.WNtc_exists = 0;
                end
                
                
                
                % go to the same position in cell_list3 (if it exists)
                while IDnumber3 ~= IDnumber & IDnumber3 < IDnumber+1
                    line=fgetl(fid3);
                    if  ~ischar(line), break, end %break at end of file
                    while isempty(line)
                        line=fgetl(fid3);
                    end
                    if strfind(line,'cell ID#')
                        ind = findstr(line,'ID#');
                        IDnumber3 = str2num(line(ind(end)+3:ind(end)+12));
                    end
                end
                if IDnumber3 == IDnumber
                    posPath3 = ftell(fid3);           pathstr3=fgetl(fid3);
                    posFN3 = ftell(fid3);             filenamestr3=fgetl(fid3);
                    posClustQual3 = ftell(fid3);      clusterqualstr3=fgetl(fid3);
                    posPV3 = ftell(fid3);             pvstr3=fgetl(fid3);
                    posProtocol3 = ftell(fid3);       protocolstr3=fgetl(fid3);
                    posNow3 = ftell(fid3);
                    
                    fprintf('\n%s\n',line);
                    fprintf('\n%s\n',pathstr3);
                    fprintf('\n%s\n',filenamestr3);
                    fprintf('\n%s\n',clusterqualstr3);
                    % load the waveforms for GAPtc
                    datadir=strsplit(pathstr3, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr3, ': ');
                    filename=filename{2}(1:end-2);
                    str = [newdatadir '\' filename '-wv.mat'];
                    load(str)
                    Y(3,:,:) = mWV';
                    S(3,:,:) = sWV';
                    
                    flag.GPtc_exists = 1;
                else
                    fprintf('\n%s\n','GPtc file does not exist');
                    flag.GPtc_exists = 0;
                end
                
                
                % go to the same position in cell_list4 (if it exists)
                while IDnumber4 ~= IDnumber & IDnumber4 < IDnumber+1
                    line=fgetl(fid4);
                    if  ~ischar(line), break, end %break at end of file
                    while isempty(line)
                        line=fgetl(fid4);
                    end
                    if strfind(line,'cell ID#')
                        ind = findstr(line,'ID#');
                        IDnumber4 = str2num(line(ind(end)+3:ind(end)+12));
                    end
                end
                if IDnumber4 == IDnumber
                    posPath4 = ftell(fid4);           pathstr4=fgetl(fid4);
                    posFN4 = ftell(fid4);             filenamestr4=fgetl(fid4);
                    posClustQual4 = ftell(fid4);      clusterqualstr4=fgetl(fid4);
                    posPV4 = ftell(fid4);             pvstr4=fgetl(fid4);
                    posProtocol4 = ftell(fid4);       protocolstr4=fgetl(fid4);
                    posNow4 = ftell(fid4);
                    
                    fprintf('\n%s\n',line);
                    fprintf('\n%s\n',pathstr4);
                    fprintf('\n%s\n',filenamestr4);
                    fprintf('\n%s\n',clusterqualstr4);
                    % load the waveforms for GAPbehav
                    datadir=strsplit(pathstr4, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr4, ': ');
                    filename=filename{2}(1:end-2);
                    str = [newdatadir '\' filename '-wv.mat'];
                    load(str)
                    Y(4,:,:) = mWV';
                    S(4,:,:) = sWV';
                    
                    flag.GPbehav_exists = 1;
                else
                    fprintf('\n%s\n','GPbehav file does not exist');
                    flag.GPbehav_exists = 0;
                end
                nTrodes = size(mWV,2);
                
                
                % plot the waveforms:
                Handles.waveforms_fig = figure('position',[10 100 600 1000]);hold on
                
                X = 0:1/30:39/30;
                Xfast = interp(X,10);
                % Y is the mean waveform as test#, trode#, timept
                for itrode = 1:nTrodes
                    % always plot the pinp waveform with limits and peak
                    Y(1,itrode,:) = Y(1,itrode,:)-Y(1,itrode,1);
                    subplot(nTrodes,1,itrode); hold on
                    plot(X,squeeze(Y(1,itrode,:)),'color',[0 .6 .6], 'linewidth',2)
                    plot(X,squeeze(Y(1,itrode,:))+squeeze(S(1,itrode,:)),'.')
                    plot(X,squeeze(Y(1,itrode,:))-squeeze(S(1,itrode,:)),'.')
                    yl(1,:) = ylim;
                    
                    if flag.WNtc_exists
                        Y(2,itrode,:) = Y(2,itrode,:)-Y(2,itrode,1);
                        plot(X,squeeze(Y(2,itrode,:)),'color',[.6 .6 0], 'linewidth',2)
                        yl(2,:) = ylim;
                    end
                    if flag.GPtc_exists
                        Y(3,itrode,:) = Y(3,itrode,:)-Y(3,itrode,1);
                        plot(X,squeeze(Y(3,itrode,:)),'color',[.6 .2 .6], 'linewidth',2)
                        yl(3,:) = ylim;
                    end
                    if flag.GPbehav_exists
                        Y(4,itrode,:) = Y(4,itrode,:)-Y(4,itrode,1);
                        plot(X,squeeze(Y(4,itrode,:)),'color',[.9 0 .9], 'linewidth',2)
                        yl(4,:) = ylim;
                    end
                    
                end
                
                Handles.waveforms_axes = get(Handles.waveforms_fig,'children');
                for itrode=1:nTrodes
                    yl(itrode,:)=Handles.waveforms_axes(itrode).YLim;
                end
                yl = [min(yl(:,1)) max(yl(:,2))];
                for itrode = 1:nTrodes
                    axes(Handles.waveforms_axes(nTrodes+1-itrode))
                    ylim(yl);
                    % try to get width at half-height
                    Yfast = interp(Y(1,itrode,:),10);
                    [M, Mind] = max(Yfast);
                    %hMind = find(Yfast>=M/2);
                    hMind = [Mind Mind+find(Yfast(Mind+1:end)>=M/exp(1))];
                    text(X(15),yl(2),['tau (repol): ' num2str(round(range(Xfast(hMind)),2)) ' ms'],'fontsize',12)
                    %text(X(15),yl(2),['spike width: ' num2str(round(range(Xfast(hMind)),2)) ' ms'],'fontsize',12)
                    
                    text(X(15),yl(2)+.1*diff(yl(1,:)),['cell# ' num2str(IDnumber) '  wire# ' num2str(itrode)],'fontsize',12)
                    if itrode == 1 
                        plot(Xfast(hMind),Yfast(hMind),'color',[.6 .6 .6],'linewidth',2)
                        text(X(35),yl(2), 'pinp','color',[0 .6 .6],'fontsize',12)
                        if flag.WNtc_exists
                            text(X(35),yl(2)*.8, 'WNtc','color',[.6 .6 0],'fontsize',12)
                        end
                        if flag.GPtc_exists
                            text(X(35),yl(2)*.6, 'GPtc','color',[.6 .2 .6],'fontsize',12)
                        end
                        if flag.GPbehav_exists
                            text(X(35),yl(2)*.4, 'GPbehav','color',[.9  0 .9],'fontsize',12)
                        end
                    end
                end
                
                % plot the ISI histogram (read the .t files)
                % bins = logspace(log10(.001),log10(1000),500);
                str = [pathstr(7:end) '\' filenamestr(11:end)];
                spiketimes=read_MClust_output(str);
                HistISI(spiketimes/10000);
                title('ISI histogram')
                hold on
                Handles.figISI = gcf;
                Handles.axesISI = gca;
                set(gcf,'position',[640 680 560 420])
                
                if flag.WNtc_exists
                    str = [pathstr2(7:end) '\' filenamestr2(11:end)];
                    spiketimes2=read_MClust_output(str);
                    figure(Handles.figISI);
                    hold on
                    [H, binsUsed] = HistISI(spiketimes2/10000,'axesHandle',Handles.axesISI,'myColor',[.6 .6 0]);
                    
                    %%%% plot tuning curve specified by cell_list2
                    datadir=strsplit(pathstr2, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr2, ': ');
                    filename=filename{2};
                    
                    PlotTC_PSTH_single(newdatadir, filename, xlimits);
                    Handles.WNtc_fig = gcf;
                    Handles.helperLine1 = annotation('line','position',[ 0.4383 0.9460 0 -0.8500]);    % just a helper line
                    set(Handles.WNtc_fig,'position',[860 90 600 1000])
                    Handles.WNtc_axes = get(Handles.WNtc_fig,'children');
                    axes(Handles.WNtc_axes(1));
                    xl = xlim; yl = ylim;
                    text(xl(2)+xl(1)*1.5,yl(1)-range(yl)/2,['WN tuning curve ID#' num2str(IDnumber)],'fontsize',14)
                    
                    PlotTC_PSTHendAligned_single(newdatadir, filename, xlimits);
                    Handles.WNtcEndAligned_fig = gcf;
                    Handles.helperLine2 = annotation('line','position',[ 0.4883 0.9460 0 -0.8500]);    % just a helper line
                    set(Handles.WNtcEndAligned_fig,'position',[1300 90 600 1000])
                    Handles.WNtcEA_axes = get(Handles.WNtcEndAligned_fig,'children');
                    axes(Handles.WNtcEA_axes(1));
                    xl = xlim; yl = ylim;
                    text(xl(2)+xl(1)*1.5,yl(1)-range(yl)/2,['WN TC - END Aligned ID#' num2str(IDnumber)],'fontsize',14)
                
                    analTemp1(newdatadir, filename, xlimits);
                end
                
                if flag.GPtc_exists
                    str = [pathstr3(7:end) '\' filenamestr3(11:end)];
                    spiketimes3=read_MClust_output(str);
                    figure(Handles.figISI);
                    hold on
                    [H, binsUsed] = HistISI(spiketimes3/10000,'axesHandle',Handles.axesISI,'myColor',[.6 .2 .6]);
                    
                    %%%% plot tuning curve specified by cell_list2
                    datadir=strsplit(pathstr3, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr3, ': ');
                    filename=filename{2};
                    
                    PlotGPIAS_PSTH_single(newdatadir, filename, xlimits);
                    Handles.GPIAStc_fig = gcf;
                    set(Handles.GPIAStc_fig,'position',[190 60 660 980])
                    Handles.GPIAStc_children = get(Handles.GPIAStc_fig,'Children');
                    axes(Handles.GPIAStc_children(1));
                    xl = xlim; yl = ylim;
                    text(xl(2)+xl(1)*1.5,yl(1)-range(yl)/2,['Gap with/out burst ID#' num2str(IDnumber)],'fontsize',14)
                    
                    analTemp1(newdatadir, filename, xlimits);                   
                end
                
                if flag.GPbehav_exists
                    str = [pathstr4(7:end) '\' filenamestr4(11:end)];
                    spiketimes4=read_MClust_output(str);
                    figure(Handles.figISI);
                    hold on
                    [H, binsUsed] = HistISI(spiketimes4/10000,'axesHandle',Handles.axesISI,'myColor',[.8 0 .8]);
                    
                    %%%% plot tuning curve (or whatever) specified by cell_list2
                    datadir=strsplit(pathstr4, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr4, ': ');
                    filename=filename{2};
                    
                    PlotGPIAS_PSTH_single(newdatadir, filename, xlimits);
                    Handles.GPIASbehav_fig = gcf;
                    set(Handles.GPIASbehav_fig,'position',[870 330 660 710])
                    Handles.GPIASbehav_children = get(Handles.GPIASbehav_fig,'Children');
                    axes(Handles.GPIASbehav_children(1));
                    xl = xlim; yl = ylim;
                    text(xl(2)+xl(1)*1.5,yl(1)-range(yl)/2,['Gap with burst  ID#' num2str(IDnumber)],'fontsize',14)
                                        
                    analTemp1(newdatadir, filename, xlimits);
                end
                
                
                fprintf('\n\n')
                input('<CR> to move on')
                close all
                cd(path_home)
            end         % ismember(IDnumber,flag.IDnumber)
            
        else
            if str2num(pvstr(end))      % process all PV cells
                datadir=strsplit(pathstr, ': ');
                newdatadir=datadir{2};
                filename=strsplit(filenamestr, ': ');
                filename=filename{2};
                
                fprintf('\n%s\n', newdatadir);
                close all
                %JustPlotIt([newdatadir '\' filename])
                %PlotTC_PSTH_single(newdatadir, filename, xlimits)
                [h,p] = PlotPINP_PSTH_single(newdatadir, filename, xlimits, []);

                if h
                    fprintf('\n%s\n',['# ' num2str(IDnumber) '   This file meets the hypothesis of a PV cell response']);
                    clear temp
                    %temp = input('\nEnter a <CR> to keep as is, or a number to change something');
                    temp =[];
                    if isempty(temp)
                        PVlist(iPV) = IDnumber;
                        iPV=iPV+1;
                    else
                        fprintf('\n%s\n',pathstr);
                        temp1 = input('\nEnter a <CR> to keep the path, or a number to change it');
                        if ~isempty(temp1)
                            str = '';
                            while length(str) ~= length(pathstr)
                                str = input(['\nEnter a new path with exactly the same length as the old path (' num2str(length(pathstr)) ')']);
                            end
                            fseek(fid1,posPath,-1);
                            fprintf(fid1,'%s',str);
                        end
                        fprintf('\n%s\n',filenamestr);
                        temp1 = input('\nEnter a <CR> to keep the FN, or a number to change it');
                        if ~isempty(temp1)
                            str = '';
                            while length(str) ~= length(filenamestr)
                                str = input(['\nEnter a new FN with exactly the same length as the old FN (' num2str(length(filenamestr)) ')']);
                            end
                            fseek(fid1,posFN,-1);
                            fprintf(fid1,'%s',str);
                        end
                        fprintf('\n%s\n',clusterqualstr);
                        temp1 = input('\nEnter a <CR> to keep the quality, or a number to change it');
                        if ~isempty(temp1)
                            str = '';
                            while length(str) ~= length(clusterqualstr)
                                str = input(['\nEnter a new quality with exactly the same length as the old quality (' num2str(length(clusterqualstr)) ')']);
                            end
                            fseek(fid1,posClustQual,-1);
                            fprintf(fid1,'%s',str);
                        end
                        fprintf('\n%s\n',pvstr);
                        temp1 = input('\nEnter a <CR> to keep the PV value, or a number to change it');
                        if ~isempty(temp1)
                            str = '';
                            while length(str) ~= length(pvstr)
                                str = input(['\nEnter a new PV value with exactly the same length as the old PV value (' num2str(length(pvstr)) ')']);
                            end
                            fseek(fid1,posPV,-1);
                            fprintf(fid1,'%s',str);
                        end
                        fprintf('\n%s\n',protocolstr);
                        temp1 = input('\nEnter a <CR> to keep the Protocol, or a number to change it');
                        if ~isempty(temp1)
                            str = '';
                            while length(str) ~= length(protocolstr)
                                str = input(['\nEnter a new Protocol with exactly the same length as the old Protocol (' num2str(length(protocolstr)) ')']);
                            end
                            fseek(fid1,posProtocol,-1);
                            fprintf(fid1,'%s',str);
                        end
                    end
                else
                    clear temp1
                    temp1 = input(['\n\n # ' num2str(IDnumber) '   WARNING - NOT A PV CELL. Enter a <CR> to make PV==0']);
                    if isempty(temp1)
                        pvstr(end)=0;
                        fseek(fid1,posPV,-1);
                        fprintf(fid1,'%s',pvstr);
                        nonPV(iNonPV) = IDnumber;
                        iNonPV=iNonPV+1;
                    end
                end         % if h
                fseek(fid1,posNow,-1);
            end             % if str2num(pvstr(end))
            close all
        end                 % if flag.findSpecificCells (else)
    end                     % if strfind(line,'cell ID#')
    
end             % while
fclose(fid1); %close the output file
fclose(fid2); %close the output file
fclose(fid3); %close the output file
fclose(fid4); %close the output file