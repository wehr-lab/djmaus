% process_GPtcCellList - based on process_all_cells_kip3, 
% which is loosely based on process_all_cells by MIke
% processes and plots some summary stats for GPtc tests from a GPtc-only cell list
% calls analTemp3 to actually do the analyses

global CELL_STATS
global ProcessedCellNum

cd('D:\lab\djmaus\Data\apw')
path_home = pwd;
binwidth=1;
cell_list1='80dbnoise.txt';

%%%%%%%%
flag.minQual = 3;
flag.findSpecificCells = 1;
flag.plot(1:5)=0;       % 1)pinp  2) WNtc 3)GPtc  4) GPbehav 5)ISI histograms (required for some stats)
xlimits=[-350 350];

fid1=fopen(cell_list1,'r+');
while 1         % processes until end of file is reached, then breaks
    line=fgetl(fid1);
    if  ~ischar(line), break, end %break at end of file
    while isempty(line)
        line=fgetl(fid1);
    end
    if strfind(line,'cell ID#')
        ind = findstr(line,'ID#');
        temp = str2num(line(ind(end)+3:ind(end)+12));
    end
end
nCells=temp;
fprintf('\nFound %d cells in pinp list\n',nCells)
fseek(fid1, 0, 1); %fastforward to end, to get file size
filesize(1)=ftell(fid1);
fseek(fid1, 0, -1); %rewind to start
flag.IDnumber = 1:nCells;

IDnumber = 0;
CELL_STATS = struct;
ProcessedCellNum = 0;

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
        posPV = ftell(fid1);             pvstr=fgetl(fid1); ispv = str2num(pvstr(end));
        posProtocol = ftell(fid1);       protocolstr=fgetl(fid1);
        posNow = ftell(fid1);
        
        if flag.findSpecificCells & str2num(clusterqualstr(end))>= flag.minQual
            
            if ismember(IDnumber,flag.IDnumber)        % process a particular cell
                fprintf('\nIDnumber: %d\n',IDnumber);
                fprintf('\n%s\n',pathstr);
                fprintf('\n%s\n',filenamestr);
                fprintf('\n%s\n',clusterqualstr);
                
                ProcessedCellNum = ProcessedCellNum+1;
                CELL_STATS(ProcessedCellNum).IDnumber = IDnumber;
                ind = strfind(pathstr,'mouse-');
                if ~isempty(ind)
                    CELL_STATS(ProcessedCellNum).mouseID = pathstr(ind+6:end);
                else
                    CELL_STATS(ProcessedCellNum).mouseID = 'unknown';
                end
                CELL_STATS(ProcessedCellNum).quality = str2num(clusterqualstr(end));
                
                
                    
                    %%%% plot tuning curve specified by cell_list3
                    datadir=strsplit(pathstr, ': ');
                    newdatadir=datadir{2};
                    filename=strsplit(filenamestr, ': ');
                    filename=filename{2};
                    try
                        if flag.plot(1)
                            PlotGPIAS_PSTH_single(newdatadir, filename, xlimits);
                            Handles.GPIAStc_fig = gcf;
                            set(Handles.GPIAStc_fig,'position',[190 60 660 980])
                            Handles.GPIAStc_children = get(Handles.GPIAStc_fig,'Children');
                            axes(Handles.GPIAStc_children(1));
                            xl = xlim; yl = ylim;
                            text(xl(2)+xl(1)*1.5,yl(1)-range(yl)/2,['Gap with/out burst ID#' num2str(IDnumber)],'fontsize',14)
                        end
                        Type = 'GPIAS w/o burst';
                        analTemp3(newdatadir, filename, xlimits,[],[],Type);
                    catch
                        fprintf('\nThere was an error with %s. Moving on..?\n',filename);
                        pause
                    end
                end         % ismember(IDnumber,flag.IDnumber)
                                
                %fprintf('\n\n')
                %input('<CR> to move on')
                close all
                cd(path_home)
            
        end                 % if flag.findSpecificCells (else)
    end                     % if strfind(line,'cell ID#')
    
end             % while
fclose(fid1); %close the input file


% plot something
nCells = length(CELL_STATS);
gapdurs = CELL_STATS(1).GPtc.gapdurs;
ngapDurs = length(gapdurs);
H_PGIminDur_nonPV = zeros(1,ngapDurs);
H_PGIbestDur_nonPV = zeros(1,ngapDurs);

for icell = 1:nCells
    if ~isempty(CELL_STATS(icell).GPtc)        
        if ~isempty(CELL_STATS(icell).GPtc.PGI_minDur)
            minDur = CELL_STATS(icell).GPtc.PGI_minDur;
            if CELL_STATS(icell).GPtc.PGI_respSIGN(minDur)
                H_PGIminDur_nonPV(minDur) = H_PGIminDur_nonPV(minDur)+1;
            end
            bestDur = CELL_STATS(icell).GPtc.PGI_bestDur;
            if CELL_STATS(icell).GPtc.PGI_respSIGN(bestDur)
                H_PGIbestDur_nonPV(bestDur) = H_PGIbestDur_nonPV(bestDur)+1;
            end
        end
    end
end

figure; hold on
bar(log10(gapdurs),H_PGIbestDur_nonPV)
set(gca,'XTick',log10(gapdurs))
set(gca,'XTicklabel',gapdurs)
title('GTR    best Dur')
xlabel('GAP duration (ms)')
ylabel('#cells')

H = get(gcf,'children');
H1 = get(H,'children');
Xdata = get(H1,'XData');
Ydata = get(H1,'YData');
nCells = sum(Ydata);
for i = 1:length(Ydata)
    text(Xdata(i)*.95,Ydata(i)+1,[num2str(round(100*Ydata(i)/nCells)) '%'])
end
