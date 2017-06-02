%this is an example script that you can use as a template to process all
%cells in a cell list that you have built using CellListBuilder.
%mike 12-2016


%it's best to send output to the Data directory, since by default the
%Data directory is not backed up to GitHub. And you don't want a bunch of
%data plots on GitHub. We only want general-purpose code to be backed up to
%GitHub.

%we use postscript because it supports appending multiple figures, and you
%can then open it with a pdf browser such as Adobe Reader or Preview.
%(matlab can print directly to pdf format but then you can't append)


%delete the output file if it already exists

%the path to the data and the cell list will depend on the machine and platform you are
%running on. If you're processing data on a remote machine you probably
%need to mount it first (outside of matlab -i.e. in windows explorer or
%finder).
%cd('/Volumes/D/lab/djmaus/Data/lab') %this is mac format
cd('D:\lab\djmaus\Data\lab') %this is windows format
psfilename=['D:\lab\djmaus\Data\lab\',sprintf('plot_all_cells_%s.ps',datestr(now,'YY-mm-DD_hh-MM'))];
cell_list='ACIC1stPass.txt';


xlimits=[-300 300];

fid=fopen(cell_list);
fseek(fid, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid);
fseek(fid, 0, -1); %rewind to start
wb=waitbar(0, 'processing data');

%here you could fast forward to specific spot in your cell list.
%for example, you could put *** on a line by itself halfway through your
%cell list, and only process the second half of the cells
% line=fgetl(fid);
% while ~strcmp(line, '***')
%     line=fgetl(fid);
% end
i=1;
while 1 %processes until end of file is reached, then breaks
    tic;
    line=fgetl(fid);
    waitbar(ftell(fid)/filesize, wb);
    if  ~ischar(line), break, end %break at end of file
    while isempty(line)
        line=fgetl(fid);
    end
    if strncmp(line, 'cell',4)
        pathstr=fgetl(fid);
        filenamestr=fgetl(fid);
        clusterqualstr=fgetl(fid);
        pvstr=fgetl(fid);
        flashstr=fgetl(fid);

        datadir=strsplit(pathstr, ': ');
        datadir=datadir{2};
        filename=strsplit(filenamestr, ': ');
        filename=filename{2};

        %here you need to reformat the datapath to match what your local machine
        %expects. The datapath in the cell list can vary depending on which
        %machine was used when you built the cell list. It might take a little
        %trial and error to get the datapath to be recognized properly.
        %newdatadir=strrep(datadir, '\\Wehrrig1b', '/Volumes');
        %newdatadir=strrep(newdatadir, '\\wehrrig1b', '/Volumes');
        %newdatadir=strrep(datadir, 'D:', '/Volumes/D');
        %newdatadir=strrep(newdatadir, '\', '/'); %format for mac

        %fprintf('\n%s', datadir);
        %here you would use the appropriate function to process and/or plot your
        %data. You could just call a PlotXXX function, and it should automatically
        %process the data if necessary. Or you could explicitly call a ProcessXXX
        %function, for example if you want to force all cells to processed with a
        %certain xlimits, etc.
        %ProcessArchPVRev2(newdatadir, filename, xlimits)
        close all
        PlotFlashtrain_PSTH(datadir, filename(3),filename(18))
        if i~=1
           print('-f1',psfilename, '-dpsc2', '-append')
        else
           figure(2)
           orient tall
           print(psfilename, '-dpsc2', '-append')
           figure(1)
           orient tall
           print(psfilename,'-dpsc2', '-append')
        end
        i=i+1;
        %make sure you are in the right directory for the output file,
        %since PlotXXX or ProcessXXX will have taken you to the individual data
        %directory
        
        close all
    end

    %update waitbar with estimated time remaining
    etime=toc;
    fpos=ftell(fid)/filesize;
    fileremain=1-fpos;
    t_remain=etime*fileremain/fpos;
    str=sprintf('processing data... about %.1f minutes remaining', t_remain/60);
    waitbar(fpos,wb,  str)

end
fclose(fid) %close the output file
close(wb) %close the waitbar window
