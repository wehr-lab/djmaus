%this is an example script that you can use as a template to process all
%cells in a cell list that you have built using CellListBuilder.
%mike 12-2016

% djhome %go to the djmaus root directory on your local machine
% cd Data
%it's best to send output to the Data directory, since by default the
%Data directory is not backed up to GitHub. And you don't want a bunch of
%data plots on GitHub. We only want general-purpose code to be backed up to
%GitHub.
psfilename='Alzheimers_behavior.ps';
%we use postscript because it supports appending multiple figures, and you
%can then open it with a pdf browser such as Adobe Reader or Preview.
%(matlab can print directly to pdf format but then you can't append)

delete(psfilename)
%delete the output file if it already exists

%the path to the data and the cell list will depend on the machine and platform you are
%running on. If you're processing data on a remote machine you probably
%need to mount it first (outside of matlab -i.e. in windows explorer or
%finder).
% cd('/Volumes/D/lab/djmaus/Data/lab') %this is mac format
% cd('\\WEHRRIG3\djmaus-data\ira') %this is windows format

    cell_list='behavior_list.m';


xlimits=[-100 150];

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

while 1 %processes until end of file is reached, then breaks
    tic;
    line=fgetl(fid);
    waitbar(ftell(fid)/filesize, wb);
    if  ~ischar(line), break, end %break at end of file
    while isempty(line)
        line=fgetl(fid);
    end
    if strcmp(line, 'cell')
        pathstr=fgetl(fid);
        filenamestr=fgetl(fid);
        groupstr=fgetl(fid);
        agestr=fgetl(fid);
        typestr=fgetl(fid);
       
        
        datadir=strsplit(pathstr, ': ');
        datadir=datadir{2};
        filename=strsplit(filenamestr, ': ');
        filename=filename{2};
        group=strsplit(groupstr, ': ');
        group=group{2};
        age=strsplit(agestr, ': ');
        age=age{2};
        type=strsplit(typestr, ': ');
        type=type{2};
        
        %here you need to reformat the datapath to match what your local machine
        %expects. The datapath in the cell list can vary depending on which
        %machine was used when you built the cell list. It might take a little
        %trial and error to get the datapath to be recognized properly.

        
        fprintf('\n%s', datadir);
        %here you would use the appropriate function to process and/or plot your
        %data. You could just call a PlotXXX function, and it should automatically
        %process the data if necessary. Or you could explicitly call a ProcessXXX
        %function, for example if you want to force all cells to processed with a
        %certain xlimits, etc.
        %ProcessArchPVRev2(newdatadir, filename, xlimits)
        close all
%          PlotSoundfile(newdatadir, filename(3), filename(18))
%         PlotTC_PSTH(newdatadir, filename(3), filename(18))
%         PlotSustainedSuppression(newdatadir, filename(3), filename(18))



GPIAS_Behavior_kip(datadir)



%make sure you are in the right directory for the output file,
%since PlotXXX or ProcessXXX will have taken you to the individual data
%directory
djhome
cd '\\WEHRRIG2b\c\lab\djmaus\Data\Kat'




figure(10);
print(psfilename,'-dpsc2', '-append');
figure(11);
print(psfilename,'-dpsc2', '-append');


cd(datadir)
outfilename=sprintf('outGPIAS_Behavior.mat');
d=dir(outfilename);
if ~isempty(d)
    load(outfilename)
else
    ProcessGPIAS_Behavior(datadir)
    load(outfilename);
end



% IL=out.IL; %whether there are any interleaved laser trials
% numpulseamps=out.numpulseamps;
% numgapdurs=out.numgapdurs;
% pulseamps=out.pulseamps;
% gapdurs=out.gapdurs;
% gapdelay=out.gapdelay;
% samprate=out.samprate; %in Hz
% mM1ON=out.mM1ON;
% mM1OFF=out.mM1OFF;
% M1ON=out.M1ON;
% M1OFF=out.M1OFF;
% nrepsON=out.nrepsON;
% nrepsOFF=out.nrepsOFF;
% soa=out.soa;
% isi=out.isi;
% soaflag=out.soaflag;
% PeakON=out.PeakON;
% PeakOFF=out.PeakOFF;
% mPeakON=out.mPeakON;
mPeakOFF=out.mPeakOFF;
% semPeakON=out.semPeakON;
% semPeakOFF=out.semPeakOFF;
% percentGPIAS_ON=out.percentGPIAS_ON;
percentGPIAS_OFF=out.percentGPIAS_OFF;
% pON = out.pON;
% pOFF = out.pOFF;
% samprate=out.samprate;
% M1ONstim=out.M1ONstim;
% M1OFFstim=out.M1OFFstim;
% mM1ONstim=out.mM1ONstim;
% mM1OFFstim=out.mM1OFFstim;
% xlimits=out.xlimits;

% %find optimal axis limits
% if ~isempty(mM1OFF)
%     ylimits(1)=min(mM1OFF(:));
%     ylimits(2)=max(mM1OFF(:));
% elseif ~isempty(mM1ON)
%     ylimits(1)=min(mM1ON(:));
%     ylimits(2)=max(mM1ON(:));
% end       
       


% datadir=strsplit(pathstr, ': ');
%         datadir=datadir{2};
%         filename=strsplit(filenamestr, ': ');
%         filename=filename{2};
%         group=strsplit(groupstr, ': ');
%         group=group{2};
%         age=strsplit(agestr, ': ');
%         age=age{2};
%         type=strsplit(typestr, ': ');
%         type=type{2};


x=length(percentGPIAS_OFF)
B=reshape(percentGPIAS_OFF, [x,1])


        cd(    '\\Wehrrig2b\c\lab\djmaus\Data\Kat\')

 
        
fid3=fopen('Alzheimers_Behavior_Data.mat', 'a');
fprintf(fid3, '\n%s %s\t %s\t %s\t %s\t', datadir, filename, group, age, type);
fprintf(fid3, '%.2f\t', percentGPIAS_OFF);
fprintf(fid3, '%.2f\t', mPeakOFF);




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
fclose(fid); %close the output file
fclose(fid3);
close(wb); %close the waitbar window