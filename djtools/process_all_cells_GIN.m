%this is an example script that you can use as a template to process all
%cells in a cell list that you have built using CellListBuilder.
%mike 12-2016

% cd ('\\wehrrig1b\D\lab\djmaus\data\apw') %why doesn't this work??
% delete('list of cells.txt')

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
% cd('/Volumes/D/lab/djmaus/Data/lab') %this is mac format
cd('\\wehrrig1b\D\lab\djmaus\Data\Emily') %this is windows format
psfilename=['Test',sprintf('Plot%s.ps',datestr(now,'YY-mm-DD_hh-MM'))];


% cell_list1='testGIN.txt';
% cell_list2='testPINP.txt';

cell_list1='VIPCHR2GIN.txt';
cell_list2='VIPCHR2PINP.txt';


xlimits=[-300 300];

fid1=fopen(cell_list1);
fseek(fid1, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid1);
fseek(fid1, 0, -1); %rewind to start

fid2=fopen(cell_list2);
fseek(fid2, 0, 1); %fastforward to end, to get file size
filesize=ftell(fid2);
fseek(fid2, 0, -1); %rewind to start
wb=waitbar(0, 'processing data');

%here you could fast forward to specific spot in your cell list.
%for example, you could put *** on a line by itself halfway through your
%cell list, and only process the second half of the cells
% line=fgetl(fid2);
% while ~strcmp(line, '***')
%     line=fgetl(fid2);
%     waitbar(ftell(fid2)/filesize, wb);
% end



i = 1;

while 1 %processes until end of file is reached, then breaks
    tic;
    line1=fgetl(fid1);
    waitbar(ftell(fid1)/filesize, wb);
    line2=fgetl(fid2);
    waitbar(ftell(fid2)/filesize, wb);
    
    if  ~ischar(line1), break, end %break at end of file
    while isempty(line1)
        line1=fgetl(fid1);
    end
    
    if  ~ischar(line2), break, end %break at end of file
    while isempty(line2)
        line2=fgetl(fid2);
    end
    
    if strncmp(line1, 'cell', 4)
        pathstr=fgetl(fid1);
        filenamestr=fgetl(fid1);
        clusterqualstr=fgetl(fid1);
        pvstr=fgetl(fid1);
        archstr=fgetl(fid1);
        
        datadir=strsplit(pathstr, ': ');
        datadir=datadir{2};
        filename=strsplit(filenamestr, ': ');
        filename=filename{2};
        
        pathstr2=fgetl(fid2);
        filenamestr2=fgetl(fid2);
        clusterqualstr2=fgetl(fid2);
        pvstr2=fgetl(fid2);
        archstr2=fgetl(fid2);
        
        datadir2=strsplit(pathstr2, ': ');
        datadir2=datadir2{2};
        filename2=strsplit(filenamestr2, ': ');
        filename2=filename2{2};
        
        %here you need to reformat the datapath to match what your local machine
        %expects. The datapath in the cell list can vary depending on which
        %machine was used when you built the cell list. It might take a little
        %trial and error to get the datapath to be recognized properly.
        %newdatadir=strrep(datadir, '\\Wehrrig1b', '/Volumes');
        %newdatadir=strrep(newdatadir, '\\wehrrig1b', '/Volumes');
        %newdatadir=strrep(datadir, 'D:', '/Volumes/D');
        %newdatadir=strrep(newdatadir, '\', '/'); %format for mac
        newdatadir1=strrep(datadir, 'D:' , '\\wehrrig1b\D')
        newdatadir2=strrep(datadir2, 'D:' , '\\wehrrig1b\D');
        
       % fprintf('\n%s', newdatadir);
       
        %here you would use the appropriate function to process and/or plot your
        %data. You could just call a PlotXXX function, and it should automatically
        %process the data if necessary. Or you could explicitly call a ProcessXXX
        %function, for example if you want to force all cells to processed with a
        %certain xlimits, etc.
        %ProcessArchPVRev2(newdatadir, filename, xlimits)
        close all
        
        %PlotGPIAS_PSTH(newdatadir, filename(3), filename(18))
        
        
        
        
        cd(newdatadir1)
        [p,f,ext]=fileparts(filename);
        split=strsplit(f, '_');
        ch=strsplit(split{1}, 'ch');
        channel=str2num(ch{2});
        clust=str2num(split{end});
        outfilename1=sprintf('outPSTH_ch%dc%d.mat',channel, clust) ;
        
        fprintf('\n%s', outfilename1')
        load(outfilename1)
        
        
        
        GTRsOFF=zeros(out.numgapdurs, max(out.nrepsOFF));
        GTRsON=zeros(out.numgapdurs, max(out.nrepsON));
        
        for gd=1:out.numgapdurs
            for rep=1:out.nrepsOFF(gd)
                stOFF=out.M1OFF(gd, 1, rep).spiketimes;
                GTRspikecountOFF=length(find(stOFF>0 & stOFF<51));
                GTRsOFF(gd, rep)=GTRspikecountOFF;
                %                 OFFSETspikecountOFF=length(find(stOFF(0-out.gapdurs) & stOFF<0));
                %                 OFFSETsOFF(gd, rep)=OFFSETspikecountOFF
            end
            for rep=1:out.nrepsON(gd)
                stON=out.M1ON(gd, 1, rep).spiketimes;
                GTRspikecountON=length(find(stON>0 & stON<51));
                GTRsON(gd, rep)=GTRspikecountON;
            end
        end
        
        gap0=GTRsOFF(1, :);
        K1=[];
        K2=[];
        K3=[];
        K4=[];
        
        for gd=2:out.numgapdurs
            [h3, p3]=ttest(GTRsOFF(gd, :), gap0, [], 'right'); %  "right" is for GTRs>Gap0
%            fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p3);
            if isnan(h3)
                h3=0;
            end
            H3(gd)=h3;
            if h3
                K1 = 'GTR';
                K2 = [K2 out.gapdurs(gd)];
            end
        end
        
        
        
        for gd=2:out.numgapdurs
            [h4, p4]=ttest(GTRsOFF(gd, :), gap0, [], 'left'); %  "left" is for GTRs<Gap0
%            fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p4);
            if isnan(h4)
                h4=0;
            end
            H4(gd)=h4;
            if h4
                K3 = 'RDPGI';
                K4 = [K4 out.gapdurs(gd)];
            end
        end

        
%         for gd=1:out.numgapdurs
%             [h, p]=ttest(GTRsON(gd, :),GTRsOFF(gd, :), [], 'left'); %"right" is for Laser RI
%                  fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p);
%             if isnan(h) h=0;end
%             H(gd)=h;
%             if h
%                 
%                 K=[K out.gapdurs(gd)];
%                
%             end
%         end
        
        cd(newdatadir2)
        [pname,f,ext]=fileparts(filename2);
        split=strsplit(f, '_');
        ch=strsplit(split{1}, 'ch');
        channel=str2num(ch{2});
        clust=str2num(split{end});
        outfilename2=sprintf('outPSTH_ch%dc%d.mat',channel, clust) ;
        
        K5 = [];
        K6 = [];
        
        fprintf('\n%s', outfilename2')
        load(outfilename2);
        
%  if strcmp(outfilename2, 'outPSTH_ch5c1.mat')
%      keyboard
%  end

        for pwindex=1:out.numpulsewidths
            
            try
                binwidth=varargin{5};
            catch
                binwidth=5;
            end
            
            laserstarts=out.laserstarts;
            numlaserstarts=out.numlaserstarts;
            if numlaserstarts>1 warning('\nusing only minimum laserstart for t-test');end
            laserstart=min(laserstarts);
            ssdindex=1;
                stop=laserstart+75;
            

            

            clear ON OFF
            %How about a t-test for effect of Laser in first 75 ms
            for rep=1:out.nrepsPulse(pwindex)
                st=out.MPulse(pwindex,rep).spiketimes;
                spiketimes=st(st>(laserstart+20) & st<stop); % spiketimes in region
                ON(rep)=length(spiketimes);
            end
            for rep=1:out.nrepsOFF(ssdindex)
                st=out.MSilentSoundOFF(rep).spiketimes;
                spiketimes=st(st>(laserstart+20) & st<stop); % spiketimes in region
                OFF(rep)=length(spiketimes);
            end
            [h1,p1]=ttest2(ON, OFF, 'tail', 'right');
           % fprintf('\n p=%.4f', p1);
            if isnan(h1) h1=0;end
            if h1
                K5 = [K5 'Delayed'];
            else
                K5 = [K5 ''];
            end
            
            [h2,p2]=ttest2(ON, OFF, 'tail', 'left');
           % fprintf('\n p=%.4f', p2);
            if isnan(h2) h2=0;end
            if h2
                K6 = [K6 'Suppd'];
            else
                K6 = [K6 'nada  '];
            end
        end
        
        if length(ON) ~=  out.nrepsPulse(pwindex)
            error('sanity check failure')
        end

        cd(    '\\wehrrig1b\D\lab\djmaus\Data\Emily')
        fid3=fopen('SummaryPINPwin20to75.txt', 'a');
        fprintf(fid3, '\n%s %s %s ', newdatadir1, outfilename1, outfilename2);
        

        fprintf(fid3, '%s ', K5);
        fprintf(fid3, '%s ', K6);
        fprintf(fid3, '%s ', K1);
        fprintf(fid3, '%d ', K2);
        fprintf(fid3, '%s ', K3);
        fprintf(fid3, '%d ', K4);
        
       fclose(fid3);

%         cd ('\\wehrrig1b\D\lab\djmaus\Data\Emily')
%         
%        
%        for f=2
%             figure(f)
%             orient tall
%             print(psfilename, '-dpsc2', '-append')
%         end
%         for f=1
%             figure(f)
%             orient tall
%             print(psfilename, '-dpsc2', '-append')
%         end
%         i=i+1;
        
        %make sure you are in the right directory for the output file,
        %since PlotXXX or ProcessXXX will have taken you to the individual data
        %directory
        
        
        close all
    end
    
    %update waitbar with estimated time remaining
    etime=toc;
    fpos=ftell(fid1)/filesize;
    fileremain=1-fpos;
    t_remain=etime*fileremain/fpos;
    str=sprintf('processing data... about %.1f minutes remaining', t_remain/60);
    waitbar(fpos,wb,  str)
    
end
fclose(fid1); %close the output file
close(wb) %close the waitbar window
