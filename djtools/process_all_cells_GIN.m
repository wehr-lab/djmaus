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
% cd('\\wehrrig1b\D\lab\djmaus\Data\Emily') %this is windows format
% psfilename=['Test',sprintf('Plot%s.ps',datestr(now,'YY-mm-DD_hh-MM'))];

% 
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
        clustqual=strsplit(clusterqualstr,':')
        clustqual=clustqual{2}
        
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
      K9=[];
      K10=[];
      GiNBL=[];
      GiNBLHZ=[];
      repGiN=0;
      
        for gd=1:out.numgapdurs
            gapdurs=1:out.gapdurs;
            for repGiN=1:out.nrepsOFF(gd)
                stOFF=out.M1OFF(gd, 1, repGiN).spiketimes;
                GTRspikecountOFF=length(find(stOFF>0 & stOFF<76));
                GTRsOFF(gd, repGiN)=GTRspikecountOFF;
%                                 OFFSETspikecountOFF=length(find(stOFF>0-gapdurs & stOFF<0));
%                                 OFFSETsOFF(gd, rep)=OFFSETspikecountOFF
            end
            
            %             test(gd)=sum(OFFSETsOFF(gd, :))
            GTRspikesOFF(gd)= sum(GTRsOFF(gd, :)); %use 'k' for black
            
             K9 = [GTRspikesOFF out.gapdurs(gd)];
            
            for repGiN=1:out.nrepsON(gd)
                stON=out.M1ON(gd, 1, repGiN).spiketimes;
                GTRspikecountON=length(find(stON>0 & stON<76));
                GTRsON(gd, repGiN)=GTRspikecountON;
            end
            
            GTRspikesON(gd)= sum(GTRsON(gd, :));  %use 'b' for blue
              K10 = [GTRspikesON out.gapdurs(gd)];
        end
        
        for gd=1
            GiNBL=sum(GTRsOFF(gd,:));
        end
               
        if GiNBL<0
            GiNBLHZ=0
        else
            GiNBLHZ = (GiNBL*(1000/75))/repGiN
        end
        
        %         plot(xaxis, y1, 'k', xaxis, y2, 'b')% define xaxis, the assign y1
%         and y2 for GTRspikesOFF and GTRspikesON with related colors
        
        gap0OFF=GTRsOFF(1, :);
        gap0ON=GTRsON(1, :);
        
        K3=[];
        K3a=[];
        K4=[];
        K4a=[];
        K5=[];
        K5a=[];
        K6=[];
        K6a=[];
        K7=[];
        K7a=[];
        K8=[];
        K8a=[];
        
        H3=[];
        H4=[];
        H5=[];
        H6=[];
        H7=[];
        H8=[];
        for gd=2:out.numgapdurs
            [h3, p3]=ttest(GTRsOFF(gd, :), gap0OFF, [], 'right'); %  "right" is for GTRs>Gap0
%                        fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p3);
            if isnan(h3)
                h3=0;
            end
            H3(gd)=h3;
%                         if h3
%                             K3 = 'GTR';
%                             K3a = [K3a out.gapdurs(gd)];
%                         end
        end
        
        for gd=2:out.numgapdurs-1
            if H3(gd) & H3(gd+1)
                K3='GTR';
            end
        end
        
        
        for gd=2:out.numgapdurs
            [h4, p4]=ttest(GTRsOFF(gd, :), gap0OFF, [], 'left'); %  "left" is for GTRs<Gap0
            %            fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p4);
            if isnan(h4)
                h4=0;
            end
            H4(gd)=h4;
            %             if h4
            %                 K4 = 'RD';
            %                 K4a = [K6 out.gapdurs(gd)];
            %             end
        end
        
        for gd=2:out.numgapdurs-1
            if H4(gd) & H4(gd+1)
                K4='RD';
            end
        end
        
        for gd=2:out.numgapdurs
            [h5, p5]=ttest(GTRsON(gd, :), gap0ON, [], 'right'); %  "right" is for GTRs>Gap0
%                        fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p3);
            if isnan(h5)
                h5=0;
            end
            H5(gd)=h5;
            %             if h5
            %                 K5 = 'GTR_ON';
            %                 K5a = [K6 out.gapdurs(gd)];
            %             end
        end
        
        for gd=2:out.numgapdurs-1
            if H5(gd) & H5(gd+1)
                K5='GTR_ON';
            end
        end
        
        
        for gd=2:out.numgapdurs
            [h6, p6]=ttest(GTRsON(gd, :), gap0ON, [], 'left'); %  "left" is for GTRs<Gap0
            %            fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p4);
            if isnan(h6)
                h6=0;
            end
            H6(gd)=h6;
            %             if h6
            %                 K6 = 'RD_ON';
            %                 K6a = [K10 out.gapdurs(gd)];
            %             end
        end
        
        for gd=2:out.numgapdurs-1
            if H6(gd) & H6(gd+1)
                K6='RD_ON';
            end
        end
        
        
        for gd=1:out.numgapdurs
            [h7, p7]=ttest(GTRsON(gd, :),GTRsOFF(gd, :), [], 'right'); %"right" is for ON > OFF
            %                  fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p);
            if isnan(h7) h7=0;end
            H7(gd)=h7;
            if h7
                K7= 'L_INC';
                K7a=[K7a out.gapdurs(gd)];
            end
        end
        
        for gd=1:out.numgapdurs
            [h8, p8]=ttest(GTRsON(gd, :),GTRsOFF(gd, :), [], 'left'); %"right" is for ON > OFF
            %                  fprintf('\ngap %d: p=%.4f', out.gapdurs(gd), p);
            if isnan(h8) h8=0;end
            H8(gd)=h8;
            if h8
                K8= 'L_DEC';
                K8a=[K8a out.gapdurs(gd)];
            
            end
        end
        
        
       
        
        cd(newdatadir2)
        [pname,f,ext]=fileparts(filename2);
        split=strsplit(f, '_');
        ch=strsplit(split{1}, 'ch');
        channel=str2num(ch{2});
        clust=str2num(split{end});
        outfilename2=sprintf('outPSTH_ch%dc%d.mat',channel, clust) ;
        
        K1 = [];
        K1a = [];
        K2 = [];
       
        nmax=0;
        xmax=0;
       
        
        fprintf('\n%s', outfilename2')
        load(outfilename2);
        
        %  if strcmp(outfilename2, 'outPSTH_ch5c1.mat')
        %      keyboard
        %  end
        
        for pwindex=1:out.numpulsewidths
            
            
            binwidth=5;
            
            
            laserstarts=out.laserstarts;
            numlaserstarts=out.numlaserstarts;
            if numlaserstarts>1 warning('\nusing only minimum laserstart for t-test');end
            laserstart=min(laserstarts);
            ssdindex=1;
            stop=laserstart+75;
            
            
           repPNP = 0;
           totBL=0;
           meanBLHZ=0;
            
            clear ON OFF
            %t-test for effect of Laser in first 75 ms
            for repPNP=1:out.nrepsPulse(pwindex)
                st=out.MPulse(pwindex,repPNP).spiketimes;
                spiketimes=st(st>laserstart & st<stop); % spiketimes in region
                ON(repPNP)=length(spiketimes);
            end
            for repPNP=1:out.nrepsOFF(ssdindex)
                st=out.MSilentSoundOFF(repPNP).spiketimes;
                spiketimes=st(st>laserstart & st<stop); % spiketimes in region
                OFF(repPNP)=length(spiketimes);
            end
            
            totBL=sum(OFF);
            
            if totBL<1
                meanBLHZ=0
            else
                meanBLHZ=(totBL*(1000/75))/repPNP
            end
            
            
            st=out.mMPulse(pwindex).spiketimes;
            X=0-binwidth:binwidth:75+binwidth;
            [n , x]=hist(st, X);
            n2=n(2:end-1);
            x2=x(2:end-1);
            nmax=find(n2==max(n2));
            xmax=x2(nmax);
            
            
            
            
            
            [h1,p1]=ttest2(ON, OFF, 'tail', 'right');
            % fprintf('\n p=%.4f', p1);
            if isnan(h1) h1=0;end
            if h1
                K1 = [K1 'D-Act'];
                K1a = nmax;
               
            else
                K1 = [K1 'nada'];
                K1a = ''
                
            end
            
           
            
            [h2,p2]=ttest2(ON, OFF, 'tail', 'left');
            % fprintf('\n p=%.4f', p2);
            if isnan(h2) h2=0;end
            if h2
                K2 = [K2 'Suppd'];
            else
                K2 = [K2 'nada'];
            end
            
            
        end
        
       
        
        if length(ON) ~=  out.nrepsPulse(pwindex)
            error('sanity check failure')
        end
        
        cd(    '\\wehrrig1b\D\lab\djmaus\Data\Emily')
        fid3=fopen('meanBL.txt', 'a');
       
%                  fprintf(fid3,'\n Delayed activators or VIP cells, n=%d', totD_Act)
%        fprintf(fid3,'\n Suppressed cells, n=%d', totSuppd)
%                  fclose(fid3);
%         
%         fid3=fopen('test.txt', 'a');      


fprintf(fid3, '\n%s %s %s\t ', newdatadir1, outfilename1, outfilename2);
%     fprintf(fid3,'%s %d %s %s %s %s %s',K1,K1a,K2,K3,K4,K5,K6)
%     fprintf(fid3, '%s ', K7);
%     fprintf(fid3, '%d ', K7a);
%     fprintf(fid3, '%s ', K8);
%     fprintf(fid3, '%d ', K8a);
% fprintf(fid3, '\n%s %s OF\t ', newdatadir1, outfilename1');
% fprintf(fid3, '%d\t ', K9);
% fprintf(fid3, '\n%s %s ON\t ', newdatadir1, outfilename1);
% fprintf(fid3, '%d\t ', K10);
    fprintf(fid3, '%s\t ', clustqual);
fprintf(fid3, '%.2f\t ', meanBLHZ);
    fprintf(fid3, '%.2f ', GiNBLHZ);


        %trying to output duration tuning curves
        
        %        fid4=fopen('GTR_On_v_OFF')
        %
        %
        % %        for f1=figure
        % %            plot(GTRspikesOFF)
        % %            print(psfilename, '-dpsc2', '-append')
        % %        end
        % %        for f2=figure
        % %           plot(GTRspikesON)
        % %           print(psfilename, '-dpsc2', '-append')
        % %        end
        %        i=i+1;
        %
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
