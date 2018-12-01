function Analyze_AD_behavior
%new version 10.31.2018
%only difference is I throw out early datafiles with partial gapdur arrays
%(I switched to full gapdur arrays early on as soon as I found the "save
%time" bug)

if ispc %assume we're on wehrrig2b
    dataroot='D:\lab\djmaus\Data\Kat';
    figs_dir=dataroot;
elseif ismac
    cd /Volumes
    mkdir wehrrig2b
    system('mount_smbfs //lab:mausA1@wehrrig2b/d /Volumes/wehrrig2b');
    dataroot='/Volumes/wehrrig2b/lab/djmaus/Data/Kat';
    figs_dir='/Users/mikewehr/Documents/Manuscripts/AD paper';

end

reprocess=0;
if reprocess
    %process group data
    fprintf('\nReprocessing data...\n')
    
    cd(dataroot);
    
    cell_list='behavior_list.txt';
    fid=fopen(cell_list);
    fseek(fid, 0, 1); %fastforward to end, to get file size
    filesize=ftell(fid);
    fseek(fid, 0, -1); %rewind to start
    wb=waitbar(0, 'processing data');
    
    XFADidx=0;
    controlidx=0;
    XFAD=zeros(2,8);
    ctrlgroup=zeros(2,8);
    gdxfad=zeros(2,8);
    gdcontrol=zeros(2,8);
    
    sessioncount=0;
    excludedcount=0;
    
    
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
            age_wks=strsplit(agestr, ': ');
            age_wks=age_wks{2};
            type_=strsplit(typestr, ': ');
            type_=type_{2};
            
            sessioncount=sessioncount+1;
            fprintf('\n%s', datadir);
            
            dirstr=strsplit(datadir, '\');
            datadir2=dirstr{end};
            
            cd(dataroot)
            
            if  strfind(datadir, 'lab\djmaus\Data\lab\')
                cd ../lab
            end
            cd(datadir2)
            outfilename=sprintf('outGPIAS_Behavior.mat');
            d=dir(outfilename);
            if ~isempty(d)
                load(outfilename)
            else
                warning('outfile missing for this data directory')
              % keyboard
               PlotGPIAS_Behavior_kip
                 load(outfilename);
               
            end
            
            
            load notebook
            sex=nb.mouseSex;
            dob=nb.mouseDOB;
            mouseID=nb.mouseID;
            genotype=nb.mouseGenotype;
            rundate=stimlog(1).timestamp;
            if strcmp(dob, 'unknown') | strcmp(dob, 'age unknown')
                age_days=nan;
                warning('notebook file incomplete')
                keyboard
            else
                age_days=age(dob, rundate);
            end
            
            fprintf('\nsex: %s %s dob: %s age: %d %s', sex, genotype, dob, age_days, age_wks)
            
            if strcmp(genotype, 'unknown')
                keyboard
            end
            
            
            IL=out.IL; %whether there are any interleaved laser trials
            numpulseamps=out.numpulseamps;
            numgapdurs=out.numgapdurs;
            pulseamps=out.pulseamps;
            gapdurs=out.gapdurs;
            gapdelay=out.gapdelay;
            samprate=out.samprate; %in Hz
            mM1ON=out.mM1ON;
            mM1OFF=out.mM1OFF;
            M1ON=out.M1ON;
            M1OFF=out.M1OFF;
            nrepsON=out.nrepsON;
            nrepsOFF=out.nrepsOFF;
            soa=out.soa;
            isi=out.isi;
            soaflag=out.soaflag;
            PeakON=out.PeakON;
            PeakOFF=out.PeakOFF;
            mPeakON=out.mPeakON;
            mPeakOFF=out.mPeakOFF;
            semPeakON=out.semPeakON;
            semPeakOFF=out.semPeakOFF;
            percentGPIAS_ON=out.percentGPIAS_ON;
            percentGPIAS_OFF=out.percentGPIAS_OFF;
            pON = out.pON;
            pOFF = out.pOFF;
            samprate=out.samprate;
            xlimits=out.xlimits;
            
            paindex=1;
            
            
            
            
            switch group
                case 'XFAD'
                    if numgapdurs==8 %skip partial gapdur arrays
                        XFADidx=XFADidx+1;
                        XFAD(XFADidx,:)=mPeakOFF;
                        gdxfad(XFADidx,:)=gapdurs;
                        GPIASxfad(XFADidx,:)=percentGPIAS_OFF;
                        age_xfad(XFADidx)=age_days;
                        age_wks_xfad{XFADidx}=age_wks;
                        sex_xfad{XFADidx}=sex;
                        mouseID_xfad{XFADidx}=mouseID;
                    else
                        fprintf('\n\texcluded')
                        excludedcount=excludedcount+1;
                    end
                case '100012'
                    
                    if numgapdurs==8
                        controlidx=controlidx+1;
                        ctrlgroup(controlidx,:)=mPeakOFF;
                        gdcontrol(controlidx,:)=gapdurs;
                        GPIAScontrol(controlidx,:)=percentGPIAS_OFF;
                        age_control(controlidx)=age_days;
                        age_wkscontrol{controlidx}=age_wks;
                        sex_control{controlidx}=sex;
                        mouseID_control{controlidx}=mouseID;
                    else
                        fprintf('\n\texcluded')
                        excludedcount=excludedcount+1;
                    end
                otherwise
                    group
                    keyboard
            end
            
            
            
            
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
    close(wb); %close the waitbar window
    cd(dataroot)
    generated_by=which(mfilename);
    generated_on=datestr(now);
    save ADgroupdata
else
    fprintf('\nloading groupdata from file...\n')
    cd(dataroot)
    load ADgroupdata
end %if process groupdata

numXFAD_sessions=size(XFAD, 1);
numcontrol_sessions=size(ctrlgroup, 1);
numXFAD_mice=length(unique(mouseID_xfad));
numcontrol_mice=length(unique(mouseID_control));
log2gapdurs=log2(gapdurs(2:8));

fs=18;
ms=40;
lw=2;


cd(figs_dir)
fname=sprintf('AD-stats-output-%s.txt', datestr(now));
if ispc
fname=strrep(fname, ':', '-');
end
fid=fopen(fname, 'w');
fprintf(fid, '\nrun on %s with reprocess=%d', datestr(now), reprocess);

fprintf(fid, '\n%d total sessions', sessioncount);
fprintf(fid, '\n%d excluded sessions', excludedcount);
fprintf(fid, '\n%d XFAD sessions', numXFAD_sessions);
fprintf(fid, '\n%d control sessions', numcontrol_sessions);
fprintf(fid, '\n%d XFAD mice', numXFAD_mice);
fprintf(fid, '\n%d control mice', numcontrol_mice);
fprintf(fid, '\n');


%plot the mean peak rectified startle
figure;
hold on
for i=1:numXFAD_sessions
    mPeakOFF=XFAD(i,:);
    gd=gdxfad(i,:);
    plot(1:length(gd), mPeakOFF, 'r-o')
    set(gca, 'xtick', 1:length(gd))
    set(gca, 'xticklabel', gd)
    xlabel('gap duration')
    ylabel('startle response +- sem')
    
end

for i=1:numcontrol_sessions
    mPeakOFF=ctrlgroup(i,:);
    gd=gdcontrol(i,:);
    plot(1:length(gd), mPeakOFF, 'k-o')
    set(gca, 'xtick', 1:length(gd))
    set(gca, 'xticklabel', gd)
    xlabel('gap duration')
    ylabel('startle response +- sem')
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fig1d - bar graph of pure startles
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
hold on
xfad_purestartle=XFAD(:,1)';
control_purestartle=ctrlgroup(:,1)';
groupnames={'5XFAD'; 'control'};
width=1;
bw_colormap=[1 0 0; 0 0 0];
barweb([mean(xfad_purestartle); mean(control_purestartle) ],...
    [std(xfad_purestartle) ; std(control_purestartle)], ...
    width, groupnames, [],[],[],bw_colormap)
set(gca, 'xtick', [.85 1.15])
set(gca, 'ytick', [0:.25:1.5])
ylabel('pure startle response')
set(gca, 'fontsize', fs)
    f=findobj('type', 'bar');
    set(f(2), 'FaceColor', 'r')
    set(f(1), 'FaceColor', [1 1 1]*.25)
    
cd(figs_dir)
print -dpdf fig1d.pdf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fig1c - all sessions GPIAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%plot GPIAS
figure;
hold on
%dummies for legend
    gpias=GPIAScontrol(1,2:8);
    plot(log2gapdurs, gpias, '.k-', 'markersize', ms, 'linewidth', lw)
    gpias=GPIASxfad(1,2:8);
    gd=gdxfad(1,:);
    plot(log2gapdurs, gpias, '.r-', 'markersize', ms, 'linewidth', lw)

for i=1:numXFAD_sessions
    gpias=GPIASxfad(i,2:8);
    gd=gdxfad(i,:);
    plot(log2gapdurs, gpias, '.r-', 'markersize', ms, 'linewidth', lw)
end
for i=1:numcontrol_sessions
    gpias=GPIAScontrol(i,2:8);
    gd=gdcontrol(i,:);
    plot(log2gapdurs, gpias, '.k-', 'markersize', ms, 'linewidth', lw)
    xlabel('gap duration')
    ylabel('% gpias')
end
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [1 2 3 4 8 32 256])
ylabel('gap detection, %')
set(gca, 'fontsize', fs)
xlim([1 numgapdurs+1])
ylim([-80 100])
xlim([-.5 8.5])
xlabel('gap duration, ms')
set(gca, 'ytick', [-50:50:100])
    legend( 'control', '5XFAD', 'location', 'southeast')
cd(figs_dir)
print -dpdf fig1c.pdf

% mean startle responses, all gapdurs all sessions
figure
hold on
plot(1:8, mean(ctrlgroup, 1), 'o-k')
plot(1:8, mean(XFAD, 1), 'o-r')
set(gca, 'xtick', 1:8)
set(gca, 'xticklabel', [0 1 2 3 4 8 32 256])
title('group means')
ylabel('startle response')
legend('control', 'XFAD', 'location', 'southeast')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fig1b - group means GPIAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hold on
plot(log2gapdurs, nanmean(GPIAScontrol(:,2:8), 1), '.-k', 'markersize', ms)
plot(log2gapdurs, nanmean(GPIASxfad(:,2:8), 1), '.-r', 'markersize', ms)
e1=errorbar(log2gapdurs, nanmean(GPIAScontrol(:,2:8), 1), nanstd(GPIAScontrol(:,2:8),[], 1)./sqrt(numcontrol_sessions));
e2=errorbar(log2gapdurs, nanmean(GPIASxfad(:,2:8), 1), nanstd(GPIASxfad(:,2:8),[], 1)./sqrt(numXFAD_sessions));
set(e1, 'color', 'k', 'markersize', ms, 'linewidth', lw)
set(e2, 'color', 'r', 'markersize', ms, 'linewidth', lw)
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [1 2 3 4 8 32 256])
ylabel('gap detection, %')
legend('control', '5XFAD', 'location', 'southeast')
set(gca, 'fontsize', fs)
xlim([-.5 8.5])
ylim([-10 100])
xlabel('gap duration, ms')
set(gca, 'ytick', [0:25:100])
cd(figs_dir)
print -dpdf fig1b.pdf



% stats on group data
%convert to long form
j=0;
for i=1:numcontrol_sessions
    for gd=1:numgapdurs
        j=j+1;
        X(j)=GPIAScontrol(i, gd);
        gdgroup(j)=gd;
        grouping{j}='control';
    end
end
for i=1:numXFAD_sessions
    for gd=1:numgapdurs
        j=j+1;
        X(j)=GPIASxfad(i, gd);
        gdgroup(j)=gd;
        grouping{j}='xfad';
    end
end
[p, anovatab, stats]=kruskalwallis(X, grouping);
df=anovatab{2, 3};
chisq=anovatab{2, 5};
fprintf(fid, '\nkruskal wallis main effect of xfad vs control:');
fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);

%post-hoc test for each gapdur by ranksum
fprintf(fid,'\npost-hoc test of xfad vs control at each gap duration:');
for gd=2:numgapdurs
    p=ranksum(GPIASxfad(:, gd),GPIAScontrol(:, gd));
    fprintf(fid,'\ngd %d ms: p=%g, ', gapdurs(gd), p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now let's break it out by sex
female_xfad_sessions=find(strcmp(sex_xfad, 'f'));
male_xfad_sessions=find(strcmp(sex_xfad, 'm'));
female_control_sessions=find(strcmp(sex_control, 'f'));
male_control_sessions=find(strcmp(sex_control, 'm'));

fprintf(fid,'\n %d female xfad sessions', length(female_xfad_sessions));
fprintf(fid,'\n %d male xfad sessions', length(male_xfad_sessions));
fprintf(fid,'\n %d female control sessions', length(female_control_sessions));
fprintf(fid,'\n %d male control sessions', length(male_control_sessions));

figure
hold on
legstr={};
i=0;
plot(log2gapdurs, nanmean(GPIAScontrol(male_control_sessions,2:8), 1), 'o-c')
i=i+1;legstr{i}='male control';
plot(log2gapdurs, nanmean(GPIAScontrol(female_control_sessions,2:8), 1), 'o-m')
i=i+1;legstr{i}='female control';
plot(log2gapdurs, nanmean(GPIAScontrol(:,2:8), 1), 'o-k')
i=i+1;legstr{i}='all control';
plot(log2gapdurs, nanmean(GPIASxfad(:,2:8), 1), 'o-r')
i=i+1;legstr{i}='all xfad';
plot(log2gapdurs, nanmean(GPIASxfad(male_xfad_sessions,2:8), 1), 'o-b')
i=i+1;legstr{i}='male xfad';
plot(log2gapdurs, nanmean(GPIASxfad(female_xfad_sessions,2:8), 1), 'o-m')
i=i+1;legstr{i}='female xfad';
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [ 1 2 3 4 8 32 256])
xlim([-.5 8.5])
title('group means')
ylabel('% gpias')
legend(legstr, 'location', 'southeast')

%same thing but prettier
figure
hold on
legstr={};
i=0;p=[];
p(i+1)=plot(log2gapdurs, nanmean(GPIAScontrol(male_control_sessions,2:8), 1), 'o-c');
i=i+1;legstr{i}='male control';
p(i+1)=plot(log2gapdurs, nanmean(GPIAScontrol(female_control_sessions,2:8), 1), 'o-m');
i=i+1;legstr{i}='female control';
p(i+1)=plot(log2gapdurs, nanmean(GPIASxfad(male_xfad_sessions,2:8), 1), '.-c');
i=i+1;legstr{i}='male 5XFAD';
p(i+1)=plot(log2gapdurs, nanmean(GPIASxfad(female_xfad_sessions,2:8), 1), '.-m');
i=i+1;legstr{i}='female 5XFAD';
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [ 1 2 3 4 8 32 256])
xlim([-.5 8.5])
legend(legstr, 'location', 'southeast')
set(gca, 'fontsize', fs)
ylabel('gap detection, %')
xlabel('gap duration, ms')
set(p, 'markersize', ms, 'linewidth', lw)
set(p(1:2), 'markersize', 12, 'MarkerFaceColor', 'w')
cd(figs_dir)
print -dpdf fig2.pdf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now let's break it out by age
%plot GPIAS

cmap=cool(1+range([age_control age_xfad]));
minage=min([age_control age_xfad]);
age_ticks=prctile([age_control age_xfad], 0:25:100);
fprintf(fid, '\n\nage breakdown (min, 25%%, median, 75%%, max):');
fprintf(fid, '\n%d', round(age_ticks));

figure;
hold on
legstr={};
for i=1:numXFAD_sessions
    gpias=GPIASxfad(i,2:8);
    gd=gdxfad(i,2:8);
    h(i)=plot(log2gapdurs, gpias, '.-');
    set(h(i), 'color', cmap(1+age_xfad(i)-minage,:))
    legstr{i}=sprintf('%d days',age_xfad(i)) ;
end
set(h, 'markersize', ms, 'linewidth', lw)

xlim([-.5 8.5])
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [1 2 3 4 8 32 256])
set(gca, 'fontsize', fs)
ylabel('gap detection, %')
xlabel('gap duration, ms')
title('5XFAD by age')
% legend(legstr, 'location', 'southeast')
cb=colorbar;
colormap(cb, cool)
set(cb, 'ticks', 0:.25:1);
set(cb, 'ticklabels', round(age_ticks));
cbl=get(cb, 'label');
set(cbl, 'string', 'age in days');
cd(figs_dir)
print -dpdf fig3a.pdf


figure; hold on
legstr={};
for i=1:numcontrol_sessions
    gpias=GPIAScontrol(i,2:8);
    gd=gdcontrol(i,:);
    h2(i)=plot(log2gapdurs, gpias, '-o');
    if ~isnan(age_control(i))
        set(h2(i), 'color', cmap(1+age_control(i)-minage,:))
    else
        set(h2(i), 'color', 'k')
    end
    legstr{i}=sprintf('%d days',age_control(i)) ;
    xlabel('gap duration')
    ylabel('% gpias')
end
set(h2, 'markersize', 12, 'linewidth', lw)

xlim([-.5 8.5])
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [1 2 3 4 8 32 256])
set(gca, 'fontsize', fs)
ylabel('gap detection, %')
xlabel('gap duration, ms')
title('control by age')
% legend(legstr, 'location', 'southeast')
cb=colorbar;
colormap(cb, cool)
set(cb, 'ticks', 0:.25:1);
set(cb, 'ticklabels', round(age_ticks));
cbl=get(cb, 'label');
set(cbl, 'string', 'age in days');

cd(figs_dir)
print -dpdf fig3b.pdf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%look at GD over time
fprintf('\n')
for gdindex=2:8
    figure; hold on
    plot(age_control, GPIAScontrol(:,gdindex), 'ko')
    plot(age_xfad, GPIASxfad(:,gdindex), 'ro')
    xlabel('age, days')
    ylabel('% gpias')
    lsline
    title(sprintf('GPIAS for  %d ms gap',  gapdurs(gdindex)))
    
    Ycontrol=GPIAScontrol(:,gdindex);
    Xcontrol=[age_control(:) ones(size(age_control(:)))];
    [B,BINT,R,RINT,STATScontrol] = regress(Ycontrol,Xcontrol);
    r2control=STATScontrol(1);
    pcontrol=STATScontrol(3);
    fprintf(fid, '\ncontrol GPIAS by age, %d ms gap, r2=%.2f, p=%.4f', gapdurs(gdindex), r2control, pcontrol);
    
    Yxfad=GPIASxfad(:,gdindex);
    Xxfad=[age_xfad(:) ones(size(age_xfad(:)))];
    [B,BINT,R,RINT,STATSxfad] = regress(Yxfad,Xxfad);
    r2xfad=STATSxfad(1);
    pxfad=STATSxfad(3);
    fprintf(fid, '\nxfad GPIAS by age, %d ms gap, r2=%.2f, p=%.4f', gapdurs(gdindex), r2xfad, pxfad);
    legend(sprintf('control p=%.4f', pcontrol), sprintf('xfad p=%.4f', pxfad))
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'\n');
fprintf(fid,'\nat what age is there first a detectable effect on GD across all gapdurs');

% stats on group data
%convert to long form
j=0;
for i=1:numcontrol_sessions
    for gd=1:numgapdurs
        j=j+1;
        X(j)=GPIAScontrol(i, gd);
        gdgroup(j)=gd;
        grouping{j}='control';
        age_long(j)=age_control(i);
    end
end
for i=1:numXFAD_sessions
    for gd=1:numgapdurs
        j=j+1;
        X(j)=GPIASxfad(i, gd);
        gdgroup(j)=gd;
        grouping{j}='xfad';
        age_long(j)=age_xfad(i);
    end
end
%bin age and do stats
age_min=min(age_long);
age_max=70;
i=find(age_long<age_max & age_long>age_min);
[p, anovatab, stats]=kruskalwallis(X(i), grouping(i));
df=anovatab{2, 3};
chisq=anovatab{2, 5};
fprintf(fid,'\nkruskal wallis main effect of xfad vs control, age range %d - %d:', age_min, age_max);
fprintf(fid,'\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
fprintf(fid,'\n%d xfad sessions', length(strmatch('xfad', grouping(i)))/numgapdurs);
fprintf(fid,'\n%d control sessions', length(strmatch('control', grouping(i)))/numgapdurs);
fprintf(fid,'\n');

fclose(fid)
type(fname)

if 0
    cd(dataroot)
    delete ADfigs.ps
    f=findobj('type', 'figure');
    for i=1:length(f)
        figure(i)
        print -dpsc2 ADfigs.ps -append
    end
end

% keyboard



