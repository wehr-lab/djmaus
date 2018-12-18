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
close all hidden
reprocess=1;


if reprocess
    %process group data
    fprintf('\nReprocessing data...\n')
    
    cd(dataroot);
    
        cell_list='behavior_list.txt';
%     cell_list='behavior_list_PUPS.txt';
    fid=fopen(cell_list);
    fid_mouselist=fopen('mouse_list.txt', 'w'); %generate an output text file to compare behavior list with notebook (mouse ID, geneotype, sex, etc)
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
    
    %fast-forward
    %fseek(fid, 10000, -1)
    
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
            fprintf('\ngroup: %s', group);
            fprintf(' %s', filename);
            fprintf(fid_mouselist, '\n%s', datadir);
            fprintf(fid_mouselist, '\ngroup: %s', group);
            
            dirstr=strsplit(datadir, '\');
            datadir2=dirstr{end};
            
            cd(dataroot)
            
            if  strfind(datadir, 'lab\djmaus\Data\lab\')
                cd ../lab
            end
            cd(datadir2)
            %             outfilename=sprintf('outGPIAS_Behavior.mat');
            outfilename=filename;
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
            
            fprintf('\nsex: %s %s dob: %s age: %d %s', sex, genotype, dob, age_days, age_wks);
            fprintf(fid_mouselist, '\nsex: %s %s dob: %s age: %d %s', sex, genotype, dob, age_days, age_wks);
            
            if strcmp(genotype, 'unknown')
                %                 keyboard
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
                case '5XFAD'
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
                        fprintf('\n\texcluded');
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
                        fprintf('\n\texcluded');
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
    fclose(fid_mouselist);
    close(wb); %close the waitbar window
    cd(dataroot)
    clear figs_dir dataroot datadir datadir2
    s=whos;
    for i=1:size(s)
        if s(i).bytes>20000 clear(s(i).name);end
    end
    
    generated_by=which(mfilename);
    generated_on=datestr(now);
    save ADgroupdata
else
    fprintf('\nloading groupdata from file...\n')
    cd(dataroot)
    load ADgroupdata
end %if process groupdata

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numXFAD_sessions=size(XFAD, 1);
numcontrol_sessions=size(ctrlgroup, 1);
numXFAD_mice=length(unique(mouseID_xfad));
numcontrol_mice=length(unique(mouseID_control));
log2gapdurs=log2(gapdurs(2:8));

fs=18;
ms=40;
lw=2;

if ispc %assume we're on wehrrig2b
    dataroot='D:\lab\djmaus\Data\Kat';
    figs_dir=dataroot;
elseif ismac
    dataroot='/Volumes/wehrrig2b/lab/djmaus/Data/Kat';
    figs_dir='/Users/mikewehr/Documents/Manuscripts/AD paper';
end

cd(figs_dir)
fname=sprintf('AD-stats-output-%s.txt', datestr(now));
if ispc
    fname=strrep(fname, ':', '-');
end
fid=fopen(fname, 'w');
fprintf(fid, '\nrun on %s with reprocess=%d', datestr(now), reprocess);

fprintf(fid, '\n%d total sessions', sessioncount);
fprintf(fid, '\n%d included sessions', numXFAD_sessions+numcontrol_sessions);
fprintf(fid, '\n%d excluded sessions', excludedcount);
fprintf(fid, '\n%d XFAD sessions', numXFAD_sessions);
fprintf(fid, '\n%d control sessions', numcontrol_sessions);
fprintf(fid, '\n%d XFAD mice', numXFAD_mice);
fprintf(fid, '\n%d control mice', numcontrol_mice);
fprintf(fid, '\n%d total mice', numXFAD_mice + numcontrol_mice);
fprintf(fid, '\n%.1f average sessions per mouse', (numXFAD_sessions+numcontrol_sessions)/(numXFAD_mice + numcontrol_mice));

%how many sessions per mouse?
[umid, ~, ic]=unique(mouseID_xfad);
[xfad_sessionspermouse, bin] = histc(ic, unique(ic));
fprintf(fid, '\n%.1f average sessions per xfad mouse', mean(xfad_sessionspermouse));
fprintf(fid, '(range: %d - %d)', min(xfad_sessionspermouse),max(xfad_sessionspermouse) );

[umidc, ~, icc]=unique(mouseID_control);
[control_sessionspermouse, bin] = histc(icc, unique(icc));
fprintf(fid, '\n%.1f average sessions per control mouse', mean(control_sessionspermouse));
fprintf(fid, '(range: %d - %d)', min(control_sessionspermouse),max(control_sessionspermouse) );


fprintf(fid, '\n');


%restrict to ages > p60
XFADp60_sessions=find(age_xfad>60);
controlp60_sessions=find(age_control>60);
fprintf(fid, '\n%d over-p60 XFAD sessions', length(XFADp60_sessions));
fprintf(fid, '\n%d over-p60 control sessions', length(controlp60_sessions));


%plot the mean peak rectified startle
figure;
hold on
for i=XFADp60_sessions
    mPeakOFF=XFAD(i,:);
    gd=gdxfad(i,:);
    plot(1:length(gd), mPeakOFF, 'r-o')
    set(gca, 'xtick', 1:length(gd))
    set(gca, 'xticklabel', gd)
    xlabel('gap duration')
    ylabel('startle response +- sem')
    
end

for i=1:controlp60_sessions
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
xfad_purestartle=XFAD(XFADp60_sessions,1)';
control_purestartle=ctrlgroup(controlp60_sessions,1)';
groupnames={'5XFAD'; 'control'};
width=1;
cmap=[1 0 0; 0 0 0];
hxfad=bar(1, mean(xfad_purestartle));
hctrl=bar(2, mean(control_purestartle));
eb=errorbar([mean(xfad_purestartle); mean(control_purestartle) ], ...
    [std(xfad_purestartle) ; std(control_purestartle)]);
set(eb, 'marker', 'none', 'linewidth', 2, 'linestyle', 'none', 'color', 'k')
set(gca, 'xtick', [1 2], 'xticklabel', {'5XFAD', 'control'})
set(gca, 'ytick', [0:.25:1.5])
ylabel('pure startle response')
set(gca, 'fontsize', fs)
set(hxfad, 'FaceColor', 'r')
set(hctrl, 'FaceColor', [1 1 1]*.25)

try
    p=ranksum(xfad_purestartle, control_purestartle);
    fprintf(fid, '\n\npure startles, ranks sum p=%.4f', p);
catch
        fprintf(fid, '\nnot enough >p60 sessions for stats');
end
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

for i=XFADp60_sessions
    gpias=GPIASxfad(i,2:8);
    gd=gdxfad(i,:);
    plot(log2gapdurs, gpias, '.r-', 'markersize', ms, 'linewidth', lw)
end
for i=controlp60_sessions
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
plot(log2gapdurs, nanmean(GPIAScontrol(controlp60_sessions,2:8), 1), '.-k', 'markersize', ms)
plot(log2gapdurs, nanmean(GPIASxfad(XFADp60_sessions,2:8), 1), '.-r', 'markersize', ms)
e1=errorbar(log2gapdurs, nanmean(GPIAScontrol(controlp60_sessions,2:8), 1), nanstd(GPIAScontrol(controlp60_sessions,2:8),[], 1)./sqrt(length(controlp60_sessions)));
e2=errorbar(log2gapdurs, nanmean(GPIASxfad(XFADp60_sessions,2:8), 1), nanstd(GPIASxfad(XFADp60_sessions,2:8),[], 1)./sqrt(length(XFADp60_sessions)));
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
try
    j=0;
    for i=controlp60_sessions
        for gd=1:numgapdurs
            j=j+1;
            X(j)=GPIAScontrol(i, gd);
            gdgroup(j)=gd;
            grouping{j}='control';
        end
    end
    for i=XFADp60_sessions
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
    fprintf(fid, '\nkruskal wallis main effect of xfad vs control over-p60:');
    fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    
    % post-hoc test for each gapdur by ranksum
    fprintf(fid,'\npost-hoc test of xfad vs control at each gap duration, over-p60:');
    for gd=2:numgapdurs
        p=ranksum(GPIASxfad(XFADp60_sessions, gd),GPIAScontrol(controlp60_sessions, gd));
        fprintf(fid,'\ngd %d ms: p=%g, ', gapdurs(gd), p);
    end
    catch
        fprintf(fid, '\nnot enough >p60 sessions for stats');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now let's break it out by sex
female_xfad_sessions=find(strcmp(sex_xfad, 'f'));
male_xfad_sessions=find(strcmp(sex_xfad, 'm'));
female_control_sessions=find(strcmp(sex_control, 'f'));
male_control_sessions=find(strcmp(sex_control, 'm'));

femalep60_xfad_sessions=intersect(female_xfad_sessions, XFADp60_sessions);
malep60_xfad_sessions=intersect(male_xfad_sessions, XFADp60_sessions);
femalep60_control_sessions=intersect(female_control_sessions, controlp60_sessions);
malep60_control_sessions=intersect(male_control_sessions, controlp60_sessions);


fprintf(fid,'\n %d female over-p60 xfad sessions', length(femalep60_xfad_sessions));
fprintf(fid,'\n %d male over-p60 xfad sessions', length(malep60_xfad_sessions));
fprintf(fid,'\n %d female over-p60 control sessions', length(femalep60_control_sessions));
fprintf(fid,'\n %d male over-p60 control sessions', length(malep60_control_sessions));

figure
hold on
legstr={};
i=0;
plot(log2gapdurs, nanmean(GPIAScontrol(malep60_control_sessions,2:8), 1), 'o-c')
i=i+1;legstr{i}='male control';
plot(log2gapdurs, nanmean(GPIAScontrol(femalep60_control_sessions,2:8), 1), 'o-m')
i=i+1;legstr{i}='female control';
plot(log2gapdurs, nanmean(GPIAScontrol(:,2:8), 1), 'o-k')
i=i+1;legstr{i}='all control';
plot(log2gapdurs, nanmean(GPIASxfad(:,2:8), 1), 'o-r')
i=i+1;legstr{i}='all xfad';
plot(log2gapdurs, nanmean(GPIASxfad(malep60_xfad_sessions,2:8), 1), 'o-b')
i=i+1;legstr{i}='male xfad';
plot(log2gapdurs, nanmean(GPIASxfad(femalep60_xfad_sessions,2:8), 1), 'o-m')
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
i=0;p=[];e=[];
mgpias=nanmean(GPIAScontrol(malep60_control_sessions,2:8), 1);
sem=nanstd(GPIAScontrol(malep60_control_sessions,2:8), 1)./sqrt(length(malep60_control_sessions));
p(i+1)=plot(log2gapdurs, mgpias, 'o-c');
e(i+1)=errorbar(log2gapdurs, mgpias, sem, 'c');
i=i+1;legstr{i}='male control';
mgpias=nanmean(GPIAScontrol(femalep60_control_sessions,2:8), 1);
sem=nanstd(GPIAScontrol(femalep60_control_sessions,2:8), 1)./sqrt(length(femalep60_control_sessions));
p(i+1)=plot(log2gapdurs, mgpias, 'o-m');
e(i+1)=errorbar(log2gapdurs, mgpias, sem, 'm');
i=i+1;legstr{i}='female control';
mgpias=nanmean(GPIASxfad(malep60_xfad_sessions,2:8), 1);
sem=nanstd(GPIASxfad(malep60_xfad_sessions,2:8), 1)./sqrt(length(malep60_xfad_sessions));
p(i+1)=plot(log2gapdurs, mgpias, '.-c');
e(i+1)=errorbar(log2gapdurs, mgpias, sem, 'c');
i=i+1;legstr{i}='male 5XFAD';
mgpias=nanmean(GPIASxfad(femalep60_xfad_sessions,2:8), 1);
sem=nanstd(GPIASxfad(femalep60_xfad_sessions,2:8), 1)./sqrt(length(femalep60_xfad_sessions));
p(i+1)=plot(log2gapdurs, mgpias, '.-m');
e(i+1)=errorbar(log2gapdurs, mgpias, sem, 'm');
i=i+1;legstr{i}='female 5XFAD';
set(gca, 'xtick', log2gapdurs)
set(gca, 'xticklabel', [ 1 2 3 4 8 32 256])
xlim([-.5 8.5])
legend(p, legstr, 'location', 'southeast')
set(gca, 'fontsize', fs)
ylabel('gap detection, %')
xlabel('gap duration, ms')
set([p], 'markersize', ms, 'linewidth', lw)
set([e],  'linewidth', lw)
set(p(1:2), 'markersize', 12, 'MarkerFaceColor', 'w')
ylim([-25 100])
cd(figs_dir)
print -dpdf fig2.pdf

%stats on sex
%convert to long form
try
    clear sex
    j=0;
    for i=malep60_control_sessions
        for gd=1:numgapdurs
            j=j+1;
            X(j)=GPIAScontrol(i, gd);
            gdgroup(j)=gd;
            grouping{j}='control';
            sex{j}='male';
            mouseid_control_male(j)=mouseID_control(i);
        end
    end
    k=0;
    for i=malep60_xfad_sessions
        for gd=1:numgapdurs
            k=k+1;
            j=j+1;
            X(j)=GPIASxfad(i, gd);
            gdgroup(j)=gd;
            grouping{j}='xfad';
            sex{j}='male';
            mouseid_xfad_male(k)=mouseID_xfad(i);
        end
    end
    [p, anovatab, stats]=kruskalwallis(X, grouping);
    df=anovatab{2, 3};
    chisq=anovatab{2, 5};
    fprintf(fid, '\n\nEffect of sex (over-p60):');
    fprintf(fid, '\nkruskal wallis effect of xfad vs control for males:');
    fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    
    
    j=0;
    clear X grouping
    for i=femalep60_control_sessions
        for gd=1:numgapdurs
            j=j+1;
            X(j)=GPIAScontrol(i, gd);
            gdgroup(j)=gd;
            grouping{j}='control';
            sex{j}='female';
            mouseid_control_female(j)=mouseID_control(i);
        end
    end
    k=0;
    for i=femalep60_xfad_sessions
        for gd=1:numgapdurs
            k=k+1;
            j=j+1;
            X(j)=GPIASxfad(i, gd);
            gdgroup(j)=gd;
            grouping{j}='xfad';
            sex{j}='female';
            mouseid_xfad_female(k)=mouseID_xfad(i);
        end
    end
    [p, anovatab, stats]=kruskalwallis(X, grouping);
    df=anovatab{2, 3};
    chisq=anovatab{2, 5};
    fprintf(fid, '\nkruskal wallis effect of xfad vs control for females:');
    fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    fprintf(fid, '\n\nall age mice:')
    fprintf(fid, '\nn=%d control males', length(unique(mouseid_control_male)));
    fprintf(fid, '\nn=%d xfad males', length(unique(mouseid_xfad_male)));
    fprintf(fid, '\nn=%d control females', length(unique(mouseid_control_female)));
    fprintf(fid, '\nn=%d xfad females', length(unique(mouseid_xfad_female)));
    
    j=0;
    clear X grouping sex
    for i=femalep60_control_sessions
        for gd=1:numgapdurs
            j=j+1;
            X(j)=GPIAScontrol(i, gd);
            gdgroup(j)=gd;
            grouping{j}='control';
            sex{j}='female';
        end
    end
    for i=malep60_control_sessions
        for gd=1:numgapdurs
            j=j+1;
            X(j)=GPIAScontrol(i, gd);
            gdgroup(j)=gd;
            grouping{j}='control';
            sex{j}='male';
        end
    end
    [p, anovatab, stats]=kruskalwallis(X, sex);
    df=anovatab{2, 3};
    chisq=anovatab{2, 5};
    fprintf(fid, '\n\nkruskal wallis effect of sex for over-p60 control mice:');
    fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    
    fprintf(fid, '\nn=%d male xfad sessions', length(male_xfad_sessions));
    fprintf(fid, '\nn=%d female control sessions', length(female_control_sessions));
    fprintf(fid, '\nn=%d female xfad sessions', length(female_xfad_sessions));
    fprintf(fid, '\nn=%d male control sessions',length( male_control_sessions));
    
catch
        fprintf(fid, '\nnot enough >p60 sessions for stats');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now let's break it out by age
%plot GPIAS

cmap=cool(1+range([age_control age_xfad]));
minage=min([age_control age_xfad]);
age_ticks=prctile([age_control age_xfad], 0:25:100);
fprintf(fid, '\n\nage breakdown (min, 25%%, median, 75%%, max):');
fprintf(fid, '\n%d', round(age_ticks));

% %dot transparency - doesn't look that great
% figure;
% hold on
% legstr={};
% for i=1:numXFAD_sessions
%     gpias=GPIASxfad(i,2:8);
%     gd=gdxfad(i,2:8);
%     h(i)=scatter(log2gapdurs, gpias, 100, cmap(1+age_xfad(i)-minage,:), 'filled');
% h2(i)=plot(log2gapdurs, gpias, 'color', cmap(1+age_xfad(i)-minage,:));
%
%     %     set(h(i), 'markeredgecolor','none',...
% %         'markerfacecolor', cmap(1+age_xfad(i)-minage,:))
% %         %'color', cmap(1+age_xfad(i)-minage,:))
%     legstr{i}=sprintf('%d days',age_xfad(i)) ;
% end
% % set(h, 'markersize', 40, 'linewidth', lw)
% alpha .5

figure; hold on
[~,x]=sort(age_xfad);
for i=x
    gpias=GPIASxfad(i,2:8);
    gd=gdxfad(i,:);
    h(i)=plot(log2gapdurs, gpias, '.-');
    set(h(i), 'color', cmap(1+age_xfad(i)-minage,:))
    %      legstr{i}=sprintf('%d days',age_xfad(i)) ;
    xlabel('gap duration')
    ylabel('% gpias')
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
[~,x]=sort(age_control);
for i=x
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
fprintf(fid, '\n');
fprintf(fid, '\nlinear regression of gap detection as a function of age:\n');
figure
subplot1(4,2);
p=0;
for gdindex=2:8
    p=p+1;
    subplot1(p)
    %     set(gca, 'fontsize', fs)
    plot(age_control, GPIAScontrol(:,gdindex), 'k.', 'markersize', ms)
    plot(age_xfad, GPIASxfad(:,gdindex), 'r.', 'markersize', ms)
    xlabel('age, days')
    ylabel('% gpias')
    h=lsline;    set(h, 'linewidth', lw)
    title(sprintf('GPIAS for  %d ms gap',  gapdurs(gdindex)))
    
    Ycontrol=GPIAScontrol(:,gdindex);
    Xcontrol=[age_control(:) ones(size(age_control(:)))];
    [B,BINT,R,RINT,STATScontrol] = regress(Ycontrol,Xcontrol);
    r2control=STATScontrol(1);
    pcontrol=STATScontrol(3);
    Pcontrol(gdindex)=pcontrol;
    fprintf(fid, '\ncontrol GPIAS by age, %d ms gap, r2=%.2f,b=%.2f, p=%.4f', gapdurs(gdindex), r2control, B(1), pcontrol);
    if pcontrol<.05
        text(105, [105 1]*B, '*', 'color', 'k', 'fontsize', 32)
    end
    
    Yxfad=GPIASxfad(:,gdindex);
    Xxfad=[age_xfad(:) ones(size(age_xfad(:)))];
    [B,BINT,R,RINT,STATSxfad] = regress(Yxfad,Xxfad);
    r2xfad=STATSxfad(1);
    pxfad=STATSxfad(3);
    Pxfad(gdindex)=pxfad;
    fprintf(fid, '\nxfad GPIAS by age, %d ms gap, r2=%.2f, b=%.2f, p=%.4f', gapdurs(gdindex), r2xfad, B(1), pxfad);
    if pxfad<.05
        text(105, [105 1]*B, '*', 'color', 'r', 'fontsize', 32)
    end
    
    % legend(sprintf('control p=%.4f', pcontrol), sprintf('xfad p=%.4f', pxfad), 'location', 'southeast')
    
end
set(gcf, 'pos', [25        1136         975        1322])

% choose 4ms and 256ms as examples for fig3c,d
for gdindex=[2 3 4 5 6 7 8]
    figure; hold on
    set(gca, 'fontsize', fs)
    plot(age_control, GPIAScontrol(:,gdindex), 'k.', 'markersize', ms)
    plot(age_xfad, GPIASxfad(:,gdindex), 'r.', 'markersize', ms)
    xlabel('age, days')
    ylabel('gap detection, %')
    h=lsline;
    set(h, 'linewidth', lw)
    title(sprintf('%d ms gap',  gapdurs(gdindex)))
    legend(sprintf('control p=%.4f', Pcontrol(gdindex)), sprintf('xfad p=%.4f', Pxfad(gdindex)), 'location', 'southeast')
    cd(figs_dir)
    if gdindex==5
        print -dpdf fig3c.pdf
    elseif gdindex==8
        legend('location', 'southwest')
        print -dpdf fig3d.pdf
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Table 1
fidtable=fopen('Table 1.txt', 'w');
fprintf(fidtable, '\t\t slope \t r^2 \t p \t n');
fprintf(fidtable, '\n5XFAD');
for gdindex=2:8
    Yxfad=GPIASxfad(:,gdindex);
    Xxfad=[age_xfad(:) ones(size(age_xfad(:)))];
    [B,BINT,R,RINT,STATSxfad] = regress(Yxfad,Xxfad);
    slope_xfad=B(1);
    Slopes_xfad(gdindex-1)=slope_xfad;
    r2xfad=STATSxfad(1);
    pxfad=STATSxfad(3);
    nsess=length(age_xfad);
    nmice=length(unique(mouseID_xfad));
    fprintf(fidtable, '\n\t%d ms \t %.2f \t %.2f \t p=%.4f \t %d sessions, %d mice',...
        gapdurs(gdindex), slope_xfad, r2xfad, pxfad, nsess, nmice);
end
fprintf(fidtable, '\ncontrol');
for gdindex=2:8
    Ycontrol=GPIAScontrol(:,gdindex);
    Xcontrol=[age_control(:) ones(size(age_control(:)))];
    [B,BINT,R,RINT,STATScontrol] = regress(Ycontrol,Xcontrol);
    slope_control=B(1);
    Slopes_control(gdindex-1)=slope_control;
    r2control=STATScontrol(1);
    pcontrol=STATScontrol(3);
    nsess=length(age_control);
    nmice=length(unique(mouseID_control));
    fprintf(fidtable, '\n\t%d ms \t %.2f \t %.2f \t p=%.4f \t %d sessions, %d mice',...
        gapdurs(gdindex), slope_control, r2control, pcontrol, nsess, nmice);
    
end
fclose(fidtable);
p=ranksum(Slopes_control, Slopes_xfad, 'tail', 'right');
fprintf(fid,'\nlinear regression slopes of gap detection vs age were significantly less for xfad mice compared to controls');
fprintf(fid,'\np=%4f, 1-tailed ranksum', p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf(fid,'\n');
fprintf(fid,'\nat what age is there first a detectable effect on GD across all gapdurs?\n');

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
        sex_long{j}=sex_control{i};
        mouseID_long{j}=mouseID_control{i};
    end
end
for i=1:numXFAD_sessions
    for gd=1:numgapdurs
        j=j+1;
        X(j)=GPIASxfad(i, gd);
        gdgroup(j)=gd;
        grouping{j}='xfad';
        age_long(j)=age_xfad(i);
        sex_long{j}=sex_xfad{i};
        mouseID_long{j}=mouseID_xfad{i};
    end
end
%bin age and do stats
age_min=min(age_long);
for age_max=45:5:80;
    i=find(age_long<age_max & age_long>age_min);
    [p, anovatab, stats]=kruskalwallis(X(i), grouping(i), 'off');
    df=anovatab{2, 3};
    chisq=anovatab{2, 5};
    fprintf(fid,'\nkruskal wallis main effect of xfad vs control, age range %d - %d:', age_min, age_max);
    fprintf(fid,'\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    fprintf(fid,'\n%d xfad sessions', length(strmatch('xfad', grouping(i)))/numgapdurs);
    k=strmatch('xfad', grouping(i));
    fprintf(fid,'\n%d xfad mice',    length(unique(mouseID_long(k))));
    
    fprintf(fid,'\n%d control sessions', length(strmatch('control', grouping(i)))/numgapdurs);
    j=strmatch('control', grouping(i));
    fprintf(fid,'\n%d control mice',    length(unique(mouseID_long(j))));
    fprintf(fid,'\n');
end



% repeat but only doing a 1-tailed rank sum on 256 ms
fprintf(fid,'\nrepeat but only doing a 1-tailed rank sum on 256 ms')
for age_max=40:10:100;
    i=find(age_long<age_max & age_long>age_min & gdgroup==8);
    % for gd=2:numgapdurs
    %     p=ranksum(GPIASxfad(:, gd),GPIAScontrol(:, gd));
    %     fprintf(fid,'\ngd %d ms: p=%g, ', gapdurs(gd), p);
    % end
    
    Xi=X(i);
    Gi=grouping(i);
    Xxfad=[];
    Xcontrol=[];
    for j=1:length(Gi)
        switch Gi{j}
            case 'xfad'
                Xxfad=[Xxfad Xi(j)];
            case 'control'
                Xcontrol=[Xcontrol Xi(j)];
        end
    end
    [p, h]=ranksum(Xcontrol,Xxfad, 'tail', 'right');
    fprintf(fid,'\nranksum of gap 256 xfad vs control, age range %d - %d:', age_min, age_max);
    fprintf(fid,'\nh=%d, p=%g ',h, p);
    fprintf(fid,'\n%d xfad sessions', length(Xxfad));
    m=strmatch('xfad', grouping(i));
    fprintf(fid,'\n%d xfad mice',    length(unique(mouseID_long(m))));
    
    fprintf(fid,'\n%d control sessions', length(Xcontrol));
    n=strmatch('control', grouping(i));
    fprintf(fid,'\n%d control mice',    length(unique(mouseID_long(n))));
    fprintf(fid,'\n');
end


% do a nested GLM to jpintly test the effects of age, sex, genotype, and mouseID
% keyboard

T=table(grouping(:), gdgroup(:), age_long(:), sex_long(:), mouseID_long(:), X(:));
T.Properties.VariableNames={'genotype', 'gapdur', 'age', 'sex', 'mouseID', 'GPIAS'};
mdl=fitglm(T, 'linear')
mdl=fitglm(T, 'interactions')



fclose(fid);
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



