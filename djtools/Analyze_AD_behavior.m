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
   
% if you get this error:    mount_smbfs: mount error: /Volumes/wehrrig2b: File exists
%  it's because the volume is already mounted 

end
close all hidden
reprocess=0;


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
    
        tic;
    while 1 %processes until end of file is reached, then breaks
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
            outfilename=deblank(filename);
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
                      
            start=-xlimits(1)*samprate/1000;
            stop = start+150*samprate/1000;
            S=squeeze(M1OFF(1,1,:,start:stop));
            pre=squeeze(M1OFF(1,1,:,start-150*samprate/1000:start));
            [h,p]=ttest(max(abs(S), [], 2), max(abs(pre), [], 2));
            asr=mean(max(abs(S), [], 2));
            sdasr=mean(std(abs(S), [], 2));
            baseline=mean(max(abs(pre), [], 2));
            sdbaseline=std(max(abs(pre), [], 2));

            
            
            switch group
                case {'5XFAD', 'XFAD'}
                    if numgapdurs==8 %skip partial gapdur arrays
                        XFADidx=XFADidx+1;
                        XFAD(XFADidx,:)=mPeakOFF;
                        gdxfad(XFADidx,:)=gapdurs;
                        GPIASxfad(XFADidx,:)=percentGPIAS_OFF;
                        age_xfad(XFADidx)=age_days;
                        age_wks_xfad{XFADidx}=age_wks;
                        sex_xfad{XFADidx}=sex;
                        mouseID_xfad{XFADidx}=mouseID;
                        
                        P_xfad(XFADidx)=p;
                        ASR_xfad(XFADidx)=asr;
                        sdASR_xfad(XFADidx)=sdasr;
                        Baseline_xfad(XFADidx)=baseline;
                        sdBaseline_xfad(XFADidx)=sdbaseline;
                        SNR_xfad(XFADidx)=asr/baseline;
                       
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
                        
                        P_control(controlidx)=p;
                        ASR_control(controlidx)=asr;
                        sdASR_control(controlidx)=sdasr;
                        Baseline_control(controlidx)=baseline;
                        sdBaseline_control(controlidx)=sdbaseline;
                        SNR_control(controlidx)=asr/baseline;

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
        str=sprintf('processing data... %d sec elapsed, about %d sec remaining', round(etime), round(t_remain));
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


% filter out sessions with crappy pure startle responses

num_all_XFAD_sessions=size(XFAD, 1);
num_all_control_sessions=size(ctrlgroup, 1);
%number of sessions that aren't excluded due to crappy startles
%(but early sessions with different stimulus arrays are already excluded)

keep=find(P_control<1e-3);
num_excluded_control_sessions=num_all_control_sessions-length(keep);
ctrlgroup=ctrlgroup(keep,:);
gdcontrol=gdcontrol(keep,:);
GPIAScontrol=GPIAScontrol(keep,:);
age_control=age_control(keep);
age_wkscontrol=age_wkscontrol(keep);
sex_control=sex_control(keep);
mouseID_control=mouseID_control(keep);
P_control=P_control(keep);
ASR_control=ASR_control(keep);
sdASR_control=sdASR_control(keep);
Baseline_control=Baseline_control(keep);
sdBaseline_control=sdBaseline_control(keep);
SNR_control=SNR_control(keep);

keep=find(P_xfad<1e-3);
num_excluded_XFAD_sessions=num_all_XFAD_sessions-length(keep);
XFAD=XFAD(keep,:);
gdxfad=gdxfad(keep,:);
GPIASxfad=GPIASxfad(keep,:);
age_xfad=age_xfad(keep);
age_wks_xfad=age_wks_xfad(keep);
sex_xfad=sex_xfad(keep);
mouseID_xfad=mouseID_xfad(keep);
P_xfad=P_xfad(keep);
ASR_xfad=ASR_xfad(keep);
sdASR_xfad=sdASR_xfad(keep);
Baseline_xfad=Baseline_xfad(keep);
sdBaseline_xfad=sdBaseline_xfad(keep);
SNR_xfad=SNR_xfad(keep);



numXFAD_sessions=size(XFAD, 1);
numcontrol_sessions=size(ctrlgroup, 1);
numXFAD_mice=length(unique(mouseID_xfad));
numcontrol_mice=length(unique(mouseID_control));
log2gapdurs=log2(gapdurs(2:8));

%plotting params
fs=18;
ms=40; %fig 1b (used for '.')
ms1=12; %fig 2, 3a-b (used for 'o'
mso=6; %fig 3c-d (scatterplots)
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
fprintf(fid, '\n\n%d XFAD sessions excluded due to crappy pure startles', num_excluded_XFAD_sessions);
fprintf(fid, '\n%d control sessions excluded due to crappy pure startles', num_excluded_control_sessions);

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
% are negative gap detection values any more prevalent for xfad mice vs controls?
fprintf(fid, '\nAre negative gap detection values any more prevalent for xfad mice vs controls? enter the following into an online chi-square contingency table');
fprintf(fid,'\nnumber of negative-gpias control sessions: %d', length(find(GPIAScontrol<0)));
fprintf(fid,'\nnumber of positive-gpias control sessions: %d', length(find(GPIAScontrol>=0)));
fprintf(fid,'\nnumber of negative-gpias xfad sessions: %d', length(find(GPIASxfad<0)));
fprintf(fid,'\nnumber of positive-gpias xfad sessions: %d', length(find(GPIASxfad>=0)));
fprintf(fid, '\n');


%restrict to ages > p60
XFADp60_sessions=find(age_xfad>60);
controlp60_sessions=find(age_control>60);
XFADunderp60_sessions=find(age_xfad<=60);
controlunderp60_sessions=find(age_control<=60);
%important: these session indices are local to the age list
%i.e. XFADp60_sessions are the indices of age_xfad


j=0;
for i=XFADp60_sessions
    j=j+1;
    mouseid_overp60_xfad(j)=mouseID_xfad(i);
end
j=0;
for i=controlp60_sessions
    j=j+1;
    mouseid_overp60_control(j)=mouseID_control(i);
end
j=0;
for i=XFADunderp60_sessions
    j=j+1;
    mouseid_underp60_xfad(j)=mouseID_xfad(i);
end
j=0;
for i=controlunderp60_sessions
    j=j+1;
    mouseid_underp60_control(j)=mouseID_control(i);
end
fprintf(fid, '\n%d unique mice under p60 xfad ' , length(unique(mouseid_underp60_xfad)));
fprintf(fid, '\n%d unique mice over p60 control ' , length(unique(mouseid_overp60_control)));
fprintf(fid, '\n%d unique mice under p60 control ' , length(unique(mouseid_underp60_control)));
fprintf(fid, '\n%d unique mice over p60 xfad ' , length(unique(mouseid_overp60_xfad)));

fprintf(fid, '\n%d over-p60 XFAD sessions', length(XFADp60_sessions));
fprintf(fid, '\n%d over-p60 control sessions', length(controlp60_sessions));
fprintf(fid, '\n%d under-p60 XFAD sessions', length(XFADunderp60_sessions));
fprintf(fid, '\n%d under-p60 control sessions', length(controlunderp60_sessions));

    
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
    fprintf(fid, '\n\npure startles, ranksum p=%.4f', p);
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
plot(log2gapdurs, gpias, 'ok-', 'markersize', ms1, 'linewidth', lw)
gpias=GPIASxfad(1,2:8);
gd=gdxfad(1,:);
plot(log2gapdurs, gpias, 'or-', 'markersize', ms1, 'linewidth', lw)

for i=XFADp60_sessions
    gpias=GPIASxfad(i,2:8);
    gd=gdxfad(i,:);
    plot(log2gapdurs, gpias, 'or-', 'markersize', ms1, 'linewidth', lw)
end
for i=controlp60_sessions
    gpias=GPIAScontrol(i,2:8);
    gd=gdcontrol(i,:);
    plot(log2gapdurs, gpias, 'ok-', 'markersize', ms1, 'linewidth', lw)
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
    [p, anovatab, stats]=kruskalwallis(X, grouping, 'off');
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

% effect size for males vs females (over p60)
for gdindex=2:8
    X=GPIAScontrol(malep60_control_sessions,gdindex);
    Y=GPIASxfad(malep60_xfad_sessions,gdindex);
    [p,h,stats]=ranksum(X,Y);
    W=stats.ranksum;
    U=W-(length(X)*(length(X)+1))/2;
    r=1-(2*U)/(length(X)*length(Y)); %rank bi-serial correlation = effect size -1:1
    Rmale(gdindex)=r;
    fprintf(fid, '\neffect size for male over-p60, gd %d ms: %.4f', gapdurs(gdindex), r);
end
for gdindex=2:8
    X=GPIAScontrol(femalep60_control_sessions,gdindex);
    Y=GPIASxfad(femalep60_xfad_sessions,gdindex);
    [p,h,stats]=ranksum(X,Y);
    W=stats.ranksum;
    U=W-(length(X)*(length(X)+1))/2;
    r=1-(2*U)/(length(X)*length(Y)); %rank bi-serial correlation = effect size -1:1
    Rfemale(gdindex)=r;
    fprintf(fid, '\neffect size for female over-p60, gd %d ms: %.4f', gapdurs(gdindex), r);
end
    fprintf(fid, '\naverage effect size for male over-p60: %.4f', mean(Rmale(2:8)));
    fprintf(fid, '\naverage effect size for female over-p60: %.4f', mean(Rfemale(2:8)));


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
p(i+1)=plot(log2gapdurs, mgpias, 'o-c');
e(i+1)=errorbar(log2gapdurs, mgpias, sem, 'c');
set(p(i+1), 'markerfacecolor', 'c')
i=i+1;legstr{i}='male 5XFAD';
mgpias=nanmean(GPIASxfad(femalep60_xfad_sessions,2:8), 1);
sem=nanstd(GPIASxfad(femalep60_xfad_sessions,2:8), 1)./sqrt(length(femalep60_xfad_sessions));
p(i+1)=plot(log2gapdurs, mgpias, 'o-m');
e(i+1)=errorbar(log2gapdurs, mgpias, sem, 'm');
set(p(i+1), 'markerfacecolor', 'm')
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
set(p(:), 'markersize', ms1) %, 'MarkerFaceColor', 'w')
ylim([-25 100])
cd(figs_dir)
print -dpdf fig2a.pdf

%stats on sex
%convert to long form
try
    clear sex gdgroup X grouping mesgroup
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
    [p, anovatab, stats]=kruskalwallis(X, grouping, 'off');
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
    [p, anovatab, stats]=kruskalwallis(X, grouping, 'off');
    df=anovatab{2, 3};
    chisq=anovatab{2, 5};
    fprintf(fid, '\nkruskal wallis effect of xfad vs control for females:');
    fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    fprintf(fid, '\n\nover-p60 mice:');
    fprintf(fid, '\nn=%d over p-60 control mice', length(unique(mouseid_control_male))+length(unique(mouseid_control_female)));
    fprintf(fid, '\nn=%d over p-60 xfad mice', length(unique(mouseid_xfad_male))+length(unique(mouseid_xfad_female)));
    fprintf(fid, '\n');
    fprintf(fid, '\nn=%d control males >p60', length(unique(mouseid_control_male)));
    fprintf(fid, '\nn=%d xfad males >p60', length(unique(mouseid_xfad_male)));
    fprintf(fid, '\nn=%d control females >p60', length(unique(mouseid_control_female)));
    fprintf(fid, '\nn=%d xfad females >p60', length(unique(mouseid_xfad_female)));
    
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
    [p, anovatab, stats]=kruskalwallis(X, sex, 'off');
    df=anovatab{2, 3};
    chisq=anovatab{2, 5};
    fprintf(fid, '\n\nkruskal wallis effect of sex for over-p60 control mice:');
    fprintf(fid, '\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
    
    fprintf(fid, '\n\nall-ages total sessions:');
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
clear h
[~,x]=sort(age_xfad);
for i=x
    gpias=GPIASxfad(i,2:8);
    gd=gdxfad(i,:);
    h(i)=plot(log2gapdurs, gpias, 'o-');
    set(h(i), 'color', cmap(1+age_xfad(i)-minage,:))
    %      legstr{i}=sprintf('%d days',age_xfad(i)) ;
    xlabel('gap duration')
    ylabel('% gpias')
end
set(h, 'markersize', ms1, 'linewidth', lw)

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
print -dpdf fig3b.pdf


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
set(h2, 'markersize', ms1, 'linewidth', lw)

xlim([-.5 8.5])
ylim([-100 100])
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
print -dpdf fig3a.pdf

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
    p1=plot(age_control, GPIAScontrol(:,gdindex), 'ko', 'markersize', mso);
    p2=plot(age_xfad, GPIASxfad(:,gdindex), 'ro', 'markersize', mso);
    set(p1, 'markerfacecolor', 'k')
    set(p2, 'markerfacecolor', 'r')
    xlabel('age, days')
    ylabel('% gpias')
%     h=lsline;    set(h, 'linewidth', lw)
    title(sprintf('GPIAS for  %d ms gap',  gapdurs(gdindex)))
    
    Ycontrol=GPIAScontrol(:,gdindex);
    Xcontrol=[age_control(:) ones(size(age_control(:)))];
    [B,BINT,R,RINT,STATScontrol] = regress(Ycontrol,Xcontrol);
    r2control=STATScontrol(1);
    pcontrol=STATScontrol(3);
    Pcontrol(gdindex)=pcontrol;
    fprintf(fid, '\ncontrol GPIAS by age, %d ms gap, r2=%.2f,b=%.2f, p=%.4f', gapdurs(gdindex), r2control, B(1), pcontrol);
    if pcontrol<.05
%         text(105, [105 1]*B, '*', 'color', 'k', 'fontsize', 32)
    end
    
    Yxfad=GPIASxfad(:,gdindex);
    Xxfad=[age_xfad(:) ones(size(age_xfad(:)))];
    [B,BINT,R,RINT,STATSxfad] = regress(Yxfad,Xxfad);
    r2xfad=STATSxfad(1);
    pxfad=STATSxfad(3);
    Pxfad(gdindex)=pxfad;
    fprintf(fid, '\nxfad GPIAS by age, %d ms gap, r2=%.2f, b=%.2f, p=%.4f', gapdurs(gdindex), r2xfad, B(1), pxfad);
    if pxfad<.05
%         text(105, [105 1]*B, '*', 'color', 'r', 'fontsize', 32)
    end
    
    % legend(sprintf('control p=%.4f', pcontrol), sprintf('xfad p=%.4f', pxfad), 'location', 'southeast')
    
end
set(gcf, 'pos', [25        1136         975        1322])

% choose 4ms and 256ms as examples for fig3c,d
for gdindex=[2 3 4 5 6 7 8]
    figure; hold on
    set(gca, 'fontsize', fs)
    p1=plot(age_control, GPIAScontrol(:,gdindex), 'ko', 'markersize', mso);
    p2=plot(age_xfad, GPIASxfad(:,gdindex), 'ro', 'markersize', mso);
    set(p1, 'markerfacecolor', 'k')
    set(p2, 'markerfacecolor', 'r')
    
    xlabel('age, days')
    ylabel('gap detection, %')
    h=lsline;
    set(h, 'linewidth', lw)
    title(sprintf('%d ms gap',  gapdurs(gdindex)))
    legend(sprintf('control p=%.4f', Pcontrol(gdindex)), sprintf('xfad p=%.4f', Pxfad(gdindex)), 'location', 'southeast')
    cd(figs_dir)
%     if gdindex==5
%         print -dpdf fig3c.pdf
%     elseif gdindex==8
%         legend('location', 'southwest')
%         print -dpdf fig3d.pdf
%     end
end

fprintf(fid, '\n\nage range for all controls: %d - %d', min(age_control), max(age_control));
fprintf(fid, '\n\nage range for all XFAD:  %d - %d', min(age_xfad), max(age_xfad));


% try some sort of log regression since it isn't linear
warning off stats:nlinfit:IllConditionedJacobian
 bigfig=figure;
 orient tall
 subplot1(4,2);ph=0;
for gdindex=[2 3 4 5 6 7 8]
    littlefig=figure;
    hold on
    xdatacontrol=age_control(:);
    xdataxfad=age_xfad(:);
    ydatacontrol=GPIAScontrol(:,gdindex);
    ydataxfad=GPIASxfad(:,gdindex);

    p1=plot(age_control, (ydatacontrol), 'ko', 'markersize', 12);
    p2=plot(age_xfad, (ydataxfad), 'ro', 'markersize', 12);
    set(p1, 'markerfacecolor', 'k')
    set(p2, 'markerfacecolor', 'r')
    xlabel('age, days')
    ylabel('gap detection, %')
    ylim([-20 100])
    title(sprintf('%d ms gap',  gapdurs(gdindex)))
    legend('control', '5XFAD')
    
    figure(bigfig)
    ph=ph+1;    
    subplot1(ph);
    p1=plot(age_control, (ydatacontrol), 'ko', 'markersize', 8);
    p2=plot(age_xfad, (ydataxfad), 'ro', 'markersize', 8);
    set(p1, 'markerfacecolor', 'k')
    set(p2, 'markerfacecolor', 'r')
    xlabel('age, days')
    ylabel('% gpias')
    ylim([-20 100])
    title(sprintf('%d ms gap',  gapdurs(gdindex)))
    legend('control', '5XFAD')
    
    %fit to myfun which is a*log(b*x+c)+d
    %also try myalphafun
    x0=[ 20 100 20 100 ];
    [xcontrol,Rcontrol,J,sigmacontrol,MSEcontrol] = nlinfit(xdatacontrol,ydatacontrol, @myalphafun, x0);
    [xxfad,Rxfad,J,sigmaxfad,MSExfad] = nlinfit(xdataxfad,ydataxfad, @myalphafun, x0);
    Betacontrol(gdindex,:)=xcontrol;
    Betaxfad(gdindex,:)=xxfad;
    
    %compare to linear regression
    xdatacontrol_2=[age_control(:) ones(size(age_control(:)))];
    [B,BINT,R,RINT,STATScontrol] = regress(ydatacontrol,xdatacontrol_2);
    MSEcontrol_linear=STATScontrol(4);
    if MSEcontrol_linear<MSEcontrol
        fprintf('\n%dms: linear regression (%d) is a better fit than a*log(b*x+c) (%d)', gapdurs(gdindex), round(MSEcontrol_linear), round(MSEcontrol));
    else
        fprintf('\n%dms:  g*(e(-t/tau1) - e(-t/tau2)) (%d) is a better fit than linear regression (%d)', gapdurs(gdindex), round(MSEcontrol), round(MSEcontrol_linear));
    end
    
    
    xfitcontrol=round(xcontrol(1)): max(xdatacontrol); %xl=xlim
    xfitxfad=round(xxfad(1)): max(xdataxfad); %xl=xlim
    yhatcontrol=myalphafun(xcontrol, xfitcontrol);
    yhatxfad=myalphafun(xxfad, xfitxfad);

    figure(bigfig)
    p=plot(xfitcontrol, yhatcontrol, 'k-', xfitxfad, yhatxfad, 'r');
    set(p, 'linewidt', 2)
    xlim([min([xdatacontrol; xdataxfad])-10 max([xdatacontrol; xdataxfad])+10])
    
    figure(littlefig)
    p=plot(xfitcontrol, yhatcontrol, 'k-', xfitxfad, yhatxfad, 'r');
    set(p, 'linewidt', 2)
    xlim([min([xdatacontrol; xdataxfad])-10 max([xdatacontrol; xdataxfad])+10])
    
    if gdindex==7
        set(gca, 'fontsize', fs)
        ylim([-10 100])
        legend('control', '5XFAD', 'location', 'southeast')
        print -dpdf fig3c.pdf
    elseif gdindex==8
        ylim([-10 100])
        set(gca, 'fontsize', fs)
        legend('control', '5XFAD', 'location', 'southeast')
        print -dpdf fig3d.pdf
    end
    
    %examine the 95% confidence intervals for the log fit a*log(b*x+c)
    CIcontrol = nlparci(xcontrol,Rcontrol,'covar',sigmacontrol);
    CIxfad = nlparci(xxfad,Rxfad,'covar',sigmaxfad);
        
end

 bigfig=figure;
 orient tall


fprintf(fid, '\nfitted betas for dual exponential:');
fprintf(fid, '\nControl ');
for gdindex=[2 3 4 5 6 7 8]
fprintf(fid, '\ngd=%dms: \tt_onset: %d \ttau1:%d \ttau2:%d \tgmax: %d', gapdurs(gdindex), round(Betacontrol(gdindex,:)));
end
fprintf(fid, '\n\n5XFAD ');
for gdindex=[2 3 4 5 6 7 8]
fprintf(fid, '\ngd=%dms: \tt_onset: %d \ttau1:%d \ttau2:%d \tgmax: %d', gapdurs(gdindex), round(Betaxfad(gdindex,:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Make Table 1
%this is a table of linear regression params,
% let's restrict it to >p60 and then we can say "significantly deteriorated
% with age over p60'
% % allp60_control_sessions=[malep60_control_sessions femalep60_control_sessions];
% % allp60_xfad_sessions=[malep60_xfad_sessions femalep60_xfad_sessions];

%I see now that these are almost duplicates of XFADp60_sessions controlp60_sessions
%you might think that they should match
% but they don't, bc XFADp60_sessions are indices of age_xfad
%whereas allp60_xfad_sessions are indices of XFADp60_sessions
%confusing, right?

%try using XFADp60_sessions instead of allp60_xfad_sessions - I think
%XFADp60_sessions is correct - for some reason it doesn't make a
%difference.

fprintf(fid, '\n%d unique mice under p60 xfad ' , length(unique(mouseid_underp60_xfad)));
fprintf(fid, '\n%d unique mice over p60 control ' , length(unique(mouseid_overp60_control)));
fprintf(fid, '\n%d unique mice under p60 control ' , length(unique(mouseid_underp60_control)));
fprintf(fid, '\n%d unique mice over p60 xfad ' , length(unique(mouseid_overp60_xfad)));

fprintf(fid, '\n%d over-p60 XFAD sessions', length(XFADp60_sessions));
fprintf(fid, '\n%d over-p60 control sessions', length(controlp60_sessions));
fprintf(fid, '\n%d under-p60 XFAD sessions', length(XFADunderp60_sessions));
fprintf(fid, '\n%d under-p60 control sessions', length(controlunderp60_sessions));

fidtable=fopen('Table 1.txt', 'w');
fprintf(fidtable,'\nThis is a table of linear regression slopes for age > p60 %s\n\n', datestr(datenum(now)));

fprintf(fidtable, '\t\t slope \t r^2 \t p \t n');
fprintf(fidtable, '\n5XFAD');
figure;subplot1(4, 2);
for gdindex=2:8
    Yxfad=GPIASxfad(XFADp60_sessions,gdindex);
    Xxfad=[age_xfad(XFADp60_sessions)' ones(size(age_xfad(XFADp60_sessions)))'];
subplot1(gdindex)
plot(Xxfad(:,1), Yxfad, 'ro');
% h=lsline;set(h, 'color', 'r')
lsline
    [B,BINT,R,RINT,STATSxfad] = regress(Yxfad,Xxfad);
    slope_xfad=B(1);
    Slopes_xfad(gdindex-1)=slope_xfad;
    r2xfad=STATSxfad(1);
    pxfad=STATSxfad(3);
    nsess=length(XFADp60_sessions);
    nmice=length(unique(mouseid_overp60_xfad));
    fprintf(fidtable, '\n\t%d ms \t %.2f \t %.2f \t p=%.4f \t %d sessions, %d mice',...
        gapdurs(gdindex), slope_xfad, r2xfad, pxfad, nsess, nmice);
end
fprintf(fidtable, '\ncontrol');
for gdindex=2:8
    Ycontrol=GPIAScontrol(controlp60_sessions,gdindex);
    Xcontrol=[age_control(controlp60_sessions)' ones(size(age_control(controlp60_sessions)))'];
   subplot1(gdindex)
plot(Xcontrol(:,1), Ycontrol, 'ko');
lsline
 [B,BINT,R,RINT,STATScontrol] = regress(Ycontrol,Xcontrol);
    slope_control=B(1);
    Slopes_control(gdindex-1)=slope_control;
    r2control=STATScontrol(1);
    pcontrol=STATScontrol(3);
    nsess= length(controlp60_sessions);
    nmice= length(unique(mouseid_overp60_control));
    fprintf(fidtable, '\n\t%d ms \t %.2f \t %.2f \t p=%.4f \t %d sessions, %d mice',...
        gapdurs(gdindex), slope_control, r2control, pcontrol, nsess, nmice);
    
end
% fclose(fidtable);
subplot1(1)
title('over p60 linear regressions')
fprintf(fid,'\n\nSee Table 2.txt for a table of linear regression slopes (for age > p60)\n')
p=ranksum(Slopes_control, Slopes_xfad, 'tail', 'right');
fprintf(fid,'\nlinear regression slopes of gap detection vs age (>p60) were significantly less for xfad mice compared to controls');
fprintf(fid,'\np=%4f, 1-tailed ranksum', p);


%%%%%%%%%%%%%%%%%%%%
% repeat for young mice (<p60)


fprintf(fidtable,'\n\n\nBelow are values of linear regression slopes for age <= p60\n\n');

fprintf(fidtable, '\t\t slope \t r^2 \t p \t n');
fprintf(fidtable, '\n5XFAD');
figure;subplot1(4, 2);
for gdindex=2:8
    Yxfad=GPIASxfad(XFADunderp60_sessions,gdindex);
    Xxfad=[age_xfad(XFADunderp60_sessions)' ones(size(age_xfad(XFADunderp60_sessions)))'];
subplot1(gdindex)
plot(Xxfad(:,1), Yxfad, 'ro');
% h=lsline;set(h, 'color', 'r')
lsline
    [B,BINT,R,RINT,STATSxfad] = regress(Yxfad,Xxfad);
    slope_xfad=B(1);
    Slopes_xfad(gdindex-1)=slope_xfad;
    r2xfad=STATSxfad(1);
    pxfad=STATSxfad(3);
    nsess=length(XFADunderp60_sessions);
    nmice=length(unique(mouseid_underp60_xfad));
    fprintf(fidtable, '\n\t%d ms \t %.2f \t %.2f \t p=%.4f \t %d sessions, %d mice',...
        gapdurs(gdindex), slope_xfad, r2xfad, pxfad, nsess, nmice);
end
fprintf(fidtable, '\ncontrol');
for gdindex=2:8
    Ycontrol=GPIAScontrol(controlunderp60_sessions,gdindex);
    Xcontrol=[age_control(controlunderp60_sessions)' ones(size(age_control(controlunderp60_sessions)))'];
   subplot1(gdindex)
plot(Xcontrol(:,1), Ycontrol, 'ko');
lsline
 [B,BINT,R,RINT,STATScontrol] = regress(Ycontrol,Xcontrol);
    slope_control=B(1);
    Slopes_control(gdindex-1)=slope_control;
    r2control=STATScontrol(1);
    pcontrol=STATScontrol(3);
    nsess=length(controlunderp60_sessions);
    nmice=length(unique(mouseid_underp60_control));
    fprintf(fidtable, '\n\t%d ms \t %.2f \t %.2f \t p=%.4f \t %d sessions, %d mice',...
        gapdurs(gdindex), slope_control, r2control, pcontrol, nsess, nmice);
    
end
subplot1(1)
title('under p60 linear regressions')

p=ranksum(Slopes_control, Slopes_xfad);
fprintf(fid,'\nlinear regression slopes of gap detection vs age (<=p60), xfad compared to controls');
fprintf(fid,'\np=%4f, 2-tailed ranksum', p);
fprintf(fid,'\nmean slopes for control <=p60: %.2f', Slopes_control);
fprintf(fid,'\nmean slopes for xfad <=p60: %.2f', Slopes_xfad);

fclose(fidtable);
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
        gdms(j)=gapdurs(gd);
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
        gdms(j)=gapdurs(gd);
        grouping{j}='xfad';
        age_long(j)=age_xfad(i);
        sex_long{j}=sex_xfad{i};
        mouseID_long{j}=mouseID_xfad{i};
    end
end

%bin age and do stats
age_min=min(age_long);
% for age_max=45:5:80;
%     i=find(age_long<age_max & age_long>age_min);
%     [p, anovatab, stats]=kruskalwallis(X(i), grouping(i), 'off');
%     df=anovatab{2, 3};
%     chisq=anovatab{2, 5};
%     fprintf(fid,'\nkruskal wallis main effect of xfad vs control, age range %d - %d:', age_min, age_max);
%     fprintf(fid,'\np=%g, d.f.=%d, chi-squared=%.2f ', p, df, chisq);
%     fprintf(fid,'\n%d xfad sessions', length(strmatch('xfad', grouping(i)))/numgapdurs);
%     k=strmatch('xfad', grouping(i));
%     fprintf(fid,'\n%d xfad mice',    length(unique(mouseID_long(k))));
%     
%     fprintf(fid,'\n%d control sessions', length(strmatch('control', grouping(i)))/numgapdurs);
%     j=strmatch('control', grouping(i));
%     fprintf(fid,'\n%d control mice',    length(unique(mouseID_long(j))));
%     fprintf(fid,'\n');
% end



% repeat but only doing a 1-tailed rank sum on 256 ms
fprintf(fid,'\n only doing a 1-tailed rank sum on 256 ms');
for age_max=40:10:100;
    age_min=age_max-10;
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
    fprintf(fid,'\n%d xfad mice',    length(unique(mouseID_long(i(m)))));
    
    fprintf(fid,'\n%d control sessions', length(Xcontrol));
    n=strmatch('control', grouping(i));
    fprintf(fid,'\n%d control mice',    length(unique(mouseID_long(i(n)))));
    fprintf(fid,'\n');
end



% do a nested (mixed-effects) GLM to jpintly test the effects of age, sex, genotype, and mouseID

%mean-center age and gapdur
% age_centered=age_long-mean(age_long);
% gd_centered=gdms-mean(gdms);

T=table(grouping(:), gdgroup(:), age_long(:), sex_long(:), mouseID_long(:), X(:));
%T=table(grouping(:), gd_centered(:), age_centered(:), sex_long(:), mouseID_long(:), X(:));
T.Properties.VariableNames={'genotype', 'gapdur', 'age', 'sex', 'mouseID', 'GPIAS'};

% does the centering matter?
%swapping between the T definitions above makes a big difference
% it changes whether genotype, genotype:gapdur are significant
%(they are significant without centering, non-significant with centering)
% the coefficients are also wildly different (changed signs, values change by
% orders of magnitude)

%'genotype', 'gapdur', 'age', 'sex', 'mouseID'
fprintf(fid,'full glme (age and gapdur are NOT mean-centered)):\n');
  
fullformula=['GPIAS ~  genotype + gapdur + age + sex    ', ... 
    '  + genotype:gapdur  ', ...
    '  + genotype:age  ', ...
    '  + genotype:sex  ', ...
    '  + gapdur:age  ', ...
    '  + gapdur:sex  ', ...
    '  + age:sex  ', ...
    '+ (1 + gapdur | mouseID) + (1 + age | mouseID)'];
fprintf('\nfitting glme...')
fullglme=fitglme(T, fullformula);
% anova(fullglme)
fprintf('done\n')




% I could specify the distributions (normal, poisson, binomial, etc)

fprintf('\n\npiecewise linear glme (age early and age late):\n');
%how to use a piecewise linear fit for age?
%I will use 60 days as the inflection point ("knot" or "break-point")
% https://www.lexjansen.com/pharmasug-cn/2015/ST/PharmaSUG-China-2015-ST08.pdf
bk=60; %breakpoint
% NO mean-centering
%k=bk-mean(age_long); %adjust for mean centering
age_early=zeros(size(age_long(:)));
age_late=zeros(size(age_long(:)));
%age_centered_minus_k=age_centered-k;
% age_early(find(age_centered<=k))=age_centered_minus_k(find(age_centered<=k));
% age_late(find(age_centered>k))=age_centered_minus_k(find(age_centered>k));
age_early(find(age_long<=bk))=age_long(find(age_long<=bk));
age_late(find(age_long>bk))=age_long(find(age_long>bk));



%T2=table(grouping(:), find(age_long<=bk)(:), age_early(:), age_late(:), sex_long(:), mouseID_long(:), X(:));
T2=table(grouping(:), gdgroup(:), age_early(:), age_late(:), sex_long(:), mouseID_long(:), X(:));
T2.Properties.VariableNames={'genotype', 'gapdur', 'age_early', 'age_late', 'sex', 'mouseID', 'GPIAS'};

piecewiseformula=['GPIAS ~  genotype + gapdur + age_early + age_late + sex    ', ... 
    '  + genotype*gapdur  ', ...
    '  + genotype*age_early  ', ...
    '  + genotype*age_late  ', ...
    '  + genotype*sex  ', ...
    '  + gapdur*age_early  ', ...
    '  + gapdur*age_late  ', ...
    '  + gapdur*sex  ', ...
    '  + age_early*sex  ', ...
    '  + age_late*sex  ', ...
    '+ (1 + gapdur | mouseID) + (1 + age_early | mouseID) + (1 + age_late | mouseID)'];
fprintf('\nfitting piecewise glme...')
piecewiseglme=fitglme(T2, piecewiseformula);
% anova(piecewiseglme)
fprintf('done\n')

figure
subplot(321)
plotPartialDependence(piecewiseglme, 1)
subplot(322)
plotPartialDependence(piecewiseglme, 2)
subplot(323)
plotPartialDependence(piecewiseglme, 3)
subplot(324)
plotPartialDependence(piecewiseglme, 4)
subplot(325)
plotPartialDependence(piecewiseglme, 5)
subplot(326)
plotPartialDependence(piecewiseglme, 6)



% mdl=fitglm(T, 'linear')
% mdl=fitglm(T, 'interactions')
% 
% mdl=stepwiseglm(T, 'linear')
% plotSlice(mdl)
% figure
% plotDiagnostics(mdl)
% 
% formula=['GPIAS ~ genotype + age + gapdur + sex +  ', ... 
%     'genotype*age + genotype*sex  + genotype*gapdur*age ']
%     
% feglme=fitglme(T, formula)
% anova(feglme)
% 
% formula=['GPIAS ~ genotype + age + gapdur + sex +  ', ... 
%     'genotype*age + genotype*sex  + genotype*gapdur*age + ', ...
%     '(1 | mouseID)']
% glme=fitglme(T, formula)
% anova(glme)
% 
% %compare 2 models
% results = compare(feglme,glme,'CheckNesting',true)
% 
% %the idea here (like the hiker ditching a gallon of water, http://www.bodowinter.com/tutorial/bw_LME_tutorial2.pdf
% %is to compare two models that ditch a single parameter, and see if it
% %makes a significant difference. If it does, keep it.
% 
% % should we try a model with random slopes in addition to random intercepts?
% formula1=['GPIAS ~  age + gapdur    ', ... 
%     '   + genotype*age + ', ...
%     '(1 + gapdur| mouseID)']
% glme1=fitglme(T, formula1)
% 
% 
% 
% results = compare(glme1,glme2,'CheckNesting',true) %the second model is supposed to have one additional param
% 
% anova(glme2)
% genotype*gapdur*age matters
% genotype*age matters
% gapdur matters
% age matters
%
% genotype doesn't matter
% sex doesn't matter
% genotype*sex doesn't matter

% adding in random slope (1 + GPIAS | mouseID) I might need to check all
% the terms again

% I don't understand how to compare models re: the interaction terms, since
% the individual terms seem redundant with the interactions (same numbers,
% p=0)

%notes: since I am using linear/normal glme, it's really just a linear mixed model
% https://stats.idre.ucla.edu/other/mult-pkg/introduction-to-linear-mixed-models/

%after talking to sanjay:
%centering and/or scaling the data is really important to be able to
%interpret the coefficients. Especially when there are interactions. In
%that case the coefficient is telling you what the value of one thing is
%when the other is zero. This is easiest to interpret when zero is a
%meaningful place in your data, like e.g. the middle. This is why you
%mean-center your data. For example, if there is an age*genotype
%interaction, the coefficient means what the effect of genotype is when age
%is zero. This isn't very helpful because p0 is outside our data range. If
%you mean-center age, then the coefficient tells you the effect of genotype
%for the average-aged mouse.    
%it's also worth thinking about how gapdur is scaled. Right now I'm using
%an index, which is sort of like a log2 transform (but not quite). It might
%make more sense to use actual gapdur in ms, or log2(gapdur).
%
%model selection is not cut and dried. There are clearly multiple
%defensible ways to approach it. One philosophy is to start minimally and
%build up the model, running a compare each time to see if it's
%significantly different and lowers the AIC/BIC. (those include a penalty
%for model complexity, but differ on weighting). This includes checking fixed and
%random effects. The 2 models being compared have to be nested (one has 1
%additional effect than the other).
%Another philosophy is to just use the full model with all the
%interactions. The resulting output may look complicated but it's not
%unlike an anova where it reports all the possible interactions. Terms that
%have coefficients close to zero can just be ignored. Or alternatively you
%can focus only on the terms that correspond to your hypotheses and ignore
%other ones. Sometime the full model won't converge, in which case you are
%fully justified in dropping terms that allow you to get the model to
%converge.
%So far this is all linear. Probably the dependence on age is nonlinear
%(development + degeneration). If there was a reason to suspect a certain
%form like quadratic, that would be fine, but usually there isn't. In that
%case (i.e. in our case) a specific form might give you a bump facing the
%right direction but would be unlikely to be a very good fit for the data.
%So a piecewise linear fit is probably the best idea. Any a priori reason
%to impose a kink (e.g. Abeta plaques start at p60) is helpful, but it's
%probably OK to pick a transition point based on the data. To do this you
%would create a new variable (call it age2) that is defined as one slope
%before p60 and another after p60. This is widely done and there may even
%be some matlab tutorial on it.
%as for random slopes, it seems like age and gapdur both might reasonbly
%show random slope effects, since a given mouse might age faster or show a
%different gap-detection curve.


figure
plotResiduals(piecewiseglme)
plotResiduals(piecewiseglme,'fitted')
plotResiduals(piecewiseglme,'probability')
plotResiduals(piecewiseglme,'histogram','ResidualType','Pearson')


fclose(fid);

diary(fname)
anova(fullglme)
fullglme

piecewiseformula
anova(piecewiseglme)
piecewiseglme

diary off

diary('Table 2.txt')
piecewiseformula
anova(piecewiseglme)
piecewiseglme
diary off

type(fname)



if 1
    cd(figs_dir)
    delete ADfigs.ps
    f=findobj('type', 'figure');
    for i=1:length(f)
        figure(i)
        print -bestfit -dpsc2 ADfigs.ps -append
    end
end

%  keyboard



