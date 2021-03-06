function SpiketrainCorrelationAll
%simple script that compares the MClust outfiles to the Kilosort outfiles
%in the current directory (all possible combinations, so it takes a while
%to run)
%plots a correlation matrix of r values and the cross-correlations 
clear textprogressbar

mcdir=dir('outPSTH*.mat');
ksdir=dir('KS_outPSTH*.mat');

nmc=length(mcdir);
nks=length(ksdir);

fprintf('\n\n%s\n%s \n%d mclust outfiles \n%d kilsort outfiles', mfilename, pwd, nmc, nks)

fprintf('\n')
textprogressbar('loading mclust outfiles ');
for i=1:nmc
    load(mcdir(i).name);
    textprogressbar(100*i/nmc);
    mcspiketimes(i).spiketimes=out.spiketimes;
end
textprogressbar(' ') %terminate
textprogressbar('loading kilosort outfiles ')
for i=1:nks
    load(ksdir(i).name);
    textprogressbar(100*i/nks);
    ksspiketimes(i).spiketimes=out.spiketimes;
end
textprogressbar(' ') %terminate

 for i=1:nmc
     min_all(i)=min(mcspiketimes(i).spiketimes);
     max_all(i)=max(mcspiketimes(i).spiketimes);
 end
 for i=1:nks
     min_all2(i)=min(ksspiketimes(i).spiketimes);
     max_all2(i)=max(ksspiketimes(i).spiketimes);
 end
 xlimits=[min([min_all, min_all2]) max([max_all, max_all2])]*1000;
 
 
textprogressbar('smoothing mclust outfiles ')
for i=1:nmc
    textprogressbar(100*i/nmc);
    [t, frmc(i,:)]=GaussSmooth(mcspiketimes(i).spiketimes*1000, 1, xlimits);
end
textprogressbar(' ') %terminate
textprogressbar('smoothing kilosort outfiles ')
for i=1:nks
    textprogressbar(100*i/nks);
    [t, frks(i,:)]=GaussSmooth(ksspiketimes(i).spiketimes*1000, 1, xlimits);
end
frmc=frmc(:,1:end-100);
frks=frks(:,1:end-100);
textprogressbar(' ') %terminate


textprogressbar('calculating correlation coefficients ');
p=0;
for i=1:nmc
    for j=1:nks
        p=p+1;
        textprogressbar(100*p/(nmc*nks));
        r=corrcoef(frmc(i,:), frks(j,:));
        r=r(2);
        C(i,j)=r;
    end
end
textprogressbar(' ') %terminate

figure
imagesc(C)
% set(gca, 'ydir', 'normal', 'xtick', 1:nks, 'yticks', 1:nmc) 
set(gca,  'xtick', 1:nks, 'ytick', 1:nmc) 
xl=xlim;yl=ylim;
t1=text( -0.075 , .5, 'mclust cells', 'units', 'normal', 'horizontalalign', 'center');
 set(t1, 'rotation', 90, 'fontsize', 18)
t2=text(.5, -.075, 'kilosort cells', 'units', 'normal', 'horizontalalign', 'center');
 set(t2, 'fontsize', 18)
colorbar
print -dpdf 'mclust-kilosort-correlation-matrix.pdf'

textprogressbar('calculating cross correlations ');
p=0;
for i=1:nmc
    for j=1:nks
        p=p+1;
        textprogressbar(100*p/(nmc*nks));
        [xc, lags]=xcorr(frmc(i,:), frks(j,:), 256);
        XC(i,j,:)=xc;
    end
end
textprogressbar(' ') %terminate

figure;
subplot1(nmc, nks)
p=0;
for i=1:nmc
    for j=1:nks
        p=p+1;
        subplot1(p);
        plot(lags, squeeze(XC(i,j,:)), 'r', 'linewidth', 2)
        yl(p,:)=ylim;
    end
end
YL=[min(yl(:,1)) max(yl(:,2))];
p=0;
for i=1:nmc
    for j=1:nks
        p=p+1;
        subplot1(p);
        ylim(YL)
        line([0 0], [ylim], 'linestyle', '--')
        line([xlim], [0 0])
        text(0, .9*YL(2), sprintf('%.3f',C(i,j)), 'horizontalalign', 'center')
    end
end

%label 
p=0;
for i=1:nmc
    for j=1:nks
        p=p+1;
        subplot1(p);
        mcfn=mcdir(i).name;
        mcfn=erase(mcfn, 'outPSTH_ch');
        mcfn=erase(mcfn, '.mat');
        mcfn=['m', mcfn];

        ksfn=ksdir(j).name;
        ksfn=erase(ksfn, 'KS_outPSTH_ch');
        ksfn=erase(ksfn, '.mat');
        ksfn=['k', ksfn];

        if i==nmc
           xh= xlabel(ksfn);
            set(xh, 'rotation', 20)

        end
        if j==1
           yh= ylabel(mcfn);
        end
    end
end
    orient tall
print -dpdf 'mclust-kilosort-cross-correlations.pdf'

for i=1:nmc
    for j=1:nks
Names{i,j}= [mcdir(i).name, ', ', ksdir(j).name];
    end
end

fid=fopen('mclust-kilosort-correlation-values.txt', 'w');
fprintf(fid, 'sorted from highest correlation to lowest:');
[a,I]=sort(C(:), 'descend');

for i=1:length(a)
fprintf(fid, '\n\n%s, %s', Names{I(i)});
fprintf(fid, '\nr=%.4f', a(i));
end

fprintf(fid, '\n\n===================================\n\nsorted numerically:');
for i=1:nmc
    for j=1:nks
fprintf(fid, '\n\n%s, %s', mcdir(i).name,  ksdir(j).name);
fprintf(fid, '\nr=%.4f, MClust: %d spikes, Kilosort: %d spikes', C(i,j), length(mcspiketimes(i).spiketimes), length(ksspiketimes(j).spiketimes));
    end
end
fclose(fid)

% save mclust-kilosort-correlation
%file is too big
