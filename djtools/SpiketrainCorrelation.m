 [filename1, pathname1] = uigetfile('out*.mat',   'Pick first outfile (MClust)');
 [filename2, pathname2] = uigetfile('KS_out*.mat',   'Pick first outfile (Kilosort)');
 cd(pathname1)
 out=load(filename1);
 cd(pathname2)
 KSout=load(filename2);
 
 
    figure; hold on
xlimits1=min([out.out.spiketimes KSout.out.spiketimes])*1000;
xlimits2=max([out.out.spiketimes KSout.out.spiketimes])*1000;
xlimits=[xlimits1 xlimits2];

[t, fr]=GaussSmooth(out.out.spiketimes*1000, 1, xlimits);
%plot(t, fr)

[tKS, frKS]=GaussSmooth(KSout.out.spiketimes*1000, 1, xlimits);
%plot(tKS, frKS)

r=corrcoef(fr, frKS);
r=r(2);

[xc, lags]=xcorr(fr, frKS, 256);
figure;
plot(lags, xc)
title(sprintf('r = %.4f', r))
text(lags(1), min(xc), [filename1, '     ',  filename2], 'interpreter', 'none')




