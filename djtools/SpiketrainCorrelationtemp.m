%plots a cross-correlation and prints out the correlation coefficient

 
          c=1;
          st1=sp.st(sp.clu == sp.cids(good_cells1(c)));
            st1 = st1(st1> start & st1 < L(silence_dir_indx));
            st1 = st1 - start;

            c=2;
          st2=sp.st(sp.clu == sp.cids(good_cells1(c)));
            st2 = st2(st2> start & st2 < L(silence_dir_indx));
            st2 = st2 - start;

    figure; hold on
xlimits1=0;
xlimits2=max([st1(:); st2(:)])*1000;
xlimits=[xlimits1 xlimits2];

[t, fr]=GaussSmooth(st1*1000, 1, xlimits);
%plot(t, fr)
 
[t2, fr2]=GaussSmooth(st2*1000, 1, xlimits);
%plot(tKS, frKS)

% r=corrcoef(fr, fr2);
% r=r(2);

[xc, lags]=xcorr(fr, fr2, 256);
figure;
plot(lags, xc)
title(sprintf('r = %.4f', r))
text(lags(1), min(xc), [filename1, '     ',  filename2], 'interpreter', 'none')




