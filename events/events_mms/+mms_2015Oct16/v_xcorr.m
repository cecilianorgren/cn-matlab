tint = irf.tint('2015-10-16T10:33:30.00Z/2015-10-16T10:33:30.60Z');
tint = irf.tint('2015-10-16T10:33:29.00Z/2015-10-16T10:33:31.00Z');
tint = irf.tint('2015-10-16T10:33:28.00Z/2015-10-16T10:33:32.00Z');

[V, dV] = c_4_v_xcorr(tint,gseVe1,gseVe2.resample(gseVe1),gseVe3.resample(gseVe1),gseVe4.resample(gseVe1),...
  gseR1.resample(gseVe1),gseR2.resample(gseVe1),gseR3.resample(gseVe1),gseR4.resample(gseVe1));
%%
[V, dV] = c_4_v_xcorr(tint,gseB1,gseB2.resample(gseB1),gseB3.resample(gseB1),gseB4.resample(gseB1),...
  gseR1.resample(gseB1),gseR2.resample(gseB1),gseR3.resample(gseB1),gseR4.resample(gseB1));


%%
figure
nRows = 4;
nCols = 1;
for ii = 1:nRows*nCols; h(ii) = subplot(nRows,nCols,ii); end

isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,data1)

hca = h(isub); isub = isub + 1;
plot(hca,data2)

hca = h(isub); isub = isub + 1;
plot(hca,corr)

hca = h(isub); isub = isub + 1;
plot(hca,corr_plus)