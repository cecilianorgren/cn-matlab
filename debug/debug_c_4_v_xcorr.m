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