N = 100000;
v = -100:100;
vt = 50;
nbin = 30;

cdfy = -0.5 + 1*rand(N,1);
cdfx = erfinv(2*cdfy)*vt;
pdfnorm = vt*randn(N,1)/sqrt(2);

nPlot = 4;
nRows = 4;
nCols = 1;
for k = 1:nPlot; h(k) = subplot(nRows,nPlot/nRows,k); end
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,x,0.5*erf(x))

hca = h(isub); isub = isub + 1;
hist(hca,cdfy,nbin)

hca = h(isub); isub = isub + 1;
hist(hca,cdfx,nbin)

hca = h(isub); isub = isub + 1;
hist(hca,[cdfx pdfnorm],nbin)
legend(hca,'cdfx','randn')