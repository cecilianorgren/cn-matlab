ic = 1;
fftTint = irf.tint('2015-11-30T00:24:00.00Z/2015-11-30T00:25:00.00Z'); 
units = irf_units;

flim = 5;
c_eval('gseB0? = gseB?.filt(0,1,[],5);',ic)

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('[gseB?scmpar,gseB?scmperp] = irf_dec_parperp(gseB?,gseB?scm); gseB?scmpar.name = ''B scm par''; gseB?scmperp.name = ''B scm perp'';',ic)

Exyzfac = irf_convert_fac(Exyz,Bxyz,[1 0 0]);
Bscmfac = irf_convert_fac(Bscm,Bxyz,[1 0 0]);


c_eval('phiB?par = gseB0?.resample(gseB?scmpar).abs*gseB?scmpar.filt(flim,0,[],5)*1e-18/(ne?.filt(0,1,[],5)*1e6*units.e*units.mu0);',ic)
c_eval('phiB?perp = gseB0?.resample(gseB?scmperp).abs*gseB?scmperp.filt(flim,0,[],5)*1e-18/(ne?.filt(0,1,[],5)*1e6*units.e*units.mu0);',ic)
c_eval('intE? = irf_integrate(gseE?perp.abs.filt(flim,0,[],5));',ic)



fsE = 1/(gseE1.time(2)-gseE1.time(1));
fsB = 1/(gseB1scm.time(2)-gseB1scm.time(1));

nfft = floor(log2(gseE1perp.length));
nfft = 2^ceil(log2((fftTint(2)-fftTint(1))*fsE));
nfftTint = fftTint(1)+[0 nfft*1.1/fsE];
indT0 = find(gseE1perp.time>fftTint(1),1,'first');

% mV/m*s*km/s = mV*km/m = Vm/m
c_eval('tic; pfftE? = irf_powerfft(gseE?perp(indT0+[0 nfft]).abs,nfft,fsE,0); toc',ic)
c_eval('tic; pfftB? = irf_powerfft(gseB?scmpar(indT0+[0 nfft]).abs,nfft,fsB,0); toc',ic)

nfft = 2^14;
c_eval('tic; pfftE? = irf_powerfft(gseE?perp.abs,nfft,fsE,0); toc',ic)
c_eval('tic; pfftB?par = irf_powerfft(gseB?scmpar.abs,nfft,fsB,0); toc',ic)
c_eval('tic; pfftB?perp = irf_powerfft(gseB?scmperp.abs,nfft,fsB,0); toc',ic)

c_eval('tic; pfftPhiB?par = irf_powerfft(phiB?par.abs,nfft,fsB,0); toc',ic)
c_eval('tic; pfftPhiB?perp = irf_powerfft(phiB?perp.abs,nfft,fsB,0); toc',ic)
c_eval('tic; pfftIntE? = irf_powerfft(intE?.abs,nfft,fsB,0); toc',ic)
%%
v = 4000;
c_eval('pfftPhiE? = pfftE?; scaleE = v*repmat(torow(1./pfftPhiE?.f),length(pfftPhiE?.t),1); pfftPhiE?.p{1}=pfftPhiE?.p{1}.*scaleE;',ic)

c_eval('pfftAvB?par = pfftB?par; pfftAvB?par.p{1} = sum(pfftAvB?par.p{1});',ic)
c_eval('pfftAvPhiB?par = pfftPhiB?par; pfftAvPhiB?par.p{1} = sum(pfftAvPhiB?par.p{1});',ic)
c_eval('pfftAvPhiB?perp = pfftPhiB?perp; pfftAvPhiB?perp.p{1} = sum(pfftAvPhiB?perp.p{1});',ic)
c_eval('pfftAvPhiE? = pfftPhiE?; pfftAvPhiE?.p{1} = sum(pfftAvPhiE?.p{1});',ic)
c_eval('pfftAvIntE? = pfftIntE?; pfftAvIntE?.p{1} = sum(pfftAvIntE?.p{1});',ic)
c_eval('pfftAvE? = pfftE?; pfftAvE?.p{1} = sum(pfftAvE?.p{1});',ic)

loglog(pfftAvPhiB1par.f, pfftAvPhiB1par.p{1},...
       pfftAvPhiB1perp.f,pfftAvPhiB1perp.p{1},...
       pfftAvPhiE1.f,pfftAvPhiE1.p{1},...
       pfftAvIntE1.f,pfftAvIntE1.p{1}*v,...
       pfftAvB1par.f,pfftAvB1par.p{1},...
       pfftAvE1.f,pfftAvE1.p{1})
legend('PhiBpar','PhiBperp','PhiE','IntE*v','BnT','avE'  )

%% Test dispersion relation
v = 4000;
c_eval('pfftAvE = pfftE?; pfftAvE.p{1} = sum(pfftAvE.p{1});',ic)
c_eval('pfftAvB = pfftB?; pfftAvB.p{1} = sum(pfftAvB.p{1});',ic)

nRows = 3;
nCols = 1;

isub = 1;
hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,pfftAvE.f,pfftAvE.p{1},...
         pfftAvB.f,pfftAvB.p{1})
hca.XScale = 'log'; hca.YScale = 'log';

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,pfftAvE.f,pfftAvE.p{1},...
         pfftAvB.f,pfftAvB.p{1})
hca.XScale = 'log'; hca.YScale = 'log';

