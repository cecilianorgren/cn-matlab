ic = 1;
fftTint = irf.tint('2015-11-30T00:22:00.00Z/2015-11-30T00:23:00.00Z'); 

c_eval('flh_av = mean(flh?.data);',ic)
units = irf_units;

flim = 1;

c_eval('ne = ne?',ic)
c_eval('gseB0 = gseB?.filt(0,flim,[],5);',ic)
c_eval('[parE,perpE] = irf_dec_parperp(gseB?,gseE?);  parE.name = ''E par''; perpE.name = ''E perp'';',ic)
c_eval('[parBscm,perpBscm] = irf_dec_parperp(gseB?,gseB?scm); parBscm.name = ''B scm par''; perpB?scm.name = ''B scm perp'';',ic)
c_eval('facE = irf_convert_fac(gseE?,gseB?,[1 0 0]); facE.name = ''fac E'';',ic)
c_eval('facBscm = irf_convert_fac(gseB?scm,gseB?,[1 0 0]); facB.name = ''fac B'';',ic)


phiB = gseB0.resample(parBscm).abs*parBscm.filt(flim,0,[],5)*1e-18/(ne*1e6*units.e*units.mu0);
%c_eval('phiB?perp = gseB0?.resample(gseB?scmperp).abs*gseB?scmperp.filt(flim,0,[],5)*1e-18/(ne?.filt(0,1,[],5)*1e6*units.e*units.mu0);',ic)
%c_eval('intE? = irf_integrate(gseE?perp.abs.filt(flim,0,[],5));',ic)
%%
fsE = 1/(perpE.time(2)-perpE.time(1));
fsB = 1/(perpBscm.time(2)-perpBscm.time(1));

nfft = floor(log2(facE.length));
nfft = 2^ceil(log2((fftTint(2)-fftTint(1))*fsE));
nfftTint = fftTint(1)+[0 nfft*1.1/fsE];
indT0 = find(gseE1perp.time>fftTint(1),1,'first');

% Make spectra
pfftE = irf_powerfft(facE,nfft,fsE,0);
pfftBscm = irf_powerfft(facBscm,nfft,fsE,0);
pfftPhiB = irf_powerfft(phiB,nfft,fsE,0);

%% Test dispersion relation
fnoise = 50;
v = 100;
%pfftEperp = pfftE; pfftEperp.p = {mean(pfftEperp.p{1}+pfftEperp.p{2},1)}; pfftEperp.t = mean(pfftEperp.t,1);
pfftEperp = pfftE; pfftEperp.p = {mean(pfftEperp.p{1},1)}; pfftEperp.t = mean(pfftEperp.t,1);
pfftBpar = pfftBscm; pfftBpar.p = {mean(pfftBpar.p{3},1)}; pfftBpar.t = mean(pfftBpar.t,1);
pfftPhiB = pfftPhiB; pfftPhiB.p = {mean(pfftPhiB.p{1},1)}; pfftPhiB.t = mean(pfftPhiB.t,1);

nRows = 2;
nCols = 2;

isub = 1;
hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,pfftEperp.f,pfftEperp.p{1},...
         pfftBpar.f,pfftBpar.p{1},...
         pfftEperp.f,pfftEperp.p{1}./(pfftEperp.f.^2*2*pi)')
hca.XScale = 'log'; hca.YScale = 'log';
legend(hca,'E_{perp} (mV/m)','B_{||} (nT)','E_{perp}v/f (V)')
hca.YLabel.String = 'Power   [(   )^2/Hz]';
hca.XLabel.String = 'f (Hz)';

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,pfftEperp.f,v^2*pfftEperp.p{1}./(pfftEperp.f.^2*2*pi)',...
         pfftPhiB.f,pfftPhiB.p{1})
hca.XScale = 'log'; hca.YScale = 'log';
legend(hca,'E_{perp}v/f','\phi_{B_{||}}=(B_0/ne\mu_0)\delta B')
hca.YLabel.String = 'Power [V^2/Hz]';
hca.XLabel.String = 'f (Hz)';
hca.Title.String = ['v = ' num2str(v) ' km/s'];

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,sqrt(pfftEperp.p{1})./sqrt(pfftPhiB.p{1}),pfftEperp.f,'.')
hca.XScale = 'log'; hca.YScale = 'log';
hca.YLabel.String = 'f (Hz)';
hca.XLabel.String = 'k (km^{-1})';
hold(hca,'on')
%plot(hca,logspace(-2,1,10),logspace(-2,1,10)*v)
hold(hca,'off')

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,sqrt(pfftEperp.p{1})./sqrt(pfftPhiB.p{1}),pfftEperp.f,'.')
%hca.XScale = 'log'; hca.YScale = 'log';
hca.YLabel.String = 'f (Hz)';
hca.XLabel.String = 'k (km^{-1})';
hca.YLim = [0 80];
hold(hca,'on')
%hv = plot(hca,logspace(-2,1,10),logspace(-2,1,10)*v);
hold(hca,'off')


hca = subplot(nRows,nCols,1);
hca.Title.String = ['Time: ' tintUTC(1,12:22) ' - ' tintUTC(2,12:22)];

for ii = 1:2
  hca = subplot(nRows,nCols,ii);
  hold(hca,'on')
  hf = plot(hca,flh_av*[1 1],hca.YLim,...
                flim*[1 1],hca.YLim);%,...
           %fnoise*[1 1],hca.YLim)  
  %irf_legend(hca,{'f_{LH}','f_{B0}'},[0.1 0.1])
  hold(hca,'off')
end

for ii = 3:4
  hca = subplot(nRows,nCols,ii);
  hold(hca,'on')
  hf = plot(hca,hca.XLim,flh_av*[1 1]);%,...
                %hca.XLim,fnoise*[1 1]);
  hv = plot(hca,logspace(-2,1,10),logspace(-2,1,10)*v);            
  hold(hca,'off')
  legend([hf(1) hv],{'f_{LH}',['v = ' num2str(v) ' km/s']})
end


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


