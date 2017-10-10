% dispersion_relation_compare.m
% makes dispersion relation using phi-matching

tint = toepoch([2008 04 22 17 37 12.0;2008 04 22 17 37 13.0])'; % em waves too in high beta
sc = 1;
csys = 'gsm';
flim = 0.1;
tool.single_event
dispersion_relation_yoon2008_run
tool.makedisprel
pfftE = irf_powerfft(E,size(E,1),450,0.5);
pfftB = irf_powerfft(B,size(B,1),450,0.5);
%%
hca = subplot(1,4,1);
plot(hca,save_k(save_c>clim)*re,f_center(save_c>clim),'*',save_k*re,f_center,...
     ky/sqrt(2),x_real_store2/omega_lh,ky/sqrt(2),x_imag_store2/omega_lh)
hca.YLabel.String = '\omega/\omega_{LH}';
hca.XLabel.String = 'k\rho_e';

hca = subplot(1,4,2); 
EBscale = B0.^2/(n*1e6*e*4*pi*1e-7);
semilogx(hca,pfftE.p{1},pfftE.f/(omega_lh/2/pi),pfftB.p{1},pfftB.f/(omega_lh/2/pi))
hca.YLabel.String = '\omega/\omega_{LH}';
hca.XLabel.String = 'Power';
hca.XLim = [1e-9 1e2];

hca = subplot(1,4,3); 
semilogx(hca,pfftE.p{1},pfftE.f/(omega_lh/2/pi),pfftE.p{3},pfftE.f/(omega_lh/2/pi),pfftE.p{3},pfftE.f/(omega_lh/2/pi))
semilogx(hca,pfftB.p{1},pfftB.f/(omega_lh/2/pi),pfftB.p{3},pfftB.f/(omega_lh/2/pi),pfftB.p{3},pfftB.f/(omega_lh/2/pi))
semilogx(hca,sqrt(pfftE.p{1}./pfftB.p{1}),pfftB.f/(omega_lh/2/pi),sqrt(pfftE.p{2}./pfftB.p{2}),pfftB.f/(omega_lh/2/pi),sqrt(pfftE.p{3}./pfftB.p{3}),pfftB.f/(omega_lh/2/pi))
hca.YLabel.String = '\omega/\omega_{LH}';
hca.XLabel.String = 'E/B [10^3 km/s]';
hca.XLim = [1e-1 1e3];

hca = subplot(1,4,4); 
semilogy(hca,pfftE.f/(omega_lh/2/pi),pfftE.p{1},pfftB.f/(omega_lh/2/pi),pfftB.p{1})
tisemilogy(hca,pfftE.f/(omega_lh/2/pi),smooth(pfftE.p{1}),pfftB.f/(omega_lh/2/pi),smooth(pfftB.p{1}))
hca.XLabel.String = '\omega/\omega_{LH}';
hca.YLabel.String = 'Power';
hca.YLim = [1e-9 1e2];



