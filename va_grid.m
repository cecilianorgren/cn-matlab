units = irf_units;
b = linspace(2,30,30);
n = logspace(-2,2,20);
[B,N] = meshgrid(b,n);
va = @(b,n) b*1e-9./sqrt(units.mu0*units.mp*n*1e6)*1e-3;
VA = va(B,N);

hca = subplot(1,1,1);
pcolor(hca,N,B,log10(VA))
hca.XScale = 'log';
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'log_{10}v_A (km/s)';
hca.XScale = 'log';
hca.XLabel.String = 'n (cm^{-3})';
hca.YLabel.String = 'B (nT)';


