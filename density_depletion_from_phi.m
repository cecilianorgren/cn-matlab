lpar = 1;
lperp = linspace(0.5,3,1000);
dn3D1D = @(lpar,lperp) 2*(lpar./lperp).^2+1;

hca = subplot(1,1,1);
plot(hca,lperp/lpar,dn3D1D(lpar,lperp),lperp/lpar,lperp./lperp,'k--')
hca.XLabel.String = 'l_{||}/l_{\perp}';
hca.XLabel.String = 'cigar <-    l_{\perp}/l_{||}     -> pancake';
hca.YLabel.String = '(n_e-n_i)_{3D}/(n_e-n_i)_{1D}';