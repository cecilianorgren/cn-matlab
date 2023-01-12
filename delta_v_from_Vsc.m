units = irf_units;

Vsc = 20;
Ue = logspace(log10(Vsc),3,100);


v = @(U) sqrt(U*units.eV*2/units.me); % m/s

hca = subplot(2,1,1);
plot(hca,Ue,v(Ue)*1e-3,Ue,v(Ue-Vsc)*1e-3)

hca = subplot(2,1,2);
plot(hca,Ue,(v(Ue)-v(Ue-Vsc))*1e-3)
