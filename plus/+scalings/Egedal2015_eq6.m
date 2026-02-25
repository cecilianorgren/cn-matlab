
b1b2 = linspace(1,6,110);
va = linspace(1000,4000,100);

[B1B2,VA] = meshgrid(b1b2,va);

UEPAR = VA.*B1B2.*(B1B2-1);

figure(74)
hca = subplot(1,1,1);
%imagesc(hca,b1b2,va,UEPAR'*1e-3); 
contourf(hca,B1B2,VA,UEPAR*1e-3,0:5:50);
hca.YLabel.String = 'v_A (km/s)';
hca.XLabel.String = 'B_1/B_2';

hcb = colorbar('peer',hca);
hcb.YLabel.String = 'u_{e,||} (10^3 km/s)';