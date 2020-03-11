% The mangetic field pressure in the lobes should balance the plasma
% thermal pressure in the centre of the current sheet.

units = irf_units;
B0 = 20e-9; % T
PB = B0/2/units.mu0; % Pa

TiTe = 5;
fun_P = @(n,T) units.eV*T*(1+1/TiTe)*n;
fun_n = @(T,B) B.^2/2/units.mu0;

nTi = 100;
minTi = 1000;
maxTi = 10000;
vecTi = linspace(minTi,maxTi,nTi);

nN = 100;
minN= 0.1*1e6;
maxN = 1*1e6;
vecN = linspace(minN,maxN,nN);

[N,T] = meshgrid(vecN,vecTi);


h = setup_subplots(2,1);
isub = 1;

hca = h(isub); isub = isub + 1;
P = fun_P(N,T);
[c,hclab] = contour(hca,N*1e-6,T,P'*1e9);
clabel(c,hclab)
hold(hca,'on')
%[c,hclab] = contour(hca,N*1e-6,T,[1 1]*PB'*1e9);
clabel(c,hclab)
hold(hca,'off')
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'P (nPa)';
hca.XLabel.String = 'n (cm^{-3})';
hca.YLabel.String = 'T_i (eV)';


hca = h(isub); isub = isub + 1;
P = fun_P(N,T);
[c,hclab] = contour(hca,N*1e-6,T,P'/PB);
clabel(c,hclab)
hold(hca,'on')
%[c,hclab] = contour(hca,N*1e-6,T,[1 1]*PB'*1e9);
clabel(c,hclab)
hold(hca,'off')
hcb = colorbar('peer',hca);
hcb.YLabel.String = sprintf('P/P_B(B_0 = %g nT)',B0*1e9);
hca.XLabel.String = 'n (cm^{-3})';
hca.YLabel.String = 'T_i (eV)';