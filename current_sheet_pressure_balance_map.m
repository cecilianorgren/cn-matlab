% The mangetic field pressure in the lobes should balance the plasma
% thermal pressure in the centre of the current sheet.

units = irf_units;
B0 = 20e-9; % T
PB = B0/2/units.mu0;

TiTe = 5;
fun_P = @(n,T) units.eV*T*(1+1/TiTe)*n;
fun_n = @(T,B) B.^2/2/units.mu0;

nTi = 100;
minTi = 1000;
maxTi = 10000;
vecTi = linspace(minTi,maxTi,nTi);

nN = 100;
minN= 0.1*1e-6;
maxN = 1*1e-6;
vecN = linspace(minN,maxN,nN);

[N,T] = meshgrid(vecN,vecTi);

hca = subplot(1,1,1);
P = fun_P(N,T);
contour(hca,N,T,P');
colorbar('peer',hca)