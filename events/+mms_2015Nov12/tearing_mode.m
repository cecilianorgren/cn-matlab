units = irf_units;
Te = 50;
Ti = 500;
BL0 = 10; % nT
BM0 = 5; % nT
L = 15; % km
% V_Te =   4.19 Mm/s (40 eV)
vte = sqrt(2*units.eV*Te/units.me);
vti = sqrt(2*units.eV*Ti/units.mp);
wce = units.e*BL0*1e-9/units.me;
wci = units.e*BL0*1e-9/units.mp;
rhoe = vte/wce*1e-3; % km
rhoi = vti/wci*1e-3; % km
  
BMstar = @(L) BL0*1e-9*(1/sqrt(2))*(rhoi/L)^0.5*(Te*units.me/Ti/units.mp)^0.25*1e9; % nT
