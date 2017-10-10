% See bounce times for the electron hole
units = irf_units;
E0 = 20e-3; % 20 mV/m
phi0 = 200; % V
lpar = 5e3;
lam = 2*lpar;
k = 2*pi/lam;

obE = sqrt(units.e*E0*k/units.me);
tbE = 2*pi/obE; % s
tbE_ms = tbE*1e3 % ms

obP = sqrt(2*units.e*phi0*k^2/units.me);
tbP = 2*pi/obP; % s
tbP_ms= tbP*1e3 % ms