units = irf_units;
a = 2*pi/100e3; % 1/m
phi0 = 100; % V

dne = @(x,a) phi0*units.eps0*a^2/units.e*sin(x*a);
vph = 450e3; % m/s
x = (ne1.time-ne1.time(1))*vph;

tsdne = irf.ts_scalar(ne1.time,dne(x,a)*1e-6); tsdne.name = 'n(\phi_0,a,x) (cc)';
tsdphi = ne1*1e6*units.e/units.eps0/a^2*(-1);

h = irf_plot({gseB1,gseE1,ne1,tsdne,tsdphi});
irf_zoom(h,'x',irf.tint('2017-07-11T22:33:19.00Z/2017-07-11T22:33:21.00Z')); 
irf_zoom(h,'y');
