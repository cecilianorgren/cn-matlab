% What is the thermal flux of a given plasma population. What is the total
% thermal flux if differnent populations intermix parallel to B?
units = irf_units;
m = units.me;
T = 50; % eV

f_vt = @(T,m) sqrt(2*units.eV*T/m);

%f_3D = @(v,T,n,vd) n/((pi)^(3/2)*w(T).^3)*exp(-(v-vd).^2./w(T)./w(T)); 
f_max1D = @(v,T,n,vd,m) n/((pi)^(1/2)*vt(T,m))*exp(-(v-vd).^2./vt(T,m)./vt(T,m));


T1 = 20;
m1 = units.me;
n1 = 25;

f1 = 




