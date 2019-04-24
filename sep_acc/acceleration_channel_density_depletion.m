units = irf_units;
e = units.e;
me = units.me;

phi = 1000;
Te = 100;
n = 0.05*1e6; % m^-3
vd = -2*vt(Te);

%U = units.me*(V-vph).^2/2 - units.e*double(phi(X));
v(phi) = 2*units.e*phi;
vt = @(T) sqrt(2*e*T/me); % m/s
vacc = @(phi) sqrt(2*e*phi/me); % m/s
%f = @(v,n,T,vd,phi) n/pi^0.5./sqrt(2*e*T/me).*exp(-(v-vd)^2/(2*e*T/me));

f = @(v,n,T,vd,phi) n./pi^0.5./vt(T).*exp(-sqrt((v-vd).^2./(vt(T).^2) - vacc(phi).^2/(vt(T).^2)));


nv = 1000;
vvec = linspace(-2*vt(1000),2*vt(1000),nv);
plot(vvec*1e-6,f(vvec,n,Te,vd,0),vvec*1e-6,0*2e4*f(vvec,n,Te,vd,1000))

