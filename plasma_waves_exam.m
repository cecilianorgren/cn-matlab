% plasma_waves_exam
%% alfven resonator
units = irf_units;
B = @(B0,r) B0./r.^2;
nO = @(z,h,nO0) nO0*exp(-z/h);
nH =@(r,p,nH0) nH0*r.^-p;


B0 = 60; % 10-6 T (microT)
nO0 = 1e5; % cc
p = 1;
h = 600; % km
nH0 = 20; % cc

RE = 6400; % km

z = linspace(RE*0,4*RE,200);
r = (z+RE)/RE;

znorm = find(z<500);

n_atm = [linspace(0,1,numel(znorm)).^2];

nhyd = nH(r,p,nH0);
nox = nO(z,h,nO0).*(1.0005-exp(-z/(10*h)));

ntot = nhyd + nox;

Bfield = B(B0,r);
mu0 = 4*pi*1e-7;
eps0 = 8.85e-12;
mp = 1.67e-17;
me = 9.1e-31;
mo = 16*mp;
e = 1.6e-19;
vA = Bfield*1e-6./sqrt(mu0*(nhyd*mp+nox*mo)*1e6); % km/s

opeoce = sqrt(ntot*1e6*me/eps0./(Bfield*1e-6).^2);
overoce = me./e./(Bfield*1e-6);
overope = 1./sqrt(e^2*ntot*1e6/me/eps0); 
ope = sqrt(e^2*ntot*1e6/me/eps0); 

nplot = 7;
subplot(nplot,1,1)
semilogy(z/RE,ntot)
title(gca,'Altitude profiles of the ionosphere')
ylabel('n_e [cc]')
xlabel('r/RE')

subplot(nplot,1,2)
plot(z/RE,Bfield)
ylabel('|B| [10^{-6} T]')
xlabel('r/RE')

subplot(nplot,1,3)
plot(z/RE,vA)
ylabel('v_A [km/s]')
xlabel('r/RE')

subplot(nplot,1,4)
plot(z/RE,opeoce)
ylabel('ope/oce ')
xlabel('r/RE')

omega = 1e6; % smaller than peak plasma frequency which is ~ 3e6
subplot(nplot,1,5)
plot(z/RE,1./overope,z/RE,ones(1,numel(z))*omega)
text(2,1.5*omega,'omega')
ylabel('ope ')
xlabel('r/RE')



subplot(nplot,1,6)
plot(z/RE,omega.*overope,z/RE,ones(1,numel(z)))
ylabel('omega/ope ')
xlabel('r/RE')
subplot(nplot,1,7)
plot(z/RE,omega.*overoce)
ylabel('omega/oce ')
xlabel('r/RE')

%% Dispersion surfaces for lower hybrid waves
olh = 1;
ol = 1;
omega = @(olh,kpar,kper,ol) olh + ol*kpar./kper;

kpar = logspace(-2,2,100);
kper = logspace(-1,2,100);

[KPAR,KPER] = meshgrid(kpar,kper);

surf(log(KPAR),log(KPER),omega(olh,KPAR,KPER,ol)); shading flat;
xlabel('kpar')
ylabel('kper')

%%
K = sqrt(KPAR.^2+KPER.^2);
TH = atan(KPER./KPAR);
omega = @(k,th) k.^2*(1+tan(th).^2./(tan(th).^3));

contour(KPAR,KPER,omega(K,TH))
xlabel('kpar')
ylabel('kper')

%%
kz = @(o,c,e,n0,eps,m,z,h) sqrt(o^2/c^2-e^2*n0*exp(-z/h)/2/eps/m/(c^2));

o = 10^5;
c =3e8; % m/s
e = 1.6e-19;
nO0 = 1e5; % cc
nO0 = nO0*1e6; % m^-3
eps0 = 8.85e-12;
me = 9.1e-31;
h = 600; % km

RE = 6400; % km
z = linspace(RE*0,4*RE,200);

nhyd = nH(r,p,nH0);
nox = nO(z,h,nO0).*(1.0005-exp(-z/(10*h)));
ntot = nhyd + nox;

KZ=sqrt(o^2-e^2*ntot/eps0/me)/c;


KZ = kz(o,c,e,nO0,eps0,me,z,h);

subplot(2,1,1); plot(z/RE,KZ);
title(gca,'k(z)')
xlabel(gca,'z/RE')

for ii = 10:numel(z)
    intKZ(ii) = trapz(z(1:ii),KZ(1:ii));
end

subplot(2,1,2); plot(z/RE,intKZ)
title(gca,'int(k(z)dz )')
xlabel(gca,'z/RE')

%%
l = 2*pi;
m = 1;
x = linspace(0,l,100);
v = linspace(-3,3,100);

[X,V] = meshgrid(x,v);

x1 = sin(2*pi*X*m/l);
v1 = cos(2*pi*X*m/l);

subplot(2,1,1); surf(X,V,x1); view([0 0 1]); shading flat;
subplot(2,1,2); surf(X,V,v1); view([0 0 1]); shading flat;









