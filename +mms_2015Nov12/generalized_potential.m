m = 1;
phi0 = 0;
q = -1;
L = 1;
B0 = 10;
Bg = 0;

Az = @(x) B0*L*log(cosh(x/L));
Ay = @(x) -Bg*x;

Gamma = @(py,pz,phi0,x,L) 0.5/m*((py-q*B0*L*log(cosh(x/L))).^2+(pz+q*Bg*x).^2) + q*phi0;
 

x = linspace(-2*L,2*L,50);

py = 0.1;
pz = 0.1;

plot(x,Gamma(py,pz,phi0,x,L))


%% Physical units
units = irf_units;
m = units.me;
phi0 = 0;
q = -units.e;
L = 10; % km

T0 = 50; % eV
v0 = sqrt(2*units.eV*T0/units.me); % m/s
E0 = 0.5*m*v0.^2;

B0 = 10e-9; % T
Bg = 5e-9;

Ay = @(x) B0*L*log(cosh(x/L));
Az = @(x) -Bg*x;

theta = 90;
pz = 0;m*v0+q*Az(0.5*L);
py = -1e-26;
pz = 0.5*py;

Gamma = @(py,pz,phi0,x,L) 0.5/m*((py-q*B0*L*log(cosh(x/L))).^2+(pz+q*Bg*x/L).^2) + q*phi0;
Gamma0 = Gamma(py,pz,phi0,0,L);

x = linspace(-2*L,2*L,50);

nRows = 3; nCols = 1; isub = 1;

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,x,q*Ay(x))

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,x,q*Az(x))

hca = subplot(nRows,nCols,isub); isub = isub + 1;
plot(hca,x,Gamma(py,pz,phi0,x,L))

%% B = B0*x
m = 1;
phi0 = 0;
q = -1;
L = 10;
B0 = 10;
Bg = 5;

Ay = @(x) 0.5*B0*x.^2;
Az = @(x) -Bg*x;

Gamma = @(py,pz,phi0,x,L) 0.5/m*((py-q*Ay(x)).^2+(pz+q*Az(x)).^2) + q*phi0;

x = linspace(-2*L,2*L,50);

py = 10;
pz = 10;

plot(x,Gamma(py,pz,phi0,x,L))
%plot(x,Gamma(py,pz,phi0,x,L))

%%
x = -8:0.1:8;
fx = @(x,a) (a-x.^2).^2;
plot(x,fx(x,9))
