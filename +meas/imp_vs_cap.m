% imp. vs cap.
Cap = @(eps0,L,a) 2*pi*eps0/(log(L/2/a)-1);
Ind = @(mu0,L,a) mu0*L*log(L/a)/2/pi;

units= irf_units;
a = 0.01; % m
L = 100; % m
mu0 = units.mu0;
eps0 = units.eps0;

om = logspace(6,8,100);

XL = @(om) om*Ind(mu0,L,a);
XC = @(om) 1./om./Cap(eps0,L,a);

plot(om*1e-6,XL(om),om*1e-6,XC(om))
legend('XL','XC')
ylabel('X')
xlabel('frequency [MHz]')


%% finite wabelength effect

x = linspace(0.1,100,1000);
L1 = 850; k1 = 2*pi/L1;
L2 = 102; k2 = 2*pi/L2;
L3 = 80; k3 = 2*pi/L3;
wave = @(k,x) cos(k*x);

plot(x,wave(k1,x),...
     x,wave(k2,x),...
     x,wave(k3,x))
set(gca,'xlim',[x(1),x(end)])

title('finite wavelength effects')
xlabel('length along antenna')
ylabel('normalized electric field')
legend('L=850 m','L=102 m','L=80 m')

%%
x = linspace(1,100,1000);
L = 60; k =2*pi/L;
E = zeros(1,1000);
for ii = 1:1000
    E(ii) = (wave(k,x(end))-wave(k,x(ii)))/x(ii);
end

plot(x,E)