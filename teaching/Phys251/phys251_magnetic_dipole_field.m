units = irf_units;
mu0 = units.mu0;
M = 8.22e22; % A m^2
RE = units.RE*1e-3;

Br = @(r,theta) -mu0*M./(4*pi*r.^2)*2.*cosd(theta);
Bt = @(r,theta) -mu0*M./(4*pi*r.^2).*sind(theta);

r = RE*linspace(1,10,100);
theta = 0:180;

[R,TH] = ndgrid(r,theta);

Z = R.*cosd(TH);
X = R.*sind(TH);


%%
f_r2 = @(r1,th1,th2) r1.*sind(th2).^2./(sind(th1).^2);

h = 0; % km
th1 = 90;
th2 = [10:1:170];
r1 = 5*RE; % km

r2 = f_r2(r1,th1,th2);

z = r2.*cosd(th2);
x = r2.*sind(th2);


hca = subplot(1,1,1);

plot(hca,x/RE,z/RE,'-')


hold(hca,'on')
plot(hca,RE*cosd(0:360)/RE,RE*sind(0:360)/RE,'k')
hold(hca,'off')

axis(hca,'equal')
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'x (km)';
hca.YLabel.String = 'z (km)';





