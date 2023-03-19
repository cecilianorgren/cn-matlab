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

h = 150;
th1 = 30;
th2 = 30:150;
r1 = RE + h; % km

r2 = f_r2(r1,th1,th2);

z = r2.*cosd(th2);
x = r2.*sind(th2);


hca = subplot(1,1,1);

plot(hca,x,z)


hold(hca,'on')
plot(hca,RE*cosd(0:360),RE*sind(0:360),'k')
hold(hca,'off')

axis(hca,'equal')
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.XLabel.String = 'x (km)';
hca.YLabel.String = 'z (km)';





