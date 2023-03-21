wE = 2*pi*24*60*60;
RE = 6371200; % m
me = 9.1094e-31;
mp = 1.6726e-27;
e = 1.6022e-19;
kB = 1.3806e-23;
mu0 = 1.2566e-06;
M = 8.22e22; % A m^2

Er = @(r) -wE*mu0*M/4/pi./r.^2;


r = RE*linspace(1,7,10);
lon = linspace(0,2*pi,100);

[R,LON] = ndgrid(r,lon);

X = R.*cos(LON);
Y = R.*sin(LON);

ER = Er(R);

EX = ER.*cos(LON);
EY = ER.*sin(LON);


E_conv = 1e7; % V/m
EX_conv = 0;
EY_conv = E_conv;

hca = subplot(1,1,1);
scale = NaN;
quiver(hca,X,Y,EX,EY)
hold(hca,'on')
quiver(hca,X,Y,EX*0+EX_conv,EY*0+EY_conv)
quiver(hca,X,Y,EX+EX*0+EX_conv,EY+EY*0+EY_conv)

surf(X,Y,X*0-1,EY+EY_conv)
shading(hca,'flat')
colormap(pic_colors('blue_red'))
hca.CLim = max(abs(hca.CLim))*[-1 1];
hcb = colorbar(hca);

hold(hca,'off')



%%
% transformation matrix between cartesian and spherical
T = [sin(colat).*cos(lon), sin(colat).*sin(lon), cos(colat);...
     cos(colat).*cos(lon), cos(colat).*sin(lon), -sin(colat);...
     -sin(lon)         , cos(lon)        ,  0];
  
T_sp2cart = T';

Ecart = T_sp2cart*Bsphere;
  
bx = Bcart(1,:);
by = Bcart(2,:);
bz = Bcart(3,:);
