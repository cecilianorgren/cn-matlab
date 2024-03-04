% mms_book_guide_field

units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

Er = 0.5e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 0*5*1e-9; % guide field, T
Bn = 1.0*1e-9; % normal field, T
phi0 = 1000*units.e;
Lphi = L;

BH = 10e-9;
LxBn = 0.01*L;

% General current sheet structure:
% Harris current sheet, incl. parameters in functions to make more
% versatile later on
fBx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
fBy = @(x,y,z) x.*z.*(-BH)/LxBn/LxBn.*exp(-x.^2/LxBn.^2 - z.^2/LxBn.^2*10) + y*0 + z*0 + Bg;
fBz = @(x,y,z) x*0 + y*0 + z*0 + -Bn*x/LxBn;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z,Er) x*0 + y*0 + z*0 + Er;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ez = @(x,y,z,phi0,Lphi) x*0 + y*0 + z*0 -phi0/Lphi;

x0 = linspace(-10,10,5)*1e3;
y0 = linspace(0,0,1)*1e3;
z0 = linspace(-5,5,4)*1e3;

x0 = linspace(-10,-10,2)*1e3;
y0 = linspace(0,0,1)*1e3;
z0 = linspace(-5,5,10)*1e3;

x0 = linspace(-30,30,3)*1e3;
y0 = linspace(0,0,1)*1e3;
z0 = linspace(-5,5,10)*1e3;

[X0,Y0,Z0] = ndgrid(x0,y0,z0);

r0 = [X0(:),Y0(:),Z0(:)];
dr = [0.01, 0.01, 0.01]*1e3;
dr = dr*1e0;
%rmax = [];

xvec = linspace(-10,10,100)*1e3;
zvec = linspace(-5,5,70)*1e3;
[X,Z] = ndgrid(xvec,zvec);
BX = fBx(X,0,Z);
BY = fBy(X,0,Z);
BZ = fBz(X,0,Z);


blines = struct([]);
for ib = 1:numel(X0)
  
  x = r0(ib,1);
  y = r0(ib,2);
  z = r0(ib,3);
  
  while numel(x) < 3000 %|| x(end) > xvec(1) || x(end) < xvec(end)
    Bx = fBx(x(end),y(end),z(end));
    By = fBy(x(end),y(end),z(end));
    Bz = fBz(x(end),y(end),z(end));
    Babs = sqrt(Bx.^2 + By.^2 + Bz.^2);
    bx = Bx/Babs;
    by = By/Babs;
    bz = Bz/Babs;
    x(end+1) = x(end) + bx*dr(1);
    y(end+1) = y(end) + by*dr(2);
    z(end+1) = z(end) + bz*dr(3);
  end
  blines(ib).x = x;
  blines(ib).y = y;
  blines(ib).z = z;
end

h = setup_subplots(3,2);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,xvec,zvec,BX')
shading(hca,'flat')
hb = colorbar(hca);

hca = h(isub); isub = isub + 1;
pcolor(hca,xvec,zvec,BY')
shading(hca,'flat')
hb = colorbar(hca);

hca = h(isub); isub = isub + 1;
pcolor(hca,xvec,zvec,BZ')
shading(hca,'flat')
hb = colorbar(hca);

hca = h(isub); isub = isub + 1;
plot(hca,0,0)
hold(hca,'on')
for ib = 1:numel(blines)  
  plot(hca,blines(ib).x(1),blines(ib).z(1),'go')
  plot(hca,blines(ib).x(end),blines(ib).z(end),'rx')
  plot(hca,blines(ib).x,blines(ib).z)
end
hold(hca,'off')

hca = h(isub); isub = isub + 1;
plot3(hca,0,0,0)
hold(hca,'on')
for ib = 1:numel(blines)
  plot3(hca,blines(ib).x(1),blines(ib).y(1),blines(ib).z(1),'go')
  plot3(hca,blines(ib).x(end),blines(ib).y(end),blines(ib).z(end),'rx')
  plot3(hca,blines(ib).x,blines(ib).y,blines(ib).z)
end
hold(hca,'off')

hca = h(isub); isub = isub + 1;
plot(hca,0,0)
hold(hca,'on')
for ib = 1:numel(blines)
  plot(hca,fBy(blines(ib).x,blines(ib).y,blines(ib).z));
end
hold(hca,'off')





























