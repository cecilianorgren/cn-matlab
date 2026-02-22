%
% Define magnetic field
a = 10*1e3;
b = 1*1e3;
xvec = 0;
zvec = b*linspace(-10,10,100);
yvec = 0;

[X,Y,Z] = ndgrid(xvec,yvec,zvec);
%dx = x(2) - x(1);
%dy = y(2) - y(1);
%dz = z(2) - z(1);
%x_xline = x;
%y_xline = x*b/a;

%Ay = @(x,y,z) (x/a).^2 - (z/b).^2;
%AY0 = Ay(X,Y,Z);

%Ay = @(x,y,z,t) B0*(b)*log(sech(z/b)) + B0*(b/a)*x; % Bz = dAy/dx, Bx = -dAy/dx
%AY0 = Ay(X,Y,Z);

% Integrate orbits
units = irf_units;
m = units.mp;
q = units.e;
B0 = 10e-9; % T
E0 = 1e-3; % V/m
lz = 100e3; % m
Bx = @(x,y,z) B0*(z/lz);
By = @(x,y,z) x*0;
Bz = @(x,y,z) 0;
Ex = @(x,y,z) x*0;
Ey = @(x,y,z) x*0 + z*0; % V/m
Ez = @(x,y,z) x*0;
Ez = @(x,y,z) 0;

options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,xvec([1 end])),...
                 'AbsTol',1e-6);
options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,2,[-400*1e3 400*1e4]),...
                 'AbsTol',1e-12);
%options = odeset();
EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 20;

nP = 1;

vz0 = E0/B0; % m/s
Te0 = 20; % eV
Te = Te0 + 0.3*Te0*randn([nP 1]); % eV
Te(Te<0) = 0;
vt = sqrt(2*Te*units.eV/m); % m/s
wce = units.e*B0/m;
rhoe = vt/wce;
z0 = 0*lz;
ph1 = rand(nP,1)*360;
ph2 = rand(nP,1)*180;
vx = vt.*cosd(ph1).*sind(ph2); vx = 0;
vy = vt.*sind(ph1).*sind(ph2);
vz = vt.*cosd(ph2);


x_init_all = [repmat([0 0 z0],nP,1) vy*0 vy vz];


clear p
xlim = [0 0];
ylim = [0 0];
zlim = [0 0];
for ip = 1:size(x_init_all,1)
  x_init = x_init_all(ip,:);
  [t,x_sol] = ode45(EoM,[tstart tstop],x_init,options); % 
  p(ip).t = t;
  p(ip).x = x_sol(:,1);
  p(ip).y = x_sol(:,2);
  p(ip).z = x_sol(:,3);
  p(ip).vx = x_sol(:,4);
  p(ip).vy = x_sol(:,5);
  p(ip).vz = x_sol(:,6);
  p(ip).v = vt(ip);
  xlim(1) = min([xlim(1); p(ip).x]);
  xlim(2) = max([xlim(2); p(ip).x]);
  ylim(1) = min([ylim(1); p(ip).y]);
  ylim(2) = max([ylim(2); p(ip).y]);
  zlim(1) = min([zlim(1); p(ip).z]);
  zlim(2) = max([zlim(2); p(ip).z]);

  
end
nP = numel(p);

wb = sqrt(units.e*B0*p.v/m/lz);
Tb = 2*pi/wb;
%%


nRows = 1;
nCols = 1;
h = setup_subplots(nRows,nCols);
isub = 1;

iPs = 1:nP; 
if 1 % (x,z)
  hca = h(isub); isub = isub + 1;
  
  holdon = 0;
  
  for ip = iPs
    plot(hca,p(ip).t,p(ip).z*1e-3,'LineWidth',2)
    if not(hold(hca,'on')), hold(hca,'on'); holdon = 1; end
  end
  hold(hca,'off')
end
    









