function [bx,by,bz] = phys251_fun_magnetic_dipole_field(x,y,z)
units = irf_units;
mu0 = units.mu0;
M = 8.22e22; % A m^2
RE = units.RE*1e-3;

pref = mu0*M/(4*pi);

Br = @(r,colat) -pref/(r.^2)*2.*cos(colat);
Bcolat = @(r,colat) -pref/(r.^2).*sin(colat);
Blon = @(r,colat) r*0;


% Get particle polar angle and radius
r = sqrt(x.^2 + y.^2 + z.^2);
colat = atan(sqrt(x.^2 + y.^2)./z);
lon = atan2(y,x);


Bsphere = [Br(r,colat); Bcolat(r,colat); Blon(r,colat)];

if 1
  % transformation matrix between cartesian and spherical
  T = [sin(colat).*cos(lon), sin(colat).*sin(lon), cos(colat);...
       cos(colat).*cos(lon), cos(colat).*sin(lon), -sin(colat);...
       -sin(lon)         , cos(lon)        ,  0];
  
  T_sp2cart = T';

  Bcart = T_sp2cart*Bsphere;
  
  bx = Bcart(1,:);
  by = Bcart(2,:);
  bz = Bcart(3,:);

else
  T(1,1,:) = sin(colat).*cos(lon);
  T(2,1,:) = sin(colat).*sin(lon);
  T(3,1,:) = cos(colat);
  T(1,2,:) = cos(colat).*cos(lon);
  T(2,2,:) = cos(colat).*sin(lon);
  T(3,2,:) = -sin(colat);
  T(1,3,:) = -sin(lon);
  T(2,3,:) = cos(lon);
  T(3,3,:) = 0;
  
  
  T_sp2cart = permute(T,[2,1,3]);  
  
  bx = squeeze(T_sp2cart(1,1,:)).*Bsphere(1,:)' + squeeze(T_sp2cart(2,1,:)).*Bsphere(2,:)' + squeeze(T_sp2cart(3,1,:)).*Bsphere(3,:)';
  by = squeeze(T_sp2cart(1,2,:)).*Bsphere(1,:)' + squeeze(T_sp2cart(2,2,:)).*Bsphere(2,:)' + squeeze(T_sp2cart(3,2,:)).*Bsphere(3,:)';
  bz = squeeze(T_sp2cart(1,3,:)).*Bsphere(1,:)' + squeeze(T_sp2cart(2,3,:)).*Bsphere(2,:)' + squeeze(T_sp2cart(3,3,:)).*Bsphere(3,:)';


end