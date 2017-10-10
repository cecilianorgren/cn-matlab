function B = dipole(m,xmax,ymax,zmax)
% CN.DIPOLE    Takes a dipole moment, and returns mesh of dipole field.
%       B = (mu_0/4/pi)*3*r*(m dot r)/r^5-m/r^3)
%       [Bx,By,Bz] = CN.DIPOLE(m,xmax,ymax,zmax)
%           m - dipole moment 1x3 vector
%           xmax,ymax,zmax - max values over which to calculate B
%
%       NOT FUNCTIONING !!!

x = linspace(-xmax,xmax,100);
y = linspace(-ymax,ymax,100);
z = linspace(-zmax,zmax,100);

[X,Y,Z] = meshgrid(x,y,z);
R = sqrt(X.^2+Y.^2+Z.^2);
M = 

units = irf_units;
coeff = (units.mu0/4/pi);

Bx = @(m,X,R)( 3*X.*(bsxfun(@times,))./(R.^5) - m(1));
