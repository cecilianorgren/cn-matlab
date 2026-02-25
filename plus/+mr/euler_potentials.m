

% Magnetic field of simple x-line geometry
Bx = @(x,y,z) x*0 + y*1 + z*0;
By = @(x,y,z) x*1 + y*0 + z*0;
Bz = @(x,y,z) x*0 + y*0 + z*0 + 1;

% Euler potentials, constant along field lines
A = @(x,y,z) x.*cosh(z) - y.*sinh(z);
B = @(x,y,z) y.*cosh(z) - x.*sinh(z);

xmax = 1;
ymax = 1;

x = linspace(-xmax,xmax,100);
y = linspace(-ymax,ymax,100);

[X,Y] = meshgrid(x,y);
Z = ones(size(X));

nrows = 3;
ncols = 2;
npanels = nrows*ncols;
for ip = 1:npanels
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

if 1 % B, magnetic field quivers
  hca = h(isub); isub = isub + 1;
  contour(hca,A(X,Y,Z))
end
if 1 % B, magnetic field amplitude
  hca = h(isub); isub = isub + 1;
  contour(hca,sqrt(X.^2+Y.^2+Z.^2))
end
if 1 % A, euler potential
  hca = h(isub); isub = isub + 1;
  contour(hca,A(X,Y,Z))
end
if 1 % B, euler potential
  hca = h(isub); isub = isub + 1;
  contour(hca,B(X,Y,Z))
end
if 1 % A^2-B^2, = zero at separatrix
  hca = h(isub); isub = isub + 1;
  contour(hca,A(X,Y,Z).^2-B(X,Y,Z).^2)
end