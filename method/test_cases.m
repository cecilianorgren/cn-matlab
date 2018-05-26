nx = 900;
ny = 1000;
xmax = 3;
ymax = 3;
x = linspace(-xmax,xmax,nx);
y = linspace(-ymax,ymax,ny);
[X,Y] = meshgrid(x,y);

% physcial constant
units = irf_units;
e = units.e;

% magnetic field
Bz = 10;
Babs = abs(Bz);

% density
n = 1;

% potential
Z = peaks(X,Y);
PHI = Z;

% electric field
nE = 100;
dx = x(2)-x(1);
dy = y(2)-y(1);
[Ex,Ey] = gradient(PHI,dx,dy);

% ExB 
ExBx = Ey*Bz/Babs/Babs;
ExBy = -Ex*Bz/Babs/Babs;

% J from ExB
Jx = -n*e*ExBx;
Jy = -n*e*ExBy;

% parallel wave B from J
% rot(B) = mu0*J
% dz = infty
% dy Bz = mu0 Jx = -mu0*n*e*ExBx
% dx Bz = -mu0 Jz = mu0*n*e*ExBy
% use pauls method for finding A from B: rot(A) = B
% Bz = -mu0*n*e*ExBx/dy



% Plot
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
ipanel = 0;
for irow = 1:nrows
  for icol = 1:ncols
    ipanel = ipanel + 1;
    h(ipanel) = subplot(nrows,ncols,ipanel);
  end
end

isub = 1;
if 1 % PHI
  hca = h(isub); isub = isub + 1;
  surf(hca,X,Y,PHI)
  shading(hca,'flat')
  view(hca,[0 0 1])
  hca.Title.String = 'PHI';
end
if 0 % Ex
  hca = h(isub); isub = isub + 1; 
  surf(hca,X,Y,Ex)
  shading(hca,'flat')
  view(hca,[0 0 1])
  hca.Title.String = 'E_x';
end
if 0 % Ey
  hca = h(isub); isub = isub + 1; 
  surf(hca,X,Y,Ey)
  shading(hca,'flat')
  view(hca,[0 0 1])
  hca.Title.String = 'E_y';
end
if 1 % Phi and Ex,Ey
  hca = h(isub); isub = isub + 1; 
  dxE = 50;
  dyE = 50;
  ix = 1:dxE:nx;
  iy = 1:dyE:ny;
  hcont = contour(hca,X,Y,PHI);
  hold(hca,'on')
  hquiv = quiver(hca,X(iy,ix),Y(iy,ix),Ex(iy,ix),Ey(iy,ix));
  hold(hca,'off')
  hca.Title.String = 'E';
end
if 1 % Phi and ExBx,ExBy
  hca = h(isub); isub = isub + 1; 
  dxJ = 50;
  dyJ = 50;
  ix = 1:dxJ:nx;
  iy = 1:dyJ:ny;
  hcont = contour(hca,X,Y,PHI);
  hold(hca,'on')
  hquiv = quiver(hca,X(iy,ix),Y(iy,ix),ExBx(iy,ix),ExBy(iy,ix));
  hold(hca,'off') 
  hca.Title.String = 'V_{ExB}'; 
end
if 1 % Phi and Jx,Jy
  hca = h(isub); isub = isub + 1; 
  dxJ = 50;
  dyJ = 50;
  ix = 1:dxJ:nx;
  iy = 1:dyJ:ny;
  hcont = contour(hca,X,Y,PHI);
  hold(hca,'on')
  hquiv = quiver(hca,X(iy,ix),Y(iy,ix),Jx(iy,ix),Jy(iy,ix));
  hold(hca,'off') 
  hca.Title.String = 'J_{ExB}'; 
end
