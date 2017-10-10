RE = 6400*1e3;
B0 = 1e-9;
K = 50000e-9*(RE)^3;
A_dipole = @(theta,r) -K*sin(theta)./r^2;
A_imf = @(theta,r) 0.5*B0*r.^1.*sin(theta);

Xmax = 20;
Ymax = Xmax;
x = linspace(-Xmax*RE,Xmax*RE,500);
y = linspace(-Ymax*RE,Ymax*RE,500);

[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);
THETA = acos(Y./R);

nRows = 3;
nCols = 1;
isub = 1;

Aphi = A_dipole(THETA,R);
Ax = Aphi.*cos(THETA);
Ay = Aphi.*sin(THETA);


hca = subplot(nRows,nCols,isub); isub = isub + 1;
contour(hca,X,Y,A_dipole(THETA,R),20)
colorbar; axis square;

hca = subplot(nRows,nCols,isub); isub = isub + 1;
contour(hca,X,Y,A_imf(THETA,R),20)
colorbar; axis square;

hca = subplot(nRows,nCols,isub); isub = isub + 1;
contour(hca,X,Y,A_dipole(THETA,R)+A_imf(THETA,R),20)
%pcolor(hca,X,Y,A_dipole(THETA,R)+A_imf(THETA,R))
colorbar; axis square; shading flat;




