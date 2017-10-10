% B has size nx*ny*nz
% rSurf has size 1*nr
% rGrid has size 1*(nr+1)
%% x, y, z coordinates
[X,Y] = meshgrid(xSurf,ySurf); % km
nth = 120; thSurf = linspace(0,2*pi*(1-1/nth),nth);
[fR,fTH] = meshgrid(rSurf,thSurf); % km
fX = fR.*cos(fTH);
fY = fR.*sin(fTH);
%sX = sX*1e3; sY = sY*1e3; sZ = sZ*1e3;    % m
%sR = sqrt(sX.^2 + sY.^2);                 % m