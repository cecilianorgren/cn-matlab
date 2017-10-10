function [R,PHI,X,Y] = makegrid(r1,r2,nr,phi1,phi2,nphi)
% make meshed grid based on input
% [R,PHI,X,Y] = makegrid(r1,r2,nr,phi1,phi2,nphi);

phi = linspace(phi1,phi2,nphi);
r = linspace(r1,r2,nr);

[R,PHI] = meshgrid(r,phi);
X = R.*cos(PHI);
Y = R.*sin(PHI);