function [Ex,Ey,Ez]  = E3(x,y,z,phi0,lr,lz)
% Electrc field profile from Tao2011, and Chen2004.
%   [Ex,Ey,Ez] = art2.E3(x,y,z,phi0,lr,lz)

Ex = (x/lr^2)*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ey = (y/lr^2)*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ez = (z/lz^2)*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);