function [Er Ez]  = E(phi0,r,z,lr,lz)
% Electrc field profile from Tao2011, and Chen2004.
%   [Er Ez] = art2.E(phi0,r,z,lr,lz)

Er = (r/lr^2).*phi0*exp(-0.5*(r/lr).^2-0.5*(z/lz).^2);
Ez = (z/lz^2).*phi0.*exp(-0.5*(r/lr).^2-0.5*(z/lz).^2);