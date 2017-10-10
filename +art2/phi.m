function out = phi(phi0,x,y,z,lr,lz)
% Potential profile from Tao2011, and Chen2004.
%   out = phi(phi0,x,y,z,lr,lz);

out = phi0*exp(-0.5*(x.^2+y.^2)/lr^2-0.5*(z/lz).^2);