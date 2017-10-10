function cB = tdB(B0,n,phi0,lr,lz)
% Calculates the theoretical center magnetic field of an electron hole.
%   cB = tdB(B0,n,phi0,lr,lz)

mu0 = 1.2566e-6; e = 1.6022e-19;
cB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT