function B = dB(phi0,n,B0,lr,lz)
% Calculates the center magnetic field, see Tao2011.
%   B = dB(phi0,n,B0,lr,lz)
%       phi - center potential, V
%       n - density, cc
%       B0 - ambient magnetic field strength, nT
%       lr - perpendicular half width, km
%       lz - parallel half width, km

e = 1.6022e-19;
mu0 = 1.257e-6;

B = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT