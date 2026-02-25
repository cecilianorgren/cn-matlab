function varargout = gyrokinetic_scaling(n,B)
% dperp/dpar = sqrt(1+wpe^2/wce^2);
%  dperp/dpar = GYROKINETIC_SCALING(n,B);
%  [dperp/dpar,wpe/wce] = GYROKINETIC_SCALING(n,B);
%   n - cc
%   B nT

units = irf_units; 
n_SI = n*1e-6;
B_SI = B*1e-9;

wpe = sqrt(units.e*n_SI^2/units.eps0/units.me);
wce = units.e*B_SI/units.me;

ratio = sqrt(1+wpe^2/wce^2);

varargout{1} = ratio;
if nargout == 2
  varargout{2} = wpe/wce;
end