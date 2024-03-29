function [F,Ffree,Ftrap] = get_f_free(V,n,vt,vd,PHI,VPH,species)
% [F,Ffree,Ftrap] = GET_F_FLAT(V,n,vt,vd,PHI,VPH)
units = irf_units;
if nargin > 6 && species == 1
  m = units.mp;
  q = units.e;
else
  m = units.me;
  q = -units.e;
end

% Collect input
nx = size(PHI,1);
v = V(1,:);
dv = v(2)-v(1);
nPop = numel(n);

% Initialize variables
F = V*0;
Ffree = V*0;
Ftrap = V*0;
V0 = V*0;

% Set up f0
f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
eval(f0_str)

% Get free and trapped indices
U = m*(V-VPH).^2/2 + q*PHI;
all_ind = 1:numel(V);
ifree = find(U > 0);
itrap = setdiff(all_ind,ifree);
iabove = find(V-VPH>0);
ibelow = find(V-VPH<0);

V0(iabove) = VPH(iabove) + ((V(iabove)-VPH(iabove)).^2 + 2*q*PHI(iabove)/m).^0.5; % same as Schamel free streaming
V0(ibelow) = VPH(ibelow) - ((V(ibelow)-VPH(ibelow)).^2 + 2*q*PHI(ibelow)/m).^0.5;
V0(itrap) = NaN;

% Get distributions
Ffree(ifree) = f0(V0(ifree));
Ftrap(itrap) = NaN;
F(ifree) = Ffree(ifree);
F(itrap) = Ftrap(itrap);
end