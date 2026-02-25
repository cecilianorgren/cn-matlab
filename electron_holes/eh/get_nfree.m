function [nfree,Ffree] = get_f_flat(V,n,vt,vd,PHI,VPH)
units = irf_units;

% Collect input
%phi = PHI(:,1);
%vph = VPH(1,1);
nx = size(PHI,1);
v = V(1,:);
dv = v(2)-v(1);
nPop = numel(n);

%f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2) + ...
%          n(2)*(1/pi./vt(2).^2)^(1/2)*exp(-(v-vd(2)).^2./vt(2).^2);

f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop)',4,1))];
f0_str = [f0_str(1:end-1) ';'];
eval(f0_str)
fsep = f0(vph);

F = V*0;
Ffree = V*0;
Ftrap = V*0;
Ftrap_flat = V*0;
V0 = V*0;

E = units.me*(V-vph).^2/2 - units.e*PHI;
all_ind = 1:numel(V);
ifree = find(E > 0);
itrap = setdiff(all_ind,ifree);
iabove = find(V-vph>0);
ibelow = find(V-vph<0);

V0(iabove) = vph + ((V(iabove)-VPH).^2 - 2*units.e*(PHI(iabove))/units.me).^0.5; % same as Schamel free streaming
V0(ibelow) = vph - ((V(ibelow)-VPH).^2 - 2*units.e*(PHI(ibelow))/units.me).^0.5;
V0(itrap) = NaN;

% free particles
Ffree(ifree) = f0(V0(ifree));
Ftrap(itrap) = fsep;
  
end