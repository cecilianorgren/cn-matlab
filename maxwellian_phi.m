function fout = maxwellian_phi(v,n,vt,vd,phi,vph,beta)

units = irf_units;
vtrap = sqrt(2*units.e*phi/units.me);

all_ind = 1:numel(v);
ifree = find((v).^2 - vtrap.^2 > 0);
itrap = setdiff(all_ind,ifree);
fout = nan(size(v));
vd = 0;
%beta = -10;

%fout(ind_grt_zero) = n*(1/pi/vt^2)^(1/2)*exp(-(sqrt((v(ind_grt_zero)-vd).^2 - (sqrt(2*units.e*phi(ind_grt_zero)/units.me) +0*vph).^2) + 0*vph).^2/vt^2);
%fout(ind_grt_zero) = n*(1/pi/vt^2)^(1/2)*exp(-(sqrt((v(ind_grt_zero)-vd).^2 - 2*units.e*phi(ind_grt_zero)/units.me) + 1*vph).^2/vt^2);
fout(ifree) = n*(1/pi/vt^2)^(1/2)*exp(-( sign(v(ifree)).* sqrt( (v(ifree)-vd).^2 - vtrap(ifree).^2) + 1*vph(ifree) ).^2/vt^2);
fout(itrap) = n*(1/pi/vt^2)^(1/2)*exp(-beta*((v(itrap)-vd).^2 - vtrap(itrap).^2)/vt^2 - vph(itrap).^2/vt^2);

% 
% if v.^2 - 2*units.e*phi/units.me > 0
%   fout = n*(1/pi/vt^2)^(1/2)*exp(-(sqrt(v.^2 - 2*units.e*phi/units.me) + vph).^2/vt^2);
% else
%   fout = n*(1/pi/vt^2)^(1/2)*exp(-(sqrt(v.^2 - 2*units.e*phi/units.me)).^2/vt^2 - 0.5*vph^2/vt^2);
% end