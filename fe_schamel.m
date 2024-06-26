function varargout = fe_schamel(v,n,vt,vd,phi,vph,beta)

units = irf_units;
vtrap = sqrt(2*units.e*phi/units.me);

v = v-vph; % move into reference frame of wave
vph = vph-vd;

all_ind = 1:numel(v);
ifree = find((v).^2 - vtrap.^2 > 0);
itrap = setdiff(all_ind,ifree);
fout = nan(size(v));
%vd = 0;

fout(ifree) = n*(1/pi./vt.^2)^(1/2)*exp(-( sign(v(ifree)).* sqrt( v(ifree).^2 - vtrap(ifree).^2) + vph(ifree) ).^2./vt.^2);
fout(itrap) = n*(1/pi./vt.^2)^(1/2)*exp(-beta*(v(itrap).^2 - vtrap(itrap).^2)./vt.^2 - vph(itrap).^2./vt.^2);

%fout(ifree) = n*(1/pi./vt.^2)^(1/2)*exp(-( sign(v(ifree)).* sqrt( v(ifree).^2 - vtrap(ifree).^2) + vph(ifree) ).^2./vt.^2);
%fout(itrap) = n*(1/pi./vt.^2)^(1/2)*exp(-beta*(v(itrap).^4 - vtrap(itrap).^4)./vt.^4 - vph(itrap).^4./vt.^4);


varargout{1} = fout;
if nargout > 0
  f.all = fout;
  f.free = fout; f.free(itrap) = NaN;
  f.trap = fout; f.trap(ifree) = NaN;
  f.ifree = ifree;
  f.itrap = itrap;
  varargout{2} = f;
end