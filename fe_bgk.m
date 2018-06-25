function fout = fe_bgk(v_orig,n,vt,vd,phi,vph_orig,return_mode)
%%
units = irf_units;
vtrap = sqrt(2*units.e*phi/units.me);

vph = vph_orig;
v = v_orig; % move into reference frame of wave

all_ind = 1:numel(v);
ifree = find((v-vph).^2 - vtrap.^2 > 0);
itrap = setdiff(all_ind,ifree);
iabove = find(v-vph>0);
ibelow = find(v-vph<0);

fout = nan(size(v));
v0 = v*0;

% set maxwellian at reference/infinity
f0 = @(v) n*(1/pi./vt.^2)^(1/2)*exp(-v.^2./vt.^2);

% find v0
v0(iabove) = vph(iabove) + ((v(iabove)-vph(iabove)).^2 - vtrap(iabove).^2).^0.5; % same as Schamel free streaming
v0(ibelow) = vph(ibelow) - ((v(ibelow)-vph(ibelow)).^2 - vtrap(ibelow).^2).^0.5;
v0(itrap) = NaN;
%v0(ibelow) = NaN;

fout(ifree) = f0(v0(ifree));

if exist('return_mode','var')
  switch return_mode
    case 'trap'
      fout = fout(itrap);
    case 'free'
      fout = fout(ifree);
    case 'all'
      fout = fout;
    otherwise
      fout = fout;
  end
end

if 0
  %%
  hca = subplot(3,1,1);
  pcolor(hca,vtrap)
  shading flat
  colorbar('peer',hca)
  hca = subplot(3,1,2);
  pcolor(hca,v0)
  shading flat
  colorbar('peer',hca)
  hca = subplot(3,1,3);
  pcolor(hca,fout)
  shading flat
  colorbar('peer',hca)
end
