function varargout = fe_bgk(v_orig,n,vt,vd,phi,vph_orig,obs_density_diff,phi_vec)
%%
units = irf_units;
vtrap = sqrt(2*units.e*phi/units.me);

vph = vph_orig;
v = v_orig; % move into reference frame of wave

E = units.me*((v-vph).^2 - vtrap.^2)/2;
all_ind = 1:numel(v);
ifree = find(E > 0);
itrap = setdiff(all_ind,ifree);
iabove = find(v-vph>0);
ibelow = find(v-vph<0);

fall = nan(size(v));
ffree = nan(size(v));
ftrap = nan(size(v));

v0 = v*0;

% set maxwellian at reference/infinity
switch numel(n)
  case 1
    f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2);
  case 2
    f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2) + ...
              n(2)*(1/pi./vt(2).^2)^(1/2)*exp(-(v-vd(2)).^2./vt(2).^2);
  case 3
    f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2) + ...
              n(2)*(1/pi./vt(2).^2)^(1/2)*exp(-(v-vd(2)).^2./vt(2).^2) + ...
              n(3)*(1/pi./vt(3).^2)^(1/2)*exp(-(v-vd(3)).^2./vt(3).^2);
  otherwise
    error(sprintf('%.0f populations not supported.',numel(n)))
end

% find v0
v0(iabove) = vph(iabove) + ((v(iabove)-vph(iabove)).^2 - vtrap(iabove).^2).^0.5; % same as Schamel free streaming
v0(ibelow) = vph(ibelow) - ((v(ibelow)-vph(ibelow)).^2 - vtrap(ibelow).^2).^0.5;
v0(itrap) = NaN;
%v0(ibelow) = NaN;

ffree(ifree) = f0(v0(ifree));
f.free = ffree;

dv = abs(v(1,1) - v(2,1));
nef = nansum(ffree*dv,1);

% ftrap
% get net from phi
% should have net as a function if e*phi
units = irf_units;
n_diff =  torow([obs_density_diff(1); obs_density_diff]);
phi_vec_diff = diff(phi_vec);
phi_diff = torow([phi_vec_diff(1); phi_vec_diff]);
net = sum(n) - nef + n_diff;
net_diff = [diff(net) net(end)-net(end-1)];
net_diff_ephi = net_diff./phi_diff/units.e;

net_diff_ephi(net_diff_ephi<0) = 0;
% from fitting
net_prime = @(phi) 1/0.0002333/phi;
 
if 0 % plot diagnostics
  %%
  figure(63)
  clear h;
  nrows = 6;
  ncols = 1;
  npanels = nrows*ncols;
  isub = 1;
  for ipanel = 1:npanels
    h(ipanel) = subplot(nrows,ncols,ipanel);    
  end
  isub = 1;
  
  hca = h(isub); isub = isub + 1;
  plot(hca,phi_vec)
  hca.YLabel.String = 'phi';
  
  hca = h(isub); isub = isub + 1;
  plot(hca,net)
  hca.YLabel.String = 'net';
  
  hca = h(isub); isub = isub + 1;
  plot(hca,net,phi_vec,'*')
  hca.XLabel.String = 'net';
  hca.YLabel.String = 'phi';
  
  hca = h(isub); isub = isub + 1;
  plot(hca,net_diff)
  hca.YLabel.String = '\Delta net';
  
  hca = h(isub); isub = isub + 1;
  plot(hca,phi_diff)
  hca.YLabel.String = '\Delta phi';
end
% loop over all x (second index)
for ix = 1:size(v,2)
  v_tmp = v(:,ix);
  E_tmp = E(:,ix);
  itrap_tmp = find(E_tmp < 0);
  itrap_center = find(E_tmp == min(E_tmp));
  
  integrand = @(E,s) sqrt(units.me/2)/pi*net_prime(phi_vec(ix))./sqrt(-E-s);  
  ns = 1000;
  
  % loop over v,E
  if ~isempty(itrap_tmp)
    % itrap_half goes from outer edge/separatrix to center
    itrap_half = itrap_tmp(1):itrap_center;
    all_s =          nan(numel(itrap_half),ns);
    all_E =          nan(numel(itrap_half),ns);
    all_integrands = nan(numel(itrap_half),ns);
    for iE = 1:numel(itrap_half)       
      E_tmp_tmp = E_tmp(itrap_half(iE));            
      s = linspace(0,-E_tmp_tmp*0.999,ns);
      all_s(iE,:) = s;
      all_E(iE,:) = E_tmp_tmp;
      all_integrands(iE,:) = integrand(E_tmp_tmp,s)/units.e;
      integral = trapz(s,integrand(E_tmp_tmp,s)/units.e);
      %integral = ffree(itrap_tmp(1)-1)-integral;
      %integral = -integral;
      ftrap(itrap_half(iE),ix) = integral;
      ftrap(itrap_center-itrap_half(iE)+itrap_center,ix) = integral;
      %all_trap = [itrap_half(iE) itrap_center-itrap_half(iE)+itrap_center];
    end    
    %ftrap(itrap_tmp,ix) = ffree(itrap_tmp(1)-1,ix)-ftrap(itrap_tmp,ix);
  end
  if phi_vec(ix)>400
    1;
  end
  if 0 % diagnostic plot
    %%
    figure(64)
    clear h;
    nrows = 8; ncols = 1; npanels = nrows*ncols;
    isub = 1;
    for ipanel = 1:npanels
      h(ipanel) = subplot(nrows,ncols,ipanel);    
    end
    isub = 1;

    hca = h(isub); isub = isub + 1;
    plot(hca,v(:,ix),E(:,ix))
    hca.XLabel.String = 'v'; hca.YLabel.String = 'E';
    hca.XLim = v([1 end],ix)
    
    hca = h(isub); isub = isub + 1;
    plot(hca,v(:,ix),E(:,ix),v(:,ix),v(:,ix)*0,v(:,ix),v(:,ix)*0-phi_vec(ix)*units.e)
    hca.XLabel.String = 'v'; hca.YLabel.String = 'E';
    hca.XLim = v([1 end],ix)
    hca.YLim = abs(min(E(:,ix)))*[-1.1 0.5];
    legend(hca,{'E','0','-e\phi'})
    %hca.YGrid = 'on';
    
    hca = h(isub); isub = isub + 1;
    plot(hca,v(:,ix),ftrap(:,ix))
    hca.XLabel.String = 'v'; hca.YLabel.String = 'f_{trap}';
    hca.XLim = v([1 end],ix)
    
    hca = h(isub); isub = isub + 1;
    plot(hca,v(:,ix),ffree(:,ix))
    hca.XLabel.String = 'v'; hca.YLabel.String = 'f_{free}';
    hca.XLim = v([1 end],ix)
    
    hca = h(isub); isub = isub + 1;
    plot(hca,v(:,ix),ffree(:,ix),v(:,ix),ftrap(:,ix))
    hca.XLabel.String = 'v'; hca.YLabel.String = 'f_{free}+f_{trap}';
    hca.XLim = v([1 end],ix)
        
    hca = h(isub); isub = isub + 1;
    [S,ETMP] = meshgrid(s,E_tmp(itrap_tmp));
    TOPLOT = -ETMP-S;
    TOPLOT(-S<ETMP) = NaN;
    surf(hca,s,E_tmp(itrap_tmp),TOPLOT)
    shading(hca,'flat')
    view(hca,[0 0 1])
    hca.XLabel.String = 's'; hca.YLabel.String = 'E';
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '-E-S';
      
    hca = h(isub); isub = isub + 1;
    [S,ETMP] = meshgrid(s,E_tmp(itrap_tmp));
    TOPLOT = -ETMP-S;
    TOPLOT(-S<ETMP) = NaN;
    TOPLOT = 1./sqrt(TOPLOT);
    surf(hca,s,E_tmp(itrap_tmp),(TOPLOT))
    shading(hca,'flat')
    view(hca,[0 0 1])
    hca.XLabel.String = 's'; hca.YLabel.String = 'E';
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = '1/(-E-S)^{1/2}';
        
    hca = h(isub); isub = isub + 1;
    surf(hca,all_s,all_E,log(all_integrands))
    shading(hca,'flat')
    view(hca,[0 0 1])
    hca.XLabel.String = 's'; hca.YLabel.String = 'E';
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log(integrand)';
      
    if 0
      hca = h(isub); isub = isub + 1;
      [S,ETMP] = meshgrid(s,E_tmp(itrap_tmp));
      INTEGRAND = (integrand(ETMP,S));
      INTEGRAND(-S>0.99*ETMP) = NaN;
      surf(hca,S,ETMP,INTEGRAND)
      shading(hca,'flat')
      view(hca,[0 0 1])
      hca.XLabel.String = 'E'; hca.YLabel.String = 's';
      hcb = colorbar('peer',hca);
      hcb.YLabel.String = 'log(integrand)';
    end
  end
end
%imagesc(ftrap)


fall(itrap) = ftrap(itrap);
fall(ifree) = ffree(ifree);
imagesc(fall)

%%
varargout{1} = fall;
if nargout > 0
  f.all = fall;
  f.free = ffree; % f.free(itrap) = NaN;
  f.trap = ftrap; % f.trap(ifree) = NaN;
  f.ifree = ifree;
  f.itrap = itrap;
  varargout{2} = f;
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
