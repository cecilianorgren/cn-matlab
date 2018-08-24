function [ntrap_flat,Ftrap_flat] = get_ntrap_flat(V,n,vt,vd,PHI,VPH)
  units = irf_units;
  phi = PHI(:,1);
  vph = VPH(1,1);
  nx = size(PHI,1);
  dv = V(1,2) - V(1,1);
  
  f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2) + ...
            n(2)*(1/pi./vt(2).^2)^(1/2)*exp(-(v-vd(2)).^2./vt(2).^2);
  fsep = f0(vph);
  E = units.me*(V-vph).^2/2 - units.e*PHI;
  all_ind = 1:numel(V);
  itrap = find(E < 0);  
  Ftrap_flat = V*0;
  Ftrap_flat(itrap) = fsep;
  ntrap_flat = sum(Ftrap_flat,2)*dv;
end