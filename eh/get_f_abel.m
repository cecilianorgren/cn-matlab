function [Fsave,Ffreesave,Ftrapsave] = get_f_abel(V,n,vt,vd,PHI,VPH,nt)  
% [F,Ffree,Ftrap] = GET_F_ABEL(V,n,vt,vd,PHI,VPH,nt);  

units = irf_units;
phi = PHI(:,1);
vph = VPH(1,1);
nx = size(PHI,1);
v = V(1,:);
dv = v(2) - v(1);
  
% Set up f0
nPop = numel(n);
f0_str = ['f0 = @(v) ' sprintf('n(%g)*(1/pi./vt(%g).^2)^(1/2)*exp(-(v-vd(%g)).^2./vt(%g).^2)+',repmat((1:nPop),4,1))];
f0_str = [f0_str(1:end-1) ';'];
eval(f0_str)
          
if numel(unique(VPH)) == 1  
  fsep = f0(VPH(1,1)); % only works for single vph  
else
  error('Only single vph supported.')
end
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

V0(iabove) = vph + ((V(iabove)-vph).^2 - 2*units.e*(PHI(iabove))/units.me).^0.5; % same as Schamel free streaming
V0(ibelow) = vph - ((V(ibelow)-vph).^2 - 2*units.e*(PHI(ibelow))/units.me).^0.5;
V0(itrap) = NaN;

% free particles
Ffree(ifree) = f0(V0(ifree));
if 0 % plot
  figure(101)
  hca = subplot(2,1,1);
  E_ = E; E_(itrap) = NaN;
  pcolor(hca,1:nx,v,E_'); shading flat;
  hca = subplot(2,1,2);
  V0_ = V0; V0_(ibelow) = NaN;
  pcolor(hca,1:nx,v,real(V0_')); shading flat;
end

% trapped flattop  
Ftrap_flat(itrap) = fsep;
ntrap_flat = sum(Ftrap_flat,2)*dv;

if 0
  t_indices = 769:797;
  disp(sprintf('Used indices %g:%g.',t_indices(1),t_indices(end)))
  [fitresult, gof, fun_net, fun_net_prime] = createFit(phi(t_indices), nt(t_indices));
else
  [fitresult, gof, fun_net, fun_net_prime] = createFit(tocolumn(phi), tocolumn(nt));
end
%[fitresult, gof, fun_net, fun_net_prime] = createFit(phi(phi>10), nt(phi>10));
%phi_fit = phi;
%max_ind = find(phi_fit == max(phi_fit));
%[fitresult, gof, fun_net, fun_net_prime] = createFit(phi(max_ind:end), nt(max_ind:end));

if 0 % plot
  figure(22)
  plot(phi,nt,phi,fun_net(phi))
end

%fun_net_prime = @(V,a,b,c,d) b*a*V.^(b-1) + d*c*V.^(d-1);  
for ix = 1:nx % something wrong with uneven numel(itrap_ngtvE)
  itrap_ngtvE = find(E(ix,:) < 0); % perhaps find closest E=0
  if itrap_ngtvE > 0
    if mod(numel(itrap_ngtvE),2) % uneven
      i_midE = ceil(numel(itrap_ngtvE)/2);
    else
      i_midE = ceil(numel(itrap_ngtvE)/2) + 1;
    end
    itrap_v = intersect(itrap_ngtvE,itrap_ngtvE(i_midE):numel(E(ix,:)));
  else
    continue;      
  end
  %fprintf('ix = %g, numel(itrap_ngtvE) = %g\n',ix,numel(itrap_ngtvE))

  Etrap = E(ix,itrap_v); 
  nu = numel(itrap_v);
  phi_tmp = phi(ix);

  fun_int = @(u,V) (1/units.e)*(units.me/2)^0.5/pi*fun_net_prime(V/units.e).*(-(units.me*(u-vph).^2/2-units.e*phi_tmp)-V).^-0.5;
  %fun_int = @(u,V) (1/units.e)*(2*units.me)^0.5/pi*fun_net_prime(V/units.e,a,b,c,d).*(-(units.me*u.^2/2-units.e*phi_tmp)-V).^-0.5;    
  for iu_ = 1:nu
    u_tmp = v(itrap_v(iu_));
    nV = 500;    
    Vmax = units.e*phi_tmp-units.me*(u_tmp-vph)^2/2;    
    Vmin = 0.000*Vmax;
    V_tmp = linspace(Vmin,Vmax*0.9995,nV);

    if 1 % integral()
      fun_int_tmp = @(V)fun_int(u_tmp,V);
      f_tmp = integral(fun_int_tmp,0,Vmax);
    else % trapz()
      INT = fun_int(u_tmp,V_tmp);
      f_tmp = trapz(V_tmp,INT);
    end
    %f_tmp_cum = cumtrapz(V_tmp,INT);
    %plot(V_tmp,f_tmp_cum)
    %hold(gca,'on')

    ivpos = itrap_v(iu_);
    if mod(numel(itrap_ngtvE),2) % uneven
      ivneg = itrap_v(1)-0-(ivpos-itrap_v(1));
    else % even
      ivneg = itrap_v(1)-1-(ivpos-itrap_v(1));
    end
    Ftrap(ix,ivpos) = f_tmp;
    Ftrap(ix,ivneg) = f_tmp;      

    if 0 % plot
      hca = subplot(2,1,1);
      plot(hca,V_tmp,INT)
      hca.Title.String = sprintf('ix=%g, iu=%g (%g)',ix,iu_,nu);
      %hca = subplot(2,1,2);
      %plot(hca,u_tmp,cumsum(INT),u_tmp,real(cumsum(INT)),u_tmp,imag(cumsum(INT)))
      %legend(hca,{'all','real','imag'})
      %plot(hca,u,Ftrap(ix,itrap_v))
      %hca.Title.String = sprintf('ix=%g, iu=%g (%g)',ix,iu_,nu);
      1;
    end     
  end
  if 0 % start from f0 at center of trapping region, then Liouville from there, in a similar way as outside
    V0_(iabove) = 0 + ((V(iabove)).^2 + 2*units.e*(phimax-PHI(iabove))/units.me).^0.5; % same as Schamel free streaming
    V0_(ibelow) = 0 - ((V(ibelow)).^2 + 2*units.e*(phimax-PHI(ibelow))/units.me).^0.5;

    V0(itrap) = V0_(itrap);
    %V0(ifree) = NaN;
    f0_inner = Ftrap(fix(nx/2),:);
    for itrap_tmp_ = 1:numel(itrap)
      itrap_tmp = itrap(itrap_tmp_);
      iv_trap = find(abs(v-V0(itrap_tmp)) == min(abs(v-V0(itrap_tmp))));
      Ftrap(itrap_tmp) = f0_inner(iv_trap);
    end
  end
end

% collect: F = Ffree + Ftrap
F(ifree) = Ffree(ifree);
%F(itrap) = Ftrap(itrap) + Ftrap_flat(itrap); decide this outside
F(itrap) = Ftrap(itrap);

Fsave = F;
Ffreesave = Ffree;
Ftrapsave = Ftrap + Ftrap_flat;% Ftrapsave(itrap) = Ftrapsave(itrap) + Ftrap_flat(itrap);
end
