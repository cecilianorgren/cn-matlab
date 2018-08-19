%% Set up grid and phi
units = irf_units;
fun_phi = @(phimax,x,lx) phimax*exp(-x.^2/2/lx.^2);

data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(obs_velocity);
c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
% charge separation assuming potential structure is single gaussian
dn = units.eps0/units.e*obs_potential_max./(obs_lpp*1e3)*1e-6; % cc


%for ih = 2:numel(obs_velocity)
phi_input = 3;
switch phi_input
  case 1
    lx = 1300; 
    nx = 207;
    x = linspace(-5*lx,5*lx,nx);
    phimax = 100;
    phi = fun_phi(phimax,x,lx);
    vph = -10000e3;
  case 2 % single one from observations
    phimax = obs_potential_max(ih);
    lx = obs_lpp(ih)/2*1e3;
    nx = 200;
    x = linspace(-5*lx,5*lx,nx);
    vph = obs_velocity(ih)*1e3;
    phi = fun_phi(phimax,x,lx);
  case 3 % all observations
    %load /Users/cecilia/Data/20170706_135303_basic_eh
    tint_phi = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.620Z');
    t0 = tint_phi(1);
    ffilt = 30; % Hz
    c_eval('Etoint? = gseE?par;')
    c_eval('Etoint? = Etoint?.filt(ffilt,0,[],3).tlim(tint_phi);');    
    c_eval('intEdt? = irf_integrate(Etoint?);');    
    minpeakdistance = 50;
    c_eval('[PKS?,LOCS?,W?] = findpeaks(intEdt?.data,''MinPeakDistance'',minpeakdistance);')
    c_eval('intEdt?_detrend = intEdt?; intEdt?_detrend.data = detrend(intEdt?_detrend.data,''linear'',LOCS?);')

    % Plotting options
    doT = 1; % otherwise plot x;
    x_vec = intEdt1.time - t0; % seconds
    dx_vec = x_vec(2)-x_vec(1);
    nx = numel(x_vec);
    x_vec_diff1 = x_vec(1:end-1)+0.5*dx_vec;
    x_vec_diff2 = x_vec(2:end-1);

    % Remove phi baselevel
    %c_eval('phi_baselevel = interp_linear_piecewise(intEdt?.data,x_vec,x_vec(LOCS?));')
    c_eval('ts_locs? = irf.ts_scalar(intEdt?.time([1; LOCS?; end]),intEdt?.data([1; LOCS?; end]));')
    %c_eval('phi_baselevel? = interp_linear_piecewise(intEdt?.data([1; LOCS?; end]),x_vec([1; LOCS?; end]),x_vec);')
    %c_eval('ts_phi_baselevel? = irf.ts_scalar(intEdt?.time,phi_baselevel?);')
    c_eval('ts_phi_baselevel? = ts_locs?.resample(intEdt?);')
    c_eval('intEdt?_detrend = intEdt?-ts_phi_baselevel?;')
end

% F0
n = [0.02 0.02]*1e6; % m-3
ntot = 0.04*1e6;
R = 0.7;
n = [R (1-R)]*ntot;
T = [150 4000]; % eV
vd = [-10000e3 4000e3]; % ms-1 
vt = sqrt(2*units.e*T./units.me); % m/s
n0 = sum(n);

% Get electric field and n
dx =  x(2) - x(1);
x_diff1 = x(1:end-1) + 0.5*dx;
x_diff2 = x(2:end-1) + dx;
Efield = -diff(phi)/dx;
dn = diff(phi,2)/dx/dx*units.eps0/units.e;
dn = [dn(1); dn'; dn(end)]; % assume phi -> at edges

vtrap = sqrt(2*units.e*phimax/units.me);
nv = 2000;
vmax = max([5*vtrap, vd + 2*vt]); 
v = linspace(-vmax,vmax,nv);
dv = v(2) - v(1);

[X,V] = meshgrid(x,v); X = permute(X,[2 1]); V = permute(V,[2 1]);
PHI = fun_phi(phimax,X,lx);
VPH = V*0 + vph;
E = units.me*(V-vph).^2/2 - units.e*PHI;


[nfree,Ffree] = get_nfree(V,n,vt,vd,PHI,VPH);
[ntrap_flat,Ftrap_flat] = get_ntrap_flat(V,n,vt,vd,PHI,VPH);
ntrap = n0 - nfree + dn;
dntrap = ntrap - ntrap_flat;
Fflat = V*0; Fflat(E<0) = Ftrap_flat(E<0); Fflat(E>0) = Ffree(E>0);

%plot(x,nfree,x,ntrap,x,n0+dn)

[F_abel,F_abel_free,F_abel_trap] = get_f_abel(V,n,vt,vd,PHI,VPH,ntrap);
%F_abel = F_abel + Fflat;

[F_scha,F_scha_free,F_scha_trap,beta] = get_f_schamel(V,n,vt,vd,PHI,VPH,ntrap);
%plot(x,ntrap,x,nansum(F_scha_trap,2)*dv,x,nansum(F_abel_trap,2)*dv)

% Plot
figure(92)
nrows = 3;
ncols = 3;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;

vlim = 30000e3;
if 1 % F flat
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V,Fflat)
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{flat} (...)';  
  hca.YLim = vlim*[-1 1];
  colormap(hca,cn.cmap('white_blue'))    
end
if 1 % F_trap - Fflat
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V,F_abel-Fflat)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{Abel}-f_{flat} (...)';  
  colormap(hca,cn.cmap('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.YLim = vlim*[-1 1];
end
if 1 % F abel
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V,real(F_abel))  
  shading(hca,'flat') 
  if 0 % E contours
    hold(hca,'on')
    levels_E = linspace(min(E(:)),0,6); 
    %levels_E = logspace(log(min(E(:))),log(max(E(:))),100);     
    levels_E = levels_E + min(abs(levels_E));
    hc = contour(hca,X,V,E,levels_E);
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{Abel} (...)';  
  hca.YLim = vlim*[-1 1];
  colormap(hca,cn.cmap('white_blue'))    
end
if 1 % F schamel
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V,F_scha)  
  shading(hca,'flat') 
  if 0 % E contours
    hold(hca,'on')
    levels_E = linspace(min(E(:)),0,6); 
    %levels_E = logspace(log(min(E(:))),log(max(E(:))),100);     
    levels_E = levels_E + min(abs(levels_E));
    hc = contour(hca,X,V,E,levels_E);
    hold(hca,'off')
  end
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{Scham} (...)';
  hca.YLim = vlim*[-1 1];
  colormap(hca,cn.cmap('white_blue'))    
end
if 1 % F_Abel - F_Scha
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V,F_abel-F_scha)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{Abel}-f_{Scha} (...)';  
  colormap(hca,cn.cmap('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.YLim = vlim*[-1 1];
end
if 1 % Trapped densities
  hca = h(isub); isub = isub + 1;
  plot(hca,x,ntrap,x,nansum(F_abel_trap,2)*dv,x,nansum(F_scha_trap,2)*dv)
  hca.YLabel.String = 'n_t';
  legend(hca,{'n_0-n_f+\delta n','n_{t,Abel}','n_{t,Scha}'},'Box','off')
  irf_legend(hca,{sprintf('B_{Scha}=%g',beta)},[0.99 0.1],'color',[0 0 0])
end
if 1 % Total densities
  hca = h(isub); isub = isub + 1;
  plot(hca,x,n0+dn,x,nansum(F_abel,2)*dv,x,nansum(F_scha,2)*dv)
  hca.YLabel.String = 'n';
  legend(hca,{'n_0+dn','n_{Abel}','n_{Scha}'},'Box','off')
  irf_legend(hca,{sprintf('B_{Scha}=%g',beta)},[0.99 0.1],'color',[0 0 0])
end
if 1 % nt(phi)
  hca = h(isub); isub = isub + 1;
  [fitreslts,gof,fun_net,fun_net_prime] = createFit(torow(phi),torow(dntrap));
%     a = fitreslts.a;
%     b = fitreslts.b; % b = 0.5;
%     c = fitreslts.c;
%     d = fitreslts.d;
  plot(hca,phi,dntrap,'.',phi,fun_net(phi))
  hca.YLabel.String = 'n_t (m^{-3})';  
  hca.XLabel.String = '\phi (V)';  
end
if 1 % F at a few x's
  hca = h(isub); isub = isub + 1; 
  colors = mms_colors('matlab');  
  set(hca,'ColorOrder',colors(1:3,:))
  indx = fix(nx./[6 4 2]);
  %set(hca,'LineStyleOrder',{'-','--'})
%   plot(hca,v,F_abel(fix(nx/6),:),v,F_abel(fix(nx/4),:),v,F_abel(fix(nx/2),:),...
%            v,F_scha(fix(nx/6),:),v,F_scha(fix(nx/4),:),v,F_scha(fix(nx/2),:))
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,v,F_abel(indx(1),:),v,F_abel(indx(2),:),v,F_abel(indx(3),:))
  hold(hca,'on'),
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,v,F_scha(indx(1),:),'--',v,F_scha(indx(2),:),'--',v,F_scha(indx(3),:),'--')
  hold(hca,'off')
  %hca.YScale = 'log';
  hca.YLabel.String = 'F';  
  hca.XLabel.String = 'v';
  hca.XLim = vph+vtrap*[-1 1]*1.5;
  irf_legend(hca,{sprintf('x=%.0f',x(indx(1))),sprintf('x=%.0f',x(indx(2))),sprintf('x=%.0f',x(indx(3)))},[0.1 0.15])
  irf_legend(hca,{'- Abel','--Scha'},[0.1 0.05],'color',[0 0 0])
end
cn.print(sprintf('AbelScha_obs_eh%g',ih))
%end
%% Plot
figure(91)
nrows = 8;
ncols = 1;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;

if 1 % phi(x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x,phi)
  hca.XLabel.String = 'x';
  hca.YLabel.String = '\phi (V)';
end
if 1 % E
  hca = h(isub); isub = isub + 1;
  levels_E = linspace(min(E(:)),max(E(:)),20); levels_E = levels_E + min(abs(levels_E));
  %levels_E = [0 linspace(min(E(:)),max(E(:)),15)];
  contourf(hca,X,V,E/units.e,levels_E/units.e)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'E/e (V)';
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colormap(hca,cn.cmap('blue_red'))  
end
if 1 % electric field (x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x_diff1,Efield)
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'E (V/m)';
end
if 1 % dn(x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x,dn*1e-6)
  hca.XLabel.String = 'x';
  hca.YLabel.String = '\delta n (cc)';  
end
if 1 % dn(phi)
  hca = h(isub); isub = isub + 1;
  hp = plot(hca,phi,dn*1e-6,phi,nt*1e-6,phi,fun_net(phi,a,b,c,d)*1e-6);
  hp(2).LineWidth = 1.5;
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = '\delta n (cc)';
  legend(hca,{'\delta n','n_{et}','fit to n_{et}'},'Box','off')  
end
if 1 % Ffree
  hca = h(isub); isub = isub + 1;
  levels_E = linspace(min(E(:)),max(E(:)),20); levels_E = levels_E + min(abs(levels_E));
  %levels_E = [0 linspace(min(E(:)),max(E(:)),15)];
  pcolor(hca,X,V,Ffree)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{free} (...)';
  %hca.CLim = max(abs(hca.CLim))*[0 1];
  colormap(hca,cn.cmap('white_blue'))  
  %colormap(hca,'jet')  
end
if 1 % n(x)
  hca = h(isub); isub = isub + 1;
  hp = plot(hca,x,x*0+n0*1e-6,x,dn*1e-6,x,nf*1e-6,x,nt*1e-6);  
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'n (cc)';
  legend(hca,{'n_0, n_i','n_e-n_i','n_{ef}','n_{et} = n_i-n_{ef}+\delta n'},'Box','off')
end
if 1 % dn/dphi
  hca = h(isub); isub = isub + 1;
  plot(hca,x,fun_net_prime(phi,a,b,c,d))
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'dn_{et}/d\phi (m^{-3}V^{-1})';
end
if 0 % -E-e*PHI
  hca = h(isub); isub = isub + 1;
  %levels_E = linspace(min(E(:)),max(E(:)),20); levels_E = levels_E + min(abs(levels_E));
  %levels_E = [0 linspace(min(E(:)),max(E(:)),15)];
  pcolor(hca,X,V,INT_A/units.e)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = '(-E-e\phi)/e (...)';
  %hcb.YLabel.String = '(-E-e\phi)^{-1/2} (...)';
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  %colormap(hca,cn.cmap('blue_red'))  
end
if 0 % 1/sqrt(-E-PHI)
  hca = h(isub); isub = isub + 1;
  %levels_E = linspace(min(E(:)),max(E(:)),20); levels_E = levels_E + min(abs(levels_E));
  %levels_E = [0 linspace(min(E(:)),max(E(:)),15)];
  contourf(hca,X,V,INT_A/units.e)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = '(-E-e\phi)/e (...)';
  %hcb.YLabel.String = '(-E-e\phi)^{-1/2} (...)';
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  %colormap(hca,cn.cmap('blue_red'))  
end

%% Aid functions
function [nfree,Ffree] = get_nfree(V,n,vt,vd,PHI,VPH)
  beta = 1; % not used for free streaming population
  nPop = numel(n);F = zeros(size(V));
  Ffree = zeros(size(V));
  Ftrap = zeros(size(V));
  dv = V(1,2) - V(1,1);
  for iPop = 1:nPop
      [Ftmp,Ftmp_all] = fe_schamel(V,n(iPop),vt(iPop),vd(iPop),PHI,VPH,beta); % distribution function, from Schamel 1986
      F = F + Ftmp; 
      Ffree = Ffree + Ftmp_all.free; 
      Ftrap = Ftrap + Ftmp_all.trap;
  end
  nfree = nansum(Ffree,2)*dv;
end
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
function [Fsave,Ffreesave,Ftrapsave] = get_f_abel(V,n,vt,vd,PHI,VPH,nt)  
  units = irf_units;
  phi = PHI(:,1);
  vph = VPH(1,1);
  nx = size(PHI,1);
  v = V(1,:);
  dv = v(2)-v(1);
  f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2) + ...
            n(2)*(1/pi./vt(2).^2)^(1/2)*exp(-(v-vd(2)).^2./vt(2).^2);
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
  
  [fitresult, gof, fun_net, fun_net_prime] = createFit(phi, nt-ntrap_flat);
  
	if 0 % plot
    figure(22)
    plot(phi,ntt-ntrap_flat,phi,fun_net(phi))
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

    fun_int = @(u,V) 0.51*(1/units.e)*(2*units.me)^0.5/pi*fun_net_prime(V/units.e).*(-(units.me*(u-vph).^2/2-units.e*phi_tmp)-V).^-0.5;
    %fun_int = @(u,V) (1/units.e)*(2*units.me)^0.5/pi*fun_net_prime(V/units.e,a,b,c,d).*(-(units.me*u.^2/2-units.e*phi_tmp)-V).^-0.5;    
    for iu_ = 1:nu
      u_tmp = v(itrap_v(iu_));
      nV = 500;    
      Vmax = units.e*phi_tmp-units.me*(u_tmp-vph)^2/2;    
      Vmin = 0.001*Vmax;
      V_tmp = linspace(Vmin,Vmax*0.999,nV);
      
      INT = fun_int(u_tmp,V_tmp);
      f_tmp = trapz(V_tmp,INT);
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
  F(itrap) = Ftrap(itrap) + Ftrap_flat(itrap);
  
  Fsave = F;
  Ffreesave = Ffree;
  Ftrapsave = Ftrap + Ftrap_flat;% Ftrapsave(itrap) = Ftrapsave(itrap) + Ftrap_flat(itrap);
end
function [Fsave,Ffreesave,Ftrapsave,beta] = get_f_schamel(V,n,vt,vd,PHI,VPH,nt)
  %beta = -0.8;
  units = irf_units;
  E = units.me*(V-VPH).^2/2 - units.e*PHI;
  dv = V(1,2) - V(1,1);
  beta_vec = linspace(-5,0,20);
  nbeta = numel(beta_vec);
  C = zeros(nbeta,1);
  nt_diff = Inf;
  loop_iter = 0;
  %ibeta_best = 1;
  while loop_iter < 3
    loop_iter = loop_iter + 1;
    for ibeta = 1:nbeta
      F = zeros(size(V));
      Ffree = zeros(size(V));
      Ftrap = zeros(size(V));
      nPop = numel(n);
      for iPop = 1:nPop
        [Ftmp,Ftmp_all] = fe_schamel(V,n(iPop),vt(iPop),vd(iPop),PHI,VPH,beta_vec(ibeta)); % distribution function, from Schamel 1986
        F = F + Ftmp; 
        Ffree = Ffree + Ftmp_all.free; 
        Ftrap = Ftrap + Ftmp_all.trap;
      end
      nt_tmp = nansum(Ftrap,2)*dv;
      nt_diff_tmp = sum(abs(nt_tmp-nt));
      if 0 % plot
        hca = subplot(2,1,1);
        imagesc(hca,F')
        hold(hca,'on')
        contour(hca,E)
        hold(hca,'off')
        hcb = colorbar('peer',hca);
        hca = subplot(2,1,2);
        plot(hca,1:200,nt,1:200,nt_tmp)
        1;
      end
      if nt_diff_tmp < nt_diff
        beta = beta_vec(ibeta);
        nt_diff = nt_diff_tmp;
        Fsave = F;
        Ffreesave = Ffree;
        Ftrapsave = Ftrap;
        ibeta_best = ibeta;
      end      
    end
    if ibeta_best == 1
      beta_vec = linspace(beta_vec(1)-1,beta_vec(2),20);
    elseif ibeta_best == nbeta
      beta_vec = linspace(beta_vec(nbeta-1)-1,-0.001,20);
    else
      beta_vec = linspace(beta_vec(ibeta_best-1),beta_vec(ibeta_best+1),20);
    end    
    nbeta = numel(beta_vec);
  end
end
function out = fun_net_prime_(phi,nt) 
  %out = 0.5147*4825*phi.^(0.5147-1) + 1.075*-295*phi.^(1.075-1);
  out = b*a*phi.^(b-1) + d*c*phi.^(d-1);
  %out(out<0) = 0;
  out = out*1;
  
	[fitresult, gof, fun_out, fun_deriv_out] = createFit(phi, nt);
  out = fun_deriv_out(phi);
end
function out = fun_net_(phi,nt) 
	[fitresult, gof, fun_out, fun_deriv_out] = createFit(phi, nt);
  out = fun_out(phi);
end
function [fitresult, gof] = createFit_power(phi, nt)
  %CREATEFIT(PHI,NT)
  %  Create a fit.
  %
  %  Data for 'untitled fit 1' fit:
  %      X Input : phi
  %      Y Output: nt
  %  Output:
  %      fitresult : a fit object representing the fit.
  %      gof : structure with goodness-of fit info.
  %
  %  See also FIT, CFIT, SFIT.

  %  Auto-generated by MATLAB on 16-Aug-2018 23:12:04

  %% Fit: 'untitled fit 1'.
  [xData, yData] = prepareCurveData( phi, nt );

  % Set up fittype and options.
  ft = fittype( 'a*x^b+c*x^d', 'independent', 'x', 'dependent', 'y' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'Off';
  opts.StartPoint = [-300 1 400 0.8];

  % Fit model to data.
  [fitresult, gof] = fit( xData, yData, ft, opts );

  if 0 % Plot fit with data.
    figure( 'Name', 'untitled fit 1' );
    h = plot( fitresult, xData, yData );
    legend( h, 'nt vs. phi', 'untitled fit 1', 'Location', 'NorthEast' );
    % Label axes
    xlabel phi
    ylabel nt
    grid on
  end
end
function [fitresult, gof, fun_out, fun_deriv_out] = createFit(phi, dntrap)
%CREATEFIT(PHI,DNTRAP)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : phi
%      Y Output: dntrap
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 17-Aug-2018 17:40:05


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( phi, dntrap );

% Set up fittype and options.
ft = fittype( 'poly8' );

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft );

fun_eval_str = sprintf(...
  'fun_out = @(x) %g*x.^8 + %g*x.^7 + %g*x.^6 + %g*x.^5 + %g*x.^4 + %g*x.^3 + %g*x.^2 + %g*x + %g;',...
  fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7,fitresult.p8,fitresult.p9);
eval(fun_eval_str)

fun_eval_str_deriv = sprintf(...
  'fun_deriv_out = @(x) 8*%g*x.^(8-1) + 7*%g*x.^(7-1) + 6*%g*x.^(6-1) + 5*%g*x.^(5-1) + 4*%g*x.^(4-1) + 3*%g*x.^(3-1) + 2*%g*x.^(2-1) + 1*%g*x.^(1-1) + 0*%g;',...
  fitresult.p1,fitresult.p2,fitresult.p3,fitresult.p4,fitresult.p5,fitresult.p6,fitresult.p7,fitresult.p8,fitresult.p9);
eval(fun_eval_str_deriv)

% Plot fit with data.
% figure( 'Name', 'untitled fit 1' );
% h = plot( fitresult, xData, yData );
% legend( h, 'dntrap vs. phi', 'untitled fit 1', 'Location', 'NorthEast' );
% % Label axes
% xlabel phi
% ylabel dntrap
% grid on
end


