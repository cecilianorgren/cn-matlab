%% Set up grid and phi
fun_phi = @(phimax,x,lx) phimax*exp(-x.^2/2/lx.^2);

lx = 500; 
nx = 200;
x = linspace(-5*lx,5*lx,nx);
phimax = 10;
phi = fun_phi(phimax,x,lx);

units = irf_units;
vtrap = sqrt(2*units.e*phimax/units.me);
nv = 1000;
v = linspace(-4*vtrap,4*vtrap,nv);
dv = v(2) - v(1);

[X,V] = meshgrid(x,v); X = permute(X,[2 1]); V = permute(V,[2 1]);
PHI = fun_phi(phimax,X,lx);
E = units.me*V.^2/2 - units.e*PHI;

%% Get electric field and n
dx =  x(2) - x(1);
x_diff1 = x(1:end-1) + 0.5*dx;
x_diff2 = x(2:end-1) + dx;
Efield = -diff(phi)/dx;
dn = diff(phi,2)/dx/dx*units.eps0/units.e;
dn = [dn(1); dn'; dn(end)]; % assume phi -> at edges

%% Get free population
n = [0.02 0.02]*1e6; % m-3
T = [50 50]; % eV
vd = [4000e3 0000e3]; % ms-1 
vt = sqrt(2*units.e*T./units.me); % m/s
n0 = sum(n);

f0 = @(v) n(1)*(1/pi./vt(1).^2)^(1/2)*exp(-(v-vd(1)).^2./vt(1).^2) + ...
          n(2)*(1/pi./vt(2).^2)^(1/2)*exp(-(v-vd(2)).^2./vt(2).^2);

all_ind = 1:numel(V);
ifree = find(E > 0);
itrap = setdiff(all_ind,ifree);
iabove = find(V>0);
ibelow = find(V<0);

fall = nan(size(v));
ffree = nan(size(v));
ftrap = nan(size(v));

V0 = V*0;

V0(iabove) = 0 + ((V(iabove)).^2 - 2*units.e*(PHI(iabove))/units.me).^0.5; % same as Schamel free streaming
V0(ibelow) = 0 - ((V(ibelow)).^2 - 2*units.e*(PHI(ibelow))/units.me).^0.5;
V0(itrap) = NaN;

Ffree = V0*0;
Ftrap = V0*NaN;

Ffree(ifree) = f0(V0(ifree));

nf = nansum(Ffree,2)*dv;
nt = n0 - nf + dn;

% Get trapped distribution
dntdphi = [(nt(2)-nt(1))/(phi(2)-phi(1)); diff(nt)./diff(phi');] ; % noisy, make functional fit instead
dntdphi = [(nt(2)-nt(1))/(phi(2)-phi(1)); (nt(3:end)-nt(1:end-2))./(phi(3:end)'-phi(1:end-2)'); (nt(end)-nt(end-1))/(phi(end)-phi(end-1))]; % noisy, make functional fit instead
dntdphi(phi<0.01*phimax) = NaN;
abcd = 0;
switch abcd
  case 0 % automated
    [fitreslts,gif] = createFit(torow(phi),torow(nt));
    a = fitreslts.a;
    b = fitreslts.b;
    c = fitreslts.c;
    d = fitreslts.d;
  case 1 % original fit
    a = 4825; 
    b = 0.5147;
    c = -295; 
    d = 1.075;
  case 2
    a = 2800;
    b = 0.53;
    c = -140;
    d = 1.00;
  case 3
    a = 9000;
    b = 0.5147;
    c = -220;
    d = 1.11;
  case 4 % 
    a = 1500; 
    b = 0.50;
    c = -320; 
    d = 1.07;
  case 5 % phimax = 100; lx = 1000;
    a = -239; 
    b = 1.085;
    c = 4797; 
    d = 0.5204;
  case 6 % dphi/dt = 1;
    a = 2000; 
    b = 0.5;
    c = 0; 
    d = 0;
  case 7 % original fit mult by 0.32
    a = 0.32*4825; 
    b = 0.5147;
    c = -0.32*295; 
    d = 1.075;
  case 8 % phimax = 100; lx = 1000; mult. by 0.32.    
    a = -239*0.37; 
    b = 1.085;
    c = 4797*0.31; 
    d = 0.5204;
end
%fun_net = @(phi) 4825*phi.^0.5147 + -295*phi.^1.075;
%fun_net = @(phi) 4825*phi.^0.5147 + -295*phi.^1.075;
%fun_net = @(phi) 4825*phi.^0.5147 + -350*phi.^1.02;
%fun_net_prime = @(phi) 0.5147*4825*phi.^(0.5147-1) + 1.075*-295*phi.^(1.075-1);
%figure(11);plot(phi,nt,phi,fun_net(phi))

for ix = 1:nx
  itrap_ngtvE = find(E(ix,:) < 0);
  itrap_grt0 = find(v>0);
  itrap_v = intersect(itrap_ngtvE,itrap_grt0);
  
  Etrap = E(ix,itrap_v); 
  nu = numel(itrap_v);
  
  if 0 % fit, changes with phi and lx
    %%
    net_prime = fun_net_prime(phi(ix),a,b,c,d); 
    net_prime = fun_net_prime(V_tmp/units.e,a,b,c,d); 
  else % approximate derivative by finite differences
    if ix == 1
      net_prime = (nt(2) - nt(1))/(phi(2)-phi(1));
    elseif ix == nx
      net_prime = (nt(end) - nt(end-1))/(phi(end)-phi(end-1));
    else
      net_prime = (nt(ix+1) - 0*2*nt(ix) - nt(ix-1))/(phi(ix+1) - 0*2*phi(ix)- phi(ix-1));
      net_prime = (nt(ix+1) + nt(ix-1))/(phi(ix+1) + phi(ix-1));
      if net_prime == Inf || net_prime == NaN
        net_prime = 0;
      end
    end
  end
  phi_tmp = phi(ix);
  %fun_int = @(u,phimax,phi) fun_net_prime(phi)*(units.e*phimax-units.me*u.^2/2-units.e*phi).^-0.5;
  if 0
    fun_int = @(u,V) (1/units.e)*(2*units.me)^0.5/pi*net_prime*(-(units.me*u.^2/2-units.e*phi_tmp)-V).^-0.5;
    fun_int_dV = @(u,V) (1/units.e)*(2*units.me)^0.5/pi*net_prime*(-(units.me*u.^2/2-units.e*phi_tmp)-V).^0.5*(-1);
  else % dnt/dphi(V) included as a function
    fun_int = @(u,V) (1/units.e)*(2*units.me)^0.5/pi*fun_net_prime(V/units.e,a,b,c,d).*(-(units.me*u.^2/2-units.e*phi_tmp)-V).^-0.5;
    fun_int_dV = @(u,V) (1/units.e)*(2*units.me)^0.5/pi*fun_net_prime(V/units.e,a,b,c,d).*(-(units.me*u.^2/2-units.e*phi_tmp)-V).^0.5*(-1);
  end
  for iu_ = 1:nu
    u_tmp = v(itrap_v(iu_));
    nV = 300;    
    Vmax = units.e*phi_tmp-units.me*u_tmp^2/2;    
    Vmin = 0.00001*Vmax;
    V_tmp = linspace(Vmin,Vmax*0.9999,nV);
    if 1 % different net_prime, must be used to include dnt/dphi
      net_prime = fun_net_prime(phi(ix),a,b,c,d); 
      net_prime = fun_net_prime(V_tmp/units.e,a,b,c,d); 
    end
    if 1 % do integration numerically
      %V_tmp = linspace(Vmin,Vmax*0.9999,nV);
      INT = fun_int(u_tmp,V_tmp);
      f_tmp = trapz(V_tmp,INT);
    else
      f_tmp = fun_int_dV(u_tmp,Vmax)-fun_int_dV(u_tmp,Vmin);
    end    
    
    if 0 % reverse index
      ii_pos = itrap_v(1)+itrap_v(nu)-itrap_v(iu_);
      ii_neg = itrap_v(1)-1-(ii_pos-itrap_v(1));
      Ftrap(ix,ii_pos) = f_tmp;
      Ftrap(ix,ii_neg) = f_tmp;
    else
      Ftrap(ix,itrap_v(iu_)) = f_tmp;
      Ftrap(ix,itrap_v(1)-1-(itrap_v(iu_)-itrap_v(1))) = f_tmp;      
    end
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
  if 0 && nu > 0 % magic 1, aligns separatrices
    v0 = 0 + ((V(iabove)).^2 - 2*units.e*PHI(iabove)/units.me).^0.5; % same as Schamel free streaming
    fsep = f0(2*units.e*phi_tmp/units.me);
    Ftrap(ix,:) = Ftrap(ix,:) - Ftrap(ix,itrap_v(nu)) + Ffree(ix,itrap_v(nu)+1);
  end
  if 0 && nu > 0 % magic 2
    v0 = 0 + ((V(iabove)).^2 - 2*units.e*PHI(iabove)/units.me).^0.5; % same as Schamel free streaming
    ftrap_mean = nanmean(Ftrap(ix,:));
    dftrap = ftrap_mean-Ftrap(ix,:);
    %fsep = f0(2*units.e*phi_tmp/units.me);
    Ftrap(ix,:) = ftrap_mean + dftrap;
  end
  if 0%nu > 0 % magic 1
    v0 = 0 + ((V(iabove)).^2 - 2*units.e*PHI(iabove)/units.me).^0.5; % same as Schamel free streaming
    fsep = f0(2*units.e*phi_tmp/units.me);
    Ftrap(ix,:) = Ftrap(ix,:) - Ftrap(ix,itrap_v(nu)) + Ffree(ix,itrap_v(nu)+1);
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

% Add free and trapped pops together
F = V0*NaN;
F(itrap) = Ftrap(itrap);
F(ifree) = Ffree(ifree);

% Schamels for comparison
beta = -0.8;
nPop = numel(n);
F_sch = zeros(size(V));
Ffree_sch = zeros(size(V));
Ftrap_sch = zeros(size(V));
for iPop = 1:nPop
[Ftmp,Ftmp_all] = fe_schamel(V,n(iPop),vt(iPop),vd(iPop),PHI,V*0,beta); % distribution function, from Schamel 1986
F_sch = F_sch + Ftmp; 
Ffree_sch = Ffree_sch + Ftmp_all.free; 
Ftrap_sch = Ftrap_sch + Ftmp_all.trap;
end

figure(92)
nrows = 3;
ncols = 2;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;

if 0 % dnet/dphi
  hca = h(isub); isub = isub + 1;
  net_prime_plot = fun_net_prime(phi,a,b,c,d);
  plot(hca,x,net_prime_plot)
  hca.YLabel.String = 'dn_{et}/d\phi';
end
if 1 % dn(phi)
  hca = h(isub); isub = isub + 1;
  hp = plot(hca,phi,nt*1e-6,phi,fun_net(phi,a,b,c,d)*1e-6);  
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = '\delta n (cc)';
  legend(hca,{'n_{et}','fit to n_{et}'},'Box','off','location','southeast')  
end
if 0 % Ffree
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
if 0 % Ftrap
  hca = h(isub); isub = isub + 1;
  levels_E = linspace(min(E(:)),max(E(:)),20); levels_E = levels_E + min(abs(levels_E));
  %levels_E = [0 linspace(min(E(:)),max(E(:)),15)];
  pcolor(hca,X,V,Ftrap)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{trap} (...)';
  %hca.CLim = max(abs(hca.CLim))*[0 1];
  colormap(hca,cn.cmap('white_blue'))  
  %colormap(hca,'jet')  
end
if 1 % F
  hca = h(isub); isub = isub + 1;
  levels_E = linspace(min(E(:)),max(E(:)),20); levels_E = levels_E + min(abs(levels_E));
  %levels_E = [0 linspace(min(E(:)),max(E(:)),15)];
  pcolor(hca,X,V,F)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{} (.all..)';
  %hca.CLim = max(abs(hca.CLim))*[0 1];
  colormap(hca,cn.cmap('white_blue'))  
  %colormap(hca,'jet')  
end
if 1 % F schamel
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V,F_sch)  
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'v';
  hcb.YLabel.String = 'f_{} (.all..)';
  %hca.CLim = max(abs(hca.CLim))*[0 1];
  colormap(hca,cn.cmap('white_blue'))  
  %colormap(hca,'jet')  
end
if 0 % add up densitites
  hca = h(isub); isub = isub + 1;
  plotyy(hca,x,nansum(Ftrap,2),x,dn)
  hca.YLabel.String = 'sum(f_{trap})'; 
end
if 1 % compare trapped densities
  hca = h(isub); isub = isub + 1;
  plot(hca,x,nt,x,nansum(Ftrap,2)*dv,x,nansum(Ftrap_sch,2)*dv)
  hca.YLabel.String = 'n';
  legend(hca,{'n_f','sum(f_t)','sum(f_{t,sch})'},'Box','off')
end
if 1 % add up densitites
  hca = h(isub); isub = isub + 1;
  plot(hca,x,nf,x,nansum(Ftrap,2)*dv,x,sum(F,2)*dv,x,nt,x,nansum(Ftrap_sch,2)*dv,x,sum(F_sch,2)*dv)
  hca.YLabel.String = 'n';
  legend(hca,{'n_f','sum(f_t)','sum(f_f+f_t)','n_t','n_{t,sch}','sum(f_f+f_t)_{sch}'},'Box','off')
end
if 1 % F at a few x's
  hca = h(isub); isub = isub + 1;  
  plot(hca,v,F(fix(nx/6),:),v,F(fix(nx/4),:),v,F(fix(nx/2),:),v,F(fix(nx/2)+5,:))
  %hca.YScale = 'log';
  hca.YLabel.String = 'F';
  hca.XLabel.String = 'v';
end
1;

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
function out = fun_net_prime(phi,a,b,c,d) 
  %out = 0.5147*4825*phi.^(0.5147-1) + 1.075*-295*phi.^(1.075-1);
  out = b*a*phi.^(b-1) + d*c*phi.^(d-1);
  %out(out<0) = 0;
  out = out*1;
end
function out = fun_net(phi,a,b,c,d) 
  out = a*phi.^b + c*phi.^d;  
end
function [fitresult, gof] = createFit(phi, nt)
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
  opts.StartPoint = [0.883105865634136 0.826556034240979 0.612634266868576 0.964888535199277];

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

