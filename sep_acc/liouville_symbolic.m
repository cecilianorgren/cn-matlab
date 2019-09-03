%% f0
units = irf_units;
%n = [0.02 0.02]*1e6; % m-3
n = 0.05*1e6;
T = 100; % eV
vd = 0000e3; % ms-1 oa
vt = sqrt(2*units.e*T./units.me); % m/s

str_info = {'unperturbed f:';...
            ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
            ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
            ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
            };
          
%% Gaussian potential
lx = 3100; 
nx = 207;
phimax = 2000; % 340
vph = -5000e3;

x_vec = linspace(-7*lx,7*lx,nx);
v_max = 40e6;
nv = 400;
v_vec = linspace(-v_max,v_max,nv);
dv = v_vec(2)-v_vec(1);
[X,V] = meshgrid(x_vec,v_vec); X = permute(X,[2 1]); V = permute(V,[2 1]);

%% symbolic expression, 120 s
tic
syms x v %phi(x) f0(v) ff(x,v) v0(x,v)
%phi(x) = phimax*exp(-x.^2/2/lx.^2); 
phi(x) = phimax*(1+tanh(x/lx));
E = -diff(phi,x);
dn = diff(units.eps0/units.e*phi,x,2); % dn = -diff(E,x); % the same, dn = ne-ni

U = units.me*(V-vph).^2/2 - units.e*double(phi(X));
PHI = repmat(tocolumn(double(phi(x_vec))),1,nv);

f0(v) = n(1)/pi^0.5/vt(1)*exp(-(v-vd(1))^2/vt(1)^2);
v0(x,v) = vph + sign(v-vph)*((v-vph)^2-2*units.e*phi(x)/units.me)^0.5;
v_sep_bot(x) = vph - (2*units.e*phi(x)/units.me)^0.5;
v_sep_top(x) = vph + (2*units.e*phi(x)/units.me)^0.5;
ff(x,v) = f0(v0(x,v));

v_sep = v0(x,vph);
f_sep = double(f0(vph));
FF = double(ff(X,V));
V0 = double(v0(X,V));
FF(U<0) = 0;

toc

%% numerical integration: nt = int(ff)dv
%nt = int(ff);
tic
for ix = 1:nx
  x_tmp = x_vec(ix);
  v_top = vph + sqrt(2*units.e*double(phi(x_tmp))/units.me);
  v_bot = vph - sqrt(2*units.e*double(phi(x_tmp))/units.me);  
  ff_tmp(v) = ff(x_tmp,v);
  ff_fun = matlabFunction(ff_tmp);
  %nt_tmp = integral(ff_tmp(v),v);
  nt_flat(ix) = 2*sqrt(2*units.e*double(phi(x_tmp))/units.me)*f_sep;
  nf_top(ix) = integral(ff_fun,v_top,Inf);
  nf_bot(ix) = integral(ff_fun,-Inf,v_bot);    
  if 0 % plot
    %%
    nrows = 2; ncols = 2; isub = 1;
    hca = subplot(nrows,ncols,isub); isub = isub + 1;
    plot(v_vec,ff_fun(v_vec),v_top,ff_fun(v_top),'*',v_bot,ff_fun(v_bot),'*')
  end
end
toc

%nf_bot = real(nf_bot);
%nf_top = real(nf_top);
nf = nf_bot + nf_top;
nt = ni + double(dn(x_vec)) - nf;
dnt = nt-nt_flat;

%%
tic;
[Fflat,Fflat_free,Fflat_trap] = get_f_flat(V,n,vt,vd,PHI,PHI*0+vph);
nfree_mod = nansum(Fflat_free,2)*dv;
ntrap_flat_mod = nansum(Fflat_trap,2)*dv;
ntrap_mod = ntot - torow(nfree_mod) + torow(dn(x_vec));
dntrap_mod = ntrap_mod - torow(ntrap_flat_mod);
[Fabel,Fabel_free,Fabel_trap] = get_f_abel(V,n,vt,vd,PHI,PHI*0+vph,dnt);
Fabel = Fabel + Fflat_trap;              
[Fabel_2,Fabel_free_2,Fabel_trap_2] = get_f_abel(V,n,vt,vd,PHI,PHI*0+vph,dntrap_mod);
Fabel_2 = Fabel_2 + Fflat_trap;

[fitresult, gof, fun_net, fun_net_prime] = createFit(double(phi(x_vec)), dnt);
[fitresult_2, gof_2, fun_net_2, fun_net_prime_2] = createFit(double(phi(x_vec)), dntrap_mod);
fun_net_str_cell = tokenize(char(fun_net_2),'+');
fun_net_prime_str_cell = tokenize(char(fun_net_prime_2),'+');
%fun_net_str = 
%for ipol = 1:numel(fun_net_str_cell)
  
%[Fscha_mod,Fscha_free_mod,Fscha_trap_mod,beta_mod] = get_f_schamel(V_mod,n,vt,vd,PHI_mod,VPH_mod,ntrap_mod,beta_range);
toc
%% Plot results
figure(95)
nrows = 3;
ncols = 4;
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
if 1 % phi(x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec*1e-3,phi(x_vec))
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = '\phi (V)';
end
if 1 % E(x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec*1e-3,E(x_vec)*1e-3)
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'E (mV/m)';
end
if 1 % dn(x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec*1e-3,dn(x_vec)*1e-6)
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = '(\epsilon_0/e)\nabla^2\phi (cm^{-3})';
end
if 0 % f0(v)
  hca = h(isub); isub = isub + 1;
  plot(hca,v_vec*1e-6,f0(v_vec))
  hca.XLabel.String = 'v (10^{3} km/s)';
  hca.YLabel.String = 'f_0 (s^1m^{-4})';
  hca.XLim = [-20 20];
end
if 1 % v0(x,v)
  hca = h(isub); isub = isub + 1;
  contourf(hca,X*1e-3,V*1e-6,real(V0),31);
  shading(hca,'flat');
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'v (10^3 km/s)';
  hca.Title.String = 'v_{0}';
  hca.YLim = [-20 20];
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  colormap(hca,cn.cmap('blue_red'));
end
if 1 % ff(x,v)
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X*1e-3,V*1e-6,FF);
  shading(hca,'flat');
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'v (10^3 km/s)';
  hca.Title.String = 'f_{f}';
  hca.YLim = [-20 20];
end
if 1 % n(x)
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec*1e-3,ni+x_vec*0,x_vec*1e-3,ni+double(dn(x_vec)),x_vec*1e-3,nf,x_vec*1e-3,nt,x_vec*1e-3,(nansum(Fabel,2)*dv))
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'n (m^{-3})';
  irf_legend(hca,{'n_i';'n_i+(\epsilon_0/e)\nabla^2\phi';'n_f';'n_t';'n_{mod}'},[0.98 0.8])  
end
if 1 % dn(x), phi(x)
  hca = h(isub); isub = isub + 1;
  ax = plotyy(hca,x_vec*1e-3,phi(x_vec),x_vec*1e-3,[nt;dnt]);  
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = '\phi (V)';
  ax(2).YLabel.String = 'n (m^{-3})';
  irf_legend(hca,{'\phi';'n_t';'n_t-n_{t,flat}'},[0.98 0.98])
end
if 1 % nt, dnt vs phi
  hca = h(isub); isub = isub + 1;
  plot(hca,phi(x_vec),nt,phi(x_vec),dnt,phi(x_vec),fun_net(phi(x_vec))) % dnt_flat = nt-nt_flat;
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = 'n_t (m^{-3})';
  irf_legend(hca,{'n_t-n_{t,flat}:';fun_net_str_cell},[0.98 0.95],'color',hca.Children(end-1).Color,'fontsize',8)
  irf_legend(hca,{'n_t'},[0.01 0.98])
end
if 1 % dnt-dntfit
  hca = h(isub); isub = isub + 1;
  plot(hca,phi(x_vec),dnt-fun_net(phi(x_vec))) % dnt_flat = nt-nt_flat;
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = 'n_t-n_t^{fit} (m^{-3})';
  %irf_legend(hca,{'n_t-n_{t,flat}:';fun_net_str_cell},[0.98 0.95],'color',hca.Children(end-1).Color,'fontsize',8)  
end
if 1 % nt, dnt vs phi
  hca = h(isub); isub = isub + 1;
  plot(hca,phi(x_vec),fun_net_prime(phi(x_vec))) % dnt_flat = nt-nt_flat;
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = 'd(n_t-n_{t,flat})/d\phi (m^{-3}/V)';
  hleg = irf_legend(hca,{'n_t-n_{t,flat}:';fun_net_prime_str_cell},[0.98 0.85],'color',hca.Children(1).Color,'fontsize',8);
  hleg(2).Interpreter = 'none';
end
if 1 % F obs, Abel
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X*1e-3,V*1e-6,Fabel)
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb.YLabel.String = 'f_{Abel} (s^1m^{-4})';  
  hca.YLim = vlim*[-1 1]*1e-6;
  colormap(hca,cn.cmap('white_blue'))
  hca.Title.String = 'Abel';    
end
if 1 % F obs, Abel
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X*1e-3,V*1e-6,Fabel_2)
  shading(hca,'flat') 
  hcb = colorbar('peer',hca);
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'v (10^3 km/s)';
  hcb.YLabel.String = 'f_{Abel} (s^1m^{-4})';  
  hca.YLim = vlim*[-1 1]*1e-6;
  colormap(hca,cn.cmap('white_blue'))
  hca.Title.String = 'Abel 2';    
end
