%% load data
doLoadData = 0;
if doLoadData
ic = 1;
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
units = irf_units;

event = 1;
sep.get_tints;

mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)

c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
tint_fred = irf.tint('2017-07-06T08:16:37.00Z/2017-07-06T08:16:41.00Z');
c_eval('eDist = ePDist?.tlim(tint_fred);',ic)

% Remove background
nSecondary = 5;
nPhoto = 1;
tic; eDist_nobg = mms.remove_edist_background(eDist,'nSecondary',nSecondary,'Nphotoe_art',nPhoto); toc;

% Reduced electron distribution
eint = [00 40000];
lowerelim = 000;
vint = [-Inf Inf];
%tint_fred = tint_fred;%tint_phi;


eDist_orig = eDist;

c_eval('scpot = scPot?.resample(eDist);',ic)
c_eval('dir_red = dmpaB?.resample(eDist).norm;',ic)
energies = sqrt(2*eDist.depend{1}(1,:)*units.eV/units.me)/1000; % km/s
vgmax = 70000;
vg = -vgmax:1000:vgmax;
vg(abs(vg)>70000) = [];

nMC = 200;
tic; ef1D_orig = eDist_orig.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
tic; ef1D_nobg = eDist_nobg.reduce('1D',dir_red,'scpot',scpot,'lowerelim',lowerelim,'nMC',nMC,'vg',vg); toc % reduced distribution along B
end

%% f0
units = irf_units;
%n = [0.02 0.02]*1e6; % m-3
n = 0.04*1e6;
T = 350; % eV
vd = 0000e3; % ms-1 oa
vt = sqrt(2*units.e*T./units.me); % m/s

str_info = {'unperturbed f:';...
            ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
            ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
            ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
            };
          
%% Gaussian potential
%lx = 1000; 
lt = 0.4; % s
t0 = 1.2; % s
%nx = 100;
%nt = 100;
phimax = 800; % 340
vph = -15000e3;

%x_vec = linspace(-3*lx,4*lx,nx);
%t_vec = linspace(-3*lt,4*lt,nt);
t_vec = eDist.time-eDist.time(1);
nt = numel(t_vec);
v_max = 60e6;
nv = 2000;
v_vec = linspace(-v_max,v_max,nv);
dv = v_vec(2)-v_vec(1);
%[X,V] = meshgrid(x_vec,v_vec); X = permute(X,[2 1]); V = permute(V,[2 1]);
[T,V] = meshgrid(t_vec,v_vec); T = permute(T,[2 1]); V = permute(V,[2 1]);

%% Set up potential
%phi = @(x,lx) phimax*(1+tanh(x/lx));
phi = @(t,lt) phimax*(1+tanh((t-t0)/lt));
U = units.me*(V-vph).^2/2 - units.e*double(phi(T,lt));
%U = units.me*(V-vph).^2/2 - units.e*double(phi(X,lx));
%PHI = repmat(tocolumn(double(phi(x_vec,lx))),1,nv);
PHI = repmat(tocolumn(double(phi(t_vec,lt))),1,nv);
f0 = @(v) n(1)./pi^0.5./vt(1).*exp(-(v-vd(1)).^2./vt(1).^2);

%v0 = @(x,v) vph + sign(v-vph).*((v-vph).^2-2*units.e*phi(x)./units.me).^0.5;
v0 = @(t,v) vph + sign(v-vph).*((v-vph).^2-2*units.e*phi(t)./units.me).^0.5;
%v_sep_bot(x) = vph - (2*units.e*phi(x)/units.me)^0.5;
%v_sep_top(x) = vph + (2*units.e*phi(x)/units.me)^0.5;
%ff(x,v) = f0(v0(x,v));

%v_sep = v0(x,vph);
%f_sep = double(f0(vph));
%FF = double(ff(X,V));
%V0 = double(v0(X,V));
%FF(U<0) = 0;

%% Peform liouville mapping
[Fflat,Fflat_free,Fflat_trap,V0] = get_f_flat(V,n,vt,vd,1*PHI,PHI*0+vph);
F_ = Fflat_free;
F_(V<vph) = NaN;
n_ = nansum(F_,2)*dv;
nfree_mod = nansum(Fflat_free,2)*dv;

%% plot results
tint_model = irf.tint('2017-07-06T08:16:37.00Z/2017-07-06T08:16:39.50Z');
ts_phi_mod = irf.ts_scalar(eDist.time,phi(t_vec,lt));
ts_v_phi_mod = irf.ts_scalar(eDist.time,sqrt(2*units.e*phi(t_vec,lt)./units.me));
ts_n_mod = irf.ts_scalar(eDist.time,n_*1e-6);
ts_f_mod = PDist(eDist.time,flipdim(F_,2),'line (reduced)',v_vec*1e-3);
ts_V0 = irf.ts_scalar(eDist.time,V0(:,1:20:end)'*1e-6);

nlim = [0.01 0.04];
vlim = 40*[-1 1];
flim = [-6 -2.5];

nrows = 3;
ncols = 1;
%h = setup_subplots(nrows,ncols);
h = irf_plot(nrows);
isub = 1;


if 0
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,ts_V0,'k')
  hca.YLim = [-40 40];
end
if 0
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,-1*ts_phi_mod,'k')
end
if 0
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,-1*ts_v_phi_mod*1e-6,'k')
end
if 1 % observed distribution
  hca = h(isub); isub = isub + 1;
  %fred_to_plot = ef1D_orig.tlim(tint_model);
  fred_to_plot = ef1D_nobg.tlim(tint_model);
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));
  irf_zoom(hca,'x',tint_model)
  colormap(hca,'jet')
  hca.YLim = vlim;
  hca.CLim = flim;
  hold(hca,'on')
  irf_plot(hca,-1*ts_v_phi_mod*1e-6,'k')
  hold(hca,'off')
end
if 1 % mapped distribution, ts
  hca = h(isub); isub = isub + 1;
  %fred_to_plot = ef1D_orig.tlim(tint_model);
  fred_to_plot = ts_f_mod;
  irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));
  irf_zoom(hca,'x',tint_model)
  colormap(hca,'jet')
  hca.YLim = vlim;
  hca.CLim = flim;
  hold(hca,'on')
  irf_plot(hca,-1*ts_v_phi_mod*1e-6,'k')
  hold(hca,'off')
end
if 0 % liouville mapped distribution
  hca = h(isub); isub = isub + 1;
  variable = F_;
  variable(U<0) = NaN;
  pcolor(hca,t_vec,-v_vec*1e-6,log10(variable'))
  shading(hca,'flat')
  colormap(hca,'jet')
  hb = colorbar('peer',hca);
  %hca.CLim = [-6 -2.2];
  hca.XLabel.String = 'distance (km)';
  hca.YLabel.String = 'v (10^3 km/s)';
  hca.YLim = vlim;
  hca.CLim = flim;
  hold(hca,'on')
  irf_plot(hca,-1*ts_v_phi_mod*1e-6,'k')
  hold(hca,'off')
end
if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,t_vec,n_*1e-6)
  hca.YLim = nlim;
end
if 1
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{ne1.tlim(tint_model),ts_n_mod},'comp')
  irf_zoom(hca,'x',tint_model)
  hca.YLim = nlim;
end
if 0
  hca = h(isub); isub = isub + 1;
end
if 0
  hca = h(isub); isub = isub + 1;
end
irf_plot_axis_align(h)
irf_zoom(h,'x',tint_model)

%%
if 0
%%
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
end
