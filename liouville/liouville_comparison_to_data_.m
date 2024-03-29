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

  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('[ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)

  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
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
n = 0.035*1e6;
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
phimax = 5000; % 340
vph = -1*17000.0001e3;

%x_vec = linspace(-3*lx,4*lx,nx);
%t_vec = linspace(-3*lt,4*lt,nt);
t_vec = eDist.time-eDist.time(1);
nt = numel(t_vec);
v_max = 100e6;
nv = 2000;
v_vec = linspace(-v_max,v_max,nv);
dv = v_vec(2)-v_vec(1);
%[X,V] = meshgrid(x_vec,v_vec); X = permute(X,[2 1]); V = permute(V,[2 1]);
[T,V] = meshgrid(t_vec,v_vec); T = permute(T,[2 1]); V = permute(V,[2 1]);

%% Set up potential
%phi = @(x,lx) phimax*(1+tanh(x/lx));
phi = @(t,lt) phimax*0.5*(1+tanh((t-t0)/lt));
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
F_(V<vph) = 0; % make it dark blue instead of white (NaN)
n_ = nansum(F_,2)*dv;
nfree_mod = nansum(Fflat_free,2)*dv;
[n_lb,n_sep] = paper_electron_acceleration.liouville_mapped_nf(n0*1e-6,T0,0,phimax,-vph*1e-3);
n_(end)*1e-6;

%% plot results
tint_model = irf.tint('2017-07-06T08:16:37.00Z/2017-07-06T08:16:39.50Z');
ts_phi_mod = irf.ts_scalar(eDist.time,phi(t_vec,lt));
ts_v_phi_mod = irf.ts_scalar(eDist.time,sqrt(2*units.e*phi(t_vec,lt)./units.me));
ts_n_mod = irf.ts_scalar(eDist.time,n_*1e-6);
ts_f_mod = PDist(eDist.time,flipdim(F_,2),'line (reduced)',v_vec*1e-3);
%ts_V0 = irf.ts_scalar(eDist.time,V0(:,1:20:end)'*1e-6);
v_start = (-15:5:20)*1e3; v_start(v_start>-vph*1e-3) = [];
ts_v0 = ts_get_v0(ts_phi_mod,-vph*1e-3,v_start);
ts_v00 = ts_get_v0(ts_phi_mod,-vph*1e-3,0);
ts_v_sep = ts_get_v0(ts_phi_mod,-vph*1e-3,-0.99*vph*1e-3);


nlim = [0.001 0.039];
vlim = 39.9*[-1 1];
flim = [-6 -2.5];

nrows = 4;
ncols = 1;
%h = setup_subplots(nrows,ncols);
h = irf_plot(nrows);
isub = 1;


if 0
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,ts_v0*1e-3,'k')
  %hca.YLim = [-40 40];
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
  hve = irf_plot(hca,gseVe1par*0,'k:');
  hve = irf_plot(hca,gseVe1par*1e-3,'k');
  hve.LineWidth = 2;
  irf_plot(hca,ts_v0*1e-3,'k')
  irf_plot(hca,ts_v_sep*1e-3,'color',[0 0 0],'linewidth',1.5)
  hold(hca,'off')
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
end
if 1 % mapped distribution, ts
  hca = h(isub); isub = isub + 1;
  %fred_to_plot = ef1D_orig.tlim(tint_model);
  fred_to_plot = ts_f_mod;
  [~, hcb] = irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));
  hcb.YLabel.String = {'PSD','s/m^4'};
  irf_zoom(hca,'x',tint_model)
  colormap(hca,'jet')
  hca.YLim = vlim;
  hca.CLim = flim;
  hold(hca,'on')
  irf_plot(hca,ts_v0*1e-3,'k')
  hold(hca,'off')
  hleg = irf_legend(hca,{['\psi = ' sprintf('%.0f V',phimax)];['v_{\psi}^{lb} = ' sprintf('%.0f km/s',round(vph/1000))]},[0.50 0.92],'w');
  hleg(1).Color = [1 1 1];
  hleg(2).Color = [1 1 1];
  hleg(1).HorizontalAlignment = 'left';
  hleg(2).HorizontalAlignment = 'left';
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
end
if 0 % f_obs-f_map (beam, same clim)
  hca = h(isub); isub = isub + 1;
  
  fred_obs = ef1D_nobg.tlim(tint_model);
  fred_mod = ts_f_mod.tlim(tint_model);
  
  % need to interpolate data to same velocities
  % could also reduce it using different vg
  [V_obs,T_obs] = meshgrid(fred_obs.depend{1}(1,:),1:fred_obs.length);
  [V_mod,T_mod] = meshgrid(fred_mod.depend{1}(1,:),1:fred_mod.length);
  data_fred_mod = fred_mod.data;  
  data_fred_mod_interp = interp2(V_mod,T_mod,data_fred_mod,V_obs,T_obs);
  
  fred_to_plot = fred_obs;
  fred_to_plot.data = -1*fred_to_plot.data + data_fred_mod_interp;
  fred_to_plot_log = fred_to_plot; fred_to_plot_log.data = log10(abs(fred_to_plot.data)).*sign(fred_to_plot.data);
  
  %%
  [~, hcb] = irf_spectrogram(hca,fred_to_plot_log.specrec('velocity_1D','10^3 km/s'),'lin');
  hcb.YLabel.String = {'PSD','s/m^4'};
  irf_zoom(hca,'x',tint_model)
  colormap(hca,'jet')
  hca.YLim = vlim;
  %hca.CLim = max(abs(flim))*[-1 1];
  %hca.CLim = 3e-3*[-1 1];
  %hca.CLim = 6*[-1 1];
  hca.CLim = flim;
  hold(hca,'on')
  irf_plot(hca,ts_v0*1e-3,'k')
  hold(hca,'off')
  hleg = irf_legend(hca,{['\psi = ' sprintf('%.0f V',phimax)];['v_{\psi}^{lb} = ' sprintf('%.0f km/s',round(vph/1000))]},[0.50 0.92],'w');
  hleg(1).Color = [1 1 1];
  hleg(2).Color = [1 1 1];
  hleg(1).HorizontalAlignment = 'left';
  hleg(2).HorizontalAlignment = 'left';
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
  %colormap(hca,pic_colors('blue_red'))
end
if 1 % f_obs-f_map ()
  hca = h(isub); isub = isub + 1;
  
  fred_obs = ef1D_nobg.tlim(tint_model);
  fred_mod = ts_f_mod.tlim(tint_model);
  
  % need to interpolate data to same velocities
  % could also reduce it using different vg
  [V_obs,T_obs] = meshgrid(fred_obs.depend{1}(1,:),1:fred_obs.length);
  [V_mod,T_mod] = meshgrid(fred_mod.depend{1}(1,:),1:fred_mod.length);
  data_fred_mod = fred_mod.data;  
  data_fred_mod_interp = interp2(V_mod,T_mod,data_fred_mod,V_obs,T_obs);
  
  fred_to_plot = fred_obs;
  fred_to_plot.data = fred_to_plot.data - data_fred_mod_interp;  
  
  
  [~, hcb] = irf_spectrogram(hca,fred_to_plot.specrec('velocity_1D','10^3 km/s'));
  hcb.YLabel.String = {'PSD','s/m^4'};
  irf_zoom(hca,'x',tint_model)
  colormap(hca,'jet')
  hca.YLim = vlim;  
  hca.CLim = 3e-3*[-1 1];
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hold(hca,'on')
  irf_plot(hca,ts_v0*1e-3,'k')
  irf_plot(hca,ts_v_sep*1e-3,'color',[0 0 0],'linewidth',1.5)
  
  hold(hca,'off')
  hleg = irf_legend(hca,{['\psi = ' sprintf('%.0f V',phimax)];['v_{\psi}^{lb} = ' sprintf('%.0f km/s',round(vph/1000))]},[0.50 0.92],'w');
  hleg(1).Color = [1 1 1];
  hleg(2).Color = [1 1 1];
  hleg(1).HorizontalAlignment = 'left';
  hleg(2).HorizontalAlignment = 'left';
  hca.YLabel.String = 'v_{||} (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
  colormap(hca,pic_colors('blue_red'))
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
if 1 % densities, including partial densities
  hca = h(isub); isub = isub + 1;
  
  fred_obs = ef1D_nobg.tlim(tint_model);
  fred_mod = ts_f_mod.tlim(tint_model);
  
  % need to interpolate data to same velocities
  % could also reduce it using different vg
  [V_obs,T_obs] = meshgrid(fred_obs.depend{1}(1,:),1:fred_obs.length);
  [V_mod,T_mod] = meshgrid(fred_mod.depend{1}(1,:),1:fred_mod.length);
  data_fred_mod = fred_mod.data;  
  data_fred_mod_interp = interp2(V_mod,T_mod,data_fred_mod,V_obs,T_obs);
  
  fred_to_plot = fred_obs;
  fred_to_plot.data = fred_to_plot.data - data_fred_mod_interp;  
  
  [V_dep,V_sep] = meshgrid(fred_obs.depend{1}(1,:),ts_v_sep.tlim(tint_model).data);
  fred_passing = fred_obs;  
  fred_trapped = fred_obs;
  
  fred_passing.data(V_dep>V_sep) = NaN;
  fred_trapped.data(V_dep<V_sep) = NaN;
  dv = fred_trapped.depend{1}(1,2) - fred_trapped.depend{1}(1,1);
  n_passing = irf.ts_scalar(fred_passing.time,nansum(fred_passing.data,2));
  n_trapped = irf.ts_scalar(fred_trapped.time,nansum(fred_trapped.data,2));
  
  
  if 0
    [~, hcb] = irf_spectrogram(hca,fred_trapped.specrec('velocity_1D','10^3 km/s'));
  else
    colors = pic_colors('matlab');
    set(hca,'ColorOrder',colors);
    irf_plot(hca,{ne1.tlim(tint_model),ts_n_mod,n_passing,n_trapped,n_passing+n_trapped},'comp')
    irf_zoom(hca,'x',tint_model)
    hca.YLim = nlim;
    hca.YLabel.String = 'n_e (cm^{-3})';
    irf_legend(hca,{'FPI';'mapped';'n_{map}';'n_{excess}';'n_{map}+n_{excess}'},[0.01 0.1])
  end
end
if 0 % densities
  hca = h(isub); isub = isub + 1;
  colors = pic_colors('matlab');
  set(hca,'ColorOrder',colors(1:2,:));
  irf_plot(hca,{ne1.tlim(tint_model),ts_n_mod},'comp')
  irf_zoom(hca,'x',tint_model)
  hca.YLim = nlim;
  hca.YLabel.String = 'n_e (cm^{-3})';
  irf_legend(hca,{'FPI';'mapped'},[0.49 0.98])
end
if 0
  hca = h(isub); isub = isub + 1;
end
if 0
  hca = h(isub); isub = isub + 1;
end
irf_plot_axis_align(h)
irf_zoom(h,'x',tint_model)

irf_legend(h(1),'(d)',[0.02 0.85],'color',[0 0 0],'fontsize',16)
irf_legend(h(2),'(e)',[0.02 0.97],'color',[1 1 1],'fontsize',16)
irf_legend(h(3),'(f)',[0.02 0.6],'color',[0 0 0],'fontsize',16)
