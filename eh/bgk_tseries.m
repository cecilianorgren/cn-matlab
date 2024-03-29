%% Set up grid and phi
% diff in flux

vph_all = [-9000*1e3 -10000*1e3];
phi_mult_all = [1.0 1.5];
iff_all = 1:17;

vph_all = [-9000*1e3];
phi_mult_all = [1.0];
iff_all = 19;

n_vph_all = numel(vph_all);
n_phi_mult_all = numel(phi_mult_all);
n_iff_all = numel(iff_all);

flux_all = cell(n_vph_all,n_phi_mult_all,n_iff_all);
psd_all = cell(n_vph_all,n_phi_mult_all,n_iff_all);
n_all = cell(n_vph_all,n_phi_mult_all,n_iff_all);

for i_vph = 1:n_vph_all
for i_phi_mult = 1:n_phi_mult_all
for i_iff = 1:n_iff_all
vph = vph_all(i_vph);
phi_mult = phi_mult_all(i_phi_mult);
iff = iff_all(i_iff);
  
units = irf_units;
mms_id = 1;
tint_all = irf.tint('2017-07-06T13:54:05.52Z/2017-07-06T13:54:05.630Z');
tint = tint_all;

% Can probably change this to tint_all;
scPot = double(mean(scPot1.tlim(tint_phi).data));

% Set velocity grid
vmax = 90000e3;
nv = 2000;
v_vec = linspace(-vmax,vmax,nv);
dv = v_vec(2) - v_vec(1);

% Get EDI data
%mms.load_data_edi;

%% EDI energy and corresponding velocity
E_edi = 500 - 0*scPot; % eV
v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
dE_edi = 25; % eV

E_edi_plus = E_edi + dE_edi;
E_edi_minus = E_edi - dE_edi;
v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
v_edi_plusminus = v_edi_plus-v_edi_minus;
dv_edi_minus = v_edi_minus - v_edi;
dv_edi_plus = v_edi_plus - v_edi;
dv_edi = dv_edi_plus - dv_edi_minus; % m/s
%%
% EDI flux
edi_nodes = 1;
% if 0
%   c_eval('ts_edi_flux180_mms? = irf.ts_scalar(flux180_mms?.tlim(tint).time,mean(flux180_mms?.tlim(tint).data(:,edi_nodes),2));')
%   c_eval('ts_edi_flux0_mms? = irf.ts_scalar(flux0_mms?.tlim(tint).time,mean(flux0_mms?.tlim(tint).data(:,edi_nodes),2));')
%   c_eval('ts_edi_flux180 = ts_edi_flux180_mms?',mms_id)
%   c_eval('ts_edi_flux0 = ts_edi_flux0_mms?',mms_id)
% else  
  c_eval('ts_edi_flux180_mms? = ePitch?_flux_edi.palim([168.75 180]);')
  c_eval('ts_edi_flux0_mms? = ePitch?_flux_edi.palim([0 11.25]);')
  c_eval('ts_edi_flux180 = ts_edi_flux180_mms?;',mms_id)
  c_eval('ts_edi_flux0 = ts_edi_flux0_mms?;',mms_id)
% end

% Indices corresponding to edi energy/velocity interval
vind_edi_0 = intersect(find(v_vec>v_edi_minus),find(v_vec<v_edi_plus));
vind_edi_180 = intersect(find(v_vec<-v_edi_minus),find(v_vec>-v_edi_plus));

% Get properties of individual eh's
data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_vtrap_all = sqrt(2*units.e*obs_potential/units.me);
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(obs_velocity);
c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
c_eval('obs_vph = obs_vph?;',mms_id)

% charge separation assuming potential structure is single gaussian
dn = units.eps0/units.e*obs_potential_max./(obs_lpp*1e3)*1e-6; % cc
dn = dn;

% Get potential from Efield
%vph = -9000*1e3;    
tsVph = irf.ts_scalar(tint_all,vph*ones(tint_all.length,1));
tint_single = tint_all;
tint_utc = tint_single.utc;
c_eval('Epar = gseE?par;',mms_id)
[phi,phi_progressive,phi_ancillary] = get_phi(Epar,vph,tint_all,tint_single);

%phi_add = 50; % bugged at the moment
phi = phi*phi_mult;%+phi_add; 


dntrap_all = cell(numel(obs_lpp),1);
phi_all = cell(numel(obs_lpp),1);
x_all = cell(numel(obs_lpp),1);
fun_fit_all = cell(numel(obs_lpp),1);


% F0
if 1
  %iff = 17;
  iff = 20;
  [f0,params] = mms_20170706_135303.get_f0(iff);
  n = params.n;
  ntot = sum(n);
  R = n(1)/ntot;
  T = params.T;
  vd = params.vd;
  vt = params.vt;
  n0 = sum(n);
  ff4 = f0(v_vec,n,vd,vt);
else
  n = [0.02 0.02]*1e6; % m-3
  ntot = sum(n);
  R = 0.7;
  n = [R (1-R)]*ntot;
  T = [150 4000]; % eV
  vd = [-11000e3 4000e3]; % ms-1 
  vt = sqrt(2*units.e*T./units.me); % m/s
  n0 = sum(n);
end



str_info = {'unperturbed f:';...
            ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
            ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
            ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
            };
               
beta_range = [-2.5 0];
mms_id = 1;

%% For entire time series
  
phi_vec = phi.data;
phi_obs = phi_vec;
x_vec = phi_ancillary.x_vec;
x_obs = x_vec;

% Get electric field and n
dx_obs = x_obs(2) - x_obs(1);
x_obs_diff1 = x_obs(1:end-1) + 0.5*dx_obs;
x_obs_diff2 = x_obs(2:end-1) + dx_obs;
EfieldFromPhi = -diff(phi_obs)/(-dx_obs);
c_eval('Efield_obs = gseE?par.tlim(tint_all).data;',mms_id)
if 0 % from phi
  dn_obs = diff(phi_obs,2)/dx_obs/dx_obs*units.eps0/units.e;
  dn_obs = [dn_obs(1); tocolumn(dn_obs); dn_obs(end)]; % assume phi -> phi at edges
else % from E  
  c_eval('Etmp = gseE?par.tlim(tint_all);',mms_id)
  dtE = Etmp.time(2) - Etmp.time(1);
  ttmp = [Etmp.time+-0.5*dtE Etmp.time(end)+0.5*dtE];
  Efield_obs_todiff = Etmp.resample(ttmp);
  dn_obs = -diff(Efield_obs_todiff.data,1)*1e-3/(sign(vph)*dx_obs)*units.eps0/units.e;
end
dn_obs = dn_obs*1;

% Divergence of E
if 1
  %%
  c_eval('R? = gseR?.resample(gseE1).tlim(tint_all)*1e3;',1:4)
  c_eval('E? = gseE?.resample(gseE1).tlim(tint_all)*1e-3;',1:4)
  [divE,avE]=c_4_grad('R?','E?','div'); divE.name = 'div E';
  dn_divE = divE*units.eps0/units.e ; dn_divE.name = 'dn from div E';
end



tsEfieldFromPhi = irf.ts_scalar(phi.time(1:end-1) + 0.5*(phi.time(2)-phi.time(1)),EfieldFromPhi*1e3);
tsEfieldFromPhi.units = 'mV/m';
tsDnFromPhi = irf.ts_scalar(phi.time,dn_obs*1e-6);
tsDnFromPhi.units = 'cm^{-3}';
%%
[X_obs,V_obs] = meshgrid(x_obs,v_vec); X_obs = permute(X_obs,[2 1]); V_obs = permute(V_obs,[2 1]);
PHI_obs = repmat(tocolumn(phi_obs),1,nv);
VPH_obs = V_obs*0 + vph;
E_obs = units.me*(V_obs-vph).^2/2 - units.e*PHI_obs;
ifree_obs = find(E_obs>0); itrap_obs = find(E_obs<=0);

[Fflat_obs,Fflat_free_obs,Fflat_trap_obs] = get_f_flat(V_obs,n,vt,vd,PHI_obs,VPH_obs);
nfree_obs = nansum(Fflat_free_obs,2)*dv;
ntrap_flat_obs = nansum(Fflat_trap_obs,2)*dv;
ntrap_obs = ntot - torow(nfree_obs) + torow(dn_obs);
dntrap_obs = ntrap_obs - torow(ntrap_flat_obs);
%[Fabel_obs,Fabel_free_obs,Fabel_trap_obs] = get_f_abel(V_obs,n,vt,vd,PHI_obs,VPH_obs,dntrap_obs);
[Fabel_obs,Fabel_free_obs,Fabel_trap_obs] = get_f_abel_split_by_phi_at_zero(V_obs,n,vt,vd,PHI_obs,VPH_obs,dntrap_obs);
Fabel_obs = Fabel_obs + Fflat_trap_obs;              
%[Fscha_obs,Fscha_free_obs,Fscha_trap_obs,beta_obs] = get_f_schamel(V_obs,n,vt,vd,PHI_obs,VPH_obs,ntrap_obs,beta_range);

%FVscha_obs = Fscha_obs.*V_obs;
FVabel_obs = Fabel_obs.*V_obs;

tsFabel_nfree_obs = irf.ts_scalar(phi.time,nansum(Fabel_free_obs,2)*dv);
tsFabel_ntrap_obs = irf.ts_scalar(phi.time,nansum(Fabel_trap_obs,2)*dv);

% Set up specrec for plotting
Fspecrec.t = phi.time.epochUnix;  
Fspecrec.p = Fabel_obs;  
Fspecrec.p_label = 'f_e (m^{-4}s^{-1})';
Fspecrec.f = v_vec*1e-6;
Fspecrec.f_label = 'v (10^{3} km/s)';

%% Flux at edi energy/velocity interval
vind_edi_0 = intersect(find(v_vec>v_edi_minus),find(v_vec<v_edi_plus));
vind_edi_180 = intersect(find(v_vec<-v_edi_minus),find(v_vec>-v_edi_plus));
FVdv_abel_obs_edi_0 = nansum(FVabel_obs(:,vind_edi_0),2)*dv;
FVdv_abel_obs_edi_180 = nansum(FVabel_obs(:,vind_edi_180),2)*dv;
%FVdv_scha_obs_edi_0 = nansum(FVscha_obs(:,vind_edi_0),2)*dv;
%FVdv_scha_obs_edi_180 = nansum(FVscha_obs(:,vind_edi_180),2)*dv;

fluxModel0 = irf.ts_scalar(phi.time,abs(FVdv_abel_obs_edi_0)*1e-4); % 1/m2 -> 1/cm2
fluxModel0.units = 'cm^{-2}s^{-1}';
fluxModel0.name = 'Model flux at pa 180 in EDI energy range';
fluxModel180 = irf.ts_scalar(phi.time,abs(FVdv_abel_obs_edi_180)*1e-4); % 1/m2 -> 1/cm2
fluxModel180.units = 'cm^{-2}s^{-1}';
fluxModel0.name = 'Model flux at pa 0 in EDI energy range';

tsDnModel = irf.ts_scalar(phi.time,sum(Fabel_obs*dv,2)*1e-6);
tsDnModel.units = 'cm^{-3}';
tsDnModel.name = 'Model electron density ';

tsDnModel_all = irf.ts_scalar(phi.time,[sum(Fabel_obs*dv,2)*1e-6, sum(Fabel_free_obs*dv,2)*1e-6,sum(Fabel_trap_obs*dv,2)*1e-6]);
tsDnModel_all.units = 'cm^{-3}';
tsDnModel_all.name = 'Model electron density ';


phi_long = phi_vec;
dntrap_long = dntrap_obs;
[fitreslts,gof,fun_net_long,fun_net_prime_long] = createFit(torow(phi_vec),torow(dntrap_obs));
  
%% Plot
tint_phi = tint_all;
timeline = Epar.time;

if 0 % plot 1
  figure(93)
  nrows = 2;
  ncols = 5;
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
    plot(hca,x_obs*1e-3,phi_obs,x_mod*1e-3,phi_mod,x_obs*1e-3,phi_obs_nodetrend)
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = '\phi (V)';
    %legend(hca,{'obs','mod, Gaussian','obs no detrend'},'Box','off')
    irf_legend(hca,{'obs';{'mod:','Gaussian'};{'obs:','no detrend'}},[0.02,0.98])
    hca.Title.String = sprintf('tint_{obs} = %s - %s',tint_utc(1,:),tint_utc(2,:));
  end
  if 1 % f0 info 
    hca = h(isub); isub = isub + 1;
    if exist('h_info_str','var'); delete(h_info_str);end
    h_info_str = irf_legend(hca,str_info,[0.01 0.99],'color',[0 0 0]);  
    hca.Visible = 'off';
  end
  if 1 % F obs, Schamel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fscha_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'f_{Scha,obs} (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'))
    hca.Title.String = sprintf('Schamel, beta_{obs}=%.2f',beta_obs);
  end
  if 1 % F obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fabel_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'f_{Abel,obs} (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'))
    hca.Title.String = 'Abel obs';    
  end
  if 1 % F mod, Schamel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_mod*1e-3,V_mod*1e-6,Fscha_mod)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'f_{Scha,mod} (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'))
    hca.Title.String = sprintf('Schamel, beta_{mod}=%.2f',beta_mod);
  end
  if 1 % F mod, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_mod*1e-3,V_mod*1e-6,Fabel_mod)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'f_{Abel,mod} (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'))
    hca.Title.String = 'Abel mod (Gaussian)';
  end
  if 1 % Trapped densities
    hca = h(isub); isub = isub + 1;
    hlines = plot(hca,x_obs*1e-3,ntrap_obs*1e-6,x_obs*1e-3,nansum(Fabel_trap_obs,2)*dv*1e-6,x_obs*1e-3,nansum(Fscha_trap_obs,2)*dv*1e-6,...
                      x_mod*1e-3,ntrap_mod*1e-6,x_mod*1e-3,nansum(Fabel_trap_mod,2)*dv*1e-6,x_mod*1e-3,nansum(Fscha_trap_mod,2)*dv*1e-6);
    %hlines(1).LineWidth = 1.5;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'n_t (cm^{-3})';
    hca.Title.String = 'n_t = n_0-n_f+(\epsilon_0/e)\nabla^2\phi';
%     legend(hca,{'n_{t,obs}','n_{t,Abel,obs}','n_{t,Scha,obs}',...
%                 'n_{t,mod}','n_{t,Abel,mod}','n_{t,Scha,mod}'},'Box','off','Location','bestoutside')
    irf_legend(hca,{'n_{t,obs}';'n_{t,Abel,obs}';'n_{t,Scha,obs}';...
                    'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'},...
                    [0.35 0.1])

    irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
  end
  if 1 % Total densities
    hca = h(isub); isub = isub + 1;
    hlines = plot(hca,x_obs*1e-3,(n0+dn_obs)*1e-6,x_obs*1e-3,nansum(Fabel_obs,2)*dv*1e-6,x_obs*1e-3,nansum(Fscha_obs,2)*dv*1e-6,...
                      x_mod*1e-3,(n0+dn_mod)*1e-6,x_mod*1e-3,nansum(Fabel_mod,2)*dv*1e-6,x_mod*1e-3,nansum(Fscha_mod,2)*dv*1e-6);
    hlines(1).LineWidth = 1.5;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'n (cm^{-3})';
    hca.Title.String = 'n = n_0+(\epsilon_0/e)\nabla^2\phi';
%     legend(hca,{'n_{obs}','n_{Abel,obs}','n_{Scha,obs}',...
%                 'n_{mod}','n_{Abel,mod}','n_{Scha,mod}'},'Box','off','Location','bestoutside')
    irf_legend(hca,{'n_{t,obs}';'n_{t,Abel,obs}';'n_{t,Scha,obs}';...
                    'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'},...
                    [0.99 0.02])
    irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.1],'color',[0 0 0])
  end
  if 1 % nt(phi)
    hca = h(isub); isub = isub + 1;
    [fitreslts,gof,fun_net_obs,fun_net_prime_obs] = createFit(torow(phi_obs),torow(dntrap_obs));
    [fitreslts,gof,fun_net_mod,fun_net_prime_mod] = createFit(torow(phi_mod),torow(dntrap_mod));
    plot(hca,phi_obs,dntrap_obs,'o',phi_obs,fun_net_obs(phi_obs),...
             phi_mod,dntrap_mod,'+',phi_mod,fun_net_mod(phi_mod))
    hca.YLabel.String = 'n_t-n_{t,flat} (m^{-3})';  
    hca.XLabel.String = '\phi (V)';  
    legend(hca,{'n_{t,obs}','n_{t,obs} fit','n_{t,mod}','n_{t,mod} fit'},'Box','off')
  end
  if 0 % dnt/dphi)
    hca = h(isub); isub = isub + 1;
    [fitreslts,gof,fun_net,fun_net_prime] = createFit(torow(phi),torow(dntrap));
    plot(hca,phi,fun_net_prime(phi))
    hca.YLabel.String = 'dn_t/d\phi (m^{-3}/V)';  
    hca.XLabel.String = '~\phi (V)';
  end
  if 1 % F at a few x's
    hca = h(isub); isub = isub + 1; 
    colors = mms_colors('matlab');  
    set(hca,'ColorOrder',colors(1:3,:))
    indx = fix(nx_obs./[4 3 2]);
    %set(hca,'LineStyleOrder',{'-','--'})
  %   plot(hca,v,F_abel(fix(nx/6),:),v,F_abel(fix(nx/4),:),v,F_abel(fix(nx/2),:),...
  %            v,F_scha(fix(nx/6),:),v,F_scha(fix(nx/4),:),v,F_scha(fix(nx/2),:))
    set(hca,'ColorOrder',colors(1:3,:))
    plot(hca,v*1e-6,Fabel_obs(indx(1),:),v*1e-6,Fabel_obs(indx(2),:),v*1e-6,Fabel_obs(indx(3),:))
    hold(hca,'on'),
    set(hca,'ColorOrder',colors(1:3,:))
    plot(hca,v*1e-6,Fscha_obs(indx(1),:),'--',v*1e-6,Fscha_obs(indx(2),:),'--',v*1e-6,Fscha_obs(indx(3),:),'--')
    hold(hca,'off')
    %hca.YScale = 'log';
    hca.YLabel.String = 'f (s^1m^{-4})';  
    hca.XLabel.String = 'v (10^3 km)';
    hca.XLim = (vph+vtrap*[-1 1]*1.5)*1e-6;
    irf_legend(hca,{sprintf('x=%.0f km',x_obs(indx(1))*1e-3),sprintf('x=%.0f km',x_obs(indx(2))*1e-3),sprintf('x=%.0f km',x_obs(indx(3))*1e-3)},[0.1 0.15])
    irf_legend(hca,{'- Abel','--Scha'},[0.1 0.05],'color',[0 0 0])
    hca.Title.String = 'based on \phi_{obs}';
  end
  %cn.print(sprintf('AbelScha_obs_eh%g_mms%g',ih,mms_id))
end
if 0 % plot 2
    figure(93)
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
      plot(hca,x_obs*1e-3,phi_obs)
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = '\phi (V)';      
      hca.Title.String = sprintf('tint_{obs} = %s - %s',tint_utc(1,:),tint_utc(2,:));
      irf_legend(hca,{sprintf('v_{ph}= %.0f km/s',vph*1e-3)},[0.99,0.99],'Color',[0 0 0])
    end
    if 0 % f0 info 
      hca = h(isub); isub = isub + 1;
      if exist('h_info_str','var'); delete(h_info_str);end
      h_info_str = irf_legend(hca,str_info,[0.01 0.99],'color',[0 0 0]);  
      hca.Visible = 'off';
    end
    if 1 % E
      hca = h(isub); isub = isub + 1;
      pcolor(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e)
      shading(hca,'flat') 
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'v (10^3 km/s)';
      hcb.YLabel.String = 'E/e (eV)';  
      hca.YLim = vlim*[-1 1]*1e-6;
      hca.CLim = max(abs(E_obs(:)/units.e))*[-1 1];
      hca.CLim = 2e3*[-1 1];
      colormap(hca,cn.cmap('blue_red'))
       if 1 % E contours
        hold(hca,'on')
        %levels_E = [linspace(min(E_obs(:)),0,6) linspace(0,max(E_obs(:)),50)]; 
        levels_E = linspace(min(E_obs(:)),0,5); 
        levels_E = [levels_E 0:levels_E(2)-levels_E(1):max(E_obs(:))];                   
        [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,levels_E/units.e,'k');
        [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,0,'k','LineWidth',2);
        hold(hca,'off')
      end
      hcb.YLim = [min(E_obs(:)) hca.CLim(2)*units.e]/units.e;
      %hca.Title.String = sprintf('Schamel, beta_{obs}=%.2f',beta_obs);
    end
    if 1 % F obs, Abel
      hca = h(isub); isub = isub + 1;
      pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fabel_obs)
      shading(hca,'flat') 
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'v (10^3 km/s)';
      hcb.YLabel.String = 'f_{Abel,obs} (s^1m^{-4})';  
      hca.YLim = vlim*[-1 1]*1e-6;
      colormap(hca,cn.cmap('white_blue'))
      hca.Title.String = 'Abel obs';    
    end
    if 1 % F obs, Schamel
      hca = h(isub); isub = isub + 1;
      pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fscha_obs)
      shading(hca,'flat') 
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'v (10^3 km/s)';
      hcb.YLabel.String = 'f_{Scha,obs} (s^1m^{-4})';  
      hca.YLim = vlim*[-1 1]*1e-6;
      colormap(hca,cn.cmap('white_blue'))
      hca.Title.String = sprintf('Schamel, beta_{obs}=%.2f',beta_obs);
      if 0 % E contours, not seen
        hold(hca,'on')      
        contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs,0,'c');      
        hold(hca,'off')
      end
    end
    if 1 % Free, trapped, and total densities
      hca = h(isub); isub = isub + 1;
      hlines = plot(hca,x_obs*1e-3,(n0+dn_obs)*1e-6,...
                        x_obs*1e-3,nansum(Fabel_free_obs,2)*dv*1e-6,...
                        x_obs*1e-3,nansum(Fabel_trap_obs,2)*dv*1e-6...
                        );
      %hlines(1).LineWidth = 1.5;
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'n (cm^{-3})';
      hca.Title.String = 'n_t = n_0-n_f+(\epsilon_0/e)\nabla^2\phi';
      irf_legend(hca,{'n_{0}+\delta n';'n_{f}';'n_{t}'},...
                      [0.02 0.3])
  %                     'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'
      %irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
    end
    if 1 % FV obs, Abel
      hca = h(isub); isub = isub + 1;
      pcolor(hca,X_obs*1e-3,V_obs*1e-6,FVabel_obs)
      shading(hca,'flat') 
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'v (10^3 km/s)';
      hcb.YLabel.String = 'flux_{Abel,obs} (m^{-3})';  
      hca.YLim = vlim*[-1 1]*1e-6;
      hca.CLim = max(abs(FVabel_obs(:)))*[-1 1];    
      colormap(hca,cn.cmap('blue_red'))    
      hca.Title.String = 'Abel obs';
      if 1 % EDI energies
        hold(hca,'on')
        hlines = plot(hca,x_obs([1 end]),v_edi*1e-6*[1 1],x_obs([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
        hlines = plot(hca,...
          x_obs([1 end]),(v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(v_edi_minus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_minus)*1e-6*[1 1],...
          'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
        irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);   
        hold(hca,'off')
      end
    end
    if 1 % FV obs, Schamel
      hca = h(isub); isub = isub + 1;    
      pcolor(hca,X_obs*1e-3,V_obs*1e-6,FVscha_obs)
      shading(hca,'flat') 
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'v (10^3 km/s)';
      hcb.YLabel.String = 'flux_{Scha,obs} (m^{-3})';  
      hca.YLim = vlim*[-1 1]*1e-6;
      hca.CLim = max(abs(FVscha_obs(:)))*[-1 1];    
      colormap(hca,cn.cmap('blue_red'))        
      hca.Title.String = sprintf('Schamel, beta_{obs}=%.2f',beta_obs);
      if 0 % E contours, not seen
        hold(hca,'on')      
        contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs,0,'c');      
        hold(hca,'off')
      end    
      if 1 % EDI energies
        hold(hca,'on')
        hlines = plot(hca,x_obs([1 end]),v_edi*1e-6*[1 1],x_obs([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
        hlines = plot(hca,...
          x_obs([1 end]),(v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(v_edi_minus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_minus)*1e-6*[1 1],...
          'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
        irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);   
        hold(hca,'off')
      end
    end
    if 0 % FV mod, Abel
      hca = h(isub); isub = isub + 1;
      pcolor(hca,X_mod*1e-3,V_mod*1e-6,FVabel_mod)
      shading(hca,'flat') 
      hcb = colorbar('peer',hca);
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'v (10^3 km/s)';
      hcb.YLabel.String = 'flux_{Abel,obs} (m^{-3})';  
      hca.YLim = vlim*[-1 1]*1e-6;
      hca.CLim = max(abs(FVabel_mod(:)))*[-1 1];    
      colormap(hca,cn.cmap('blue_red'))    
      hca.Title.String = 'Abel obs';
      if 1 % EDI energies
        hold(hca,'on')
        hlines = plot(hca,x_obs([1 end]),v_edi*1e-6*[1 1],x_obs([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
        hlines = plot(hca,...
          x_obs([1 end]),(v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(v_edi_minus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_minus)*1e-6*[1 1],...
          'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
        irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);   
        hold(hca,'off')
      end
    end
    if 0 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0
      hca = h(isub); isub = isub + 1;    
      nodes = 1;
      units_scale = 1e-4; % m^-2 > cm^-2 
      units_scale_2 = 1e6;
      plot_EDI = mean(flux0.data(:,nodes),2)/units_scale_2;
      plot_abel = abs(FVdv_abel_obs_edi_0)*units_scale/units_scale_2;
      plot_scha = abs(FVdv_scha_obs_edi_0)*units_scale/units_scale_2;

      plot(hca,x_obs,plot_abel,x_obs,plot_scha,x_edi_0,plot_EDI);
      hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
      irf_legend(hca,{'Abel';'Scha';'EDI'},[0.02 0.99])
      text(hca,0,hca.YLim(2),'0^o','verticalalignment','top')
      hca.YLim(1) = 0;
    end
    if 0 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
      hca = h(isub); isub = isub + 1;    
      nodes = 1;
      units_scale = 1e-4; % m^-2 > cm^-2 
      units_scale_2 = 1e6;
      plot_EDI = mean(flux180.data(:,nodes),2)/units_scale_2;
      plot_abel = abs(FVdv_abel_obs_edi_180)*units_scale/units_scale_2;
      plot_scha = abs(FVdv_scha_obs_edi_180)*units_scale/units_scale_2;

      plot(hca,x_obs,plot_abel,x_obs,plot_scha,x_edi_180,plot_EDI);
      hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
      irf_legend(hca,{'Abel';'Scha';'EDI'},[0.02 0.99])
      text(hca,0,hca.YLim(2),'180^o','verticalalignment','top')
      hca.YLim(1) = 0;
    end
    if 0 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0 and 180
      hca = h(isub); isub = isub + 1;    
      nodes = 1;
      units_scale = 1e-4; % m^-2 > cm^-2 
      units_scale_2 = 1e6;
      plot_EDI_0 = mean(flux0.data(:,nodes),2)/units_scale_2;
      plot_abel_0 = abs(FVdv_abel_obs_edi_0)*units_scale/units_scale_2;
      plot_abel_mod_0 = abs(FVdv_abel_mod_edi_0)*units_scale/units_scale_2;
      plot_scha_0 = abs(FVdv_scha_obs_edi_0)*units_scale/units_scale_2;    
      plot_EDI_180 = mean(flux180.data(:,nodes),2)/units_scale_2;
      plot_abel_180 = abs(FVdv_abel_obs_edi_180)*units_scale/units_scale_2;
      plot_abel_mod_180 = abs(FVdv_abel_mod_edi_180)*units_scale/units_scale_2;
      plot_scha_180 = abs(FVdv_scha_obs_edi_180)*units_scale/units_scale_2;    

      hlines = plot(hca,x_obs*1e-3,plot_abel_0,x_obs*1e-3,plot_scha_0,x_obs*1e-3,plot_abel_mod_0,x_edi_0*1e-3,plot_EDI_0,...
                        x_obs*1e-3,plot_abel_180,x_obs*1e-3,plot_scha_180,x_obs*1e-3,plot_abel_mod_180,x_edi_180*1e-3,plot_EDI_180);
      colors = mms_colors('matlab');
      nlines = 4;
      c_eval('hlines(?).Color = colors(?,:);',1:nlines)
      c_eval('hlines(nlines+?).Color = colors(?,:);',1:nlines)

      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
      irf_legend(hca,{'Abel';'Scha';'Abel Gauss';'EDI'},[0.02 0.6])
      text(hca,0.6*hca.XLim(2),0.2*hca.YLim(2),'0^o','verticalalignment','bottom')
      text(hca,0.6*hca.XLim(2),1.0*hca.YLim(2),'180^o','verticalalignment','top')
      hca.YLim(1) = 0;
    end
    if 0 % Trapped densities
      hca = h(isub); isub = isub + 1;
      hlines = plot(hca,x_obs*1e-3,ntrap_obs*1e-6,x_obs*1e-3,nansum(Fabel_trap_obs,2)*dv*1e-6,x_obs*1e-3,nansum(Fscha_trap_obs,2)*dv*1e-6,...
                        x_mod*1e-3,ntrap_mod*1e-6,x_mod*1e-3,nansum(Fabel_trap_mod,2)*dv*1e-6,x_mod*1e-3,nansum(Fscha_trap_mod,2)*dv*1e-6);
      %hlines(1).LineWidth = 1.5;
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'n_t (cm^{-3})';
      hca.Title.String = 'n_t = n_0-n_f+(\epsilon_0/e)\nabla^2\phi';
  %     legend(hca,{'n_{t,obs}','n_{t,Abel,obs}','n_{t,Scha,obs}',...
  %                 'n_{t,mod}','n_{t,Abel,mod}','n_{t,Scha,mod}'},'Box','off','Location','bestoutside')
      irf_legend(hca,{'n_{t,obs}';'n_{t,Abel,obs}';'n_{t,Scha,obs}';...
                      'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'},...
                      [0.35 0.1])

      irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
    end
    if 0 % Total densities
      hca = h(isub); isub = isub + 1;
      hlines = plot(hca,x_obs*1e-3,(n0+dn_obs)*1e-6,x_obs*1e-3,nansum(Fabel_obs,2)*dv*1e-6,x_obs*1e-3,nansum(Fscha_obs,2)*dv*1e-6,...
                        x_mod*1e-3,(n0+dn_mod)*1e-6,x_mod*1e-3,nansum(Fabel_mod,2)*dv*1e-6,x_mod*1e-3,nansum(Fscha_mod,2)*dv*1e-6);
      hlines(1).LineWidth = 1.5;
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'n (cm^{-3})';
      hca.Title.String = 'n = n_0+(\epsilon_0/e)\nabla^2\phi';
  %     legend(hca,{'n_{obs}','n_{Abel,obs}','n_{Scha,obs}',...
  %                 'n_{mod}','n_{Abel,mod}','n_{Scha,mod}'},'Box','off','Location','bestoutside')
      irf_legend(hca,{'n_{t,obs}';'n_{t,Abel,obs}';'n_{t,Scha,obs}';...
                      'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'},...
                      [0.99 0.02])
      irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.1],'color',[0 0 0])
    end
    if 1 % Total densities, not alll
      hca = h(isub); isub = isub + 1;
      hlines = plot(hca,x_obs*1e-3,(n0+dn_obs)*1e-6,x_obs*1e-3,nansum(Fabel_obs,2)*dv*1e-6,x_obs*1e-3,nansum(Fscha_obs,2)*dv*1e-6);
      hlines(1).LineWidth = 1.5;
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = 'n (cm^{-3})';
      %hca.Title.String = 'n = n_i+(\epsilon_0/e)\nabla^2\phi';
  %     legend(hca,{'n_{obs}','n_{Abel,obs}','n_{Scha,obs}',...
  %                 'n_{mod}','n_{Abel,mod}','n_{Scha,mod}'},'Box','off','Location','bestoutside')
      irf_legend(hca,{'n=n_i+(\epsilon_0/e)\nabla^2\phi';'n_{Abel}';'n_{Scha}'},...
                      [0.99 0.06])
      irf_legend(hca,{sprintf('beta_{Scha}=%.2f',beta_obs)},[0.02 0.1],'color',[0 0 0])
    end
    if 1 % nt(phi)
      hca = h(isub); isub = isub + 1;
      [fitreslts,gof,fun_net_obs,fun_net_prime_obs] = createFit(torow(phi_obs),torow(dntrap_obs));
      %[fitreslts,gof,fun_net_mod,fun_net_prime_mod] = createFit(torow(phi_mod),torow(dntrap_mod));
      plot(hca,phi_obs,dntrap_obs,'o',phi_obs,fun_net_obs(phi_obs))
      hca.YLabel.String = 'n_t-n_{t,flat} (m^{-3})';  
      hca.XLabel.String = '\phi (V)';  
      legend(hca,{'n_{t,obs}','n_{t,obs} fit'},'Box','off')
    end
    if 1 % dnt/dphi)
      hca = h(isub); isub = isub + 1;
      [fitreslts,gof,fun_net_obs,fun_net_prime_obs] = createFit(torow(phi_obs),torow(dntrap_obs));      
      plot(hca,phi_obs,fun_net_prime_obs(phi_obs))
      hca.YLabel.String = 'dn_t/d\phi (m^{-3}/V)';  
      hca.XLabel.String = '\phi (V)';
    end
    if 1 % F at a few x's
      hca = h(isub); isub = isub + 1; 
      colors = mms_colors('matlab');  
      set(hca,'ColorOrder',colors(1:3,:))
      indx = fix(nx_obs./[4 3 2]);
      %set(hca,'LineStyleOrder',{'-','--'})
    %   plot(hca,v,F_abel(fix(nx/6),:),v,F_abel(fix(nx/4),:),v,F_abel(fix(nx/2),:),...
    %            v,F_scha(fix(nx/6),:),v,F_scha(fix(nx/4),:),v,F_scha(fix(nx/2),:))
      set(hca,'ColorOrder',colors(1:3,:))
      plot(hca,v*1e-6,Fabel_obs(indx(1),:),v*1e-6,Fabel_obs(indx(2),:),v*1e-6,Fabel_obs(indx(3),:))
      hold(hca,'on'),
      set(hca,'ColorOrder',colors(1:3,:))
      plot(hca,v*1e-6,Fscha_obs(indx(1),:),'--',v*1e-6,Fscha_obs(indx(2),:),'--',v*1e-6,Fscha_obs(indx(3),:),'--')
      hold(hca,'off')
      %hca.YScale = 'log';
      hca.YLabel.String = 'f (s^1m^{-4})';  
      hca.XLabel.String = 'v (10^3 km)';
      hca.XLim = (vph+vtrap*[-1 1]*1.5)*1e-6;
      irf_legend(hca,{sprintf('x=%.0f km',x_obs(indx(1))*1e-3),sprintf('x=%.0f km',x_obs(indx(2))*1e-3),sprintf('x=%.0f km',x_obs(indx(3))*1e-3)},[0.1 0.15])
      irf_legend(hca,{'- Abel','--Scha'},[0.1 0.05],'color',[0 0 0])
      hca.Title.String = 'based on \phi_{obs}';
    end
    %cn.print(sprintf('AbelScha_obs_eh%g_flux_mms%g',ih,mms_id))
    1;
  end
if 0 % plot, timeseries
fig = figure(39);
h = irf_plot(8);
isub = 1;

vlim_f = 30;
if 1 % Epar
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{phi_progressive.Etoint,Epar.tlim(tint_all),tsEfieldFromPhi},'comp');  
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  hca.YLabel.Interpreter = 'tex';
end
if 0 % Efrom phi
  hca = h(isub); isub = isub + 1;  
  irf_plot(hca,tsEfieldFromPhi)
  hca.YLabel.String = {'E','(mV/m)'};
end
if 1 % PHI, TSeries, plot
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,phi);  
  irf_legend(hca,{sprintf('v_{ph}= %g km/s',vph*1e-3)},[0.1 0.99],'color',[0 0 0]);
  hca.YLabel.String = {'\phi','(V)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 1 % diff E (Poisson) (density)
  hca = h(isub); isub = isub + 1;  
  nscale = 1e-3;
  irf_plot(hca,{tsDnFromPhi/nscale,(tsDnModel-ntot*1e-6)/nscale},'comp')
  hca.YLabel.String = {'\delta n',sprintf('(10^{%.0f} cm^{-3})',log10(nscale))};
  hca.YLabel.Interpreter = 'tex';
end
if 1 % F
  hca = h(isub); isub = isub + 1;
  irf_spectrogram(hca,Fspecrec,'lin');
  if 1 % EDI energies
    hold(hca,'on')
    hlines = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'k');
    hlines = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'k');
    irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 0 % model phase velocity
    hold(hca,'on')
    line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
    hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
    irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 0 % observed phase velocity
    hold(hca,'on')
    hlines = plot(hca,obs_eh_xvec,obs_velocity*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
    irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end  
  if 0 % str info
  str_info = {'unperturbed f:';...
    ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
    ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
    ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
    sprintf('beta_{Schamel}=%g',beta);...
    };
  irf_legend(hca,str_info,[1.01 1.4],'color',hlines(1).Color);    
  end
end
if 1 % sum F (density)
  hca = h(isub); isub = isub + 1;    
  irf_plot(hca,tsDnModel);
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YLabel.Interpreter = 'tex';
end
if 0 % F free
  hca = h(isub); isub = isub + 1;
  irf_spectrogram(hca,Fspecrec,'lin');
  if 1 % EDI energies
    hold(hca,'on')
    hlines = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'k');
    hlines = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'k');
    irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 0 % model phase velocity
    hold(hca,'on')
    line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
    hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
    irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 0 % observed phase velocity
    hold(hca,'on')
    hlines = plot(hca,obs_eh_xvec,obs_velocity*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
    irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end  
  if 0 % str info
  str_info = {'unperturbed f:';...
    ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
    ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
    ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
    sprintf('beta_{Schamel}=%g',beta);...
    };
  irf_legend(hca,str_info,[1.01 1.4],'color',hlines(1).Color);    
  end
end
if 0 % sum F (density,free,trap,all)
  hca = h(isub); isub = isub + 1;                  
  plot(hca,x_vec([1 end]),ntot*1e-6*[1 1],x_vec,mod_density*1e-6,x_vec,mod_density_free*1e-6,x_vec,mod_density_trap*1e-6,x_vec,net*1e-6); % 1e-6 from 1/m3 > 1/cm3
  legend(hca,{'n_{mod,\phi=0}','n_{mod}','n_{mod,free}','n_{mod,trap}','n_{mod,\phi=0} - n_{mod,free}+\delta n'},'location','eastoutside');
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % sum F (density)
  hca = h(isub); isub = isub + 1;  
  if 1
    plot(hca,x_vec,mod_density*1e-6,x_vec_diff1,obs_density*1e-6); % 1e-6 from 1/m3 > 1/cm3
    legend(hca,{'n_{mod}','n_{0} + (\epsilon_0/e)\nabla^2\phi'},'location','southwest','box','off');
    %irf_legend(hca,{sprintf('Q_X=%.2f, Q_A=%g',Q_Density_xcorr,Q_Density_area)},[0.01 0.1],'color',[0 0 0])
  else
    plot(hca,x_vec([1 end]),ntot*1e-6*[1 1],x_vec,mod_density*1e-6,x_vec([1 end]),mean(mod_density*1e-6)*[1 1],x_vec_diff1,obs_density*1e-6); % 1e-6 from 1/m3 > 1/cm3
    legend(hca,{'n_{mod,\phi=0}','n_{mod}','<n_{mod}>','n_{obs,Poisson,1D}'},'location','eastoutside');
  end
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % flux: F*v
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,FV);
  shading(hca,'flat')
  hca.YLabel.String = {'v','(10^3 km/s)'};
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f*v (1/m^3)'; 
  hca.CLim = hca.CLim(2)*[-1 1]*1.5; 
  hcb.YLim = hca.CLim;
  hca.YLim = vlim_f*[-1 1];
  if 1 % EDI energies
    hold(hca,'on')
    hlines = plot(hca,x_vec([1 end]),v_edi*1e-6*[1 1],x_vec([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
    for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
    hlines = plot(hca,...
      x_vec([1 end]),(v_edi_plus)*1e-6*[1 1],...
      x_vec([1 end]),(v_edi_minus)*1e-6*[1 1],...
      x_vec([1 end]),(-v_edi_plus)*1e-6*[1 1],...
      x_vec([1 end]),(-v_edi_minus)*1e-6*[1 1],...
      'LineWidth',1.5);
    for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
    irf_legend(hca,{'-- EDI'},[0.01 0.99],'color',hlines(1).Color);   
    hold(hca,'off')
  end
  if 1 % model phase velocity
    hold(hca,'on')
    line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
    hlines = plot(hca,x_vec,vph_vec*1e-6,'LineWidth',1.5,'Color',line_color(1,:),'LineStyle','-.');     
    irf_legend(hca,{'-. v_{mod}'},[0.2 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end
  if 1 % observed phase velocity
    hold(hca,'on')
    hlines = plot(hca,obs_eh_xvec,obs_velocity*1e-3,'*k','LineWidth',1.5,'Color',[0 0 0]);
    irf_legend(hca,{'* v_{obs}'},[0.1 0.99],'color',hlines(1).Color);  
    hold(hca,'off')
  end  
  colormap(hca,cn.cmap('blue_red'))
end
if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
  %%
  hca = h(isub); isub = isub + 1;
  %irf_plot({fluxModel180,ts_edi_flux180},'comp')  
  irf_plot(hca,{fluxModel180},'comp')  
  hca.YLabel.String = {'flux 180^o',sprintf('(10^%g cm^{-2}s^{-1})',log10(1e6))};
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'model''EDI'},[0.05 0.98])
  %text(hca,0.002,0.99*hca.YLim(2),'180^o','verticalalignment','top')
  hca.YLim(1) = 0;  
end
if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0
  %%
  hca = h(isub); isub = isub + 1;
  %irf_plot({fluxModel0,ts_edi_flux0},'comp')  
  irf_plot(hca,{fluxModel0},'comp')  
  hca.YLabel.String = {'flux 0^o',sprintf('(10^%g cm^{-2}s^{-1})',log10(1e6))};
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'model''EDI'},[0.05 0.98])
  %text(hca,0.002,0.99*hca.YLim(2),'180^o','verticalalignment','top')
  hca.YLim(1) = 0;  
end
if 0 % F tot
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,V*1e-6,Ftot);
  shading(hca,'flat')
  hca.YLabel.String = 'v (km/s)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'F';
  hca.CLim = hca.CLim(2)*[-1 1]; 
  colormap(hca,cn.cmap('blue_red'))
end
if 0 % sum FVdv, bulk velocity
  hca = h(isub); isub = isub + 1;
  plot(hca,x_vec,sumFVdv./sumFdv*1e-3);  
  hca.YLabel.String = {'v','(km/s)'};
end          

irf_plot_axis_align

if 0 % average over time, comaprison to FPI, add to right side of figure
  %%
  figure(49)
  %h = subplot(nrows,3,[3 6 9 12]);
  %h = subplot(nrows,2,[2 4 6 8]+0);
  h = subplot(1,1,1);
  %h.Position = [0.65    0.65    0.3    0.3];
  isub = 1;
  if 1 % F, 
    hca = h(isub); isub = isub + 1;
    v_scale = 1e-3;
    hlines = plot(hca,v_vec*v_scale*1e-3,mod_f_average,v_vec*v_scale*1e-3,mod_fmax_average,v_fpi*v_scale,f_fpi,'--');
    hca.YLabel.String = {'f (s^1/m^4)'};
    hca.XLabel.String = {'v (10^3 km/s)'};
    hca.XLim = [-40 40];
    fpi_utc = ts_f_fpi.time.utc;
    str_lines = {'<f_{mod}>';'f_{mod,\phi=0}';...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};
    %legend(hlines,str_lines)
    irf_legend(hca,str_lines,[0.99 0.99])
    str_info = {['T_{in}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      sprintf('beta_{Schamel}=%g',beta);...
      };
    set(hca,'ColorOrder',zeros(10,3))
    irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];
       if 1
         plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
         irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
       end
      hold(hca,'off')
    end
  end
end
if 0 % nt(phi)        
  %%
  h = subplot(10,10,100);
  h.Position = [0.65    0.35    0.3    0.2];
  isub = 1;                
  hca = h(isub); isub = isub + 1;
  [fitresult, gof, fun_net, fun_net_prime] = createFit(phi_vec, dntrap);
  [fitresult, gof, fun_net_noscaling, fun_net_prime_noscaling] = createFit(phi_vec, dntrap/dntrap_scaling);
  plot(hca,phi_vec,dntrap,'.',phi_vec,fun_net(phi_vec),phi_vec,dntrap/dntrap_scaling,'.',phi_vec,fun_net_noscaling(phi_vec))  
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = '\delta n - n_{t,flat} (m^{-3})';                  
  irf_legend(hca,{'data used';'fit with scaling';'data original';'fit no scaling/original'},[0.98 0.98])
end
if 0 % dnt(phi) function
  %%
  h = subplot(10,10,99);
  h.Position = [0.65    0.10    0.3    0.2];
  isub = 1;                
  hca = h(isub); isub = isub + 1;
  [fitresult, gof, fun_net, fun_net_prime] = createFit(phi_vec, dntrap);
  plot(hca,phi_vec,fun_net_prime(phi_vec))  
  hca.XLabel.String = '\phi (V)';
  hca.YLabel.String = 'd(\delta n - n_{t,flat})/d\phi (m^{-3}V^{-1})';
end
if 0 % flat skymap, to check how normalization is done
  %h = subplot(nrows,3,[3 6 9 12]);
  h = subplot(nrows,2,[2 4]+8);
  h.Position = [0.65    0.35    0.3    0.2];
  isub = 1;
  if 1 % F, 
    hca = h(isub); isub = isub + 1;
    [ax,hcb] = mms.plot_skymap(hca,ePDist1,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
    %hca.YLabel.String = {'f','(s^1/m^4)'};
    %hca.XLabel.String = {'v','(10^3 km/s)'};
    str_lines = {'<f_{mod}>';'f_{mod,\phi=0}';...
      sprintf('-- fpi: %s',fpi_utc(1,12:23));...
      sprintf('-- fpi: %s',fpi_utc(2,12:23));...
      sprintf('-- fpi: %s',fpi_utc(3,12:23));...
      sprintf('-- fpi: %s',fpi_utc(4,12:23))};
    %legend(hlines,str_lines)
    %irf_legend(hca,str_lines,[0.99 0.99])                  
  end
end
if 0 % flat skymap, units of flux, to check how normalization is done
  %%
  %h = subplot(nrows,3,[3 6 9 12]);
  h = subplot(nrows,2,[2 4]+14);
  h.Position = [0.65    0.05    0.3    0.2];
  isub = 1;                
  hca = h(isub); isub = isub + 1;
  if 1
  echannel = 17;
  E_fpi = mean(ePDist1.tlim(tint_phi).depend{1}(:,echannel),1);
  f_485 = squeeze(mean(ePDist1.convertto('s^3/m^6').tlim(tint_phi).data(:,echannel,:,:)));
  ang_polar = ePDist1.tlim(tint_phi).depend{3};
  ang_azim = mean(ePDist1.tlim(tint_phi).depend{2},1);
  d_azim = pi/16;
  d_pol = 1 - cosd(180/16);
  d_vel_vec = linspace(v_edi_minus,v_edi_plus,100);
  d_vel = trapz(d_vel_vec,d_vel_vec.^2); % [v^2dv] = m3/s3
  v_fpi = sqrt(2*units.e*E_fpi./units.me); % m/s
  f_vol = d_vel*d_azim*d_pol;
  flux_fpi_485 = f_485*v_fpi*f_vol; % m-2s-1, (to make into cm-2s-1, multiply with 1e-4)
  flux_scale = 1e6;
  surf(hca,ang_azim,ang_polar,flux_fpi_485'*1e-4/flux_scale); % cm-2s-1
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(flux_scale));
  hca.YLabel.String = 'Polar angle';
  hca.XLabel.String = 'Azimuthal angle';
  hca.Title.String = 'Where electrons are going';
  view(hca,[0 0 1])
  hca.XLim = [0 360];
  hca.YLim = [0 180];
  %[ax,hcb] = mms.plot_skymap(hca,ePDist1.convertto('s^3/m^6'),'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
  else
    %%
    [ax,hcb] = mms.plot_skymap(hca,ePDist1.vd3v,'flat','tint',tint_phi,'energy',480,'vectors',{mean(dmpaB1.tlim(tint_phi).norm.data,1),'B'});
  end
  %hca.YLabel.String = {'f','(s^1/m^4)'};
  %hca.XLabel.String = {'v','(10^3 km/s)'};
  str_lines = {'<f_{mod}>';'f_{mod,\phi=0}';...
    sprintf('-- fpi: %s',fpi_utc(1,12:23));...
    sprintf('-- fpi: %s',fpi_utc(2,12:23));...
    sprintf('-- fpi: %s',fpi_utc(3,12:23));...
    sprintf('-- fpi: %s',fpi_utc(4,12:23))};
end
if 0%doPrint
  cn.print(sprintf('schamel_2F_vph%g_ntot%g_R%g_T1%g_T2%g_vd1%g_vd2%g_beta%g_phishift%g',vph,ntot*1e-6,R,T(1),T(2),vd(1)*1e-3,vd(2)*1e-3,beta,phi_shift))
end

if 1 % average over time, comaprison to FPI
  %%
  figure(40)
  hca = subplot(1,1,1);
  clear h; 
  isub = 1;
  mod_f_average = mean(Fabel_obs,1);
  mod_f0 = f0(v_vec,n,vd,vt);
  
  if 1 % F, 
    hca = h(isub); isub = isub + 1;
    v_scale = 1e-3;
    hlines = plot(hca,v_vec*v_scale*1e-3,mod_f_average,v_vec*v_scale*1e-3,mod_fmax_average,v_fpi*v_scale,f_fpi,'--');
    hca.YLabel.String = {'f','(s^1/m^4)'};
    hca.XLabel.String = {'v','(10^3 km/s)'};
    hca.XLim = [-40 40];
    str_lines = {'f_{mod}';'f_{mod,\phi=0}';'-- fpi';'-- fpi';'-- fpi';'-- fpi'};
    %legend(hlines,str_lines)
    irf_legend(hca,str_lines,[0.99 0.99])
    str_info = {['T_{in}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      sprintf('beta_{Schamel}=%g',beta);...
      };
    set(hca,'ColorOrder',zeros(10,3))
    irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];
       if 1
         plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
         irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
       end
      hold(hca,'off')
    end
  end
  if doPrint 
  cn.print(sprintf('schamel_2F_average_vph%g_ntot%g_R%g_T1%g_T2%g_vd1%g_vd2%g_beta%g',ntot*1e-6,R,T(1),T(2),vd(1)*1e-3,vd(2)*1e-3,beta))
end
end
if 0 % average over time, comaprison to FPI
  %%
  figure(34)
  hca = subplot(nrows,3,[3 6 9]);
  clear h;
  nrows = 1;
  ncols = 1;
  npanels = nrows*ncols;
  for ip = 1:npanels
    h(ip) = subplot(nrows,ncols,ip);
  end
  isub = 1;
  if 1 % F, 
    hca = h(isub); isub = isub + 1;
    v_scale = 1e-3;
    hlines = plot(hca,v_vec*v_scale*1e-3,mod_f_average,v_vec*v_scale*1e-3,mod_fmax_average,v_fpi*v_scale,f_fpi,'--');
    hca.YLabel.String = {'f','(s^1/m^4)'};
    hca.XLabel.String = {'v','(10^3 km/s)'};
    hca.XLim = [-40 40];
    str_lines = {'f_{mod}';'f_{mod,\phi=0}';'-- fpi';'-- fpi';'-- fpi';'-- fpi'};
    %legend(hlines,str_lines)
    irf_legend(hca,str_lines,[0.99 0.99])
    str_info = {['T_{in}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      sprintf('beta_{Schamel}=%g',beta);...
      };
    set(hca,'ColorOrder',zeros(10,3))
    irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];
       if 1
         plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
         irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
       end
      hold(hca,'off')
    end
  end
  if doPrint 
  cn.print(sprintf('schamel_2F_average_vph%g_ntot%g_R%g_T1%g_T2%g_vd1%g_vd2%g_beta%g',ntot*1e-6,R,T(1),T(2),vd(1)*1e-3,vd(2)*1e-3,beta))
end
end
end
if 1 % 1 % plot, timeseries, for paper
  %%
  fig = figure(39);
  npanels = 5;
  h = irf_plot(npanels);
  clear h_all;
  h_all = h;
  isub = 1;

  vlim_f = 30;
  if 1 % Epar
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,{Epar.tlim(tint_all)},'comp');  
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    hca.YLabel.Interpreter = 'tex';
  end
  if 1 % PHI, TSeries, plot
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,phi);  
    irf_legend(hca,{sprintf('v_{ph,av}= %g km/s',vph*1e-3)},[0.08 0.99],'color',[0 0 0]);
    hca.YLabel.String = {'\phi','(V)'};  
    hca.YLabel.Interpreter = 'tex';
    if 1 % add locations where dphi/dntrap is recalculated every time
      hold(hca,'on')
      irf_plot(hca,phi_progressive.ts_detrend_locs,'*','color',[0,0.4470,0.7410])
      hold(hca,'off')
    end
  end
  if 1 % F
    hca = h(isub); isub = isub + 1;
    irf_spectrogram(hca,Fspecrec,'lin');
    hca.YLabel.String = {'v_{||}','(10^3 km/s)'};
    edi_color = [1 1 1];
    vph_color = [1 1 1];
    h_all_markings = [];
    if 1 % EDI energies
      hold(hca,'on')
      %lines_EDI_plus = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'color',edi_color);
      lines_EDI_minus = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'color',edi_color);
      %hleg_EDI = irf_legend(hca,{'- EDI'},[0.08 0.99],'color',lines_EDI_minus(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_EDI_plus; lines_EDI_minus; hleg_EDI];
      %h_all_markings = [h_all_markings; lines_EDI_minus; hleg_EDI]; 
    end
    if 1 % model phase velocity
      hold(hca,'on')
      line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
      lines_vphav = irf_plot(hca,tsVph*1e-6,'LineWidth',1.0,'Color',vph_color,'LineStyle','--');
      %lines_vphav = irf_plot(hca,tsVph*1e-6,'--k');
      %hleg_vphav = irf_legend(hca,{'-- v_{ph,mod}'},[0.32 0.99],'color',lines_vphav(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_vphav; hleg_vphav];
    end
    if 1 % observed phase velocity
      hold(hca,'on')
      lines_vphobs = irf_plot(hca,obs_vph*1e-3,'*k','LineWidth',1.5,'Color',vph_color);
      lines_vphobs.MarkerSize = 4;
      %hleg_vphobs = irf_legend(hca,{'* v_{ph,obs}'},[0.18 0.99],'color',lines_vphobs(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_vphobs; hleg_vphobs];
    end  
    colormap(hca,cn.cmap('white_blue'))
    if 0 % str info
    str_info = {'unperturbed f:';...
      ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      sprintf('beta_{Schamel}=%g',beta);...
      };
    irf_legend(hca,str_info,[1.01 1.4],'color',hlines(1).Color);    
    end
    hca.YLim = vlim_f*[-1 1];
  end
  if 1 % diff E (Poisson) (density)
    hca = h(isub); isub = isub + 1;  
    nscale = 1e-3;
    irf_plot(hca,{-1*tsDnFromPhi/nscale,-1*(tsDnModel-ntot*1e-6)/nscale},'comp')
    hca.YLabel.String = {'\delta n',sprintf('(10^{%.0f} cm^{-3})',log10(nscale))};
    hca.YLabel.Interpreter = 'tex';
    %irf_legend(hca,{'n_i - (\epsilon_0/e)\partial_{||}E_{||}  ','  \int f_{model}dv_{||}'},[0.01 0.10]);
    %irf_legend(hca,{'n_{et}^{obs}  ','  n_{et}^{mod}'},[0.01 0.10]);
    irf_legend(hca,{'obs.  ','model'},[0.08 0.98]);
    doDoubleAxis = 1; % dn
    if doDoubleAxis  
      ax1 = hca;
      ax2 = axes('Position',get(ax1,'Position'));
      ax2.YLim = ax1.YLim/nscale/n0;    
      set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none','box','off'); % color of axis      
      %set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
      %set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
      irf_timeaxis(ax2,'nolabels')
      ax2.XLabel.String = [];
      ax2.YLabel.String = {'\delta n/n'};
      ax2.YLabel.Interpreter = 'tex';    
      ax2.YTick = hca.YTick*1e-3/n0;  
      h_all = [h_all,ax2];
      ax2_dn = ax2;
      ax1_dn = ax1;
      ax2.YTick = ax1.YTick;
    end      
  end  
  if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
    hca = h(isub); isub = isub + 1;
    %%        
    f_scale = 1e6;
    colors = mms_colors('matlab');
    %colors = colors([2 1],:)
    %hca = subplot(1,1,1);
    
    irf_plot(hca,ts_edi_flux180/f_scale,'color',colors(1,:));
    ax1 = hca;
    ax2 = axes('Position',get(ax1,'Position'));
    irf_plot(ax2,fluxModel180.resample(ts_edi_flux180)/f_scale,'color',colors(2,:))      
    hold(ax2,'on')
    %irf_plot(ax2,fluxModel180/f_scale,'color',colors(3,:))      
    hold(ax2,'off')
    %irf_plot(ax2,fluxModel180/f_scale,'color',colors(2,:))      
    set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
    set(ax2,'YAxisLocation','right');
    set(ax2,'Color','none','box','off'); % color of axis      
    set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
    set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
    irf_timeaxis(ax2,'nolabels')
    
    ax2.YLabel.String = {'j^{mod}',sprintf('(10^%g cm^{-2}s^{-1})',log10(f_scale))};
    ax2.YLabel.Interpreter = 'tex';
    ax1.YLabel.String = {'j^{EDI}',sprintf('(10^%g cm^{-2}s^{-1}sr^{-1})',log10(f_scale))};
    ax1.YLabel.Interpreter = 'tex';
    irf_legend(hca,{'EDI','  model'},[0.08 0.98])
    %irf_legend(hca,{'EDI: \theta = [168.75 180]^o','model: v = [-13600 -12900] km/s'},[0.08 0.98])
    %text(hca,0.002,0.99*hca.YLim(2),'180^o','verticalalignment','top')    
    ax1.YLim = [0 3.7];
    ax2.YLim = [0 3.7];
    h_all = [h_all,ax2];
    ax2_flux = ax2;
    ax1_flux = ax1;    
  end
  if 0 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0
    hca = h(isub); isub = isub + 1;
    %%        
    f_scale = 1e6;
    colors = mms_colors('matlab');
    %colors = colors([2 1],:)
    %hca = subplot(1,1,1);
    if 0
      %irf_plot({fluxModel180,ts_edi_flux180},'comp')  
    else
      irf_plot(hca,ts_edi_flux0/f_scale,'color',colors(1,:));
      ax1 = hca;
      ax2 = axes('Position',get(ax1,'Position'));
      irf_plot(ax2,fluxModel0/f_scale,'color',colors(2,:))      
      set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required      
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none','box','off'); % color of axis      
      set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
      set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
      irf_timeaxis(ax2,'nolabels')
    end
    ax2.YLabel.String = {'flux 0^o, Model',sprintf('(10^%g cm^{-2}s^{-1})',log10(f_scale))};
    ax2.YLabel.Interpreter = 'tex';
    ax1.YLabel.String = {'flux 0^o, EDI',sprintf('(10^%g cm^{-2}s^{-1}sr^{-1})',log10(f_scale))};
    ax1.YLabel.Interpreter = 'tex';
    irf_legend(hca,{'EDI','model'},[0.02 0.98])
    %text(hca,0.002,0.99*hca.YLim(2),'180^o','verticalalignment','top')    
    ax1.YLim = [0 4];
    ax2.YLim = [0 4];
    h_all = [h_all,ax2];
  end

  if 0 % sum FVdv, bulk velocity
    hca = h(isub); isub = isub + 1;
    plot(hca,x_vec,sumFVdv./sumFdv*1e-3);  
    hca.YLabel.String = {'v','(km/s)'};
  end   
    
  irf_plot_axis_align(h)
  
  h(3).YLabel.String = {'v_{||}','(10^3 km/s)'};
  h(3).YLabel.Interpreter = 'tex';
  
  iref = 1;
  ijmp = 1;
  %h_all(iref+4+ijmp).Position = h_all(iref+4).Position;  
  
  
  %h_all(iref+7).Position = h_all(iref+5).Position;
  
  %irf_zoom(h_all,'x',[tint_phi(1) ts_edi_flux180.time(end)])  % tint_phi
  irf_zoom(h_all,'x',[tint_phi(1) tint_phi(end)])  % tint_phi
  
  h_all(iref+4+ijmp).XLabel = [];
  
  h(1).YLim = [-70 70];
  h(4).YLim = [-0.9 1.7];
  
  ax2_flux.Position = ax1_flux.Position;
  ax2_dn.Position = ax1_dn.Position;
  ax2_dn.YTick = ax1_dn.YTick/1e-3/n0;
  ax2_dn.YLim = ax1_dn.YLim/1e-3/n0;

  
  h(3).YLabel.Position(1) = h(1).YLabel.Position(1);
  
  %h_all(iref+7).XLabel = [];
  %h_all(end-2).YAxisLocation = 'right';
  %h_all(end-3).YAxisLocation = 'right';
  legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
  legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
  for ipanel = 1:npanels
    irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  end
  c_eval('h_all(?).XGrid = ''off''; h_all(?).YGrid = ''off'';',1:numel(h_all))
  
  % annotation
  if 0  % for h(3).YLim = [-29 29]
    h(3).YLim = [-29.9 29.9]
    annotation('textarrow',[0.71 0.71],[0.47 0.485],'string',{'EDI range'},'fontsize',12,'horizontalalignment','center');
    annotation('textarrow',[0.76 0.76],[0.550 0.500],'string',{'v_{ph,av}'},'fontsize',12,'horizontalalignment','center');
    annotation('textarrow',[0.30 0.325],[0.550 0.510],'string',{'v_{ph,ind}'},'fontsize',12,'horizontalalignment','center');
  elseif 1 
    %%
    h(3).YLim = [-29 15]
    annotation('textarrow',[0.685 0.685],[0.47 0.495],'string',{'EDI range'},'fontsize',12,'horizontalalignment','center');
    annotation('textarrow',[0.76 0.76],[0.550 0.500]+0.02,'string',{'v_{ph,av}'},'fontsize',12,'horizontalalignment','center');
    annotation('textarrow',[0.30 0.325],[0.550 0.510]+0.02,'string',{'v_{ph,ind}'},'fontsize',12,'horizontalalignment','center');
  end
end
if 0 % 0 % plot, timeseries and averaged, for diagnostics
  %%
  fig = figure(41);
  npanels = 5;
  [h,h2] = initialize_combined_plot(npanels,2,1,0.6,'vertical'); % horizontal
  clear h_all;
  h_all = h;
  isub = 1;

  vlim_f = 30;
  if 1 % Epar
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,{Epar.tlim(tint_all)},'comp');  
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    hca.YLabel.Interpreter = 'tex';
  end
  if 1 % PHI, TSeries, plot
    hca = h(isub); isub = isub + 1;
    irf_plot(hca,phi);  
    irf_legend(hca,{sprintf('v_{ph,av}= %g km/s, phi multiplier = %.1f',vph*1e-3,phi_mult)},[0.08 0.99],'color',[0 0 0]);
    hca.YLabel.String = {'\phi','(V)'};  
    hca.YLabel.Interpreter = 'tex';
  end
  if 1 % F
    hca = h(isub); isub = isub + 1;
    irf_spectrogram(hca,Fspecrec,'lin');
    hca.YLabel.String = {'v_{||}','(10^3 km/s)'};
    edi_color = [1 1 1];
    vph_color = [1 1 1];
    h_all_markings = [];
    if 1 % EDI energies
      hold(hca,'on')
      %lines_EDI_plus = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'color',edi_color);
      lines_EDI_minus = irf_plot(hca,irf.ts_scalar(phi.time([1 end]),-[v_edi_plus v_edi_minus;v_edi_plus, v_edi_minus]*1e-6),'color',edi_color);
      hleg_EDI = irf_legend(hca,{'- EDI'},[0.08 0.99],'color',lines_EDI_minus(1).Color);  
      hold(hca,'off')
      %h_all_markings = [h_all_markings; lines_EDI_plus; lines_EDI_minus; hleg_EDI];
      h_all_markings = [h_all_markings; lines_EDI_minus; hleg_EDI]; 
    end
    if 1 % model phase velocity
      hold(hca,'on')
      line_color = [0.5 0.5 0.5]; %line_color = mms_colors('matlab');
      lines_vphav = irf_plot(hca,tsVph*1e-6,'LineWidth',1.0,'Color',vph_color,'LineStyle','--');
      %lines_vphav = irf_plot(hca,tsVph*1e-6,'--k');
      hleg_vphav = irf_legend(hca,{'-- v_{ph,mod}'},[0.32 0.99],'color',lines_vphav(1).Color);  
      hold(hca,'off')
      h_all_markings = [h_all_markings; lines_vphav; hleg_vphav];
    end
    if 1 % observed phase velocity
      hold(hca,'on')
      lines_vphobs = irf_plot(hca,obs_vph*1e-3,'*k','LineWidth',1.5,'Color',vph_color);
      lines_vphobs.MarkerSize = 4;
      hleg_vphobs = irf_legend(hca,{'* v_{ph,obs}'},[0.18 0.99],'color',lines_vphobs(1).Color);  
      hold(hca,'off')
      h_all_markings = [h_all_markings; lines_vphobs; hleg_vphobs];
    end  
    colormap(hca,cn.cmap('white_blue'))
    if 0 % str info
    str_info = {'unperturbed f:';...
      ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      sprintf('beta_{Schamel}=%g',beta);...
      };
    irf_legend(hca,str_info,[1.01 1.4],'color',hlines(1).Color);    
    end
    hca.YLim = vlim_f*[-1 1];
  end
  if 1 % diff E (Poisson) (density)
    hca = h(isub); isub = isub + 1;  
    nscale = 1e-3;
    irf_plot(hca,{-1*tsDnFromPhi/nscale,-1*(tsDnModel-ntot*1e-6)/nscale},'comp')
    hca.YLabel.String = {'\delta n',sprintf('(10^{%.0f} cm^{-3})',log10(nscale))};
    hca.YLabel.Interpreter = 'tex';
    %irf_legend(hca,{'n_i - (\epsilon_0/e)\partial_{||}E_{||}  ','  \int f_{model}dv_{||}'},[0.01 0.10]);
    %irf_legend(hca,{'n_{et}^{obs}  ','  n_{et}^{mod}'},[0.01 0.10]);
    irf_legend(hca,{'obs.  ','model'},[0.08 0.98]);
    doDoubleAxis = 1; % dn
    if doDoubleAxis  
      ax1 = hca;
      ax2 = axes('Position',get(ax1,'Position'));
      ax2.YLim = ax1.YLim/nscale/n0;    
      set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none','box','off'); % color of axis      
      %set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
      %set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
      irf_timeaxis(ax2,'nolabels')
      ax2.XLabel.String = [];
      ax2.YLabel.String = {'\delta n/n'};
      ax2.YLabel.Interpreter = 'tex';    
      ax2.YTick = hca.YTick*1e-3/n0;  
      h_all = [h_all,ax2];
      ax2_dn = ax2;
      ax1_dn = ax1;
      ax2.YTick = ax1.YTick;
    end      
  end  
  if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 180
    hca = h(isub); isub = isub + 1;
    %%        
    f_scale = 1e6;
    colors = mms_colors('matlab');
    %colors = colors([2 1],:)
    %hca = subplot(1,1,1);
    
    irf_plot(hca,ts_edi_flux180/f_scale,'color',colors(1,:));
    ax1 = hca;
    ax2 = axes('Position',get(ax1,'Position'));
    irf_plot(ax2,fluxModel180.resample(ts_edi_flux180)/f_scale,'color',colors(2,:))      
    %irf_plot(ax2,fluxModel180/f_scale,'color',colors(2,:))      
    set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required
    set(ax2,'YAxisLocation','right');
    set(ax2,'Color','none','box','off'); % color of axis      
    set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
    set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
    irf_timeaxis(ax2,'nolabels')
    
    ax2.YLabel.String = {'j^{mod}',sprintf('(10^%g cm^{-2}s^{-1})',log10(f_scale))};
    ax2.YLabel.Interpreter = 'tex';
    ax1.YLabel.String = {'j^{EDI}',sprintf('(10^%g cm^{-2}s^{-1}sr^{-1})',log10(f_scale))};
    ax1.YLabel.Interpreter = 'tex';
    irf_legend(hca,{'EDI','model'},[0.08 0.98])
    %irf_legend(hca,{'EDI: \theta = [168.75 180]^o','model: v = [-13600 -12900] km/s'},[0.08 0.98])
    %text(hca,0.002,0.99*hca.YLim(2),'180^o','verticalalignment','top')    
    ax1.YLim = [0 3.7];
    ax2.YLim = [0 3.7];
    h_all = [h_all,ax2];
    ax2_flux = ax2;
    ax1_flux = ax1;    
  end
  if 0 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0
    hca = h(isub); isub = isub + 1;
    %%        
    f_scale = 1e6;
    colors = mms_colors('matlab');
    %colors = colors([2 1],:)
    %hca = subplot(1,1,1);
    if 0
      %irf_plot({fluxModel180,ts_edi_flux180},'comp')  
    else
      irf_plot(hca,ts_edi_flux0/f_scale,'color',colors(1,:));
      ax1 = hca;
      ax2 = axes('Position',get(ax1,'Position'));
      irf_plot(ax2,fluxModel0/f_scale,'color',colors(2,:))      
      set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required      
      set(ax2,'YAxisLocation','right');
      set(ax2,'Color','none','box','off'); % color of axis      
      set(ax2,'XColor','k','YColor',colors(2,:)); % color of axis lines and numbers
      set(ax1,'XColor','k','YColor',colors(1,:)); % color of axis lines and numbers
      irf_timeaxis(ax2,'nolabels')
    end
    ax2.YLabel.String = {'flux 0^o, Model',sprintf('(10^%g cm^{-2}s^{-1})',log10(f_scale))};
    ax2.YLabel.Interpreter = 'tex';
    ax1.YLabel.String = {'flux 0^o, EDI',sprintf('(10^%g cm^{-2}s^{-1}sr^{-1})',log10(f_scale))};
    ax1.YLabel.Interpreter = 'tex';
    irf_legend(hca,{'EDI','model'},[0.02 0.98])
    %text(hca,0.002,0.99*hca.YLim(2),'180^o','verticalalignment','top')    
    ax1.YLim = [0 4];
    ax2.YLim = [0 4];
    h_all = [h_all,ax2];
  end

  if 0 % sum FVdv, bulk velocity
    hca = h(isub); isub = isub + 1;
    plot(hca,x_vec,sumFVdv./sumFdv*1e-3);  
    hca.YLabel.String = {'v','(km/s)'};
  end   
    
  irf_plot_axis_align(h)
  
  h(3).YLabel.String = {'v_{||}','(10^3 km/s)'};
  h(3).YLabel.Interpreter = 'tex';
  
  iref = 1;
  ijmp = 1;
  %h_all(iref+4+ijmp).Position = h_all(iref+4).Position;  
  
  
  %h_all(iref+7).Position = h_all(iref+5).Position;
  
  %irf_zoom(h_all,'x',[tint_phi(1) ts_edi_flux180.time(end)])  % tint_phi
  irf_zoom(h_all,'x',[tint_phi(1) tint_phi(end)])  % tint_phi
  
  h_all(iref+4+ijmp).XLabel = [];
  
  h(1).YLim = [-70 70];
  h(4).YLim = [-0.9 1.7];
  
  ax2_flux.Position = ax1_flux.Position;
  ax2_dn.Position = ax1_dn.Position;
  ax2_dn.YTick = ax1_dn.YTick/1e-3/n0;
  ax2_dn.YLim = ax1_dn.YLim/1e-3/n0;

  
  h(3).YLabel.Position(1) = h(1).YLabel.Position(1);
  
  %h_all(iref+7).XLabel = [];
  %h_all(end-2).YAxisLocation = 'right';
  %h_all(end-3).YAxisLocation = 'right';
  legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
  legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
  for ipanel = 1:npanels
    irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  end
  c_eval('h_all(?).XGrid = ''off''; h_all(?).YGrid = ''off'';',1:numel(h_all))
  
  % averaged distributions
  isub = 1;
  
  mod_f_average = mean(Fabel_obs,1);
  mod_f0 = f0(v_vec,n,vd,vt);
  edist = ePDist1.tlim(tint_phi);
  vg = -40e3:1000:40e3; % km/s
  lowerelim = 50;
  ef1D = edist.reduce('1D',dmpaB1,'vg',vg,'nMC',1000,'lowerelim',lowerelim);
  v_fpi = ef1D.depend{1}(1,:);
  f_fpi = mean(ef1D.data,1);
  f_fpi = ef1D.data;
  if 1 % F, 
    hca = h2(isub); isub = isub + 1;
    v_scale = 1e-3;
    hlines = plot(hca,v_vec*v_scale*1e-3,mod_f0,v_vec*v_scale*1e-3,mod_f_average,'linewidth',1.5);
    hold(hca,'on')
    hlines = plot(hca,v_fpi*v_scale,1*f_fpi*1e0,'--');
    hlines = plot(hca,v_fpi*v_scale,1*mean(f_fpi,1)*1e0,'-','linewidth',1.5);
    hold(hca,'off')
    hca.YLabel.String = {'f','(s^1/m^4)'};
    hca.XLabel.String = {'v','(10^3 km/s)'};
    hca.XLim = [-40 40];
    str_lines = {'f_{mod,\phi=0}';'f_{mod}';'-- fpi';'-- fpi';'-- fpi';'-- fpi';'- mean fpi'};
    %legend(hlines,str_lines)
    irf_legend(hca,str_lines,[0.99 0.99])
    str_info = {['T_{0}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{0}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,0}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      ...sprintf('beta_{Schamel}=%g',beta);...
      };
    set(hca,'ColorOrder',zeros(10,3))
    irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];
       if 1
         plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
         irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
       end
      hold(hca,'off')
    end
    if 1 % vph
      hold(hca,'on')
      plot(hca,vph*1e-6*[1 1],hca.YLim,'linewidth',1.0)
      hleg = irf_legend(hca,'vph',[0.51+vph*1e-6/hca.XLim(2)/2 0.02],[0 0 0]);     
      hold(hca,'off')
    end
  end  
  delete(h2(2))
  %cn.print(sprintf('%g_bgk_tsav_phi_mult_%.1f_scPot_corr_vph%.0f_n%g',iff,phi_mult,vph*1e-3,n0))
end
if 0 % 0 % plot, time-averaged distributions, for paper
  %%
  clear h;
  figure(40)
  h = subplot(1,1,1);   
  isub = 1;
  
  mod_f_average = mean(Fabel_obs,1);
  mod_f0 = f0(v_vec,n,vd,vt);
  %edist = ePDist1.tlim(tint_phi);
  edist = eDist_bgremoved.tlim(tint_phi);
  
  vg = -40e3:1000:40e3; % km/s
  lowerelim = 50;
  ef1D = edist.reduce('1D',dmpaB1,'vg',vg,'nMC',1000,'lowerelim',lowerelim);
  v_fpi = ef1D.depend{1}(1,:);
  f_fpi = mean(ef1D.data,1);
  f_fpi = ef1D.data;
  if 1 % F, 
    hca = h(isub); isub = isub + 1;
    v_scale = 1e-3;
    hlines = plot(hca,v_vec*v_scale*1e-3,mod_f0,v_vec*v_scale*1e-3,mod_f_average,'linewidth',1.5);
    hold(hca,'on')
    hlines = plot(hca,v_fpi*v_scale,1*f_fpi*1e0,'--');
    hlines = plot(hca,v_fpi*v_scale,1*mean(f_fpi,1)*1e0,'-','linewidth',1.5);
    hold(hca,'off')
    hca.YLabel.String = {'f','(s^1/m^4)'};
    hca.XLabel.String = {'v','(10^3 km/s)'};
    hca.XLim = [-40 40];
    str_lines = {'f_{mod,\phi=0}';'f_{mod}';'-- fpi';'-- fpi';'-- fpi';'-- fpi';'- mean fpi'};
    %legend(hlines,str_lines)
    irf_legend(hca,str_lines,[0.99 0.99])
    str_info = {['T_{0}= [' sprintf('%g  ',T) '] eV'];...
      ['n_{0}= [' sprintf('%g  ',n*1e-6) '] cc'];...
      ['v_{d,0}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...
      ...sprintf('beta_{Schamel}=%g',beta);...
      };
    set(hca,'ColorOrder',zeros(10,3))
    irf_legend(hca,str_info,[0.01 0.99],[0 0 0]);   
    if 1 % EDI velocities                
      hold(hca,'on')
      all_edi_plusminus = [v_edi_minus;  v_edi_plus;...
                 -v_edi_minus; -v_edi_plus]*[1 1];
       if 1
         plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
         irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])
       end
      hold(hca,'off')
    end
    if 1 % vph
      hold(hca,'on')
      plot(hca,vph*1e-6*[1 1],hca.YLim,'linewidth',1.0)
      hleg = irf_legend(hca,'vph',[0.51+vph*1e-6/hca.XLim(2)/2 0.02],[0 0 0]);     
      hold(hca,'off')
    end
  end  
end

%% Find f at center of EH.
ff = Fspecrec.p;
vv = Fspecrec.f;
vv_ind = find(abs(vv--9)==min(abs(vv-9)));
ff_ehcenter = ff(:,vv_ind);

%% Plot, compare dntrap/dphi
if 0
  figure(98)
  nrows = 2;
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

  ih_vec = 1:numel(phi_all);
  ih_vec(obs_potential(:,mms_id) > 200) = [];

  if 1 % vdntrap vs phi
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    for ih = ih_vec
      plot(hca,phi_all{ih},dntrap_all{ih})
    end
    hold(hca,'off')
  end
  if 1 % vdntrap vs phi, fit
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    for ih = ih_vec
      plot(hca,phi_all{ih},fun_fit_all{ih}(phi_all{ih}))
    end
    hold(hca,'off')
  end
  if 1 % vdntrap vs phi, fit
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    for ih = ih_vec
      plot(hca,phi_all{ih},torow(fun_fit_all{ih}(phi_all{ih}))-torow(dntrap_all{ih}));
    end
    hold(hca,'off')
  end
  if 1 % vdntrap vs phi
    hca = h(isub); isub = isub + 1;
    plot(hca,phi_long,dntrap_long)  
  end
  if 1 % vdntrap vs phi
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    for ih = ih_vec
      plot(hca,phi_all{ih},dntrap_all{ih},'.')
    end
    plot(hca,phi_long,dntrap_long) 
    hold(hca,'off')
  end
  if 1 % vdntrap vs phi
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    for ih = ih_vec
      plot(hca,phi_all{ih},fun_fit_all{ih}(phi_all{ih}),'.')
    end
    plot(hca,phi_long,fun_net_long(phi_long)) 
    hold(hca,'off')
  end
end

%% Collect data in cell-matrices, for comparison later on
if 0
load('/Users/cno062/MATLAB/cn-matlab/liouville/liouville_eh_abel_parameter_study.mat')
flux_all{i_vph,i_phi_mult,i_iff} = fluxModel180/f_scale;
psd_all{i_vph,i_phi_mult,i_iff} = {{v_vec*v_scale*1e-3,f0(v_vec,n,vd,vt)},{v_vec*v_scale*1e-3,mean(Fabel_obs,1)},{v_fpi*v_scale,1*f_fpi*1e0}};
n_all{i_vph,i_phi_mult,i_iff} = {-1*tsDnFromPhi/nscale,-1*(tsDnModel-ntot*1e-6)/nscale};
save('/Users/cno062/MATLAB/cn-matlab/liouville/liouville_eh_abel_parameter_study.mat','flux_all','psd_all','n_all')
%pause
end
end
end
end
