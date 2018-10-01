%% Set up grid and phi
units = irf_units;
fun_phi = @(phimax,x,lx) phimax*exp(-x.^2/2/lx.^2);

% EDI energy and corresponding velocity
E_edi = 500; % eV
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

dntrap_all = cell(numel(obs_velocity),1);
phi_all = cell(numel(obs_velocity),1);
x_all = cell(numel(obs_velocity),1);
fun_fit_all = cell(numel(obs_velocity),1);


% F0
n = [0.02 0.02]*1e6; % m-3
ntot = 0.04*1e6;
R = 0.7;
n = [R (1-R)]*ntot;
T = [150 4000]; % eV
vd = [-11000e3 4000e3]; % ms-1 
vt = sqrt(2*units.e*T./units.me); % m/s
n0 = sum(n);

str_info = {'unperturbed f:';...
            ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
            ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
            ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
            };
               
beta_range = [-2.5 0];
for ih = 10%:numel(obs_velocity)
phi_input = 3;
switch phi_input
  case 1 % Gaussian of choice
    lx = 1300; 
    nx = 207;
    x = linspace(-5*lx,5*lx,nx);
    phimax = 100;
    phi = fun_phi(phimax,x,lx);
    vph = -10000e3;
    x_obs = x;
    phi_obs = phi;
  case 2 % single one from observations, Gaussian from lx and phimax
    mms_id = 1;
        
    phimax = obs_potential_max(ih);
    lx = obs_lpp(ih)/2*1e3;
    nx = 10001;
    
    vph = obs_velocity(ih)*1e3;        
    c_eval('tint = obs_t0_epoch_mms?(ih) + 5*lx/abs(vph)*[-1 1];',mms_id)
    t0 = tint(1);
    tcenter = mean(tint-tint(1));
    c_eval('Epar_obs = gseE?par.tlim(tint);',mms_id)
    if isempty(Epar_obs); continue; end
    phi_obs = irf_integrate(Epar_obs);
    t_obs = phi_obs.time-t0;
    
    % edi
    c_eval('flux180 = flux180_mms?.tlim(tint+[-0.001 0.001]);',mms_id)
    t_edi_180 = flux180.time-t0;        
    t_edi_180 = t_edi_180 - tcenter;
    c_eval('flux0 = flux0_mms?.tlim(tint+[-0.001 0.001]);',mms_id)
    t_edi_0 = flux0.time-t0;        
    t_edi_0 = t_edi_0 - tcenter;
    
    
    t_obs = t_obs - tcenter;
    
    tint_utc = tint.utc;
    tint_utc = tint_utc(:,12:24);
    
    x_edi_180 = (t_edi_180-mean(t_obs))'*vph;    
    x_edi_0 = (t_edi_0-mean(t_obs))'*vph;    
    x_obs = t_obs'*vph;    
    x_mod = x_obs; % x_obs = linspace(x_obs(1),x_obs(end),nx); %x = linspace(-5*lx,5*lx,nx);
    
    phi_obs_nodetrend = phi_obs.data*vph*1e-3;
    phi_mod = fun_phi(obs_potential(ih,mms_id),x_mod,lx);
    
    phi_obs_detrend = detrend(phi_obs_nodetrend)';
    phi_obs_detrend = phi_obs_detrend - min(phi_obs_detrend);
    
    phi_obs = phi_obs_detrend;
          
    nx_obs = numel(x_obs);
    nx_mod = numel(x_mod);
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
    x = x_vec;

    % Remove phi baselevel
    %c_eval('phi_baselevel = interp_linear_piecewise(intEdt?.data,x_vec,x_vec(LOCS?));')
    c_eval('ts_locs? = irf.ts_scalar(intEdt?.time([1; LOCS?; end]),intEdt?.data([1; LOCS?; end]));')
    %c_eval('phi_baselevel? = interp_linear_piecewise(intEdt?.data([1; LOCS?; end]),x_vec([1; LOCS?; end]),x_vec);')
    %c_eval('ts_phi_baselevel? = irf.ts_scalar(intEdt?.time,phi_baselevel?);')
    c_eval('ts_phi_baselevel? = ts_locs?.resample(intEdt?);')
    c_eval('intEdt?_detrend = intEdt?-ts_phi_baselevel?;')

    vph = -9000e3; % m/s, representative phase velocity
    phi_shift = 00; % to keep potential > 0
    phi_scaling = 1.0; % if we underestimate phi    

    % Potential from observed E
    mms_id = 1;
    c_eval('phi_timeline = intEdt?_detrend.time;',mms_id)
    c_eval('phi?_detrend = intEdt?_detrend*vph*1e-3*phi_scaling;',mms_id)
    c_eval('phi?_detrend_shift = phi?_detrend + phi_shift;',mms_id)
    c_eval('phi?_detrend_shift.data(phi?_detrend_shift.data<0) = 0;',mms_id)                
    c_eval('phi_vec = phi?_detrend_shift.data;',mms_id)
    phi = phi_vec;
    
    c_eval('epar_vec = Etoint?.data;',mms_id)
    
    % charge density from observed phi        
    obs_density_diff = -diff(epar_vec*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
    %obs_density_diff_nofilt = -diff(epar_vec_nofilt*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
    %obs_density = mod_density_average + obs_density_diff; % ni + ne - ni = ne,  assuming ion density is unperturbed  
    ts_obs_density = irf.ts_scalar(phi_timeline(1:end-1) + 0.5*(phi_timeline(2)-phi_timeline(1)),obs_density); 
  case 4 % single from observations, but taken from the timeseries obtained for all, using time and rescaling vph from konrad
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

    % Remove phi baselevel
    c_eval('ts_locs? = irf.ts_scalar(intEdt?.time([1; LOCS?; end]),intEdt?.data([1; LOCS?; end]));')
    c_eval('ts_phi_baselevel? = ts_locs?.resample(intEdt?);')
    c_eval('intEdt?_detrend = intEdt?-ts_phi_baselevel?;')

    vph_avg = -9000e3; % m/s, representative phase velocity
    phi_shift = 00; % to keep potential > 0
    phi_scaling = 1.0; % if we underestimate phi    

    % Potential from observed E
    mms_id = 1;
    c_eval('phi_timeline = intEdt?_detrend.time;',mms_id)
    c_eval('phi?_detrend = intEdt?_detrend*vph_avg*1e-3*phi_scaling;',mms_id)
    c_eval('phi?_detrend_shift = phi?_detrend + phi_shift;',mms_id)
    c_eval('phi?_detrend_shift.data(phi?_detrend_shift.data<0) = 0;',mms_id) 
    
    % Get shorter time interval
    phimax = obs_potential_max(ih);
    lx = obs_lpp(ih)/2*1e3;
    nx = 10001;
    
    vph = obs_velocity(ih)*1e3;        
    c_eval('tint_phi_one = obs_t0_epoch_mms?(ih) + 5*lx/abs(vph)*[-1 1];',mms_id)
    t0 = tint_phi_one(1);
    tcenter = mean(tint_phi_one-tint_phi_one(1));
    
    c_eval('phi_vec = phi?_detrend_shift.tlim(tint_phi_one).data*vph/vph_avg;',mms_id)
    phi = phi_vec;
    
    c_eval('epar_vec = Etoint?.tlim(tint_phi_one).data;',mms_id)
    
    % Plotting options
    doT = 1; % otherwise plot x;
    x_vec = intEdt1.tlim(tint_phi_one).time - t0; % seconds
    dx_vec = x_vec(2)-x_vec(1);
    dx = dx_vec; 
    nx = numel(x_vec);
    x_vec_diff1 = x_vec(1:end-1)+0.5*dx_vec;
    x_vec_diff2 = x_vec(2:end-1);
    x = x_vec;
    x_obs = x_vec;
    
    % charge density from observed phi        
    obs_density_diff = -diff(epar_vec*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
    %obs_density_diff_nofilt = -diff(epar_vec_nofilt*1e-3,1)*units.eps0/units.e/(-sign(vph)*dx); % ne-ni
    %obs_density = mod_density_average + obs_density_diff; % ni + ne - ni = ne,  assuming ion density is unperturbed  
    %ts_obs_density = irf.ts_scalar(phi_timeline(1:end-1) + 0.5*(phi_timeline(2)-phi_timeline(1)),obs_density);
    
end

% Get electric field and n
dx_obs = x_obs(2) - x_obs(1);
x_obs_diff1 = x_obs(1:end-1) + 0.5*dx_obs;
x_obs_diff2 = x_obs(2:end-1) + dx_obs;
Efield_obs = -diff(phi_obs)/dx_obs;
c_eval('Efield_obs = gseE?par.tlim(tint).data;',mms_id)
dn_obs = diff(phi_obs,2)/dx_obs/dx_obs*units.eps0/units.e;
dn_obs = [dn_obs(1); tocolumn(dn_obs); dn_obs(end)]; % assume phi -> at edges
%dn_obs = dn_obs.*phi_obs'/max(phi_obs)*2;

dx_mod = x_mod(2) - x_mod(1);
x_mod_diff1 = x_mod(1:end-1) + 0.5*dx_mod;
x_mod_diff2 = x_mod(2:end-1) + dx_mod;
Efield_mod = -diff(phi_mod)/dx_mod;
dn_mod = diff(phi_mod,2)/dx_mod/dx_mod*units.eps0/units.e;
dn_mod = [dn_mod(1); tocolumn(dn_mod); dn_mod(end)]; % assume phi -> at edges


vtrap = sqrt(2*units.e*phimax/units.me);
nv = 2000;
vmax = max([5*vtrap, vd + 2*vt]); 
v = linspace(-vmax,vmax,nv);
dv = v(2) - v(1);

[X_obs,V_obs] = meshgrid(x_obs,v); X_obs = permute(X_obs,[2 1]); V_obs = permute(V_obs,[2 1]);
[X_mod,V_mod] = meshgrid(x_mod,v); X_mod = permute(X_mod,[2 1]); V_mod = permute(V_mod,[2 1]);
PHI_obs = repmat(tocolumn(phi_obs),1,nv);
PHI_mod = repmat(tocolumn(phi_mod),1,nv);
VPH_obs = V_obs*0 + vph;
VPH_mod = V_mod*0 + vph;
E_obs = units.me*(V_obs-vph).^2/2 - units.e*PHI_obs;
E_mod = units.me*(V_mod-vph).^2/2 - units.e*PHI_mod;
ifree_obs = find(E_obs>0); itrap_obs = find(E_obs<=0);
ifree_mod = find(E_mod>0); itrap_mod = find(E_mod<=0);

[Fflat_obs,Fflat_free_obs,Fflat_trap_obs] = get_f_flat(V_obs,n,vt,vd,PHI_obs,VPH_obs);
nfree_obs = nansum(Fflat_free_obs,2)*dv;
ntrap_flat_obs = nansum(Fflat_trap_obs,2)*dv;
ntrap_obs = ntot - torow(nfree_obs) + torow(dn_obs);
dntrap_obs = ntrap_obs - torow(ntrap_flat_obs);
[Fabel_obs,Fabel_free_obs,Fabel_trap_obs] = get_f_abel(V_obs,n,vt,vd,PHI_obs,VPH_obs,dntrap_obs);
Fabel_obs = Fabel_obs + Fflat_trap_obs;              
[Fscha_obs,Fscha_free_obs,Fscha_trap_obs,beta_obs] = get_f_schamel(V_obs,n,vt,vd,PHI_obs,VPH_obs,ntrap_obs,beta_range);

[Fflat_mod,Fflat_free_mod,Fflat_trap_mod] = get_f_flat(V_mod,n,vt,vd,PHI_mod,VPH_mod);
nfree_mod = nansum(Fflat_free_mod,2)*dv;
ntrap_flat_mod = nansum(Fflat_trap_mod,2)*dv;
ntrap_mod = ntot - torow(nfree_mod) + torow(dn_mod);
dntrap_mod = ntrap_mod - torow(ntrap_flat_mod);
[Fabel_mod,Fabel_free_mod,Fabel_trap_mod] = get_f_abel(V_mod,n,vt,vd,PHI_mod,VPH_mod,dntrap_mod);
Fabel_mod = Fabel_mod + Fflat_trap_mod;              
[Fscha_mod,Fscha_free_mod,Fscha_trap_mod,beta_mod] = get_f_schamel(V_mod,n,vt,vd,PHI_mod,VPH_mod,ntrap_mod,beta_range);

FVscha_obs = Fscha_obs.*V_obs;
FVabel_obs = Fabel_obs.*V_obs;
FVabel_mod = Fabel_mod.*V_mod;
    
% flux at edi energy/velocity interval
vind_edi_0 = intersect(find(v>v_edi_minus),find(v<v_edi_plus));
vind_edi_180 = intersect(find(v<-v_edi_minus),find(v>-v_edi_plus));
FVdv_abel_obs_edi_0 = nansum(FVabel_obs(:,vind_edi_0),2)*dv;
FVdv_abel_obs_edi_180 = nansum(FVabel_obs(:,vind_edi_180),2)*dv;
FVdv_scha_obs_edi_0 = nansum(FVscha_obs(:,vind_edi_0),2)*dv;
FVdv_scha_obs_edi_180 = nansum(FVscha_obs(:,vind_edi_180),2)*dv;
FVdv_abel_mod_edi_0 = nansum(FVabel_mod(:,vind_edi_0),2)*dv;
FVdv_abel_mod_edi_180 = nansum(FVabel_mod(:,vind_edi_180),2)*dv;


%% Plot
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
    plot(hca,x_obs*1e-3,phi_obs,x_mod*1e-3,phi_mod,x_obs*1e-3,phi_obs_nodetrend)
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = '\phi (V)';
    %legend(hca,{'obs','mod, Gaussian','obs no detrend'},'Box','off')
    irf_legend(hca,{'obs';{'mod:','Gaussian'};{'obs:','no detrend'}},[0.02,0.98])
    hca.Title.String = sprintf('tint_{obs} = %s - %s',tint_utc(1,:),tint_utc(2,:));
    irf_legend(hca,{sprintf('v_{ph}= %.0f km/s',vph*1e-3)},[0.99,0.99],'Color',[0 0 0])
  end
  if 1 % f0 info 
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
  if 1 % 10^6 cm^{-2}s^{-1}, comparing model flux with flux measured by EDI, at 0 and 180
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
    [fitreslts,gof,fun_net_mod,fun_net_prime_mod] = createFit(torow(phi_mod),torow(dntrap_mod));
    plot(hca,phi_obs,dntrap_obs,'o',phi_obs,fun_net_obs(phi_obs),...
             phi_mod,dntrap_mod,'+',phi_mod,fun_net_mod(phi_mod))
    hca.YLabel.String = 'n_t-n_{t,flat} (m^{-3})';  
    hca.XLabel.String = '\phi (V)';  
    legend(hca,{'n_{t,obs}','n_{t,obs} fit','n_{t,mod}','n_{t,mod} fit'},'Box','off')
  end
  if 0 % dnt/dphi)
    hca = h(isub); isub = isub + 1;
    [fitreslts,gof,fun_net_obs,fun_net_prime_obs] = createFit(torow(phi_obs),torow(dntrap_obs));
    [fitreslts,gof,fun_net_mod,fun_net_prime_mod] = createFit(torow(phi_mod),torow(dntrap_mod));
    plot(hca,phi_obs,fun_net_prime_obs(phi_obs),phi_mod,fun_net_prime_mod(phi_mod))
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
  %cn.print(sprintf('AbelScha_obs_eh%g_flux_mms%g',ih,mms_id))
end
if 1 % plot 3
  figure(93)
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

  all_hcb = {};
  ihcb = 1;
  vlim = 30000e3;
  shiftcb = 0.01;
  xlim = [-15 15];
  if 1 % phi(x)
    hca = h(isub); isub = isub + 1;
    plot(hca,x_obs*1e-3,phi_obs)
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = '\phi (V)';
    %legend(hca,{'obs','mod, Gaussian','obs no detrend'},'Box','off')
    %irf_legend(hca,{'obs';{'mod:','Gaussian'};{'obs:','no detrend'}},[0.02,0.98])
    %hca.Title.String = sprintf('tint_{obs} = %s - %s',tint_utc(1,:),tint_utc(2,:));
    irf_legend(hca,{sprintf('v_{ph}= %.0f km/s',vph*1e-3)},[0.99,0.99],'Color',[0 0 0])
    hca.XLim = xlim;
  end
  if 1 % E
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e)
    %contourf(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'U/e (eV)';  
    hca.YLim = vlim*[-1 1]*1e-6;
    hca.CLim = max(abs(E_obs(:)/units.e))*[-1 1];
    hca.CLim = 0.5e3*[-1 1];
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
    hca.XLim = xlim;
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    %hca.Title.String = sprintf('Schamel, beta_{obs}=%.2f',beta_obs);
  end
  if 1 % Ff obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fabel_free_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'f_e (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'));
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    hcb.Position(1) = hcb.Position(1)-shiftcb;
    %hca.Title.String = 'Abel obs';    
    hca.XLim = xlim;
  end
  if 0 % F obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fabel_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'f_e (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'));
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    hcb.Position(1) = hcb.Position(1)-shiftcb;
    %hca.Title.String = 'Abel obs';    
    hca.XLim = xlim;
  end  
  if 1 % free,, and total densities
    hca = h(isub); isub = isub + 1;
    hlines = plot(hca,x_obs*1e-3,(n0+dn_obs)*1e-6,...
                      x_obs*1e-3,nansum(Fabel_free_obs,2)*dv*1e-6...
                      );
    %hlines(1).LineWidth = 1.5;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'n (cm^{-3})';
    %hca.Title.String = 'n_t = n_0-n_f+(\epsilon_0/e)\nabla^2\phi';
    irf_legend(hca,{'n_{0}+(\epsilon_0/e)\nabla^2\phi';'n_{f}';''},...
                    [0.02 0.3])
%                     'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'
    %irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
    hca.YLim = [0 0.042];
    hca.XLim = xlim;
  end
  if 0 % free, trapped, and total densities
    hca = h(isub); isub = isub + 1;
    hlines = plot(hca,x_obs*1e-3,(n0+dn_obs)*1e-6,...
                      x_obs*1e-3,nansum(Fabel_free_obs,2)*dv*1e-6,...
                      x_obs*1e-3,nansum(Fabel_trap_obs,2)*dv*1e-6...
                      );
    %hlines(1).LineWidth = 1.5;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'n (cm^{-3})';
    hca.Title.String = 'n_t = n_0+(\epsilon_0/e)\nabla^2\phi-n_f';
    irf_legend(hca,{'n_{0}+(\epsilon_0/e)\nabla^2\phi';'n_{f}';'n_{t}'},...
                    [0.02 0.3])    
%                     'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'
    %irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
    hca.YLim = [0 0.042];
    hca.XLim = xlim;
  end
  if 1 % FV obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,FVabel_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v (10^3 km/s)';
    hcb.YLabel.String = 'flux/\Delta v (m^{-3})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    hca.CLim = max(abs(FVabel_obs(:)))*[-1 1];    
    colormap(hca,cn.cmap('blue_red'))  
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
    hca.XLim = xlim;
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    hcb.Position(1) = hcb.Position(1)-shiftcb;
  end
  if 1 % plotyy 10^6 cm^{-2}s^{-1}, 10^6 cm^{-2}s^{-1}sr{-1}, comparing model flux with flux measured by EDI, at 0 and 180
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
    
    ax = plotyy(hca,x_obs*1e-3,[plot_abel_0 plot_abel_180]',x_edi_0*1e-3,[plot_EDI_0 plot_EDI_180]');
    colors = mms_colors('matlab');    
    nlines = 2;
    c_eval('ax(1).Children(?).Color = colors(1,:);',1:nlines)
    c_eval('ax(2).Children(?).Color = colors(2,:);',1:nlines)
    
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
    ax(2).YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1}sr^{-1})',log10(units_scale_2));
    irf_legend(hca,{'Model';'EDI'},[0.02 0.6])
    text(hca,0.5*hca.XLim(2),0.2*hca.YLim(2),'0^o','verticalalignment','bottom')
    text(hca,0.5*hca.XLim(2),1.0*hca.YLim(2),'180^o','verticalalignment','top')
    hca.YLim(1) = 0;    
    hca.XLim = xlim;
  end  
  %cn.print(sprintf('AbelScha_obs_eh%g_flux_mms%g',ih,mms_id))
  width = h(1).Position(3)*0.8;
  for ip = 1:npanels
    h(ip).Position(3) = width;
    h(ip).FontSize = 10;
  end
  all_hcb{1}.Position(1) = h(2).Position(1) + h(2).Position(3) - all_hcb{1}.Position(3);
  all_hcb{2}.Position(1) = h(3).Position(1) + h(3).Position(3) - all_hcb{2}.Position(3);
  all_hcb{3}.Position(1) = h(5).Position(1) + h(5).Position(3) - all_hcb{3}.Position(3);
  
  %all_hcb{1}.FontSize = 10;
  if 0
  h(1).Position(1) = h(1).Position(1) - 0.06;
  h(2).Position(1) = h(2).Position(1) - 0.06;
  h(3).Position(1) = h(3).Position(1) - 0.04;
  h(4).Position(1) = h(4).Position(1) - 0.04;
  h(5).Position(1) = h(5).Position(1) - 0.02;
  h(6).Position(1) = h(6).Position(1) - 0.02;
  width = h(1).Position(3)*0.9;
  for ip = 1:npanels
    h(ip).Position(3) = width;
  end
  end
end
%% Collect for all eh
%[fitreslts,gof,fun_net,fun_net_prime] = createFit(torow(phi),torow(dntrap));
%dntrap_all{ih} = dntrap;
%x_all{ih} = x;
%phi_all{ih} = phi;
%fun_fit_all{ih} = fun_net;

end
