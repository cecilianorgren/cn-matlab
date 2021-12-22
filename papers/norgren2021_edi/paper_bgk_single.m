% old version is eh_models_external_functions
% load data
mms_id = 1;
ih = 11;
% load ESW data
manual = edi_event_manual_dt;
vph_orig = manual(ih).vpar;
vph = 100*round(vph_orig/100)*1e3; % m/s, rounded of
vph = vph_orig*1e3;
vph = -8600e3;
lx = abs(manual(ih).lpp(mms_id)/2*1e3);
nx = 10001;0
tref = EpochTT(manual(ih).t_ref); % mms 1
tref = tref + manual(ih).dt(mms_id); % adjust in case any other than mms1 is chosen

%c_eval('tint = obs_t0_epoch_mms?(ih) + 5*lx/abs(vph)*[-1 1];',mms_id)

tint = tref + 6*lx/abs(vph)*[-1 1]; % time interval for figure
% Manual
tint = irf.tint('2017-07-06T13:54:05.5836Z/2017-07-06T13:54:05.5873Z');
tint = irf.tint('2017-07-06T13:54:05.5837Z/2017-07-06T13:54:05.5873Z');
T = tint(2)-tint(1);
tref = tint(1)+0.5*T; % center of time interval

% time string for title
tint_utc = tint.utc;
tint_utc = tint_utc(:,12:24);

%t0 = tint(1); % beginning of time interval
%tcenter = mean(tint-tint(1)); % (relative) center of time interval

% Set up grid and phi
units = irf_units;

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

% Electron distribution at phi = 0, f0
%iff = 20; % must be the same as for the bgk timeseries

[f0,params] = mms_20170706_135303.get_f0;
iff = params.inp;
n = params.n;  
ntot = sum(n);
n0 = ntot;
R = n(1)/ntot;
T = params.T;
vd = params.vd;
vt = params.vt;
str_info = {'unperturbed f:';...
            ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
            ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
            ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
            };

% Integrate electric field to get potential
c_eval('Epar_obs = gseE?par.tlim(tint);',mms_id)

phi_obs = irf_integrate(Epar_obs);
t_obs = phi_obs.time - tref; % relative time vector, 0 at center of time interval

% edi, extract flux at 0 and 180 deg pitch angle
c_eval('flux180 = ePitch?_flux_edi.palim([168.75 180]).tlim(tint+[-0.001 0.001]);',mms_id)
c_eval('flux0 = ePitch?_flux_edi.palim([0 11.25]).tlim(tint+[-0.001 0.001]);',mms_id)

% Time vector
t_edi_180 = flux180.time - tref;
t_edi_0   = flux0.time - tref;

% Length vector
x_edi_180 = t_edi_180'*vph;    
x_edi_0 = t_edi_0'*vph;    
x_obs = t_obs'*vph; % from E/phi   
x_mod = x_obs;

% Potentially adjust potential
phi_obs_nodetrend = phi_obs.data*vph*1e-3;
phi_mod = fun_phi(obs_potential(ih,mms_id),x_mod,lx);

phi_obs_detrend = detrend(phi_obs_nodetrend)';
phi_obs_detrend = phi_obs_detrend - min(phi_obs_detrend);

phi_obs = phi_obs_detrend;

nx_obs = numel(x_obs);
nx_mod = numel(x_mod);
 

% Get electric field and n
dx_obs = x_obs(2) - x_obs(1);
x_obs_diff1 = x_obs(1:end-1) + 0.5*dx_obs;
x_obs_diff2 = x_obs(2:end-1) + dx_obs;
Efield_obs = -diff(phi_obs)/dx_obs;
c_eval('Efield_obs = gseE?par.tlim(tint).data;',mms_id)
dn_obs = diff(phi_obs,2)/dx_obs/dx_obs*units.eps0/units.e;
dn_obs = [dn_obs(1); tocolumn(dn_obs); dn_obs(end)]; % assume phi -> at edges
dn_obs = dn_obs*1;

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

% Reconstruct phase space
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

vv_ind = find(abs(V_obs(1,:)-vph)==min(abs(V_obs(1,:)-vph)));
F_center = Fabel_obs(:,vv_ind); % to check f at center of EH


% Plot
if 1 % plot 3
  %%
  fig = figure(93);
  nrows = 3;
  ncols = 2;
  npanels = nrows*ncols;
  isub = 0;
  if 0
    for icol = 1:ncols
      for irow = 1:nrows  
        isub = isub + 1;         
        h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
      end
    end
  else
    for ipanel = 1:npanels
      h(ipanel) = subplot(nrows,ncols,ipanel);
    end
  end
  isub = 1;

  all_hcb = {};
  ihcb = 1;
  vlim = 30000e3;
  shiftcb = 0.01;
  xlim = [-15 15];
  xlim = max(x_obs)*1e-3*[-1 1];
  if 1 % phi(x)
    hca = h(isub); isub = isub + 1;
    plot(hca,x_obs*1e-3,phi_obs)
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = '\phi (V)';
    %legend(hca,{'obs','mod, Gaussian','obs no detrend'},'Box','off')
    %irf_legend(hca,{'obs';{'mod:','Gaussian'};{'obs:','no detrend'}},[0.02,0.98])
    
    %hca.Title.String = sprintf('%s - %s',tint_utc(1,:),tint_utc(2,:));   
    hvph = irf_legend(hca,{sprintf('v_{ph}= %.0f km/s',round(vph*1e-5)*1e5*1e-3)},[0.5,1.01],'Color',[0 0 0],'fontsize',10);
    hvph.HorizontalAlignment = 'center';
    hvph.VerticalAlignment = 'bottom';
    %hca.Title.String = sprintf('v_{ph}= %.0f km/s',vph*1e-3);    
    
    hca.XLim = xlim;
   % hca.Title.Position(1)=24;
  end
  if 1 % E
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e)
    %contourf(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    hcb.YLabel.String = 'U/e (eV)';  
    hca.YLim = vlim*[-1 1]*1e-6;
    hca.CLim = max(abs(E_obs(:)/units.e))*[-1 1];
    hca.CLim = 0.9e3*[-1 1];
    cmap = cn.cmap('blue_red'); 
    cmap_ = cmap([1:96 (end-95):end],:);
    colormap(hca,cmap_)
     if 1 % E contours
      hold(hca,'on')
      %levels_E = [linspace(min(E_obs(:)),0,6) linspace(0,max(E_obs(:)),50)]; 
      levels_E = linspace(min(E_obs(:)),0,5); 
      levels_E = [levels_E 0:levels_E(2)-levels_E(1):max(E_obs(:))];                   
      [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,levels_E/units.e,'k');
      [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,[0 0],'k','LineWidth',1.0);
      hold(hca,'off')
    end
    hcb.YLim = [min(E_obs(:)) hca.CLim(2)*units.e]/units.e;
    hca.XLim = xlim;
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    %hca.Title.String = sprintf('Schamel, beta_{obs}=%.2f',beta_obs);
  end
  if 0 % Ff obs, Abel
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
  if 1 % F obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,Fabel_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    hcb.YLabel.String = 'f_e^{mod} (sm^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'));
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    hcb.Position(1) = hcb.Position(1)-shiftcb;
    %hca.Title.String = 'Abel obs';    
    hca.XLim = xlim;
    if 1 % E contours
      hold(hca,'on')   
      %[hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,levels_E/units.e,'k');
      [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,[0 0],'color',0*[0.8 .8 .8],'linewidth',1.0,'linestyle','-');
      hold(hca,'off')
      hca.CLim = [0 0.0023];
    end
  end  
  if 0 % free,, and total densities
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
                      x_obs*1e-3,nansum(Fabel_obs,2)*dv*1e-6,...
                      x_obs*1e-3,nansum(Fabel_free_obs,2)*dv*1e-6,...
                      x_obs*1e-3,nansum(Fabel_trap_obs,2)*dv*1e-6...
                      );
    %hlines(1).LineWidth = 1.5;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'n_e (cm^{-3})';
    
    %hca.Title.String = 'n_t = n_0+(\epsilon_0/e)\nabla^2\phi-n_f';
    irf_legend(hca,{...
      'n_{e}^{obs}';...%'n_{0}+(\epsilon_0/e)\nabla^2\phi';...
      'n_{e}^{mod}';...
      'n_{ef}';...
      'n_{et}'},...
      [0.02 0.25])    
%                     'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'
    %irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
    hca.YLim = [0 0.042];
    hca.XLim = xlim;
  end
  if 1 % free, trapped, and total densities, 2 yscales
    hca = h(isub); isub = isub + 1;
    colors = pic_colors('matlab');
    ploty1 = [(n0+dn_obs)*1e-6, nansum(Fabel_obs,2)*dv*1e-6];
    ploty2 = [nansum(Fabel_free_obs,2)*dv*1e-6, nansum(Fabel_trap_obs,2)*dv*1e-6];
    [AX,H1,H2] = plotyy(hca,x_obs*1e-3,ploty2',x_obs*1e-3,ploty1');
    H1(1).Color = colors(3,:);
    H1(2).Color = colors(4,:);
    H2(1).Color = colors(1,:);
    H2(2).Color = colors(2,:);
    %hlines(1).LineWidth = 1.5;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'n_{ef}, n_{ep} (cm^{-3})';
    AX(2).YLabel.String = 'n_{e}^{obs}, n_{e}^{mod} (cm^{-3})';
    %hca.Title.String = 'n_t = n_0+(\epsilon_0/e)\nabla^2\phi-n_f';
    irf_legend(hca,{...
      'n_{e}^{obs}';...%'n_{0}+(\epsilon_0/e)\nabla^2\phi';...
      'n_{e}^{mod}';...
      'n_{ep}';...
      'n_{et}'},...
      [0.02 0.25])    
    AX(2).YLim = [0.038 0.042];
    %ax(1).YTick = 0:0.01:0.04;
    %ax(2).YTick = 0:0.01:0.04;
    AX(1).YTick = AX(1).YTick(1:2:end);
    AX(2).YTick = 0.038:0.001:0.042;
%                     'n_{t,mod}';'n_{t,Abel,mod}';'n_{t,Scha,mod}'
    %irf_legend(hca,{sprintf('beta_{obs}=%.2f',beta_obs);sprintf('beta_{mod}=%.2f',beta_mod)},[0.02 0.98],'color',[0 0 0])
    %hca.YLim = [0 0.042];
    hca.XLim = xlim;
    AX(2).XLim = xlim;
  end
  if 1 % FV obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,FVabel_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    hcb.YLabel.String = 'vf_e^{mod} (m^{-3})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    hca.CLim = max(abs(FVabel_obs(:)))*[-1 1];    
    colormap(hca,cn.cmap('blue_red'))  
    if 1 % EDI energies
      hold(hca,'on')
      if 1 % solid lines                
        hlines = plot(hca,...
          x_obs([1 end]),(v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(v_edi_minus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_minus)*1e-6*[1 1],...
          'LineWidth',0.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = '-'; hlines(iline).Color = [0 0 0]; end
       % irf_legend(hca,{'- EDI'},[0.01 0.85],'color',hlines(1).Color);   
      else % dashed lines
        hlines = plot(hca,x_obs([1 end]),v_edi*1e-6*[1 1],x_obs([1 end]),-v_edi*1e-6*[1 1],'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = '--'; hlines(iline).Color = [0 0 0]; end  
        hlines = plot(hca,...
          x_obs([1 end]),(v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(v_edi_minus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_plus)*1e-6*[1 1],...
          x_obs([1 end]),(-v_edi_minus)*1e-6*[1 1],...
          'LineWidth',1.5);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = ':'; hlines(iline).Color = [0 0 0]; end
        irf_legend(hca,{'-- EDI'},[0.01 0.85],'color',hlines(1).Color);   
      end
      hold(hca,'off')
    end
    hca.XLim = xlim;
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    hcb.Position(1) = hcb.Position(1)-shiftcb;
    if 1 % E contours
      hold(hca,'on')   
      %[hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,levels_E/units.e,'k');
      [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,[0 0],'color',0*[0.8 .8 .8],'linewidth',1.0,'linestyle','-');
      hold(hca,'off')
      %hca.CLim = [0 0.0023];
    end
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
    colors = [colors(2,:); colors(1,:)];
    nlines = 2;
    c_eval('ax(1).Children(?).Color = colors(1,:);',1:nlines)
    c_eval('ax(2).Children(?).Color = colors(2,:);',1:nlines)    
    c_eval('ax(?).YColor = colors(?,:);',1:2)
    c_eval('ax(2).Children(?).Marker = ''*'';',1:nlines)    
    
        
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = sprintf('j_e^{mod} (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
    ax(2).YLabel.String = sprintf('j_e^{EDI} (10^%g cm^{-2}s^{-1}sr^{-1})',log10(units_scale_2));
    %irf_legend(hca,{'Model';'* EDI'},[0.02 0.3])
    text(hca,0.5*hca.XLim(2),0.2*hca.YLim(2),'0^o','verticalalignment','bottom')
    text(hca,0.5*hca.XLim(2),1.0*hca.YLim(2),'180^o','verticalalignment','top')
    hca.YLim(1) = 0;    
    %hca.XLim = xlim;
    c_eval('ax(?).YLim = [0 2.5];',1:2)    
    c_eval('ax(?).XLim = xlim;',1:2)    
    
    ax(1).YAxisLocation = 'right',
    ax(2).YAxisLocation = 'left';
    ax(1).YTick = 0:1:100;
    ax(2).YTick = 0:1:100;
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
  
  legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
  legends_color = {'k','w','k','k','k','k','k','k','k','k','k','k'};
  for ipanel = 1:npanels
    irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  end

  %all_hcb{1}.FontSize = 10;
  h(1).YLim = [0 320];
  h(1).Title.String = sprintf('Time interval: %s - %s',tint_utc(1,:),tint_utc(2,7:end));   
  h(1).Title.Position = [23, 377, 0];
  %hca.Title.Position;
  annotation('textarrow',[0.27 0.27],[0.55 0.525],'string',{'trapped/passing','boundary'},'fontsize',11,'horizontalalignment','center');
  annotation('textarrow',[0.25 0.25],[0.29 0.27],'string',{'EDI range'},'fontsize',11,'horizontalalignment','center');
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
  if 0
    %%
    AX(2).YLabel.String = 'n_e^{obs}, n_e^{mod} (cm^{-3})';
    AX(1).YLabel.String = 'n_{ep}, n_{et} (cm^{-3})';
    h(1).YLim = [0 370];
    annotation('textarrow',[0.27 0.27],[0.55 0.525],'string',{'trapped/passing','boundary'},'fontsize',12,'horizontalalignment','center');
    annotation('textarrow',[0.25 0.25],[0.29 0.27],'string',{'EDI range'},'fontsize',12,'horizontalalignment','center');
  end
end

