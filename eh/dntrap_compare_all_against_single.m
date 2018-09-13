%% Set up grid and phi
units = irf_units;
tint_all = irf.tint('2017-07-06T13:54:05.490Z/2017-07-06T13:54:05.680Z');

data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
obs_lpp = obs_eh_properties.Lpp; % peak to peak length

    
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
mms_id = 1;
for ih = 1:numel(obs_velocity)
  nlx = 5;
  lx = obs_lpp(ih)/2*1e3;
  vph = obs_velocity(ih)*1e3;        
  c_eval('tint_single = obs_t0_epoch_mms?(ih) + nlx*lx/abs(vph)*[-1 1];',mms_id)
  tint_utc = tint_single.utc;
  c_eval('Epar = gseE?par;',mms_id)
  [phi,phi_progressive,phi_ancillary] = get_phi(Epar,vph,tint_all,tint_single);
  
  phi_vec = phi.data;
  phi_obs = phi_vec;
  x_vec = phi_ancillary.x_vec;
  x_obs = x_vec;
  
  nx_obs = numel(x_obs);
  
  % Get electric field and n
  dx_obs = x_obs(2) - x_obs(1);
  x_obs_diff1 = x_obs(1:end-1) + 0.5*dx_obs;
  x_obs_diff2 = x_obs(2:end-1) + dx_obs;
  Efield_obs = -diff(phi_obs)/dx_obs;
  c_eval('Efield_obs = gseE?par.tlim(tint).data;',mms_id)
  dn_obs = diff(phi_obs,2)/dx_obs/dx_obs*units.eps0/units.e;
  dn_obs = [dn_obs(1); tocolumn(dn_obs); dn_obs(end)]; % assume phi -> phi at edges
    
  nv = 2000;
  vmax = max([5*sqrt(2*units.e*max(phi_vec)/units.me), vd + 2*vt]); 
  v = linspace(-vmax,vmax,nv);
  dv = v(2) - v(1);

  [X_obs,V_obs] = meshgrid(x_obs,v); X_obs = permute(X_obs,[2 1]); V_obs = permute(V_obs,[2 1]);
  PHI_obs = repmat(tocolumn(phi_obs),1,nv);
  VPH_obs = V_obs*0 + vph;
  E_obs = units.me*(V_obs-vph).^2/2 - units.e*PHI_obs;
  ifree_obs = find(E_obs>0); itrap_obs = find(E_obs<=0);

  [Fflat_obs,Fflat_free_obs,Fflat_trap_obs] = get_f_flat(V_obs,n,vt,vd,PHI_obs,VPH_obs);
  nfree_obs = nansum(Fflat_free_obs,2)*dv;
  ntrap_flat_obs = nansum(Fflat_trap_obs,2)*dv;
  ntrap_obs = ntot - torow(nfree_obs) + torow(dn_obs);
  dntrap_obs = ntrap_obs - torow(ntrap_flat_obs);
  [Fabel_obs,Fabel_free_obs,Fabel_trap_obs] = get_f_abel(V_obs,n,vt,vd,PHI_obs,VPH_obs,dntrap_obs);
  Fabel_obs = Fabel_obs + Fflat_trap_obs;              
  [Fscha_obs,Fscha_free_obs,Fscha_trap_obs,beta_obs] = get_f_schamel(V_obs,n,vt,vd,PHI_obs,VPH_obs,ntrap_obs,beta_range);
  
  FVscha_obs = Fscha_obs.*V_obs;
  FVabel_obs = Fabel_obs.*V_obs;  
    
  % flux at edi energy/velocity interval
  vind_edi_0 = intersect(find(v>v_edi_minus),find(v<v_edi_plus));
  vind_edi_180 = intersect(find(v<-v_edi_minus),find(v>-v_edi_plus));
  FVdv_abel_obs_edi_0 = nansum(FVabel_obs(:,vind_edi_0),2)*dv;
  FVdv_abel_obs_edi_180 = nansum(FVabel_obs(:,vind_edi_180),2)*dv;
  FVdv_scha_obs_edi_0 = nansum(FVscha_obs(:,vind_edi_0),2)*dv;
  FVdv_scha_obs_edi_180 = nansum(FVscha_obs(:,vind_edi_180),2)*dv;

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
  if 1 % plot 2
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
      clear hts
      hts = irf_plot(hca,{phi_progressive.phi_detrend_shift_only_pos_vals,phi},'comp');
      irf_zoom(hca,'x',tint_single+0.003*[-1 1])
      hca.XLabel.String = 'x (km)';
      hca.YLabel.String = '\phi (V)';      
      hca.Title.String = sprintf('tint_{obs} = %s - %s',tint_utc(1,:),tint_utc(2,:));
      irf_legend(hca,{sprintf('v_{ph}= %.0f km/s',vph*1e-3)},[0.99,0.99],'Color',[0 0 0])
    end
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
      hca.XLim = [-30 10];
      irf_legend(hca,{sprintf('x=%.0f km',x_obs(indx(1))*1e-3),sprintf('x=%.0f km',x_obs(indx(2))*1e-3),sprintf('x=%.0f km',x_obs(indx(3))*1e-3)},[0.1 0.15])
      irf_legend(hca,{'- Abel','--Scha'},[0.1 0.05],'color',[0 0 0])
      hca.Title.String = 'based on \phi_{obs}';
    end
    cn.print(sprintf('AbelScha_obs_eh%g_flux_mms%g_',ih,mms_id))
  end  
  %% Collect for all eh
  [fitreslts,gof,fun_net,fun_net_prime] = createFit(torow(phi_vec),torow(dntrap_obs));
  dntrap_all{ih} = dntrap_obs;
  x_all{ih} = x_vec;
  phi_all{ih} = phi_vec;
  fun_fit_all{ih} = fun_net;
  fun_prime_fit_all{ih} = fun_net_prime;
end

%% For entire time series
  
  vph = -9000*1e3;        
  tint_single = tint_all;
  tint_utc = tint_single.utc;
  c_eval('Epar = gseE?par;',mms_id)
  [phi,phi_progressive,phi_ancillary] = get_phi(Epar,vph,tint_all,tint_single);
  
  phi_vec = phi.data;
  phi_obs = phi_vec;
  x_vec = phi_ancillary.x_vec;
  x_obs = x_vec;
  
  % Get electric field and n
  dx_obs = x_obs(2) - x_obs(1);
  x_obs_diff1 = x_obs(1:end-1) + 0.5*dx_obs;
  x_obs_diff2 = x_obs(2:end-1) + dx_obs;
  Efield_obs = -diff(phi_obs)/dx_obs;
  c_eval('Efield_obs = gseE?par.tlim(tint).data;',mms_id)
  dn_obs = diff(phi_obs,2)/dx_obs/dx_obs*units.eps0/units.e;
  dn_obs = [dn_obs(1); tocolumn(dn_obs); dn_obs(end)]; % assume phi -> phi at edges
    
  nv = 2000;
  vmax = max([5*sqrt(2*units.e*max(phi_vec)/units.me), vd + 2*vt]); 
  v = linspace(-vmax,vmax,nv);
  dv = v(2) - v(1);

  [X_obs,V_obs] = meshgrid(x_obs,v); X_obs = permute(X_obs,[2 1]); V_obs = permute(V_obs,[2 1]);
  PHI_obs = repmat(tocolumn(phi_obs),1,nv);
  VPH_obs = V_obs*0 + vph;
  E_obs = units.me*(V_obs-vph).^2/2 - units.e*PHI_obs;
  ifree_obs = find(E_obs>0); itrap_obs = find(E_obs<=0);

  [Fflat_obs,Fflat_free_obs,Fflat_trap_obs] = get_f_flat(V_obs,n,vt,vd,PHI_obs,VPH_obs);
  nfree_obs = nansum(Fflat_free_obs,2)*dv;
  ntrap_flat_obs = nansum(Fflat_trap_obs,2)*dv;
  ntrap_obs = ntot - torow(nfree_obs) + torow(dn_obs);
  dntrap_obs = ntrap_obs - torow(ntrap_flat_obs);
  [Fabel_obs,Fabel_free_obs,Fabel_trap_obs] = get_f_abel(V_obs,n,vt,vd,PHI_obs,VPH_obs,dntrap_obs);
  Fabel_obs = Fabel_obs + Fflat_trap_obs;              
  [Fscha_obs,Fscha_free_obs,Fscha_trap_obs,beta_obs] = get_f_schamel(V_obs,n,vt,vd,PHI_obs,VPH_obs,ntrap_obs,beta_range);
  
  FVscha_obs = Fscha_obs.*V_obs;
  FVabel_obs = Fabel_obs.*V_obs;  
    
  % flux at edi energy/velocity interval
  vind_edi_0 = intersect(find(v>v_edi_minus),find(v<v_edi_plus));
  vind_edi_180 = intersect(find(v<-v_edi_minus),find(v>-v_edi_plus));
  FVdv_abel_obs_edi_0 = nansum(FVabel_obs(:,vind_edi_0),2)*dv;
  FVdv_abel_obs_edi_180 = nansum(FVabel_obs(:,vind_edi_180),2)*dv;
  FVdv_scha_obs_edi_0 = nansum(FVscha_obs(:,vind_edi_0),2)*dv;
  FVdv_scha_obs_edi_180 = nansum(FVscha_obs(:,vind_edi_180),2)*dv;

  phi_long = phi_vec;
  dntrap_long = dntrap_obs;
  [fitreslts,gof,fun_net_long,fun_net_prime_long] = createFit(torow(phi_vec),torow(dntrap_obs));
  
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
  if 1 % plot 2
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
  

%% Plot, compare dntrap/dphi
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






