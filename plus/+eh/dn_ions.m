%% Set up grid and phi
units = irf_units;
fun_phi = @(phimax,x,lx) phimax*exp(-x.^2/2/lx.^2);


% Fi0

ntot = 0.04*1e6;
R = 1;
n = [R (1-R)]*ntot;
T = [6000 0]; % eV
vd = [2000e3 0]; % ms-1 
vt = sqrt(2*units.e*T./units.mp); % m/s
n0 = sum(n);

str_info = {'unperturbed f:';...
            ['T_{in}= [' sprintf('%g  ',T) '] eV'];...
            ['n_{in}= [' sprintf('%g  ',n*1e-6) '] cc'];...
            ['v_{d,in}= [' sprintf('%g  ',vd*1e-3) '] km/s'];...            
            };
               
lx = 1300; 
nx = 207;
x = linspace(-5*lx,5*lx,nx);
phimax = 50;
phi = fun_phi(phimax,x,lx);
vph = -9000e3;

% Get electric field and n
dx = x(2) - x(1);
x_diff1 = x(1:end-1) + 0.5*dx;
x_diff2 = x(2:end-1) + dx;
E = -diff(phi)/dx;

dn = diff(phi,2)/dx/dx*units.eps0/units.e;
dn = [dn(1); tocolumn(dn); dn(end)]; % assume phi -> at edges

vtrap = sqrt(2*units.e*phimax/units.mp);
nv = 4000;
vmax = max([abs(5*vtrap), abs(vd + 2*vt)]); 
vmax = vph;
v = linspace(-vmax,vmax,nv);
dv = v(2) - v(1);

[X,V] = meshgrid(x,v); X = permute(X,[2 1]); V = permute(V,[2 1]);
PHI = repmat(tocolumn(phi),1,nv);
VPH = V*0 + vph;
U = units.me*(V-vph).^2/2 - units.e*PHI;
ifree = find(U>0); itrap_obs = find(U<=0);

[Ffree,Ffree_free,Ffree_trap] = get_f_free(V,n,vt,vd,PHI,VPH,1);
pcolor(X,V*1e-6,U); shading flat;
%%
[Fflat,Fflat_free,Fflat_trap] = get_f_flat(V,n,vt,vd,PHI,VPH,1);
nfree = nansum(Fflat_free,2)*dv;
ntrap_flat = nansum(Fflat_trap,2)*dv;
ntrap = ntot - torow(nfree) + torow(dn);
dntrap = ntrap - torow(ntrap_flat);

[Fabel,Fabel_free,Fabel_trap] = get_f_abel(V,n,vt,vd,PHI,VPH,dntrap);
Fabel = Fabel + Fflat_trap;              

FVscha_obs = Fscha_obs.*V_obs;
FVabel_obs = Fabel_obs.*V_obs;
FVabel_mod = Fabel_mod.*V_mod;
    
%% Plot

if 1 % plot 3
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
  if 1 % phi(x)
    hca = h(isub); isub = isub + 1;
    plot(hca,x_obs*1e-3,phi_obs)
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = '\phi (V)';
    %legend(hca,{'obs','mod, Gaussian','obs no detrend'},'Box','off')
    %irf_legend(hca,{'obs';{'mod:','Gaussian'};{'obs:','no detrend'}},[0.02,0.98])
    
    %hca.Title.String = sprintf('%s - %s',tint_utc(1,:),tint_utc(2,:));   
    hca.Title.String = sprintf('%s - %s',tint_utc(1,:),tint_utc(2,7:end));   
    hca.Title.Position
    irf_legend(hca,{sprintf('v_{ph}= %.0f km/s',round(vph*1e-5)*1e5*1e-3)},[0.13,0.98],'Color',[0 0 0],'fontsize',12)
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
    hcb.YLabel.String = 'f_e (s^1m^{-4})';  
    hca.YLim = vlim*[-1 1]*1e-6;
    colormap(hca,cn.cmap('white_blue'));
    hcb.Position(1) = hcb.Position(1)+shiftcb;
    hcb.Position(1) = hcb.Position(1)-shiftcb;
    %hca.Title.String = 'Abel obs';    
    hca.XLim = xlim;
    if 1 % E contours
      hold(hca,'on')      
      [hc_data,hc] = contour(hca,X_obs*1e-3,V_obs*1e-6,E_obs/units.e,[0 0],'color',[0.8 .8 .8],'linewidth',1.5);
      hold(hca,'off')
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
  if 1 % free, trapped, and total densities
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
  if 1 % FV obs, Abel
    hca = h(isub); isub = isub + 1;
    pcolor(hca,X_obs*1e-3,V_obs*1e-6,FVabel_obs)
    shading(hca,'flat') 
    hcb = colorbar('peer',hca);
    all_hcb{ihcb} = hcb; ihcb = ihcb + 1;
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = 'v_{||} (10^3 km/s)';
    hcb.YLabel.String = 'flux/\Delta v (m^{-3})';  
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
          'LineWidth',1.0);
        for iline = 1:numel(hlines), hlines(iline).LineStyle = '-'; hlines(iline).Color = [0 0 0]; end
        irf_legend(hca,{'- EDI'},[0.01 0.85],'color',hlines(1).Color);   
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
    c_eval('ax(?).YColor = colors(?,:);',1:2)
    c_eval('ax(2).Children(?).Marker = ''*'';',1:nlines)    
    
        
    hca.XLabel.String = 'x (km)';
    hca.YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1})',log10(units_scale_2));
    ax(2).YLabel.String = sprintf('flux (10^%g cm^{-2}s^{-1}sr^{-1})',log10(units_scale_2));
    irf_legend(hca,{'Model';'* EDI'},[0.02 0.6])
    text(hca,0.5*hca.XLim(2),0.2*hca.YLim(2),'0^o','verticalalignment','bottom')
    text(hca,0.5*hca.XLim(2),1.0*hca.YLim(2),'180^o','verticalalignment','top')
    hca.YLim(1) = 0;    
    hca.XLim = xlim;
    c_eval('ax(?).YLim = [0 2.5];',1:2)    
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


