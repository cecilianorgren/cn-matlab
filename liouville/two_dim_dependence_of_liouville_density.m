units = irf_units;
n0 = 0.1;
T0 = 100; 
vt0 = sqrt(2*units.e*T0./units.me); % m/s
v0 = 0;
n_psi = 89; min_psi = 0.01; max_psi = 100*T0;
n_vpsi = 69; min_vpsi = 0.01; max_vpsi = 10*vt0*1e-3; % km/s
doLogspace = 1;
if doLogspace
  psi = logspace(log10(min_psi),log10(max_psi),n_psi);
  vpsi = logspace(log10(min_vpsi),log10(max_vpsi),n_vpsi);
else  
  psi = linspace(min_psi,max_psi,n_psi);
  vpsi = linspace(min_vpsi,max_vpsi,n_vpsi);
end

[PSI,VPSI] = ndgrid(psi,vpsi);

n_lb_map = zeros(n_psi,n_vpsi);
n_sep_map = zeros(n_psi,n_vpsi);
v_lb_map = zeros(n_psi,n_vpsi);
v_sep_map = zeros(n_psi,n_vpsi);

tic;
for i_psi = 1:n_psi
  fprintf('%g,',i_psi)
  for i_vpsi = 1:n_vpsi
    tmp_psi = psi(i_psi);
    tmp_vpsi = vpsi(i_vpsi);
    [n_lb_tmp,n_sep_tmp,v_lb_tmp,v_sep_tmp] = paper_electron_acceleration.liouville_mapped_nf(n0,T0,v0,tmp_psi,tmp_vpsi);
    n_lb_map(i_psi,i_vpsi) = n_lb_tmp;
    n_sep_map(i_psi,i_vpsi) = n_sep_tmp;
    v_lb_map(i_psi,i_vpsi) = v_lb_tmp;
    v_sep_map(i_psi,i_vpsi) = v_sep_tmp;
  end
end
fprintf('\n')
toc

%% Calculate intersection between n and v lines
% first load phi statistics: obtained from running paper_electron_acceleration.phi_acc_statistics
% intersecting points are acc_pot./te_lobe and ne_sep./ne_lobe or ne_sep_min./ne_lobe
xmap = PSI/T0; 
ymap = VPSI*1e-3/(vt0*1e-6);
nmap = n_sep_map*1e-6/n0;
vmap = abs(v_sep_map)/vt0;

nvals = double(ne_sep./ne_lobe);
nvals_min = double(ne_sep_min./ne_lobe);
vacc = sqrt(2*units.eV*acc_pot./units.me);
vt = sqrt(2*units.eV*te_lobe./units.me);
vvals = vacc./vt; % based on acc_pot

nevents = numel(nvals);

%contour(xmap,ymap,nmap,nvals);
%hold on
%contour(xmap,ymap,vmap,vvals);

% Are these changing order in some way? Perhaps do within loop?
Cn = contourcs(xmap(:,1),ymap(1,:),nmap',nvals);
Cv = contourcs(xmap(:,1),ymap(1,:),vmap',vvals);


% Use polyxpoly to find intersections
doPlot = 1;
if doPlot, hca = subplot(1,1,1); end
xy_intersect = nan(nevents,2);
xy_intersect_min = nan(nevents,2);

for ic = 1:nevents
    Cn = contourcs(xmap(:,1),ymap(1,:),nmap',nvals(ic)*[1 1]);
    Cn_min = contourcs(xmap(:,1),ymap(1,:),nmap',nvals_min(ic)*[1 1]);
    Cv = contourcs(xmap(:,1),ymap(1,:),vmap',vvals(ic)*[1 1]);    
    
    [xint_tmp, yint_tmp] = polyxpoly(Cn(1).X, Cn(1).Y, Cv(1).X, Cv(1).Y);      
    [xint_tmp_min, yint_tmp_min] = polyxpoly(Cn_min(1).X, Cn_min(1).Y, Cv(1).X, Cv(1).Y);
    
    if doPlot
      if 1
        plot(hca,Cn(1).X, Cn(1).Y, Cv(1).X, Cv(1).Y,...
                 Cn_min(1).X, Cn_min(1).Y)
        hca.XLim = [0 50];
        hca.YLim = [0 10];
        hold(hca,'on')
        plot(hca,xint_tmp,yint_tmp,'k*',xint_tmp_min,yint_tmp_min,'ro')
        hold(hca,'off')
        pause(0.1)
      end
    end
    
    if not(isempty(xint_tmp))
      xy_intersect(ic,1:2) = [xint_tmp,yint_tmp];
      xy_intersect_min(ic,1:2) = [xint_tmp_min,yint_tmp_min];
    end
end

%% Plot map with all intersection points
h = setup_subplots(2,3,1);
isub = 1;
if 1 % contour plot n and v, overlaid, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6/n0,0:0.1:1,'r');
  %clabel(C,hc,'LabelSpacing',200);
  
  hold(hca,'on')
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),abs(v_sep_map)/vt0,0:0.5:10,'k');  
  %clabel(C,hc,'LabelSpacing',300);
  hold(hca,'off')
  
  %hb = colorbar('peer',hca);
  %hb.YLabel.String = 'v_e^{map}/v_{t0}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
  hca.Position(3) = h(1).Position(3);
  
  hold(hca,'on')
  colors = pic_colors('matlab');
  sub_events = 10:13;
  plot(hca,xy_intersect(:,1),xy_intersect(:,2),'k+',...
           xy_intersect_min(:,1),xy_intersect_min(:,2),'ko',...
           [xy_intersect(:,1) xy_intersect_min(:,1)]',[xy_intersect(:,2) xy_intersect_min(:,2)]','color',colors(1,:),'linewidth',1.5)
  %plot(hca,xy_intersect_min(:,1),xy_intersect_min(:,2),'ko')
  hold(hca,'off')
  
  
  legend(hca,{'n/n_0 - const.',' v/v_{t0} - const.'},'location','northwest')
  hca.YLim = [0 7];
end

if 1 % scatter plot of Blobe vs n_sheet
  hca = h(isub); isub = isub + 1;
  plot(hca,te_sheet,ne_sheet,'o')
  hca.XLabel.String = 'T_e^{sh} (eV)';
  hca.YLabel.String = 'n^{sh} (cm^{-3})';
  hold(hca,'on')
  hold(hca,'off')
  %legend(hca,{'mean(n_{sep})','min(n_{sep})','v_\perp= 0.5v_{A0}','v_\perp= 0.4v_{A0}','v_\perp= 0.3v_{A0}','v_\perp= 0.2v_{A0}','v_\perp= 0.1v_{A0}'},...
  %  'Location','best','Box','off')
end

% relate derived velocities to other velocities, for example the
% reconnection rate
parpperpratio = sqrt(units.mp/units.me);
vparvt0 = xy_intersect(:,2);
vparvt0_min = xy_intersect_min(:,2);
vpar = vparvt0.*vt;
vpar_min = vparvt0_min.*vt;
vperpvt0 = vparvt0/parpperpratio;
vperpvt0_min = vparvt0_min/parpperpratio;
vperp = vperpvt0.*vt; % vt is m/s
vperp_min = vperpvt0_min.*vt; % vt is m/s

vA0 = B_lobe_*1e-9./sqrt(units.mu0*ne_sheet*1e6*units.mp)*1e-3;


if 1 % scatter plot of Blobe vs n_sheet
  hca = h(isub); isub = isub + 1;
  plot(hca,B_lobe_,ne_sheet,'o')
  hca.XLabel.String = 'B^{lb} (nT)';
  hca.YLabel.String = 'n^{sh} (cm^{-3})';
  hold(hca,'on')
  hold(hca,'off')
  %legend(hca,{'mean(n_{sep})','min(n_{sep})','v_\perp= 0.5v_{A0}','v_\perp= 0.4v_{A0}','v_\perp= 0.3v_{A0}','v_\perp= 0.2v_{A0}','v_\perp= 0.1v_{A0}'},...
  %  'Location','best','Box','off')
end

if 1 % scatter plot of PBlobe vs P_sheet
  hca = h(isub); isub = isub + 1;
  PB_lb = B_lobe_.^2*1e-18/units.mu0/2*1e9; % nPa
  PT_sh = units.eV*te_sheet.*ne_sheet*1e6*1e9; % 
  plot(hca,PB_lb,PT_sh,'o')
  hca.XLabel.String = 'B_{lb}^2/2\mu_0 (nPa)';
  hca.YLabel.String = 'P_{sh} (nPa)';
  %hca.XLim = xlim;
  %hca.YLim = ylim;
  hca.XLim = [0 max([hca.XLim hca.YLim])];
  hca.YLim = hca.XLim;  
  
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,xlim,xlim,'color',[0.7 0.7 0.7],'linestyle','--')
  plot(hca,xlim,0.8*xlim,'color',[0.7 0.7 0.7],'linestyle','--')
  plot(hca,xlim,0.6*xlim,'color',[0.7 0.7 0.7],'linestyle','--')
  plot(hca,xlim,0.4*xlim,'color',[0.7 0.7 0.7],'linestyle','--')
  plot(hca,xlim,0.2*xlim,'color',[0.7 0.7 0.7],'linestyle','--')
  hold(hca,'off')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %legend(hca,{'mean(n_{sep})','min(n_{sep})','v_\perp= 0.5v_{A0}','v_\perp= 0.4v_{A0}','v_\perp= 0.3v_{A0}','v_\perp= 0.2v_{A0}','v_\perp= 0.1v_{A0}'},...
  %  'Location','best','Box','off')
end

if 1 % scatter plot of vperp vs vA
  hca = h(isub); isub = isub + 1;
  plot(hca,vA0,vperp*1e-3,'+',vA0,vperp_min*1e-3,'o')
  hca.XLabel.String = 'v_{A0} (km/s)';
  hca.YLabel.String = 'v_\perp = v_{\psi}(m_e/m_p)^{1/2} (km/s)';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  hlines = plot(hca,xlim,0.5*xlim,'--',xlim,0.4*xlim,'--',xlim,0.3*xlim,'--',xlim,0.2*xlim,'--',xlim,0.1*xlim,'--');
  hold(hca,'off')
  legend(hca,{'mean(n_{sep})','min(n_{sep})','v_\perp= 0.5v_{A0}','v_\perp= 0.4v_{A0}','v_\perp= 0.3v_{A0}','v_\perp= 0.2v_{A0}','v_\perp= 0.1v_{A0}'},...
    'Location','best','Box','off')
end

%% Plot, first
h = setup_subplots(2,2);
isub = 1;

if 0 % colormap
  hca = h(isub); isub = isub + 1;
  pcolor(hca,PSI,VPSI*1e-3,n_sep_map*1e-6)
  shading(hca,'flat');
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi (V)';
  hca.YLabel.String = 'v_{\psi} (10^3 km/s)';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
end
if 1 % contour plot
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI,VPSI*1e-3,n_sep_map*1e-6,0.002:0.002:0.04);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi (V)';
  hca.YLabel.String = 'v_{\psi} (10^3 km/s)';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
  plot(hca,2400,17,'r+')
  plot(hca,800,0,'go')
  hold(hca,'off')  
end

if 1 % contour plot, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6/n0,0:0.1:1);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  %hb.YLabel.String = 'n_e^{sep,map}/n^{lb}';
  hb.YLabel.String = 'n^{map}/n_{0}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
  %plot(hca,2400/T0,17/(vt0*1e-6),'r+')
  %plot(hca,800/T0,0,'go')
  %hold(hca,'off')  
end

%% Plot, contours of constant nmap/n0 and v, and n*v
h = setup_subplots(2,3);
isub = 1;


if 1 % contour plot, nsep/n0
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6/n0,0:0.1:1);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  %hb.YLabel.String = 'n_e^{sep,map}/n^{lb}';
  hb.YLabel.String = 'n^{map}/n_{0}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
end
if 0 % contour plot v, without normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),v_sep_map*1e-6,-100:5:100);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{sep,map} (10^3 km/s)';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 1 % contour plot v, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),abs(v_sep_map)/vt0,0:0.5:10);  
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{map}/v_{t0}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 0 % contour plot n and v, overlaid, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6/n0,0:0.1:1,'r');
  %clabel(C,hc,'LabelSpacing',200);
  
  hold(hca,'on')
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),abs(v_sep_map)/vt0,0:0.5:10,'k');  
  %clabel(C,hc,'LabelSpacing',300);
  hold(hca,'off')
  
  %hb = colorbar('peer',hca);
  %hb.YLabel.String = 'v_e^{map}/v_{t0}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
  hca.Position(3) = h(1).Position(3);
  
  try
  hold(hca,'on')
  plot(hca,xy_intersect(:,1),xy_intersect(:,2),'k*',...
           xy_intersect_min(:,1),xy_intersect_min(:,2),'ko',...
           [xy_intersect(:,1) xy_intersect_min(:,1)],[xy_intersect(:,2) xy_intersect_min(:,2)])
  %plot(hca,xy_intersect_min(:,1),xy_intersect_min(:,2),'ko')
  hold(hca,'off')
  catch
  end
  
  legend(hca,{'n/n_0 - const.',' v/v_{t0} - const.'},'location','northwest')
end

if 1 % contour plot n*v, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map.*abs(v_sep_map)*1e-6/n0/vt0,0:0.2:3);  
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{map}v_e^{map}/n_0v_{t0}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 1 % contour plot, nsep/n0 and flux
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6/n0,0:0.1:1);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  %hb.YLabel.String = 'n_e^{sep,map}/n^{lb}';
  hb.YLabel.String = 'n^{map}/n_{0}';
  hca.XLabel.String = 'v_{beam}/v_{t0} = (\psi/T_0)^{1/2}';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
end
if 0 % contour plot n*v, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),n_sep_map.*abs(v_sep_map)*1e-6/n0/vt0,0:0.2:3);  
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{map}v_e^{map}/n_0v_{t0}';
  hca.XLabel.String = 'v_{beam}/v_{t0} = (\psi/T_0)^{1/2}';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 1 % contour plot v, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),abs(v_sep_map)/vt0,0:0.5:10);  
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{map}/v_{t0}';
  hca.XLabel.String = 'v_{beam}/v_{t0} = (\psi/T_0)^{1/2}';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 1 % contour plot, nsep/n0 and flux
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6/n0,0:0.1:1,'k');
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hold(hca,'on')
  [C,hc] = contour(hca,sqrt(PSI/T0),VPSI*1e-3/(vt0*1e-6),n_sep_map.*abs(v_sep_map)*1e-6/n0/vt0,0:0.2:3);  
  %clabel(C,hc,'LabelSpacing',300);
  hold(hca,'off')
   
  %hb.YLabel.String = 'n_e^{sep,map}/n^{lb}';
  hb.YLabel.String = 'n_e^{map}v_e^{map}/n_0v_{t0}';
   hca.XLabel.String = 'v_{beam}/v_{t0} = (\psi/T_0)^{1/2}';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  %hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
end
%% Plot, including phi/Tlb statistics
% phi statistics obtained from running paper_electron_acceleration.phi_acc_statistics
% available possibly relevant data
% acc_pot
% tepar_lobe
% te_lobe
% te
% ne_sep
% ne_sep_min

psi_t0 = acc_pot./te_lobe;
nevents = numel(psi_t0);

h = setup_subplots(2,2);
isub = 1;

if 1 % contour plot of n, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6,0.002:0.002:0.04);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map} (cm^{-3})';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  %hold(hca,'on')
  %plot(hca,2400/T0,17/(vt0*1e-6),'r+')
  %plot(hca,800/T0,0,'go')
  %hold(hca,'off')  
end
if 1 % contour plot of n, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map/n0,0.0:0.1:1);
  %clabel(C,hc,'LabelSpacing',100);  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'n_e^{sep,map}/n^{lb}';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);
  hold(hca,'on')
  for ievent = 1:nevents        
    plotx = PSI/T0;            
    ploty = VPSI*1e-3/(vt0*1e-6);    
    plotc = n_sep_map*1e-6;
    
    % use interp1 to get new plotc
    % boundaries for Xq and Yq, new grid points
    minx = psi_t0(ievent);
    maxx = max_psi/T0;
    miny = 0;
    maxy = max_vpsi/vt0;
    
    % interpolate
    %plotc_interp = interp2(plotx,ploty,plotc,Xq,Yq)
    
    [remrow,remcol] = find(plotx<psi_t0(ievent));
    remrow = unique(remrow)';
    remcol = unique(remcol)';
    
    plotx(remrow,:) = [];
    ploty(remrow,:) = [];
    plotc(remrow,:) = [];
    [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep_min(ievent)*[1 1]);
    %clabel(C,hc,'LabelSpacing',100);  
    hc.LineWidth = 1;    
    hc.LineColor = [0 0 0];  
    
%     [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep(ievent)*[1 1]);    
%     hc.LineWidth = 1;
%     hc.LineColor = [0.8 0 0.3];  
  end
  hold(hca,'off')  
  
%   hold(hca,'on')
%   plot(hca,2400/T0,17/(vt0*1e-6),'r+')
%   plot(hca,800/T0,0,'go')
%   hold(hca,'off')  
end
if 1 % contour plot v, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),v_sep_map*1e-6,-100:5:100);
  clabel(C,hc,'LabelSpacing',100);
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{sep,map} (10^3 km/s)';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
end
if 1 % contour plot v and n, with normalization
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),v_sep_map*1e-6,-100:5:100,'k');
  %clabel(C,hc,'LabelSpacing',100);
  
  hold(hca,'on')
  [C,hc] = contour(hca,PSI/T0,VPSI*1e-3/(vt0*1e-6),n_sep_map*1e-6,0.002:0.002:0.04,'r');
  %clabel(C,hc,'LabelSpacing',100);
  clabel(C,hc,'Manual');
  hold(hca,'off')
  hold(hca,'on')
  for ievent = 1
    plotx = PSI/T0;            
    ploty = VPSI*1e-3/(vt0*1e-6);    
    plotc = n_sep_map*1e-6;
    
    % use interp1 to get new plotc
    % boundaries for Xq and Yq, new grid points
    minx = psi_t0(ievent);
    maxx = max_psi/T0;
    miny = 0;
    maxy = max_vpsi/vt0;
    
    % interpolate
    %plotc_interp = interp2(plotx,ploty,plotc,Xq,Yq)
    
    [remrow,remcol] = find(plotx<psi_t0(ievent));
    remrow = unique(remrow)';
    remcol = unique(remcol)';
    
    plotx(remrow,:) = [];
    ploty(remrow,:) = [];
    plotc(remrow,:) = [];
    [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep_min(ievent)*[1 1]);
    %clabel(C,hc,'LabelSpacing',100);  
    hc.LineWidth = 1;    
    hc.LineColor = [0 0 0];  
    
%     [C,hc] = contour(hca,plotx,ploty,plotc,ne_sep(ievent)*[1 1]);    
%     hc.LineWidth = 1;
%     hc.LineColor = [0.8 0 0.3];  
  end
  hold(hca,'off')  
  
  
  hb = colorbar('peer',hca);
  hb.YLabel.String = 'v_e^{sep,map} (10^3 km/s)';
  hca.XLabel.String = '\psi/T_0';
  hca.YLabel.String = 'v_{\psi}/v_{t0}';
  hca.Title.String = sprintf('n_0 = %.3f cm^{-3}, T0 = %.0f eV',n0,T0);   
  
  hold(hca,'on')
  plot(hca,2400/T0,17/(vt0*1e-6),'r+')
  plot(hca,800/T0,0,'go')
  hold(hca,'off')  
end
