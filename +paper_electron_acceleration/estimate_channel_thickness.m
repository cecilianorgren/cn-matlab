fun_phi = @(phimax,x,l) phimax.*exp(-x.^2./(2*l.^2));
fun_E  = @(phimax,x,l) (x./l.^2).*phimax.*exp(-x.^2./(2*l.^2));
phimax = 1800;

lx = 40e3; % m
x1 = -200*1e3; % m
x2 = 200*1e3; % m
nx = 500;
x = linspace(x1,x2,nx);

%%
h = setup_subplots(2,1);
isub = 1;

hca = h(isub); isub = isub + 1;
ax = plotyy(hca,x*1e-3,fun_phi(phimax,x,lx),x*1e-3,fun_E(phimax,x,lx)*1e3);
ax(1).YLabel.String = '\phi (V)';
ax(2).YLabel.String = 'E (mV/m)';

ax(1).Title.String = sprintf('lx = %g km',lx*1e-3);
ax(2).XGrid = 'on';
ax(2).YGrid = 'on';

hca = h(isub); isub = isub + 1;

plot(hca,l_vec,fun_E(phimax,l_vec,l_vec))
hca.YLabel.String = 'E^{max} (mV/m)';
hca.XLabel.String = 'lx (km)';
hca.XGrid = 'on';
hca.YGrid = 'on';

%%
nl = 101;
l_vec = linspace(20,90,nl);
nphimax = 100;
phimax_vec = linspace(500,2500,nphimax);
[L,PSI] = meshgrid(l_vec,phimax_vec);
EMAX = fun_E(PSI,L,L);

h = setup_subplots(1,1);
isub = 1;

if 1 % EMAX(L,PSI)
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,L,PSI,EMAX);
  htext = clabel(C,hc,'LabelSpacing',1000,'Color','k','FontWeight','bold');
  hca.XLabel.String = 'l_z (km)';
  hca.YLabel.String = '\psi (V)';
  colormap(hca,mms_colors('matlab'))
end
if 0
  hca = h(isub); isub = isub + 1;
  [C,hc] = contourf(hca,PSI,EMAX,L);
  clabel(C,hc,'LabelSpacing',72,'Color','b','FontWeight','bold')
end

%% Plot for paper
fontsize = 14;
tint_zoom = irf.tint('2017-07-06T00:54:11.00Z/2017-07-06T00:54:19.00Z');
time_eperp_reverse = irf_time('2017-07-06T00:54:15.18Z','utc>EpochTT');
npanels = 2;
[h,h2] = initialize_combined_plot(npanels,1,1,0.6,'vertical'); % horizontal

if 1 % ve
  hca = irf_panel('ve');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseVe1.x*1e-3,gseVe1.y*1e-3,gseVe1.z*1e-3},'comp');
  hca.YLabel.String = 'v_e (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'x','y','z'},[0.02 0.98],'fontsize',fontsize)
end
if 0 % Eperp y
  hca = irf_panel('E perp y');
  hlines = irf_plot(hca,{gseE1perp.y,gseE1perp.y.resample(gseE1perp.time(1):0.5:gseE1perp.time(end))},'comp');
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = 'E_{\perp,y} (mV/m)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'original','resampled to 2 Hz'},[0.98 0.98],'fontsize',fontsize)
end
if 1 % Eperp z
  hca = irf_panel('E perp z');
  %set(hca,'ColorOrder',mms_colors('x1'))
  hlines = irf_plot(hca,{gseE1perp.z,gseE1perp.z.resample(gseE1perp.time(1):0.5:gseE1perp.time(end))},'comp');
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = 'E_{\perp,z} (mV/m)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'original';'resampled to 2 Hz'},[0.02 0.98],'fontsize',fontsize)
end
set(irf_panel('E perp z'),'ylim',[-59 59],'ytick',[-100:20:100])

irf_zoom(h,'x',tint_zoom)


for ipanel = 1:npanels
  h(ipanel).FontSize = fontsize;
  hmark =  irf_pl_mark(h(ipanel),time_eperp_reverse.epochUnix,'k');
  hmark.LineWidth = 1.5;
end

  
  hca = h2;
  [C,hc] = contour(hca,L,PSI,EMAX,[0:5:30 40:10:100]);
  htext = clabel(C,hc,'LabelSpacing',1000,'Color','k','FontWeight','bold');
  hca.XLabel.String = 'l_z (km)';
  hca.YLabel.String = '\psi (V)';
  colormap(hca,mms_colors('matlab'))
  hca.FontSize = fontsize;
  hold(hca,'on')
  plot(l_vec,l_vec*0+1800,'--','color',[0.4 0.4 0.4])
  hold(hca,'off')
  h2.XGrid = 'on';
  h2.YGrid = 'on';
  
  h(1).Position(2) = h(1).Position(2)+0.04;
  h(2).Position(2) = h(2).Position(2)+0.04;
  h2.Position(2) = h2.Position(2)+0.04;
%% Plot for paper
fontsize = 14;
tint_zoom = irf.tint('2017-07-06T00:54:11.00Z/2017-07-06T00:54:19.00Z');
time_eperp_reverse = irf_time('2017-07-06T00:54:15.18Z','utc>EpochTT');
npanels = 4;
h = irf_plot(npanels);

if 1 % ve
  hca = irf_panel('ve');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseVe1.x*1e-3,gseVe1.y*1e-3,gseVe1.z*1e-3},'comp');
  hca.YLabel.String = 'v_e (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'x','y','z'},[0.02 0.98],'fontsize',fontsize)
end
if 0 % Eperp y
  hca = irf_panel('E perp y');
  hlines = irf_plot(hca,{gseE1perp.y,gseE1perp.y.resample(gseE1perp.time(1):0.5:gseE1perp.time(end))},'comp');
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = 'E_{\perp,y} (mV/m)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'original','resampled to 2 Hz'},[0.98 0.98],'fontsize',fontsize)
end
if 1 % Eperp z
  hca = irf_panel('E perp z');
  hlines = irf_plot(hca,{gseE1perp.z,gseE1perp.z.resample(gseE1perp.time(1):0.5:gseE1perp.time(end))},'comp');
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = 'E_{\perp,z} (mV/m)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'original';'resampled to 2 Hz'},[0.02 0.98],'fontsize',fontsize)
end
set(irf_panel('E perp z'),'ylim',[-59 59],'ytick',[-100:20:100])

irf_zoom(h(1:2),'x',tint_zoom)


for ipanel = 1:2
  h(ipanel).FontSize = fontsize;
  hmark =  irf_pl_mark(h(ipanel),time_eperp_reverse.epochUnix,'k');
  hmark.LineWidth = 1.5;
end

  isub = 3;
  hca = h(isub); isub = isub + 1;
  delete(hca)
  
  hca = h(isub); isub = isub + 1;
  [C,hc] = contour(hca,L,PSI,EMAX);
  htext = clabel(C,hc,'LabelSpacing',1000,'Color','k','FontWeight','bold');
  hca.XLabel.String = 'l_z (km)';
  hca.YLabel.String = '\psi (V)';
  colormap(hca,mms_colors('matlab'))
  hca.FontSize = fontsize;
  hold(hca,'on')
  plot(l_vec,l_vec*0+1800,'--','color',[0.4 0.4 0.4])
  hold(hca,'off')
