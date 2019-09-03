%% Load data needed for figure paper
ic = 1;
tint = irf.tint('2017-07-06T00:54:03.00Z/2017-07-06T00:56:03.00Z');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
pathLocalUser = ['/Users/' localuser '/'];

c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];
event = 3;
sep.get_tints;
load(sprintf('%s/acc_pot_data_event_%g',saveAccPotPath,event),'acc_pot_data')
tsAccPot_nobg_abs = acc_pot_data.tsAccPot_nobg_abs;
tsAccPot = tsAccPot_nobg_abs;  

%% Get acceleration potential from sep_acc.figure2_get_acc_pot
ic = 1;
itint  = 1;
tint_acc_pot = tints_sep{itint};
c_eval('ts_acc_pot = tsAccPot?_nobg_abs.tlim(tint_acc_pot);',ic)
events = [3];
nevents = numel(events);
for ievent = 1
  event = events(ievent);
  sep.get_tints;
  tints_phi{ievent} = tint_phi;
  tints_sep{ievent} = tint_sep;
  tints_lobe{ievent} = tint_lobe;
  tints_sheet{ievent} = tint_sheet;
end

%%
fun_phi = @(phimax,x,l) phimax.*exp(-x.^2./(2*l.^2));
fun_E  = @(phimax,x,l) (x./l.^2).*phimax.*exp(-x.^2./(2*l.^2));
%fun_l = @(phi,E);
phimax = 1800;

lx = 40e3; % m
x1 = -200*1e3; % m
x2 = 200*1e3; % m
nx = 500;
x = linspace(x1,x2,nx);

%%
nl = 101;
l_vec = linspace(20,90,nl);

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
  
%% Plot for paper, including potential timeseries
fontsize = 14;
tint_zoom = irf.tint('2017-07-06T00:54:11.00Z/2017-07-06T00:54:19.00Z');
time_eperp_reverse = irf_time('2017-07-06T00:54:15.18Z','utc>EpochTT');
npanels = 3;
[h,h2] = initialize_combined_plot(npanels,1,1,0.6,'vertical'); % horizontal

% resample time line, same for Eperp and potential
f_resamp = 3;
t_diff = gseE1perp.time(1)-ePDist1.time(1);
%t0 = ePDist1.time(1)-t_diff;
timeline = ePDist1.time(1):1/f_resamp:ePDist1.time(end);
timeline = timeline +- t_diff;

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
  %hlines = irf_plot(hca,{gseE1perp.z,gseE1perp.z.resample(gseE1perp.time(1):0.5:gseE1perp.time(end))},'comp');
  hlines = irf_plot(hca,{gseE1perp.z,gseE1perp.z.resample(timeline)},'comp');
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = 'E_{\perp,z} (mV/m)';
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,{'original';sprintf('resampled to %g Hz',f_resamp)},[0.02 0.98],'fontsize',fontsize)
  irf_legend(hca,{'original';sprintf('downsampled',f_resamp)},[0.02 0.98],'fontsize',fontsize)
end
if 1 % acceleration potential
  hca = irf_panel('accelersation potential');
  %set(hca,'ColorOrder',mms_colors('x1'))
  %hlines = irf_plot(hca,{ts_acc_pot,ts_acc_pot.resample((ts_acc_pot.time(1)+-0.2):0.5:ts_acc_pot.time(end))},'comp');  
  hlines = irf_plot(hca,{ts_acc_pot,ts_acc_pot.resample(timeline)},'comp');  
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = '\psi (V)';
  hca.YLabel.Interpreter = 'tex';  
  hca.YLim = [0 1999];
  %irf_legend(hca,{'original';sprintf('resampled to %g Hz',f_resamp)},[0.02 0.98],'fontsize',fontsize)
  irf_legend(hca,{'original';sprintf('downsampled',f_resamp)},[0.02 0.98],'fontsize',fontsize)
  if 1
    hold(hca,'on')    
    varplot = ts_acc_pot.resample(timeline);
    ind0 = find(~isnan(varplot.data));
    hlines = irf_plot(hca,{varplot(ind0(3:4))},'comp','ks');  
    hlines
    hold(hca,'off')
  end
end

set(irf_panel('E perp z'),'ylim',[-59 59],'ytick',[-100:20:100])

irf_zoom(h,'x',tint_zoom)


for ipanel = 1:npanels
  h(ipanel).FontSize = fontsize;
  hmark =  irf_pl_mark(h(ipanel),time_eperp_reverse.epochUnix,'k');
  hmark.LineWidth = 1.5;
end

  
  hca = h2;
  PSIscale = exp(-0.5); 
  [C,hc] = contour(hca,L,PSI*PSIscale,EMAX,[0:5:30 40:10:100]);
  htext = clabel(C,hc,'LabelSpacing',1000,'Color','k','FontWeight','bold');
  hca.XLabel.String = 'l_z (km)';
  hca.YLabel.String = '\psi exp(-1/2) (V)';
  colormap(hca,mms_colors('matlab'))
  hca.FontSize = fontsize;
  hold(hca,'on')
  %plot(l_vec,l_vec*0+1800,'--','color',[0.4 0.4 0.4])
  plot(25,1000,'s','color',[0.4 0.4 0.4])
  hold(hca,'off')
  h2.XGrid = 'on';
  h2.YGrid = 'on';
  
  hold(h2,'on')
  psi_data = ts_acc_pot.resample(timeline).tlim(tint_acc_pot).data;  
  %plot(h2,ts_acc_pot.tlim(tint_acc_pot).data)
  hold(h2,'off')
  
  
  h(1).Position(2) = h(1).Position(2)+0.04;
  h(2).Position(2) = h(2).Position(2)+0.04;
  h(3).Position(2) = h(3).Position(2)+0.04;
  h2.Position(2) = h2.Position(2)+0.04;
  
  
% labels (a), (b)...
fontsize_leg = 16;
irf_legend(h(1),{'(a)'},[-0.05 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])
irf_legend(h(2),{'(b)'},[-0.05 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])
irf_legend(h(3),{'(c)'},[-0.05 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])
irf_legend(h2(1),{'(d)'},[-0.1 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])

%% Plot for paper, including potential timeseries, only timeseries
fontsize = 14;
tint_zoom = irf.tint('2017-07-06T00:54:11.00Z/2017-07-06T00:54:19.00Z');
time_eperp_reverse = irf_time('2017-07-06T00:54:15.18Z','utc>EpochTT');
npanels = 3;
h = irf_plot(npanels);

% resample time line, same for Eperp and potential
f_resamp = 3;
t_diff = gseE1perp.time(1)-ePDist1.time(1);
%t0 = ePDist1.time(1)-t_diff;
timeline = ePDist1.time(1):1/f_resamp:ePDist1.time(end);
timeline = timeline +- t_diff;


varplot = ts_acc_pot.resample(timeline);
ind0 = find(~isnan(varplot.data));
times_lz = varplot.time(ind0(3:4));

colors = pic_colors('matlab');
if 1 % ve
  hca = irf_panel('ve');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseVe1.x*1e-3,gseVe1.y*1e-3,gseVe1.z*1e-3},'comp');
  hca.YLabel.String = 'v_e (10^3 km/s)';
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,{'x','y','z'},[0.02 0.98],'fontsize',fontsize)
end
if 1 % Eperp z
  hca = irf_panel('E perp z');
  set(hca,'ColorOrder',colors(1:2,:))
  %set(hca,'ColorOrder',mms_colors('x1'))
  %hlines = irf_plot(hca,{gseE1perp.z,gseE1perp.z.resample(gseE1perp.time(1):0.5:gseE1perp.time(end))},'comp');
  hlines = irf_plot(hca,{gseE1perp.z,gseE1perp.z.resample(timeline)},'comp');
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = 'E_{\perp,z} (mV/m)';
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,{'original';sprintf('resampled to %g Hz',f_resamp)},[0.02 0.98],'fontsize',fontsize)
  irf_legend(hca,{'original';sprintf('downsampled',f_resamp)},[0.02 0.98],'fontsize',fontsize)
    if 1
    hold(hca,'on')    
    varplot = gseE1perp.z.resample(timeline).resample(times_lz);
    %ind0 = find(~isnan(varplot.data));
    hlines = irf_plot(hca,{varplot},'comp');  
    hlines.Children(1).Color = 0*colors(3,:);
    hlines.Children(1).MarkerFaceColor = colors(3,:);
    hlines.Children(1).LineStyle = 'none';
    hlines.Children(1).Marker = 's';
    hlines.Children(1).MarkerSize = 10;
    
    hold(hca,'off')
  end
end
if 1 % acceleration potential
  hca = irf_panel('accelersation potential');
  set(hca,'ColorOrder',colors(1:2,:))
  %set(hca,'ColorOrder',mms_colors('x1'))
  %hlines = irf_plot(hca,{ts_acc_pot,ts_acc_pot.resample((ts_acc_pot.time(1)+-0.2):0.5:ts_acc_pot.time(end))},'comp');  
  hlines = irf_plot(hca,{ts_acc_pot,ts_acc_pot.resample(timeline)},'comp');  
  hlines.Children(1).LineWidth = 1.5;
  hca.YLabel.String = '\psi (V)';
  hca.YLabel.Interpreter = 'tex';  
  hca.YLim = [0 1999];
  %irf_legend(hca,{'original';sprintf('resampled to %g Hz',f_resamp)},[0.02 0.98],'fontsize',fontsize)
  irf_legend(hca,{'original';sprintf('downsampled',f_resamp)},[0.02 0.98],'fontsize',fontsize)
  if 1
    hold(hca,'on')    
    varplot = ts_acc_pot.resample(timeline);
    ind0 = find(~isnan(varplot.data));
    hlines = irf_plot(hca,{varplot(ind0(3:4))},'comp');  
    hlines.Children(1).Color = 0*colors(3,:);
    hlines.Children(1).MarkerFaceColor = colors(3,:);
    hlines.Children(1).LineStyle = 'none';
    hlines.Children(1).Marker = 's';
    hlines.Children(1).MarkerSize = 10;
    
    hold(hca,'off')
  end
end

set(irf_panel('E perp z'),'ylim',[-59 59],'ytick',[-100:20:100])

irf_zoom(h,'x',tint_zoom)


for ipanel = 1:npanels
  h(ipanel).FontSize = fontsize;
  hmark =  irf_pl_mark(h(ipanel),time_eperp_reverse.epochUnix,'k');
  hmark.LineWidth = 1.5;
  hmark.LineStyle = '--';
end

  
  
  
  h(1).Position(2) = h(1).Position(2)+0.04;
  h(2).Position(2) = h(2).Position(2)+0.04;
  h(3).Position(2) = h(3).Position(2)+0.04;  
  
  
% labels (a), (b)...
fontsize_leg = 16;
irf_legend(h(1),{'(a)'},[-0.05 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])
irf_legend(h(2),{'(b)'},[-0.05 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])
irf_legend(h(3),{'(c)'},[-0.05 0.99],'FontSize',fontsize_leg,'Color',[0 0 0])
