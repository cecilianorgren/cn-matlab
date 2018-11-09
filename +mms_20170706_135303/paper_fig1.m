%% Prepare data
tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');

eDist = ePDist1.tlim(tint).elim(eint);

eint = [000 40000];
vint = [-Inf Inf];
scpot = scPot1.resample(eDist);
lowerelim = scpot*0 + 0;

tic; ef1D = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B

tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
c_eval('[phi?,phi_progressive?,phi_ancillary?] = get_phi(gseE?par,vph,tint_zoom,tint_zoom);',1:4)

%% Plot, local plasma properties, wave properties from observations combined, larger time + zoomin
ic = 1;
npanels = 9;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.40Z/2017-07-06T13:54:06.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
tint_zoom = phi1.time([1 end]);

% load eh data
data_tmp = load(sprintf('/Users/%s/GoogleDrive/Data/Events/2017-07-06_081603/EH_properties.mat',localuser));
obs_eh_properties = data_tmp.EH_properties;
obs_lpp = obs_eh_properties.Lpp; % peak to peak length
obs_potential = obs_eh_properties.potential;
obs_vtrap = sqrt(2*units.e*obs_potential/units.me)*1e-3; % km/s
obs_potential_max = obs_eh_properties.max_potential;
obs_velocity = obs_eh_properties.vel;
obs_neh = numel(obs_velocity);
c_eval('obs_t0_epoch_mms? = obs_eh_properties.time_mms?;')
c_eval('obs_phi? = irf.ts_scalar(obs_t0_epoch_mms?,obs_potential(:,?));')
c_eval('obs_vph? = irf.ts_scalar(obs_t0_epoch_mms?,obs_velocity);')
c_eval('obs_vtrap? = irf.ts_scalar(obs_t0_epoch_mms?,obs_vtrap(:,?));')
c_eval('obs_vph?.data(isnan(obs_potential(:,?))) = NaN;')

if 0 % B abs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B abs');
  irf_plot(hca,gseB1.abs);
  hca.YLabel.String = 'B (nT)';
end
if 1 % B GSE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 0 % iDEF omni
  isub = isub + 1;
  hca = irf_panel('iDEF');  
  [hout,hcb] = irf_spectrogram(hca,iDist.convertto('s^3/m^6').omni.specrec,'log');
%   hold(hca,'on')
%   lineScpot = irf_plot(hca,scPot1,'k');
%   lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
%   hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % iPDist pa 64
  isub = isub + 1;
  hca = irf_panel('i PA e64 deflux lowe');  
  eint = [100 40000];  
  try
    c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
end
if 0 % i psd vpar
  isub = isub + 1;
  hca = irf_panel('iLine');
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVi},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = if1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_i (km/s)'; 
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi,Ve x
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vex,Vix');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gseVi?.x*10,gseVe?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_x','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'V_i','V_e'},[0.98 0.9],'fontsize',12);
end
if 0 % ePSD omni
  isub = isub + 1;
  hca = irf_panel('ePSD');
  [hout,hcb] = irf_spectrogram(hca,eDist.convertto('s^3/m^6').omni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scpot,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  lineElow = irf_plot(hca,lowerelim,'k');
  lineElow.Color = [0 0 0]; lineElow.LineWidth = 1.5; lineElow.LineStyle = '--';
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};   
  irf_legend(hca,'V_{sc}',[0.99 0.1],'color',0*[1 1 1])
end
if 1 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.deflux.omni.specrec,'log');
  if 1 % scpot
    hold(hca,'on')
    lineScpot = irf_plot(hca,scpot,'k-');
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.0;
    irf_legend(hca,'V_{sc}',[0.06 0.3],'color',0*[1 1 1])
    hold(hca,'off')
  end
  if 1 % temperature
     hold(hca,'on')
    lineTemp = irf_plot(hca,gseTe1.trace/3,'k-');
    lineTemp.Color = [0 0 0]; lineTemp.LineWidth = 1.0;
    
    irf_legend(hca,'T_{e}',[0.06 0.92],'color',0*[1 1 1])
    hold(hca,'off')
  end
  if 0 % lowerelim
    hold(hca,'on')
    lineElow = irf_plot(hca,lowerelim,'k');
    lineElow.Color = 1+[0 0 0]; lineElow.LineWidth = 1.5; lineElow.LineStyle = '--';
    irf_legend(hca,'E_{lim}',[0.01 0.32],'color',1*[1 1 1])
    hold(hca,'off')
  end
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('fred');
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  if 0 % electron moment along projection line
    hold(hca,'on')
    irf_plot(hca,{lineVe*1e-3},'comp')
    %irf_plot(hca,gseVi1)  
    hold(hca,'off')
  end
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.1],'color',1*[1 1 1])
  if unique(ef1D.ancillary.lowerelim) == 1
    irf_legend(hca,['E_{e} >' num2str(unique(ef1D.ancillary.lowerelim)) 'V_{sc}'],[0.99 0.1],'color',1*[1 1 1])
  end
end
if 0 % Te par perp, Ti
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
end
if 0 % Te par perp, no Ti
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))  
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);  
end
if 0 % E perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
end

% zoom in
isub_long = 1:isub;
isub_short = (isub+1):npanels;

% common time shifts
dt = [0.0000  -0.0012  -0.0009  -0.0012];
dt0 = 0.0008;
dt = dt + dt0;
  
if 0 % E par, 4 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % E par, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],num2str(dt(2)*1e3,format_ms),num2str(dt(3)*1e3,format_ms),num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{obs_phi1,obs_phi2,obs_phi3,obs_phi4},'comp');
  %c_eval('hh.Children(?).Marker = ''.'';',1:4)
  hca.YLabel.String = {'\Phi_{||}','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 1 % Phi, use eh_model_optimization_abel to get phi?
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt);  
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 1 % charge perturbation
  sub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('density perturbation');
  set(hca,'ColorOrder',mms_colors('1234b'))
  % cm scale 
  units_scaling = 1e-3;
  hh = irf_plot(hca,{dn_E1par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E2par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E3par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_E4par/v_for_density_scaling*1e-6/units_scaling,...
                     dn_divE/units_scaling,...
                     },'comp','dt',[dt 0]); 
  hh(1).Children(1).LineWidth = 2;
  hca.YLabel.String = {'\delta n',sprintf('(10^{%g} cm^{-3})',log10(units_scaling))};
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  if 1 % legend mms1 mms2 mms43 mms4, 4sc
    set(hca,'ColorOrder',mms_colors('1234b'))
    irf_legend(hca,{'mms1';'mms2';'mms3';'mms4';'4 sc'},[1.02 0.9],'fontsize',12);
  else
    set(hca,'ColorOrder',mms_colors('b'))
    irf_legend(hca,{'4 sc'},[1.02 0.9],'fontsize',12);
  end
end
if 0 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fred vph vtrap');
  
  [hsurf,hcbar] = irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  hold(hca,'on')
  c_eval('vmin? = obs_vph? - obs_vtrap?;',1:4)
  c_eval('vmax? = obs_vph? + obs_vtrap?;',1:4)
  
  %set(hca,'ColorOrder',mms_colors('111223344'))
  set(hca,'ColorOrder',mms_colors('111223344'))
  vscale = 1e-3;
  hlines = irf_plot(hca,{obs_vph1*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  
  htrap = hlines.Children(1:end-2);
  hvph = hlines.Children(end-1);
  c_eval('htrap(?).Marker = ''.'';',1:numel(htrap))
  hvph.LineStyle = '--';
  
  %irf_patch(hca,{vmin,vmax})
  %hca.YLim = sort(real([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]));
  hold(hca,'off')
  hca.YLabel.String = {'v_{||}','(10^3 km/s)'};  
  if 0 % label vtrap vph vtrap
    set(hca,'ColorOrder',mms_colors('122'))
    irf_legend(hca,{'v_{ph}'},[0.4 0.7],'fontsize',12);
    irf_legend(hca,{'v_{trap}'},[0.4 0.99],'fontsize',12);
    irf_legend(hca,{'v_{trap}'},[0.4 0.3],'fontsize',12);
  end
  if 1 % label -- vph
    set(hca,'ColorOrder',mms_colors('1'))
    irf_legend(hca,{'-- v_{ph}'},[0.01 0.98],'fontsize',12);    
  end
  hsurf = findobj(hca.Children,'Type','Surface');
  hsurf.FaceAlpha = 1;
  hline = findobj(hca.Children,'Type','Line');
  c_eval('hline(?).LineWidth = 1.5;',1:numel(hline))
end
if 0 % edi flux 0 180 1 sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{flux0_mms?,flux180_mms?},''comp'',''dt'',dt);',ic)
  hca.YLabel.String = {'flux 180^o','10^6 s^{-1}m^{-2})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'0^o','180^o'},[0.98 0.9],'fontsize',12);
end
if 1 % edi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  hca.YLabel.String = {'flux','(10^6 s^{-1}cm^{-2})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'EDI 180^o'},[0.05 0.99],'fontsize',12);
end

c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
c_eval('h(?).Position(2) = h(?).Position(2)+0.02;',isub_long)


irf_zoom(h(isub_long),'x',tint)
irf_zoom(h(isub_short),'x',tint_zoom)
irf_zoom(h(isub_short),'x',phi1.time)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(isub_long(end)).XLabel.String = [];

[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(isub_long(end)),h(isub_short(1))); 
if 0 % two black lines
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom(1),[0 0 0]);',isub_long)
  c_eval('h_mark2(?) = irf_pl_mark(h(?),tint_zoom(2),[0 0 0]);',isub_long)
else % colored region
  mark_color = [1 0.3 0];
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom,mark_color); ',isub_long)
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_long(1:3))
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','w','w','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
end

colormap(cn.cmap('blue_white'))
%colormap(cn.cmap('white_blue'))
hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
hca = irf_panel('Vi'); hca.YLim = [-799 399];
hca = irf_panel('edi flux'); hca.YLim = [0 8];
%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];

  
