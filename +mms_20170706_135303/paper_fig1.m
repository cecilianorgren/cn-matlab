%% Prepare data
% First run
% mms_20170706_135303.load_data
% mms_20170706_135303.prepare_data_single_sc
% mms_20170706_135303.prepare_data_multi_sc
tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar

eDist = ePDist1.tlim(tint);

% remove background
nSecondary = [5];
nPhoto = 0;
%[eDist_nobg] = mms.remove_edist_background(eDist_orig);
c_eval('[eDist_nobg?] = mms.remove_edist_background(eDist,''nSecondary'',nSecondary(?),''Nphotoe_art'',nPhoto,''ZeroNaN'',0);',1:numel(nSecondary))


ne1_mean = mean(ne1.tlim(tint_zoom).data);
n0 = 0.04;

eint = [000 40000];
vint = [-Inf Inf];
scpot = scPot1.resample(eDist);
lowerelim = scpot*0 + 00;

%tic; ef1D = eDist.reduce('1D',dmpaB1.resample(eDist).norm,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B
tic; ef1D = eDist_nobg1.reduce('1D',dmpaB1.resample(eDist).norm,'vint',vint,'scpot',scpot,'lowerelim',lowerelim,'nMC',500); toc % reduced distribution along B


vph = -9000e3;
c_eval('[phi?,phi_progressive?,phi_ancillary?] = get_phi(gseE?par,vph,tint_zoom,tint_zoom);',1:4)

% Prepare flux in closest FPI energy channel.
if 0
  %%
  c_eval('ePitch?_fpi = ePDist?.pitchangles(dmpaB?,16);',1:4); 
  c_eval('ePitch?_flux_fpi = ePitch?_fpi.flux;',1:4);
  c_eval('ePitch?_flux_fpi_apar500 = ePitch?_flux_fpi.elim(500).palim(180);',1:4);
  c_eval('ePitch?_flux_fpi_par500 = ePitch?_flux_fpi.elim(500).palim(0);',1:4);
  
  c_eval('ePitch?_fpi2 = ePDist?.pitchangles(dmpaB?,180-[22.5 0]);',1:4);
  c_eval('ePitch?_flux_fpi2 = ePitch?_fpi2.flux;',1:4);
  c_eval('ePitch?_flux_fpi2_apar500 = ePitch?_flux_fpi2.elim(500).palim(180);',1:4);
  %% Make the data stepfunction-like, so that it shows the accumulation time
  dt = 0.015;
  for iic = 1:4    
    c_eval('pitch_tmp = ePitch?_flux_fpi2_apar500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi2_apar500_step = irf.ts_scalar(newtime,newdata);',iic)
    
    c_eval('pitch_tmp = ePitch?_flux_fpi_apar500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi_apar500_step = irf.ts_scalar(newtime,newdata);',iic)
    
    c_eval('pitch_tmp = ePitch?_flux_fpi_par500;',iic)
    newtime = [pitch_tmp.time+-dt pitch_tmp.time+dt];
    [newtime,sortind] = newtime.sort;
    newdata = zeros(newtime.length,1); 
    newdata(1:2:end) = pitch_tmp.data;
    newdata(2:2:end) = pitch_tmp.data;  
    c_eval('ePitch?_flux_fpi_par500_step = irf.ts_scalar(newtime,newdata);',iic)
  end
end

%% Plot, local plasma properties, wave properties from observations combined, larger time + zoomin
ic = 1;
npanels = 10;
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
  palim = [168 180];  
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
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

%% Plot, (less overview panels) local plasma properties, wave properties from observations combined, larger time + zoomin
ic = 1;
npanels = 7;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 12;

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
if 0 % Ve
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
if 0 % eDEF omni
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
  [~,hcb_fred] = irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
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
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
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
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 1 % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
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
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [168 180];  
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end

c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
hcb_fred.Position(2) = hcb_fred.Position(2)+0.02;
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
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_long(1:2))
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','w','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

colormap(cn.cmap('blue_white'))
%colormap(cn.cmap('white_blue'))
%hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
hca = irf_panel('Vi'); hca.YLim = [-799 399];
hca = irf_panel('edi flux'); hca.YLim = [0 8];
%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];
hca = irf_panel('E par dt'); hca.YLim = [-70 60];

doDoubleAxis = 1; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
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
  ax2.FontSize = fontsize;
end  
  
irf_plot_axis_align

%% Plot, (less overview panels) local plasma properties (incl. n), wave properties from observations combined, larger time + zoomin
ic = 1;
npanels = 9;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 10;

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z')+0*[-8 8];
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
if 0 % Ve
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
if 0 % eDEF omni
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
  [~,hcb_fred] = irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
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
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  
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
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
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
if 0 % E perp, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp z dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).z,gseE2perp.filt(flow,0,[],5).z,gseE3perp.filt(flow,0,[],5).z,gseE4perp.filt(flow,0,[],5).z},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{\perp,z}','(mV/m)'};  
end
if 0 % E perp, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp y dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).y,gseE2perp.filt(flow,0,[],5).y,gseE3perp.filt(flow,0,[],5).y,gseE4perp.filt(flow,0,[],5).y},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{\perp,y}','(mV/m)'};  
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
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 1 % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
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
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [168 180];  
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end

c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
hcb_fred.Position(2) = hcb_fred.Position(2)+0.02;
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
  mark_color = [0.9290    0.6940    0.1250];
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom,mark_color); ',isub_long)
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_long(1:2))
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','w','k','k','k','k','k','k','k','k','k'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',12,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

%colormap(cn.cmap('blue_white'))
colormap(pic_colors('candy4'))
%colormap(cn.cmap('white_blue'))
%hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
hca = irf_panel('Vi'); hca.YLim = [-799 399];
hca = irf_panel('edi flux'); hca.YLim = [0 8.999];
%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];
hca = irf_panel('E par dt'); hca.YLim = [-70 60];

hca = irf_panel('E par'); hca.YLim = [-60 60];

doDoubleAxis = 1; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
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
  ax2.FontSize = fontsize;
end  
  
irf_plot_axis_align

hcl = findobj(gcf,'type','line');
c_eval('hcl(?).LineWidth = 0.5;',1:numel(hcl))
hh(1).Children(3).LineWidth = 2;
%% Plot, overview only
ic = 1;
npanels = 5;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.40Z/2017-07-06T13:54:06.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z'); % if showing 4 sc epar
tint_zoom = phi1.time([1 end]);
tint_zoom_epar = irf.tint('2017-07-06T13:54:05.00Z/2017-07-06T13:54:06.00Z'); % if showing 4 sc epar

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
if 0 % Ve
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
if 0 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.deflux.omni.specrec,'log');
  if 1 % scpot
    hold(hca,'on')
    lineScpot = irf_plot(hca,scpot,'k-');
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.0;
    irf_legend(hca,'V_{sc}',[0.03 0.3],'color',0*[1 1 1])
    hold(hca,'off')
  end
  if 0 % temperature
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
  [~,hcb] = irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  if 0 % electron moment along projection line
    hold(hca,'on')
    irf_plot(hca,{lineVe*1e-3},'comp')
    %irf_plot(hca,gseVi1)  
    hold(hca,'off')
  end
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
  hcb.YLabel.String = {'log_{10}f(v_{||})','(s/m^4)'};
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
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
end
if 1 % E par, mms1
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par},'comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
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
if 0 % E par, 4 sc, time shifted for visibility
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
  set(hca,'ColorOrder',mms_colors('a'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  hca.XGrid = 'on';
  hca.YGrid = 'on';
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
if 0 % Phi, use eh_model_optimization_abel to get phi?
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
if 0 % charge perturbation
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
if 0 % edi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('1234'))
  palim = [168 180];  
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
  hca.YLabel.String = {'flux','(10^6 s^{-1}cm^{-2})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'EDI 180^o'},[0.05 0.99],'fontsize',12);
end

c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
c_eval('h(?).Position(2) = h(?).Position(2)+0.02;',isub_long)

hcb = findobj(gcf,'type','colorbar');
c_eval('hcb(?).Position(2) = hcb(?).Position(2)+0.02;',1:numel(hcb))

irf_zoom(h(isub_long),'x',tint)
%irf_zoom(h(isub_short),'x',tint_zoom_epar)
%irf_zoom(h(isub_short),'x',phi1.time)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(isub_long(end)).XLabel.String = [];

%[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(isub_long(end)),h(isub_short(1))); 


hca = irf_panel('E par'); hca.YLim = [-60 60];
% Time interval markings
if 0 % two black lines
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom(1),[0 0 0]);',isub_long)
  c_eval('h_mark2(?) = irf_pl_mark(h(?),tint_zoom(2),[0 0 0]);',isub_long)
elseif 1 % colored region
  mark_color = [1 0.3 0];
  mark_color = [0.9290    0.6940    0.1250];
  clear h_mark1;
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom,mark_color); ',isub_long)
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_long([1 2 4]))
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom,mark_color); ',isub_short)
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_short)
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','w','w','k','k','k','k','k','k','k'};
legends_color = {'k','k','w','w','k','k','k','k','k','k','k'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
end

colormap(cn.cmap('blue_white'))
colormap(pic_colors('candy4'))
%colormap(cn.cmap('white_blue'))
%hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
hca = irf_panel('Vi'); hca.YLim = [-799 399];
%hca = irf_panel('edi flux'); hca.YLim = [0 8];
%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];

%% Plot, zoom only
ic = 1;
npanels = 6;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
fontsize = 12;

v_for_density_scaling = 9000e3; % m/s

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z')+0*[-8 8];
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
  irf_legend(hca,{'MMS 1';'MMS 2';'MMS 3';'MMS 4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
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
if 0 % E perp, abs, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp abs dt');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseE1perp.abs,gseE2perp.abs,gseE3perp.abs,gseE4perp.abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  
  hca.YLabel.String = {'|E_{\perp}|','(mV/m)'};  
end
if 1 % E perp, abs, filt, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp abs dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).abs,gseE2perp.filt(flow,0,[],5).abs,gseE3perp.filt(flow,0,[],5).abs,gseE4perp.filt(flow,0,[],5).abs},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  %irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  irf_legend(hca,{sprintf('f > %g Hz',flow)},[0.05 0.98],'fontsize',12);
  hca.YLabel.String = {'|E_{\perp}|','(mV/m)'};  
end
if 0 % E perp, z, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp z dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).z,gseE2perp.filt(flow,0,[],5).z,gseE3perp.filt(flow,0,[],5).z,gseE4perp.filt(flow,0,[],5).z},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  irf_legend(hca,{sprintf('f>%g Hz',flow)},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{\perp,z}','(mV/m)'};  
end
if 0 % E perp, y, 4 sc, time shifted for visibility
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp y dt filt');
  set(hca,'ColorOrder',mms_colors('1234'))
  flow = 3;
  irf_plot(hca,{gseE1perp.filt(flow,0,[],5).y,gseE2perp.filt(flow,0,[],5).y,gseE3perp.filt(flow,0,[],5).y,gseE4perp.filt(flow,0,[],5).y},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  format_ms = '%.1f';
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],['  ',num2str(dt(2)*1e3,format_ms)],['  ',num2str(dt(3)*1e3,format_ms)],['  ',num2str(dt(4)*1e3,format_ms)],'] ms'},[0.01 0.1],'fontsize',12);
  %irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],' ',num2str(dt(2)*1e3,format_ms),' ',num2str(dt(3)*1e3,format_ms),' ',num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
  irf_legend(hca,{sprintf('f>%g Hz',flow)},[0.98 0.1],'fontsize',12);
  hca.YLabel.String = {'E_{\perp,y}','(mV/m)'};  
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
  hca.YLabel.FontSize = fontsize;
  ylabel(hca,hca.YLabel.String,'interpreter','tex')

  if 1 % 4sc
    set(hca,'ColorOrder',mms_colors('b'))
    hleg_4sc = irf_legend(hca,{'4 sc'},[0.99 0.98],'fontsize',12);
    hleg_4sc.FontWeight = 'bold';
  elseif 1 % legend mms1 mms2 mms43 mms4, 4sc
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
if 0 % edi flux 0 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux par');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [0 12];  
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6*0.4,ePitch3_flux_edi.palim(palim)*1e-6*0.5,ePitch4_flux_edi.palim(palim)*1e-6*0.5},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','parallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [0 11.25]^o'},[0.05 0.99],'fontsize',12);
end
if 1 % edi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi flux');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{flux180_mms1*1e-6,flux180_mms2*1e-6,flux180_mms3*1e-6,flux180_mms4*1e-6},'comp','dt',dt);
  palim = [168 180];  
  irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6,ePitch3_flux_edi.palim(palim)*1e-6,ePitch4_flux_edi.palim(palim)*1e-6},'comp','dt',dt);
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6,ePitch2_flux_edi.palim(palim)*1e-6*0.4,ePitch3_flux_edi.palim(palim)*1e-6*0.5,ePitch4_flux_edi.palim(palim)*1e-6*0.5},'comp','dt',dt);  
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'\theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end

if 0 % fpi flux 0 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi_par500_step*1e-6,ePitch2_flux_fpi_par500_step*1e-6,ePitch3_flux_fpi_par500_step*1e-6,ePitch4_flux_fpi_par500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','parallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end
if 0 % fpi flux 180 4sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi_apar500_step*1e-6,ePitch2_flux_fpi_apar500_step*1e-6,ePitch3_flux_fpi_apar500_step*1e-6,ePitch4_flux_fpi_apar500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
end
if 1 % fpi flux 180 4sc, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar 2');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ePitch1_flux_fpi2_apar500_step*1e-6,ePitch2_flux_fpi2_apar500_step*1e-6,ePitch3_flux_fpi2_apar500_step*1e-6,ePitch4_flux_fpi2_apar500_step*1e-6},'comp')  
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'\theta = [157.50, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
end
if 0 % fpi flux 180 4sc, 11.25, 22.50
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('fpi flux apar 12 comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  hlines1 = irf_plot(hca,{ePitch1_flux_fpi_apar500_step*1e-6,ePitch2_flux_fpi_apar500_step*1e-6,ePitch3_flux_fpi_apar500_step*1e-6,ePitch4_flux_fpi_apar500_step*1e-6},'comp');
  hcl = hca.Children;  
  c_eval('hcl(?).LineWidth = 1;',1:4)
  hold(hca,'on')
  hlines2 = irf_plot(hca,{ePitch1_flux_fpi2_apar500_step*1e-6,ePitch2_flux_fpi2_apar500_step*1e-6,ePitch3_flux_fpi2_apar500_step*1e-6,ePitch4_flux_fpi2_apar500_step*1e-6},'comp');
  hcl = hca.Children;  
  c_eval('hcl(?).LineStyle = ''--'';',5:8)
  c_eval('hcl(?).LineWidth = 0.5;',5:8)
  c_eval('hcl(?).LineWidth = 1;',1:4)
  hold(hca,'off')
  %hca.YLabel.String = {'j_e^{FPI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.String = {'j_e^{FPI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12);
  irf_legend(hca,{'dashed - \theta = [168.75, 180]^o','    solid - \theta = [157.50, 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'no time shift'},[0.05 0.8],'fontsize',12,'color',[0 0 0]);
end

if 0 % edi, phi comparison for one spacecraft
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi phi comp');
  set(hca,'ColorOrder',mms_colors('a1'))
  c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6},''dt'',dt(1),''color'',[0 0 0]);',ic)
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6});
  hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'0^o','180^o'},[0.98 0.9],'fontsize',12);
  %ax2 = axis('position',hca.Position);
  hca.YLim = [0 3.9999];
  yyaxis('right')
  ax = gca;
  irf_plot(ax,phi1,'dt',dt(1))
  ax.YLabel.String = {'\phi (V)'};
  ax.YLabel.Interpreter = 'tex';
  %ax.YLabel.Position(1) = 1.07;
  irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.05 0.99],'fontsize',12,'color',[0 0 0]);
  irf_legend(hca,{'mms1'},[0.98 0.99],'fontsize',12,'color',[0 0 0]);
end
%c_eval('h(?).Position(2) = h(?).Position(2)-0.03;',isub_short)
%hcb_fred.Position(2) = hcb_fred.Position(2)+0.02;
%c_eval('h(?).Position(2) = h(?).Position(2)+0.02;',isub_long)


%irf_zoom(h(isub_long),'x',tint)
%irf_zoom(h(isub_short),'x',tint_zoom)
irf_zoom(h,'x',phi1.time)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(isub_long(end)).XLabel.String = [];

%[hline1,hline2] = irf_plot_zoomin_lines_between_panels(h(isub_long(end)),h(isub_short(1))); 
if 0 % two black lines
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom(1),[0 0 0]);',isub_long)
  c_eval('h_mark2(?) = irf_pl_mark(h(?),tint_zoom(2),[0 0 0]);',isub_long)
elseif 0 % colored region
  mark_color = [1 0.3 0];
  c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom,mark_color); ',isub_long)
  c_eval('h_mark1(?).FaceAlpha = 0.9;',isub_long(1:2))
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)'};
legends_color = {'k','k','k','k','k','k','k','k','k','k','k','k'};
for ipanel = 1:npanels
  irf_legend(h(ipanel),legends{ipanel},[0.01 0.99],'fontsize',14,'color',legends_color{ipanel});
  h(ipanel).YLabel.FontSize = fontsize;
  h(ipanel).FontSize = fontsize;
end

%colormap(cn.cmap('blue_white'))
%colormap(cn.cmap('white_blue'))
%hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
%hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
%hca = irf_panel('Vi'); hca.YLim = [-799 399];
hca = irf_panel('edi flux'); hca.YLim = [0 7.999];
hca = irf_panel('edi flux par'); hca.YLim = [0 7.999];
hca = irf_panel('fpi flux par'); hca.YLim = [0 7.999];
hca = irf_panel('fpi flux apar'); hca.YLim = [0 7.999];
hca = irf_panel('fpi flux apar 2'); hca.YLim = [0 7.999];
%hca = irf_panel('fpi flux apar 12 comp'); hca.YLim = [0 3.999];

%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];
hca = irf_panel('E par dt'); hca.YLim = [-70 60];

doDoubleAxis = 1; % dn
if doDoubleAxis  
  hca = irf_panel('density perturbation');
  ax1 = hca;
  ax2 = axes('Position',get(ax1,'Position'));
  ax2.YLim = ax1.YLim*1e-3/n0;    
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
  ax2.FontSize = fontsize;
end  
  
irf_plot_axis_align(h)
%ax.YLabel.Position(1) = 1.07;

set(ax2,'xtick',[],'xticklabel',[]); % remove 'xtick' if xticks required

%% Plot, only 1 panel of EDI phi comparison
h = irf_plot(1); 
if 1 % edi, phi comparison for one spacecraft
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('edi phi comp');
  set(hca,'ColorOrder',mms_colors('a1'))
  c_eval('irf_plot(hca,{ePitch?_flux_edi.palim(palim)*1e-6},''color'',[0 0 0]);',ic)
  %irf_plot(hca,{ePitch1_flux_edi.palim(palim)*1e-6});
  %hca.YLabel.String = {'j_e^{EDI}','(10^6 s^{-1}cm^{-2}sr^{-1})'};
  %hca.YLabel.String = {'j_e^{EDI}','\theta = [168.75 180]^o','(10^6 s^{-1}cm^{-2}sr^{-1})'};  
  hca.YLabel.String = {'j_e^{EDI}','antiparallel','(10^6 s^{-1}cm^{-2}sr^{-1})'};  
  hca.YLabel.Interpreter = 'tex';
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'0^o','180^o'},[0.98 0.9],'fontsize',12);
  %ax2 = axis('position',hca.Position);
  hca.YLim = [0 3.999];
  %hca.YTick = [0 1 2 3 4];
  yyaxis('right')
  ax = gca;
  irf_plot(ax,phi1)
  ax.YLabel.String = {'\phi (V)'};
  ax.YLabel.Interpreter = 'tex';
  %ax.YLabel.Position(1) = 1.07;
  %irf_legend(hca,{'EDI: \theta = [168.75 180]^o'},[0.02 0.99],'fontsize',12,'color',[0 0 0]);
  hca.Title.String = 'mms 1';
end
hca.Position(2) = 0.22;
hca.Position(4) = 0.65;
irf_zoom(h,'x',phi1.time)