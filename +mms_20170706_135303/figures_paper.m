%% Separatrix streaming
% Make reduced distribution
tintZoom = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');

strTintZoom = [irf_time(tintZoom(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tintZoom(2),'epochtt>utc_HHMMSS')];

eint = [000 40000];
vint = [-Inf Inf];

eDist = ePDist1.tlim(tintZoom).elim(eint);
%iDist = iPDist1.tlim(tintZoom).elim(eint);

%% Remove background electrons
[eDist1_bgremoved, eDist1_bg, ephoto_scale] = ...
             mms.remove_edist_background(eDist, 'tint', tintZoom);
 
%%
tintZoom = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:53:50.00Z');
[eDist1_bgremoved, eDist1_bg, ephoto_scale] = ...
             mms.remove_edist_background(ePDist1, 'tint', tintZoom, ...
             'Nphotoe_art', 0.1, 'nSecondary', 0.4, 'ZeroNaN', 0);
h=irf_plot({eDist1_bgremoved.omni.specrec,ePDist1.omni.specrec});
h(1).YScale = 'log';
h(2).YScale = 'log';
h(2).CLim = [-35 -25];
h(1).CLim = [-35 -25];
irf_zoom(h,'x',tintZoom)
%% electrons
%ve = gseVe1.tlim(eDist.time).resample(eDist);
scpot = scPot1.resample(eDist);
scpot_margin = 1.0; % keep in mind that this also affects the velocity at lower energies
lowerelim = scpot*0 + 40;
eLine = dmpaB1.resample(eDist).norm;
%tic; ef1D_ = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B
tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot,'lowerelim',lowerelim); toc % reduced distribution along B
lineVe = ve.dot(eLine); % projection of Vi on B

%% Plot, local plasma properties, wave properties from observations
ic = 1;
npanels = 12;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

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
if 1 % iDEF omni
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
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ePSD omni
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
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_e (km/s)'; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % fe*v vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced * v');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('v_f1D*v','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
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
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % E perp
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
  hca.YLabel.String = {'n','(cm^{-3})'};
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(5).CLim = [-35 -28]+12
colormap('jet');

%% Plot, including f proj and v phi and vtrap
ic = 1;
npanels = 7;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
tint_zoom = tintZoom;

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
if 0 % B GSE
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
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.convertto('s^3/m^6').omni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scpot_lim,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  irf_legend(hca,[num2str(scpot_margin) 'V_{sc}'],[0.99 0.1],'color',0*[1 1 1])
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('fred');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %set(hca,'ColorOrder',mms_colors('122'))
  %irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLim = [-80 80];
  hca.YLabel.String = 'v_e (km/s)'; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
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
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Phi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{obs_phi1,obs_phi2,obs_phi3,obs_phi4},'comp');
  c_eval('hh(1).Children(?).Color = mms_colors(''?'')',1:4)
  hca.YLabel.String = {'\Phi_{||}','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end
if 1 % v phase + trap
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('phase velocity');
  
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  hold(hca,'on')
  c_eval('vmin? = obs_vph? - obs_vtrap?;',1:4)
  c_eval('vmax? = obs_vph? + obs_vtrap?;',1:4)
  
  %set(hca,'ColorOrder',mms_colors('111223344'))
  set(hca,'ColorOrder',mms_colors('111223344'))
  vscale = 1e-3;
  hlines = irf_plot(hca,{obs_vph1*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  
  %irf_patch(hca,{vmin,vmax})
  %hca.YLim = sort(real([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]));
  hold(hca,'off')
  hca.YLabel.String = {'v_{||}','(10^3 km/s)'};  
  set(hca,'ColorOrder',mms_colors('122'))
  irf_legend(hca,{'v_{ph}'},[0.55 0.7],'fontsize',12);
  irf_legend(hca,{'v_{trap}'},[0.55 0.99],'fontsize',12);
  irf_legend(hca,{'v_{trap}'},[0.55 0.3],'fontsize',12);
  
  hsurf = findobj(hca.Children,'Type','Surface');
  hsurf.FaceAlpha = 1;
  hline = findobj(hca.Children,'Type','Line');
  c_eval('hline(?).LineWidth = 1.5;',1:numel(hline))
end
if 1 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);  
end

%h(5).CLim = [-7 -3];
hca = irf_panel('phase velocity');
hca.CLim = [-5 -2.5];
%colormap(hca,'jet')

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(4).CLim = [-35 -28]+12
colormap(cn.cmap('blue_white'));

%% Plot, local plasma properties, wave properties from observations combined, larger time + zoomin
ic = 1;
npanels = 8;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

tint = irf.tint('2017-07-06T13:53:40.00Z/2017-07-06T13:54:15.00Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.30Z/2017-07-06T13:54:05.80Z');
tint_zoom = irf.tint('2017-07-06T13:54:05.45Z/2017-07-06T13:54:05.70Z'); % if showing 4 sc epar

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
    irf_legend(hca,'V_{sc}',[0.01 0.3],'color',0*[1 1 1])
    hold(hca,'off')
  end
  if 1 % temperature
     hold(hca,'on')
    lineTemp = irf_plot(hca,gseTe1.trace/3,'k-');
    lineTemp.Color = [0 0 0]; lineTemp.LineWidth = 1.0;
    
    irf_legend(hca,'T_{e}',[0.01 0.98],'color',0*[1 1 1])
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
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  %irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.1],'color',1*[1 1 1])
  if unique(ef1D.ancillary.lowerelim) == 1
    irf_legend(hca,['E_{e} >' num2str(unique(ef1D.ancillary.lowerelim)) 'V_{sc}'],[0.99 0.1],'color',1*[1 1 1])
  end
end
if 0 % fe*v vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced * v');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('v_f1D*v','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
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
dt = dt+0.0008;
  
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
  irf_legend(hca,{'EDI 180^o'},[0.01 0.99],'fontsize',12);
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
c_eval('h_mark1(?) = irf_pl_mark(h(?),tint_zoom(1),[0 0 0]);',isub_long)
c_eval('h_mark2(?) = irf_pl_mark(h(?),tint_zoom(2),[0 0 0]);',isub_long)

colormap(cn.cmap('blue_white'))
%colormap(cn.cmap('white_blue'))
hca = irf_panel('eDEF'); hca.XGrid = 'off'; hca.YGrid = 'off'; hca.CLim = [4 7.5];
hca = irf_panel('fred'); hca.CLim = [-6.5 -2]; hca.YLim = [-70 70];
%hca = irf_panel('fred vph vtrap'); hcbar.Position(2) = hca.Position(2); hca.YLim = [-35 15]; hca.CLim = [-6.5 -2];
hca = irf_panel('Vi'); hca.YLim = [-799 399];
%hca = irf_panel('n'); hca.YLim = [0 0.199]; hca.YTick = [0 0.05 0.1 0.15];


%% Plot fred, electrons
ic = 1;
npanels = 9;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

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
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ePSD omni
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
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % fe*v vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced * v');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('v_f1D*v','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % fe*v^2 vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced * v^2');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('v_f1D*v^2','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
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
if 1 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  vscale = 1e-3;
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseVe?par*vscale},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};  
end
if 1 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);  
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
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
  hca.YLabel.String = {'n','(cm^{-3})'};
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(5).CLim = [-35 -28]+12;
colormap('jet');
colormap(irf_panel('fe reduced * v'),cn.cmap('blue_red'))
hca = irf_panel('phase velocity');
hca.CLim = [-5 -2];

%hca = irf_panel('ePSD');
%hca.CLim = [ePSDomni_min ePSDomni_max];
hca = irf_panel('fe reduced * v');
hca.CLim = 20*[-1 1];
%colormap(irf_panel('fe reduced * v^2'),cn.cmap('white_blue'))

%% Plot fred, and scaling parameters
ic = 1;
npanels = 11;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

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
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('ePSD');
  ePSDomni = eDist.convertto('s^3/m^6').omni;
  ePSDomni_elim = ePSDomni.elim([100 10000]);
  ePSDomni_elim.data(ePSDomni_elim.data==0) = NaN;
  ePSDomni_max = log10(max(max(ePSDomni_elim.data)));
  ePSDomni_min = log10(min(min(ePSDomni_elim.data)));
  [hout,hcb] = irf_spectrogram(hca,ePSDomni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scpot_lim,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  irf_legend(hca,[num2str(scpot_margin) 'V_{sc}'],[0.99 0.1],'color',0*[1 1 1])
end
if 1 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 1 % fe*v vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced * v');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('v_f1D*v','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
end
if 0 % fe*v^2 vpar
  isub = isub + 1;
  hca = irf_panel('fe reduced * v^2');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('v_f1D*v^2','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  hca.YLim = ef1D.depend{1}(1,[1 end])*1e-3;
  hca.YLabel.String = {'v_e','(10^3 km/s)'}; 
  irf_legend(hca,[num2str(vint(1),'%.0f') '<v_\perp<' num2str(vint(2),'%.0f')],[0.99 0.99],'color',1*[1 1 1])
  irf_legend(hca,['E_{e} >' num2str(scpot_margin) 'V_{sc}'],[0.01 0.99],'color',1*[1 1 1])
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
if 1 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  vscale = 1e-3;
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseVe?par*vscale},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'};  
end
if 1 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);  
end
if 0 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
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
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % beta
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{beta?e},''comp'');',ic)
  hca.YLabel.String = {'\beta',''};
end
if 1 % Egedal2015 eq. 6
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('uepar Egedal2015 eq. 6');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{uepar?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'u_{e,||}','(10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('1'))
  irf_legend(hca,{'Egedal 2015, eq. 6'},[0.98 0.9],'fontsize',12);
end
if 1 % vte
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('vte');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{vte?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'v_{te}','(10^3 km/s)'};
end
if 1 % vA
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('vA');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{vA?},''comp'');',ic)
  hca.YLabel.String = {'v_{A}','(km/s)'};
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(5).CLim = [-35 -28]+12;
colormap('jet');
colormap(irf_panel('fe reduced * v'),cn.cmap('blue_red'))
hca = irf_panel('phase velocity');
hca.CLim = [-5 -2];

hca = irf_panel('ePSD');
hca.CLim = [ePSDomni_min ePSDomni_max];
hca = irf_panel('fe reduced * v');
hca.CLim = 20*[-1 1];
%colormap(irf_panel('fe reduced * v^2'),cn.cmap('white_blue'))

%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];

%% 4 fred, vph/trapping range, EDI range, average flux from model, and f0 of model
% average over time, comaprison to FPI, add to right side of figure
hca = subplot(1,1,1);
if 1 % reduced f fpi
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
end
if 1 % trapping range
  hold(hca,'on')
  hold(hca,'off')
end
if 1 % EDI range
  hold(hca,'on')
  all_edi_plusminus = [v_edi_minus; v_edi_plus; -v_edi_minus; -v_edi_plus]*[1 1];   
  plot(hca,all_edi_plusminus*1e-6,hca.YLim,'k-.')
  %irf_legend(hca,'EDI',[0.55 + 0.5*v_edi_plus*1e-6/hca.XLim(2) 0.5],[0 0 0])   
  hold(hca,'off')
end

%% phi and flux

