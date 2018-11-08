%% Separatrix streaming
% Make reduced distribution
tintZoom = irf.tint('2017-07-06T13:53:50.00Z',5);
eint = [000 40000];
vint = [-Inf Inf];

eDist = ePDist1.tlim(tintZoom).elim(eint);
iDist = iPDist1.tlim(tintZoom).elim(eint);
ve = gseVe1.tlim(eDist.time).resample(eDist);
vi = gseVi1.tlim(iDist.time).resample(iDist);
scpot_margin = 1.5;
scpot_lim = scPot1.resample(eDist)*scpot_margin;
eLine = dmpaB1.resample(eDist).norm;
ePlane1 = eLine.cross(irf.ts_vec_xyz(eLine.time,repmat([1 0 0],eLine.length,1)));
ePlane2 = eLine.cross(ePlane1);

iLine = dmpaB1.resample(iDist).norm;
iPlane1 = iLine.cross(irf.ts_vec_xyz(iLine.time,repmat([1 0 0],iLine.length,1)));
iPlane2 = iLine.cross(iPlane1);

if2D = iDist.reduce('2D',iPlane1,iPlane2,'vint',vint); % reduced distribution perp to B
ef2D = eDist.reduce('2D',ePlane1,ePlane2,'vint',vint); % reduced distribution perp to B
%%
tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot_lim); toc % reduced distribution along B
tic; if1D = iDist.reduce('1D',iLine,'vint',vint); toc % reduced distribution along B
lineVe = ve.dot(eLine); % projection of Vi on B
lineVi = vi.dot(iLine); % projection of Vi on B


%% Plot
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
if 1 % eDEF omni
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
  irf_legend(hca,[num2str(scpot_margin) 'V_{sc}'],[0.99 0.1],'color',0*[1 1 1])
end
if 0 % e psd vpar
  isub = isub + 1;
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  hold(hca,'on')
  irf_plot(hca,{lineVe*1e-3},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_e (km/s)'; 
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
h(5).CLim = [-35 -28]+12;
colormap('jet');

%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];

%% Density perturbation
npanels = 8;
h = irf_plot(npanels);
dt = [0.0000  -0.0012  -0.0009  -0.0012];
dt = dt+0.0008;
  
if 1 % E par, 4 sc
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 1 % E par, 4 sc, time shifted for visibility
  hca = irf_panel('E par dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1par,gseE2par,gseE3par,gseE4par},'comp','dt',dt);
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{'mms1';'mms2';'mms3';'mms4'},[1.02 0.9],'fontsize',12);
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
  hca = irf_panel('Phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  hh = irf_plot(hca,{phi1,phi2,phi3,phi4},'comp','dt',dt);  
  hca.YLabel.String = {'\phi','(V)'};  
  ylabel(hca,hca.YLabel.String,'interpreter','tex')
  %set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['dt = [' num2str(dt(1)*1e3,format_ms)],num2str(dt(2)*1e3,format_ms),num2str(dt(3)*1e3,format_ms),num2str(dt(4)*1e3,format_ms),'] ms'},[0.01 0.1],'fontsize',12);
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