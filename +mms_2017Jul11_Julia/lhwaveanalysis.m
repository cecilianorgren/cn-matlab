% E/B/phi lower hybrid analysis
tintZoom = irf.tint('2017-07-11T22:34:01.00Z/2017-07-11T22:34:03.00Z'); %20151112071854
%tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
%tintZoomThin = irf.tint('2015-11-30T00:24:27.10Z/2015-11-30T00:24:27.15Z');
%tintZoomThin = irf.tint('2015-11-30T00:24:26.68Z/2015-11-30T00:24:26.75Z');
%tintZoomThin = irf.tint('2015-11-30T00:22:36.50Z/2015-11-30T00:22:37.50Z');
%tintZoomThin = irf.tint('2015-11-30T00:23:16.00Z/2015-11-30T00:23:17.00Z');
tintZoom = irf.tint('2017-07-11T22:33:28.50Z/2017-07-11T22:33:30.00Z'); %20151112071854
tintLH = tintZoom;
tintUTC = tintLH.utc;
ffilt = 3;
mmsid = 1;
c_eval('[phiEB,vbest,dirbest,thetas,corrs] = mms.lhwaveanalysis(tintLH,gseE?,gseB?scm,gseB?,ne?,''plot'',1,''lhfilt'',ffilt);',mmsid)

%% Separatrix LH waves
ic = 1:4;
tint = irf.tint('2017-07-11T22:33:28.00Z/2017-07-11T22:33:31.00Z'); % LH separatrix
t_center = irf_time('2017-07-11T22:33:29.00','utc>epochtt'); % LH separatrix center time
v_timing = 1.34e+03*[-0.85 -0.48  0.24];
c_eval('tsV?_timing = irf.ts_vec_xyz(gseE?.time,repmat(v_timing,gseE?.length,1));',ic)
v_direction = irf_norm(v_timing);
v_amplitude = sqrt(sum(v_timing.^2));


c_eval('gseEdt? = irf_integrate(gseE?perp,t_center); gseEdt? = irf.ts_vec_xyz(gseEdt?.time,gseEdt?.data);',ic)
phi_filt = 3;
c_eval('gsePhi? = gseEdt?.dot(tsV?_timing); gsePhi?_filt = gsePhi?.filt(phi_filt,0,[],3);',ic)

%% New coordinate system
tintLH = irf.tint('2017-07-11T22:33:28.60Z/2017-07-11T22:33:29.30Z'); % LH wave packet

newz = irf_norm(mean(gseB1.tlim(tintLH).data,1));
newx = cross(newz,cross(v_direction,newz));
newy = cross(newz,newx);
newxyz = [newx;newy;newz];

ic_tmp=ic;
ic = 1:4;
c_eval('bdryB? = gseB?*newxyz''',ic);
c_eval('bdryE? = gseE?*newxyz''',ic);
c_eval('bdryE?perp = gseE?perp*newxyz''',ic);
c_eval('bdryVe? = gseVe?*newxyz''',ic);
c_eval('bdryVi? = gseVi?*newxyz''',ic);
ic = ic_tmp;
%% Plot potential
npanels = 7;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % ne ni
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 0 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z,gseVe?par},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{-1*gseVexB?.x,-1*gseVexB?.y,-1*gseVexB?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx,(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E par
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,gseE?par);',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % Phi
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,gsePhi?_filt);',ic)
  hca.YLabel.String = {'Phi','(V)'};
  set(hca,'ColorOrder',mms_colors('12'))     
  irf_legend(hca,sprintf('v_{timing} = %.0f x [%.2f %.2f %.2f] km/s',v_amplitude, v_direction),[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f_{filt} = %g Hz',phi_filt),[0.98 0.10],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,ne?);',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};   
end
if 1 % sc Pot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,-1*scPot?);',ic)
  hca.YLabel.String = {'-scPot','(V)'};  
end
if 1 % sc Pot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,-1*scPot?);',ic)
  hca.YLabel.String = {'-scPot','(V)'};  
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',irf.tint('2017-07-11T22:33:28.00Z/2017-07-11T22:33:31.00Z')); % LH)
irf_zoom(h(:),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h(ii).FontSize = 12;
end
%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
%% Plot fields and potential in new coordinate system
npanels = 6;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{bdryB?.x,bdryB?.y,bdryB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);
end
if 0 % J  
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);     
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryVe?.x.tlim(tint),bdryVe?.y.tlim(tint),bdryVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{-1*gseVexB?.x,-1*gseVexB?.y,-1*gseVexB?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx,(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{bdryE?.x,bdryE?.y,bdryE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','norm','||'},[0.98 0.9],'fontsize',12);     
  irf_zoom(hca,'y')
end
if 0 % E par
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,gseE?par);',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % Phi
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,gsePhi?_filt);',ic)
  hca.YLabel.String = {'Phi','(V)'};
  set(hca,'ColorOrder',mms_colors('12'))     
  irf_legend(hca,sprintf('v_{timing} = %.0f x [%.2f %.2f %.2f] km/s',v_amplitude, v_direction),[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f_{filt} = %g Hz',phi_filt),[0.98 0.10],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))    
end
if 1 % sc Pot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,-1*scPot?);',ic)
  hca.YLabel.String = {'-scPot','(V)'};  
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',tint); % LH)
irf_zoom(h(:),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h(ii).FontSize = 12;
end
%% Compare 4 sc in new coordinate system
%% Figure 2: 4 sc 
npanels = 7;
h = irf_plot(npanels);

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 0 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{bdryB1.x.tlim(tint),bdryB2.x.tlim(tint),bdryB3.x.tlim(tint),bdryB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{v}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{bdryB1.y.tlim(tint),bdryB2.y.tlim(tint),bdryB3.y.tlim(tint),bdryB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{bdryB1.z.tlim(tint),bdryB2.z.tlim(tint),bdryB3.z.tlim(tint),bdryB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{||}','(nT)'};
end
if 1 % EX
  hca = irf_panel('EX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{bdryE1.x.tlim(tint),bdryE2.x.tlim(tint),bdryE3.x.tlim(tint),bdryE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{v}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % EY
  hca = irf_panel('EY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{bdryE1.y.tlim(tint),bdryE2.y.tlim(tint),bdryE3.y.tlim(tint),bdryE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{N}','(nT)'};
end
if 1 % EZ
  hca = irf_panel('EZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{bdryE1.z.tlim(tint),bdryE2.z.tlim(tint),bdryE3.z.tlim(tint),bdryE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{||}','(nT)'};
end
if 1 % B wave par
  hca = irf_panel('B scm par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1scmpar.filt(phi_filt,0,[],3),gseB2scmpar.filt(phi_filt,0,[],3),gseB3scmpar.filt(phi_filt,0,[],3),gseB4scmpar.filt(phi_filt,0,[],3)},'comp');
  irf_legend(hca,sprintf('f_{filt} = %g Hz',phi_filt),[0.98 0.10],'fontsize',12);
  hca.YLabel.String = {'B_{SCM}','(nT)'};  
end
if 1 % Phi
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gsePhi1_filt,gsePhi2_filt,gsePhi3_filt,gsePhi4_filt},'comp');
  hca.YLabel.String = {'Phi','(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))       
  irf_legend(hca,sprintf('v_{timing} = %.0f x [%.2f %.2f %.2f] km/s',v_amplitude, v_direction),[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f_{filt} = %g Hz',phi_filt),[0.98 0.10],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1,ne2,ne3,ne4},'comp');
  hca.YLabel.String = {'n_e','(cm^{-3})'};   
end
if 1 % sc Pot
  hca = irf_panel('scPot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-1*scPot1,-1*scPot2,-1*scPot3,-1*scPot4},'comp');
  hca.YLabel.String = {'-scPot','(V)'};  
end

irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Quiver plot
%tint = irf.tint('2017-07-11T22:33:28.00Z/2017-07-11T22:33:31.00Z');
nrows = 1;
ncols = 1;
npanels = nrows*ncols;
clear h;
isub = 0;
for irow = 1:nrows
  for icol = 1:ncols
    isub = isub + 1;
    h(isub) = subplot(nrows,ncols,isub);
  end
end


isub = 1;

tintQuivers = irf.tint('2017-07-11T22:33:28.70Z/2017-07-11T22:33:29.20Z');

if 0 % E perp
  hca = h(isub); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?perp.x,gseE?perp.y,gseE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_\perp','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
  irf_zoom(hca,'x',tintQuivers)
end
if 1 % GSE, 1 velocity
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')  
  times = gseE1.tlim(tintQuivers).time(1:50:end);
  gseVmatch = v_timing; % GSE  
  gseV = gseVmatch;
  c_eval('posR? = repmat(gseRR?,times.length,1)-(times-times(1))*gseV;')
  c_eval('posV?perp = gseVe?perp.resample(times).data;')
  c_eval('posV? = gseVe?.resample(times).data;')
  c_eval('posB? = gseB?.resample(times).data;')
  c_eval('posE?perp = gseE?perp.resample(times).data;')
  %c_eval('plot_quivers(hca,[posV?perp(:,1) posV?perp(:,2) posV?perp(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  %c_eval('plot_quivers(hca,[posV?(:,1) posV?(:,2) posV?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posE?perp(:,1) posE?perp(:,2) posE?perp(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  %c_eval('plot_quivers(hca,[posB?(:,1) posB?(:,2) posB?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''b''))')
  hold(hca,'off')

  hca.XLabel.String = 'x_{GSE}';
  hca.YLabel.String = 'y_{GSE}';
  hca.ZLabel.String = 'z_{GSE}';
  %hca.YDir = 'normal';
  %hca.XDir = 'reverse';
  axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
end
%%
if 0 % 3D
   
  times = mvaVe1.tlim(tintQuivers).time;
  
  gseVouteredge = 65.3*[0.37  0.32 -0.87];  
  gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
  gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
  gseMSP = 32*[-0.88 -0.42 0.24];
  
  lmnVmsh = gseVinneredge*[L' M' N'];
  lmnVmsp = gseMSP*[L' M' N'];
  
  vel_selection = 5;
  clear timesVUTC gseVdata
  gseV
  switch vel_selection
    case 1 % do velocities manually    
      tanV = irf_norm(cross(gseV,M));
      gseV = 55*[-0.90 -0.28 -0.33]-30*tanV;

      gseV = repmat(gseV,times.length,1);
      gseV(1:27,:) = repmat(-50*L,27,1);   
      gseV(28:31,:) = repmat(-50*N,4,1);   
      gseV(44:60,:) = repmat(gseVouteredge,17,1);  
    case 2 % define velocities at certain times, and then interpolate to other times      
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-25*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 3 % more 'vertical'
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 4
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 5
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*dot(-cross(irf_norm(gseMSP),M),gseVinneredge)...
                                                                                 -cross(irf_norm(gseMSP),M)*33*0.75;      
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge*0.75; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 6
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*dot(-cross(irf_norm(gseMSP),M),gseVinneredge); 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge;-40*tanVinneredge; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge;+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
  end
  
  dt = times(2)-times(1);
  
  lmnV = gseV.data*[L' M' N'];
  
  
  c_eval('posR? = repmat(mvaRR?,times.length,1)-dt*cumsum(lmnV,1);')
  c_eval('posV? = mvaVe?.resample(times).data*0.6;')
  c_eval('posJ? = mvaJ?.resample(times).data;')
  c_eval('posB? = mvaB?.resample(times).data;')
  c_eval('posE? = mvaE?.resample(times).data;')
  c_eval('posRe? = mvaEVexB?.resample(times).data;')
  
  hca = h2(isub); isub = isub +1;
  sclist = 1:4;
  hold(hca,'on')
  %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posV?(:,3) -posV?(:,2) posV?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
  c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
  %c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
  %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
  hold(hca,'off')

  hca.XLabel.String = 'N (km)';
  hca.YLabel.String = 'M (km)';
  hca.ZLabel.String = 'L (km)';
  hca.ZDir = 'normal';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  %axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  view(hca,[0 -1 0])
  %hca.YLim = [-20 20];
  %axis(hca,'equal')
  axis(hca,'equal')
  
  hca.Box = 'on';
  hca.ZLim = [-10 75];
  hca.XLim = [-8 27];
  tintUTCstart = tintQuivers(1).utc;
  tintUTCstop = tintQuivers(2).utc;
  hca.Title.String = ['V_e: ' tintUTCstop(12:22) ' - ' tintUTCstop(12:22)];
  hca.Title.String = ['V_e'];
  hleg=legend(hca,'MMS 1','MMS 2','MMS 3','MMS 4','location','northeast');
  fontsize = 17;
  hca.XLabel.FontSize = fontsize;
  hca.YLabel.FontSize = fontsize;
  hca.ZLabel.FontSize = fontsize;
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  hleg.FontSize = 10;
  if 0 % plot electric field arrows
    hca = h2(isub); isub = isub +1;
    hold(hca,'on')
    %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
    c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
    c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
    c_eval('plot_quivers(hca,[posJ?(:,3) -posJ?(:,2) posJ?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
    %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',1:4)
    %c_eval('plot_quivers(hca,[posRe?(:,3) -posRe?(:,2) posRe?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
    hold(hca,'off')

    hca.XLabel.String = 'N';
    hca.YLabel.String = 'M';
    hca.ZLabel.String = 'L';
    hca.ZDir = 'normal';
    hca.YDir = 'normal';
    hca.XDir = 'reverse';
    %axis(hca,'square')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.ZGrid = 'on';
    view(hca,[0 -1 0])
    %hca.YLim = [-20 20];
    %axis(hca,'equal')
    %axis(hca,'equal')
    %hca.XLim = 20*[-1 1];
    hca.Title.String = 'E, J (faint)';
  end
  if 0 % Interpolate BM in LN-plane to get Hall field color surface
    posL = double([posR1(:,1); posR2(:,1); posR3(:,1); posR4(:,1)]);
    posM = double([posR1(:,2); posR2(:,2); posR3(:,2); posR4(:,2)]);
    posN = double([posR1(:,3); posR2(:,3); posR3(:,3); posR4(:,3)]);
    posB = double([posB1(:,2); posB2(:,2); posB3(:,2); posB4(:,2)]);
    dN = 2; dL = 2;
    [NN,LL] = meshgrid(min(posN):dN:max(posN),min(posL):dL:max(posL));
    fBM = griddata(posN,posL,posB,NN,LL);
    hold(hca,'on');
    mesh(hca,NN,NN*0+abs(min(posM)),LL,fBM);  
    hmcb = colorbar('peer',hca); 
    %plot3(hca,posN,posM,posL,posB,'o');
    hold(hca,'off');
  end
  
  
  if 1 % plot magnetosheath and magnetosphere boundary planes    
    gseVouteredge = 65.3*[0.37  0.32 -0.87];  
    gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
    gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
    gseMSP = 32*[-0.88 -0.42 0.24];

    lmnVmsh = gseVinneredge*[L' M' N'];
    lmnVmsp = gseMSP*[L' M' N'];
     
    mspN = irf_norm(lmnVmsp);
    mshN = irf_norm(lmnVmsh);

    x = 90*[-1 1];
    y = 90*[-1 1];
    z = 90*[-1 1];
    
    funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
    funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
    funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);


    if exist('hmshN1'); delete(hmshN1); end
    if exist('hmshN2'); delete(hmshN2); end
    if exist('hmspN1'); delete(hmspN1); end
    if exist('hmspN2'); delete(hmspN2); end
    if exist('ht1'); delete(ht1); end
    if exist('ht2'); delete(ht2); end
    if exist('ht3'); delete(ht3); end
    if exist('ht4'); delete(ht4); end
   
    hold(hca,'on')
    if 1
      hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+70,'k-.');
      hmspN1 = plot3(hca,funZ(x,y,mspN)-2.8,x*0,x,'k-');
    else
      hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+72,'k-.');
      hmspN1 = plot3(hca,funZ(x,y,mspN)+20,x*0,x,'k-');
      hca.XLim = [0 50];
      hca.ZLim = [-10 40];
    end
    hold(hca,'off')

    ht1 = text(21.5,00,43,'MSH'); ht1.HorizontalAlignment = 'center'; ht1.FontSize = 13; ht1.Rotation = 55; 
    ht2 = text(2,0,50,'MSP'); ht2.HorizontalAlignment = 'center'; ht2.FontSize = 13; ht2.Rotation = -80;
    
    ht3 = text(16,0,-8,tintUTCstart(12:22)); ht3.HorizontalAlignment = 'center'; ht3.FontSize = 13;
    ht4 = text(16,0,72,tintUTCstop(12:22)); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
    
    hold(hca,'on') 
    quiver3(hca,17,0,4,0,0,5,2,'k')
    hold(hca,'off')
    ht4 = text(17,0,0,'time'); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
     
  end
end


%%
tintZoom = tintZoomThin + 0.05*[-1 1];
x = dirbest;
z = mean(gseB1.tlim(tintZoomThin).data,1); z = z/norm(z);
y = cross(z,x);

vfactor = 3;
dcvPot1_ = irf.ts_scalar(dcvPot1.time,dcvPot1.data(:,1:4));

lmnE1 = gseE1*[x;y;z]'; lmnE1.name = 'fac E';

h=irf_plot({gseB1.tlim(tintZoom),...
            gseB1scm.tlim(tintZoom),...
            gseE1perp.tlim(tintZoom),...
            gseE1par.tlim(tintZoom),...
            lmnE1.tlim(tintZoom),...
            phiEB.*irf.ts_scalar(phiEB.time,repmat([vfactor 1],phiEB.length,1)),...
            -gseVexB1.tlim(tintZoom),...
            gseVe1.tlim(tintZoom),...
            gseVe1.tlim(tintZoom)+-dirbest*vbest*vfactor,...
            gseVe1.tlim(tintZoom)-gseVe1.resample(irf_time('2015-11-30T00:24:26.80Z','utc>epochtt')).data,...
            scPot1.tlim(tintZoom),...
            dcvPot1.tlim(tintZoom),...
            ne1.tlim(tintZoom)});
h(1).YLabel.String ='B';          
h(2).YLabel.String ='scm B';          
irf_legend(h(5),{'Ek','E norm','Epar'},[0.95 0.05])
h(6).YLabel.String ='Phi'; irf_legend(h(6),{'phi E','phi B'},[0.95 0.05])
h(7).YLabel.String ='-vexB';
h(8).YLabel.String ='ve';
h(9).YLabel.String ='ve-vph';
h(10).YLabel.String ='ve-veref';
h(11).YLabel.String ='scPot';

add_length_on_top(h(1),vbest*vfactor,1)