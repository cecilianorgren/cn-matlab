%% Overview
%tintZoom = irf.tint('2017-07-11T22:32:50.00Z/2017-07-11T22:35:20.00Z'); %20151112071854
tint = irf.tint('2017-07-09T04:41:43.00Z/2017-07-09T04:43:32.00Z');
tintZoom = tint;
tint = tintZoom;

npanels = 9;

cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');
cmap = irf_colormap('waterfall');
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 0 % B gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % J 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x,gseJ?.y,gseJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
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
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?perp.x,gseVi?perp.y,gseVi?perp.z,gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};
  %hca.YLim = [-1100 1100];  
end
if 0 % Ve x B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % V ExB
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
end
if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
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
  irf_zoom(hca,'y')
end
if 1 % E par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};  
end
if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % e DEF omni
  isub = isub + 1;
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
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 1 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % ePDist pa 64
  isub = isub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 1 % Ti par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTi?.xx.tlim(tint),(facTi?.yy+facTi?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{i,||}','T_{i,\perp}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Plasma parameters, beta, pressure etc
ic = 1;
npanels = 9;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:38:00.00Z'); %20151112071854
tint = tintZoom;

zoomy = [];
isub = 0;
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
end
if 1 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 1 % Vi 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % iPDist deflux omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTi?.trace/3,''k'');',ic)
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
end
if 1 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet') 
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTe?.trace/3,''k'');',ic)
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 0 % ne ni
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % beta
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{beta1.tlim(tint),beta1i.tlim(tint),beta1e.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('\beta',ic),'(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'i+e','i','e'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log';
  hca.YLim = [1e-2 1e3];
  hca.YTick = 10.^[-1:2];
end
if 1 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, stacked
  hca = irf_panel('Pressure stacked');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_patch(hca,{gsePe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'gsePi?.trace/3,PB?,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3});'],ic)
  irf_patch(hca,{PB1,gsePi1.trace/3,gsePe1.trace/3})
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B','P_{tot}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures, total, 4sc
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3+PB1.resample(gsePe1)+gsePi1.resample(gsePe1).trace/3,...    
                gsePe2.trace/3+PB2.resample(gsePe2)+gsePi1.resample(gsePe2).trace/3,...    
                gsePe3.trace/3+PB3.resample(gsePe3)+gsePi1.resample(gsePe3).trace/3,...    
                gsePe4.trace/3+PB4.resample(gsePe4)+gsePi1.resample(gsePe4).trace/3},'comp');
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Ion pressures, total, 4sc
  hca = irf_panel('Ion pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePi1.trace/3,gsePi2.trace/3,gsePi3.trace/3,gsePi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron pressures, total, 4sc
  hca = irf_panel('Pe, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3,gsePe2.trace/3,gsePe3.trace/3,gsePe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron temperature, total, 4sc
  hca = irf_panel('Te, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseTe1.trace/3,gseTe2.trace/3,gseTe3.trace/3,gseTe4.trace/3},'comp');
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPi
  hca = irf_panel('gradPi');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPi.x*1e3,gseGradPi.y*1e3,gseGradPi.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_i','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v','N','B'},[0.98 0.9],'fontsize',12);  
end
if 0 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x,mvaE?perp.y,mvaE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  hca = irf_panel('E lowf');
  ffiltE = 0.5;
  c_eval('gseE?highf = gseE?.filt(3,0,[],3);',ic)
  c_eval('gseE?lowf = gseE?-gseE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?lowf.x,gseE?lowf.y,gseE?lowf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f<%g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % Lp from Eperp
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Lp?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'L_p','(km)'};  
  hca.YScale = 'log';
  hca.YLim = [1e2 5e3];
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 0 % J 
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x,mvaJ?.y,mvaJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi perp
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?perp.x,mvaVi?perp.y,mvaVi?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);   
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
for ii = 1:npanels
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
  h(ii).FontSize = 14;  
  h(ii).YLabel.FontSize = 13;
end

%% Overview, not so many panels
ic = 1;
npanels = 7;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854
tint = tintZoom;

zoomy = [];
isub = 0;
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.1],'fontsize',12);  
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'color',[0 0 0],'fontsize',12);    
end
if 1 % Vi 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % iPDist deflux omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTi?.trace/3,''k'');',ic)
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
  ylabel(hca,{'E_i','(eV)'},'interpreter','tex');  
end
if 1 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet') 
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTe?.trace/3,''k'');',ic)    
  hold(hca,'off')
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 0 % ne ni
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 0 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % beta
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{beta1.tlim(tint),beta1i.tlim(tint),beta1e.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('\beta',ic),'(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'i+e','i','e'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log';
  hca.YLim = [1e-2 1e3];
  hca.YTick = 10.^[-1:2];
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, stacked
  hca = irf_panel('Pressure stacked');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_patch(hca,{gsePe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'gsePi?.trace/3,PB?,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3});'],ic)
  irf_patch(hca,{PB1,gsePi1.trace/3,gsePe1.trace/3})
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B','P_{tot}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures, total, 4sc
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3+PB1.resample(gsePe1)+gsePi1.resample(gsePe1).trace/3,...    
                gsePe2.trace/3+PB2.resample(gsePe2)+gsePi1.resample(gsePe2).trace/3,...    
                gsePe3.trace/3+PB3.resample(gsePe3)+gsePi1.resample(gsePe3).trace/3,...    
                gsePe4.trace/3+PB4.resample(gsePe4)+gsePi1.resample(gsePe4).trace/3},'comp');
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Ion pressures, total, 4sc
  hca = irf_panel('Ion pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePi1.trace/3,gsePi2.trace/3,gsePi3.trace/3,gsePi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron pressures, total, 4sc
  hca = irf_panel('Pe, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3,gsePe2.trace/3,gsePe3.trace/3,gsePe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron temperature, total, 4sc
  hca = irf_panel('Te, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseTe1.trace/3,gseTe2.trace/3,gseTe3.trace/3,gseTe4.trace/3},'comp');
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPi
  hca = irf_panel('gradPi');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPi.x*1e3,gseGradPi.y*1e3,gseGradPi.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_i','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','N','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x,mvaE?perp.y,mvaE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E low f
  hca = irf_panel('E lowf');
  ffiltE = 2;
  c_eval('gseE?highf = gseE?.filt(3,0,[],3);',ic)
  c_eval('gseE?lowf = gseE?-gseE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?lowf.x,gseE?lowf.y,gseE?lowf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f<%g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % E high f
  hca = irf_panel('E highf');
  ffiltE = 2;
  c_eval('gseE?highf = gseE?.filt(3,0,[],3);',ic)
  c_eval('gseE?lowf = gseE?-gseE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?highf.x,gseE?highf.y,gseE?highf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f<%g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 1 % B scm filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B scm');
  ffilt = 1;
  c_eval('plotB = mvaB?scm.filt(ffilt,0,[],3);',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{plotB.x,plotB.y,plotB.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,sprintf('f > %g Hz',ffilt),[0.98 0.1],'fontsize',12);
end
if 0 % Lp from Eperp
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Lp?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'L_p','(km)'};  
  hca.YScale = 'log';
  hca.YLim = [1e2 5e3];
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 0 % J 
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x,mvaJ?.y,mvaJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi perp
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?perp.x,mvaVi?perp.y,mvaVi?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);   
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
for ii = 1:npanels
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
  h(ii).FontSize = 14;  
  h(ii).YLabel.FontSize = 13;
end

%% Boundary properties, modulation, thickness, current
ic = 1;
npanels = 8;
h = irf_plot(npanels);
tintZoom = irf.tint('2017-07-11T22:32:00.00Z/2017-07-11T22:33:40.00Z'); %20151112071854
tint = tintZoom;

zoomy = [];
isub = 0;
cmap = 'jet';
% B, n, J, eDEF, ePitchDEF, Eperp lowf, 
if 1 % B  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.1],'fontsize',12);  
end
if 1 % ne
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'color',[0 0 0],'fontsize',12);    
end
if 0 % J mom 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % iPDist deflux omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet')   
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTi?.trace/3,''k'');',ic)
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};  
  ylabel(hca,{'E_i','(eV)'},'interpreter','tex');  
end
if 1 % ePDist deflux omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,'jet') 
  hold(hca,'on')
  c_eval('irf_plot(hca,gseTe?.trace/3,''k'');',ic)    
  hold(hca,'off')
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'};  
end
if 1 % ePDist pa 64
  isub = isub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 10000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % J mom  par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J mom par');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{mvaJ?par},''comp'');',ic)  
  hca.YLabel.String = {'J_{||}','(nA/m^2)'};
end
if 0 % ne ni
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
end
if 1 % J mom perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J mom perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?perp.x.tlim(tint),mvaJ?perp.y.tlim(tint),mvaJ?perp.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J_\perp','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % J perp low f
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Jperp lowf');
  ffiltJ = 0.5;
  c_eval('mvaJ?highf = mvaJ?perp.filt(ffiltJ,0,[],3);',ic)
  c_eval('mvaJ?lowf = mvaJ?perp-mvaJ?highf;',ic)
  %c_eval('mvaJ?lowf = mvaJ?perp.filt(0,ffiltJ,[],3);',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?lowf.x,mvaJ?lowf.y,mvaJ?lowf.z},''comp'');',ic)
  hca.YLabel.String = {'J_{\perp}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f < %g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % beta
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{beta1.tlim(tint),beta1i.tlim(tint),beta1e.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('\beta',ic),'(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'i+e','i','e'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log';
  hca.YLim = [1e-2 1e3];
  hca.YTick = 10.^[-1:2];
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, stacked
  hca = irf_panel('Pressure stacked');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_patch(hca,{gsePe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'gsePi?.trace/3,PB?,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3});'],ic)
  irf_patch(hca,{PB1,gsePi1.trace/3,gsePe1.trace/3})
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B','P_{tot}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures, total, 4sc
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3+PB1.resample(gsePe1)+gsePi1.resample(gsePe1).trace/3,...    
                gsePe2.trace/3+PB2.resample(gsePe2)+gsePi1.resample(gsePe2).trace/3,...    
                gsePe3.trace/3+PB3.resample(gsePe3)+gsePi1.resample(gsePe3).trace/3,...    
                gsePe4.trace/3+PB4.resample(gsePe4)+gsePi1.resample(gsePe4).trace/3},'comp');
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Ion pressures, total, 4sc
  hca = irf_panel('Ion pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePi1.trace/3,gsePi2.trace/3,gsePi3.trace/3,gsePi4.trace/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron pressures, total, 4sc
  hca = irf_panel('Pe, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gsePe1.trace/3,gsePe2.trace/3,gsePe3.trace/3,gsePe4.trace/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Electron temperature, total, 4sc
  hca = irf_panel('Te, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_plot(hca,{gseTe1.trace/3,gseTe2.trace/3,gseTe3.trace/3,gseTe4.trace/3},'comp');
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % gradPi
  hca = irf_panel('gradPi');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPi.x*1e3,gseGradPi.y*1e3,gseGradPi.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_i','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','N','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % E perp
  hca = irf_panel('E perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x,mvaE?perp.y,mvaE?perp.z},''comp'');',ic)
  hca.YLabel.String = {'E_{\perp}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E low f
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E lowf');
  ffiltE = 1;
  c_eval('mvaE?highf = mvaE?.filt(ffiltE,0,[],3);',ic)
  c_eval('mvaE?lowf = mvaE?-mvaE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?lowf.x,mvaE?lowf.y,mvaE?lowf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f < %g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % E high f 
  hca = irf_panel('E highf par');
  ffiltE = 50;
  c_eval('mvaE?highf = mvaE?par.filt(ffiltE,0,[],3);',ic)  
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?highf},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f > %g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % E high f
  hca = irf_panel('E highf');
  ffiltE = 50;
  c_eval('mvaE?highf = mvaE?.filt(ffiltE,0,[],3);',ic)
  c_eval('mvaE?lowf = mvaE?-mvaE?highf;',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?highf.x,mvaE?highf.y,mvaE?highf.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,sprintf('f > %g Hz',ffiltE),[0.98 0.1],'fontsize',12);
end
if 0 % B scm filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B scm');
  ffilt = 1;
  c_eval('plotB = mvaB?scm.filt(ffilt,0,[],3);',ic)
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{plotB.x,plotB.y,plotB.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,sprintf('f > %g Hz',ffilt),[0.98 0.1],'fontsize',12);
end
if 0 % Lp from Eperp
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Lp?*1e-3},''comp'');',ic)
  hca.YLabel.String = {'L_p','(km)'};  
  hca.YScale = 'log';
  hca.YLim = [1e2 5e3];
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 0 % J 
  hca = irf_panel('J fpi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x,mvaJ?.y,mvaJ?.z},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 0 % Vi perp
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?perp.x,mvaVi?perp.y,mvaVi?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve x B
  hca = irf_panel('VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'-v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x,mvaVExB?.y,mvaVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);   
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
for ii = 1:npanels
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0])
  h(ii).FontSize = 14;  
  h(ii).YLabel.FontSize = 13;
end

