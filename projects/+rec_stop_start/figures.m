tint = irf.tint('2017-07-25T22:04:44.822634526Z/2017-07-25T22:11:31.820198151Z');
tint_b = irf.tint('2017-07-25T22:09:00.82Z/2017-07-25T22:11:31.82Z');

tint_c = irf.tint('2017-07-25T22:09:00.82Z/2017-07-25T22:11:00.00Z');
tint_d = irf.tint('2017-07-25T22:09:30.82Z/2017-07-25T22:11:00.00Z');

tint_action = irf.tint('2017-07-25T22:09:40.00Z/2017-07-25T22:10:50.00Z');
tint_cavity = irf.tint('2017-07-25T22:09:46.00Z/2017-07-25T22:09:56.00Z');

%% Figure: Overview 1
ic = 1;

npanels = 5;
h = irf_plot(npanels);
iisub = 0;
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
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
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
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % E gse downsampled
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gseE?.resample(gseVi1).x,gseE?.resample(gseVi1).y,gseE?.resample(gseVi1).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Ve perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
 
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 0.5;
  %c_eval('irf_plot(hca,{gseVe?perp.x.filt(0,ffilt,[],5),gseVe?perp.y.filt(0,ffilt,[],5),gseVe?perp.z.filt(0,ffilt,[],5),gseVe?par.filt(0,ffilt,[],5)},''comp'');',ic)  
  c_eval('irf_plot(hca,{gseVe?perp.x.resample(iPDist?),gseVe?perp.y.resample(iPDist?),gseVe?perp.z.resample(iPDist?),gseVe?par.resample(iPDist?)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
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
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
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
if 0 % V ExB resample
  timeline = tint(1):0.15:tint(2);
  timeline = gseVi1;
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB resample');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.resample(timeline),gseVExB?.y.resample(timeline),gseVExB?.z.resample(timeline)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);   
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
  irf_zoom(hca,'y')
end
if 0 % E par
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
if 0 % Te/i par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('1234'))
  refTi = 1;
  c_eval('tepar_plot = Te?par;',ic)
  tepar_plot.data(tepar_plot.data<0) = NaN;
  c_eval('tepar_plot = tepar_plot.resample(gseTi?);',ic)
  
  c_eval('teper_plot = Te?perp;',ic)
  teper_plot.data(teper_plot.data<0) = NaN;
  c_eval('teper_plot = teper_plot.resample(gseTi?);',ic)
  
  c_eval('irf_plot(hca,{tepar_plot,teper_plot,Ti?par.tlim(tint),Ti?perp.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}','T_{i,||}','T_{i,\perp}'}',[1.01 0.90],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
  hca.YLim = [0 1e4];
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % ve_par
    hold(hca,'on')
    c_eval('vexb = gseVe?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  elim = [700 40000];
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  
  if 1 % elim
    %%
    hold(hca,'on')
    c_eval('timeline = iPDist?.time;',ic)    
    irf_plot(hca,irf.ts_scalar(timeline,10000*ones([timeline.length,1])),'k')
    hold(hca,'off')
  end
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 1 % ePDist pa 64
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
  set(hca,'ColorOrder',mms_colors('1234'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 0 % Ti par perp
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
if 0 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  c_eval('if1D = if1Dx?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z');
  c_eval('if1D = if1Dz?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x high');
  c_eval('if1D = if1Dx?_high;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y high');
  c_eval('if1D = if1Dy?_high;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z high');
  c_eval('if1D = if1Dz?_high;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x low');
  c_eval('if1D = if1Dx?_low;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y low');
  c_eval('if1D = if1Dy?_low;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z low');
  c_eval('if1D = if1Dz?_low;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % iPitch
  isub = isub + 1;
  hca = irf_panel('iPitchj');
  c_eval('iPitch = ePitch?.elim([5000 10000]);',ic)
  irf_spectrogram(hca,iPitch.specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure: V ExB
ic = 1;

Etop_fpi = iPDist1.ancillary.energy(1,end)+iPDist1.ancillary.delta_energy_plus(1,end);

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
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
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 1 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x,gseVExB?.y,gseVExB?.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
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
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
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
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
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

if 0 % e DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('e DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ele_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_e^{FEEPS}','(eV)'};   
end
if 0 % EIS Pitch 1 lowE
  isub = isub + 1;
  hca = irf_panel('eis Pitch 1');  
  elim = [0 32000];
  irf_spectrogram(hca,eis_pa2.elim(elim).specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 0 % EIS Pitch 2 highE
  isub = isub + 1;
  hca = irf_panel('eis Pitch 2');  
  elim = [32000 100000];
  irf_spectrogram(hca,eis_pa2.elim(elim).specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end

if 0 % EIS Pitch 0 allE oxygen
  isub = isub + 1;
  hca = irf_panel('eis Pitch ox');  
  %elim = [0 32000];
  irf_spectrogram(hca,eis_pa2_op.specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 1 % EIS Pitch 0 allE
  isub = isub + 1;
  hca = irf_panel('eis Pitch all E');  
  %elim = [0 32000];
  irf_spectrogram(hca,eis_pa2.specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 1 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure: momentum equation / ohms law
ic = 2;

npanels = 8;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
timeline = tint(1):0.5:tint(2);
doResample = 1;

isub = 0;
zoomy = [];

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, 4sc average
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_plot(hca,{PBav,Piav,Peav,PDiav,Piav + Peav.resample(Piav) + PBav.resample(Piav)},'comp')
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i}','P_{e}','P_{i,dyn,x}','P_{B}+P_i+P_e'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 1 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  if doResample
    irf_plot(hca,{gseJcurl.x.resample(timeline),gseJcurl.y.resample(timeline),gseJcurl.z.resample(timeline)},'comp');  
  else
    irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  end
  hca.YLabel.String = {'J^{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % JxB, gradB, curvB, 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    scale = 1e0*1e0; units = 'T Am^{-2}';
    scale = 1e9*1e9; units = 'nT nAm^{-2}';
    
    hca = irf_panel(['JxB ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))   
    irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale},'comp');  
    %irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale,-gseDivPb.(comp)*scale-gseCurvB.(comp)*scale},'comp');  
    hca.YLabel.String = {['JxB_' comp],['(' units ')']};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  end
end
if 1 % JxB, 1 panels: x,y,z  
  isub = isub + 1;
  zoomy = [zoomy isub];      
  hca = irf_panel('JxB xyz');
  scale = 1e9*1e9; units = 'nT nAm^{-2}';
  set(hca,'ColorOrder',pic_colors('matlab'))   
  irf_plot(hca,{gseJxB.x*scale,gseJxB.y*scale,gseJxB.z*scale},'comp');      
  hca.YLabel.String = {['JxB'],['(' units ')']};
  set(hca,'ColorOrder',pic_colors('matlab'))
  %irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  irf_legend(hca,{'x','y','z'}',[1.01 0.9],'fontsize',12);    
end
if 1 % JxB, 1 panels: LMN 
  isub = isub + 1;
  zoomy = [zoomy isub];      
  hca = irf_panel('JxB lmn');
  scale = 1e9*1e9; units = 'nT nAm^{-2}';
  set(hca,'ColorOrder',pic_colors('matlab'))   
  irf_plot(hca,{lmnJxB.x*scale,lmnJxB.y*scale,lmnJxB.z*scale},'comp');      
  hca.YLabel.String = {['JxB'],['(' units ')']};
  set(hca,'ColorOrder',pic_colors('matlab'))
  %irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'}',[1.01 0.9],'fontsize',12);    
end
if 0 % electron mom eq. 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',mms_colors('xyza'))
    to_plot_gseVexBav = gseVexBav;
    to_plot_gseVexBav.data(abs(to_plot_gseVexBav.data)>100) = NaN;
    irf_plot(hca,{gseEav.(comp),-1*to_plot_gseVexBav.(comp),-1*gseGradPene.resample(gseVi1).(comp),-1*(gseEav.(comp).resample(to_plot_gseVexBav)+to_plot_gseVexBav.(comp))},'comp');  
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'E','-v_exB','-divPe/ne','-(E+v_exB)'}',[1.02 0.9],'fontsize',12);  
  end
end
if 0 % ohm's law, 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))
    if doResample
      to_plot_gseVixBav = gseVixBav.resample(timeline);
    else
      to_plot_gseVixBav = gseVixBav;
    end
    to_plot_gseVixBav.data(abs(to_plot_gseVixBav.data)>100) = NaN;
    if doResample
      irf_plot(hca,{...
        gseEav.(comp).resample(timeline),...
        to_plot_gseVixBav.(comp).resample(timeline),...
        -1*gseGradPene.(comp).resample(timeline),...
        gseJxBne_mVm.(comp).resample(timeline),...
        +1*(gseEav.(comp).resample(timeline)+to_plot_gseVixBav.(comp).resample(timeline))},'comp');  
    else
      irf_plot(hca,{gseEav.(comp),to_plot_gseVixBav.(comp),-1*gseGradPene.(comp).resample(gseVi1),gseJxBne_mVm.(comp),+1*(gseEav.(comp).resample(to_plot_gseVixBav)+to_plot_gseVixBav.(comp))},'comp');  
    end
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'E','v_ixB','-divPe/ne','JxB/ne','(E+v_ixB)'}',[1.01 0.9],'fontsize',12);  
  end
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
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
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
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
if 0 % Ve gse
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve par
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
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z,gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {['v_{e} mms' num2str(ic)],'(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp par av
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseVeperpav.x,gseVeperpav.y,gseVeperpav.z,gseVeparav},'comp');
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end


if 1 % ePDist pa 64
  isub = isub + 1;
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [100 30000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%tmark_eis = tint(1):20:tint(2);

%% Figure: V ExB
ic = 1;

c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];

if 1 % B gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 0 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
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
cs = 'lmn';
if 0 % vExB.x , Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).x,%sVi?perp.x},''comp'');',cs,cs,cs),ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {sprintf('v_{\\perp,z}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x , Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).y,%sVi?perp.y},''comp'');',cs,cs,cs),ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {sprintf('v_{\\perp,y}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x , Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).z,%sVi?perp.z},''comp'');',cs,cs,cs),ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {sprintf('v_{\\perp,z}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 1 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 1 % vExB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('irf_plot(hca,{%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
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
if 1 % J 
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
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
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
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',ic)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure: Current sheet structure
ic = 1;

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];
comps = ['z'];
timeline = tint(1):0.5:tint(2);

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z,gseB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, 4sc average
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_plot(hca,{PBav,Piav,Peav,PDiav,Piav + Peav.resample(Piav) + PBav.resample(Piav)},'comp')
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i}','P_{e}','P_{i,dyn,x}','P_{B}+P_i+P_e'},[0.98 0.9],'fontsize',12);  
end
if 1 % Pressures, single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...                         
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{dyn,i}','P_i','P_e','P_{tot}'}',[1.01 0.9],'fontsize',12);  
end
if 0 % J xyz curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  hca.YLabel.String = {'J^{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % J xyz 
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
if 1 % Jy fpi curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Jy fpi curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.y.resample(timeline),gseJcurl.y.resample(timeline)},''comp'');',ic)  
  hca.YLabel.String = {'J_y','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'FPI','Curl'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-1100 1100];  
end
if 1 % JxB, 1-3 panels: x,y,z
  for comp = comps
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    scale = 1e0*1e0; units_ = 'T Am^{-2}';
    scale = 1e9*1e9; units_ = 'nT nAm^{-2}';
    
    %scale = 1e0*1e0; units_ = 'nT nAm^{-2}';
    
    hca = irf_panel(['JxB ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))   
    irf_plot(hca,{gseJxB.resample(timeline).(comp)*scale,-gseDivPb.resample(timeline).(comp)*scale,-gseCurvB.resample(timeline).(comp)*scale},'comp');  
    %irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale,-gseDivPb.(comp)*scale-gseCurvB.(comp)*scale},'comp');  
    hca.YLabel.String = {['JxB_' comp],['(' units_ ')']};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  end
end
if 0 % electron mom eq. 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',mms_colors('xyza'))
    to_plot_gseVexBav = gseVexBav;
    to_plot_gseVexBav.data(abs(to_plot_gseVexBav.data)>100) = NaN;
    irf_plot(hca,{gseEav.(comp),-1*to_plot_gseVexBav.(comp),-1*gseGradPene.resample(gseVi1).(comp),-1*(gseEav.(comp).resample(to_plot_gseVexBav)+to_plot_gseVexBav.(comp))},'comp');  
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'E','-v_exB','-divPe/ne','-(E+v_exB)'}',[1.02 0.9],'fontsize',12);  
  end
end
if 1 % ohm's law, 1-3 panels: x,y,z
  for comp = comps
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))
    to_plot_gseVixBav = gseVixBav;
    to_plot_gseVixBav.data(abs(to_plot_gseVixBav.data)>100) = NaN;
    irf_plot(hca,{gseEav.(comp),to_plot_gseVixBav.(comp),-1*gseGradPene.resample(timeline).(comp).resample(gseVi1),gseJxBne_mVm.resample(timeline).(comp),+1*(gseEav.(comp).resample(to_plot_gseVixBav)+to_plot_gseVixBav.(comp))},'comp');  
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'E','v_ixB','-divPe/ne','JxB/ne','(E+v_ixB)'}',[1.01 0.9],'fontsize',12);  
    hca.YLim = [-50 50];
  end
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
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
if 0 % Ve gse
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
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
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z,gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {['v_{e} mms' num2str(ic)],'(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp par av
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseVeperpav.x,gseVeperpav.y,gseVeperpav.z,gseVeparav},'comp');
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%tmark_eis = tint(1):20:tint(2);

%% Oxygen presence
% Load data: rec_stop_start.load_data_fast
%h = irf_plot({gseB1_fast,nOp1_srvy,nHp1_srvy,gseVi1_fast,iPDist1_fast.deflux.omni.specrec});
npanels = 7;
h = irf_plot(npanels);

if 1 % B gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?_srvy.x,gseB?_srvy.y,gseB?_srvy.z,gseB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % E gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseE?_fast.x,gseE?_fast.y,gseE?_fast.z},''comp'');',ic)
  hca.YLabel.String = {'E_{GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % n
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('1234'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.8])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 1 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?_fast.x,gseVi?_fast.y,gseVi?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % VExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_fast.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
%    irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end


irf_zoom(h,'x',tint_fast)
irf_plot_axis_align

%% MVA fast
% Load data: rec_stop_start.load_data_fast
%h = irf_plot({gseB1_fast,nOp1_srvy,nHp1_srvy,gseVi1_fast,iPDist1_fast.deflux.omni.specrec});
npanels = 8;
h = irf_plot(npanels);
doFilt = 1;
ffilt = 1;

if 1 % B gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?_srvy.x,gseB?_srvy.y,gseB?_srvy.z,gseB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 0 % B gsm fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?_srvy.x,gsmB?_srvy.y,gsmB?_srvy.z,gsmB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % B lmn fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?_srvy.x,lmnB?_srvy.y,lmnB?_srvy.z,lmnB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{LMN}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'  [L; M; N] = ',...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(1,1),R(1,2),R(1,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(2,1),R(2,2),R(2,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(3,1),R(3,2),R(3,3))}',[1.002 0.9],'fontsize',12,'color','k');  
end
if 1 % E gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  if doFilt
    c_eval('irf_plot(hca,{gseE?_fast.filt(0,ffilt,[],5).x,gseE?_fast.filt(0,ffilt,[],5).y,gseE?_fast.filt(0,ffilt,[],5).z},''comp'');',ic)
    irf_legend(hca,sprintf('f < %g Hz',ffilt),[0.02 0.07])
  else
    c_eval('irf_plot(hca,{gseE?_fast.x,gseE?_fast.y,gseE?_fast.z},''comp'');',ic)
  end
    
  hca.YLabel.String = {'E_{GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 1 % E lmn fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  if doFilt
    c_eval('irf_plot(hca,{lmnE?_fast.filt(0,ffilt,[],5).x,lmnE?_fast.filt(0,ffilt,[],5).y,lmnE?_fast.filt(0,ffilt,[],5).z},''comp'');',ic)    
    irf_legend(hca,sprintf('f < %g Hz',ffilt),[0.02 0.07])
  else
    c_eval('irf_plot(hca,{lmnE?_fast.x,lmnE?_fast.y,lmnE?_fast.z},''comp'');',ic)    
  end
  hca.YLabel.String = {'E_{LMN}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end 
if 1 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?_fast.x,gseVi?_fast.y,gseVi?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVi?_fast.x,lmnVi?_fast.y,lmnVi?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?_fast.x,gseVe?_fast.y,gseVe?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('vplot = lmnVe?_fast;',ic)
  vplot.data(vplot.abs.data>1e4,:) = NaN;
  c_eval('irf_plot(hca,{vplot.x,vplot.y,vplot.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVExB?_srvy.x,lmnVExB?_srvy.y,lmnVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % VExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % n log scale
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('1234'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.8])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_fast.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
    irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end


irf_zoom(h,'x',tint_fast)
irf_plot_axis_align

%% Fields and scatter plots
tint_scatter = tint_action;
tint_fpr1 = irf.tint('2017-07-25T22:10:00.000Z/2017-07-25T22:10:22.000Z');
tint_fpr2 = irf.tint('2017-07-25T22:10:25.000Z/2017-07-25T22:10:32.000Z');
tint_scatter = irf.tint('2017-07-25T22:09:45.000Z/2017-07-25T22:10:40.000Z');
%tint_scatter = tint_fpr2;

cs = 'gse'; comps = 'xyz';
cs = 'lmn'; comps = 'LMN';

npanels = 4;
nrows = 3;
ncols = 2;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.4,'vertical'); % horizontal

% TSeries plots
if 1 % B srvy
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel(['B ' cs]);
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?_srvy.x,%sB?_srvy.y,%sB?_srvy.z,%sB?_srvy.abs},''comp'');',cs,cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B^{%s}',cs),'(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3),'|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % E fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel(['E ' cs]);
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sE?_fast.x,%sE?_fast.y,%sE?_fast.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('E^{%s}',cs),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},[0.98 0.9],'fontsize',12);  
end 
if 1 % E +vixB fast+brst
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel(['E + vixB' cs]);
  set(hca,'ColorOrder',mms_colors('xyza'))  
  ts_yy = eval(sprintf('%sE%g_fast.resample(%sVi%g) + %sVi%g.cross(%sB%g_srvy.resample(%sVi%g))*1e3*1e-9*1e3;',cs,ic,cs,ic,cs,ic,cs,ic,cs,ic));
  irf_plot(hca,{ts_yy.x,ts_yy.y,ts_yy.z},'comp');
  hca.YLabel.String = {sprintf('(E + v_ixB)^{%s}',cs),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},[0.98 0.9],'fontsize',12);  
end 
if 1 % Vi fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel(['Vi ' cs]);
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sVi?_fast.x,%sVi?_fast.y,%sVi?_fast.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('v_{i}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},[0.98 0.9],'fontsize',12);  
end 

%irf_zoom(h1,'x',tint_fast)
irf_zoom(h1,'x',tint_scatter + [-10 10])
irf_zoom(h1,'y')
irf_plot_axis_align

hca = irf_panel(['E ' cs]); hca.YLim = [-25 49];

drawnow

% Scatter plots
hb = gobjects(0);
ihb = 1;
isub = 1;

if 1 % (Bx,By)
  hca = h2(isub); isub = isub + 1;
  
  ts_xx = eval(sprintf('%sB%g_srvy.tlim(tint_scatter).x;',cs,ic));
  ts_yy = eval(sprintf('%sB%g_srvy.y;',cs,ic));
  ts_yy = ts_yy.resample(ts_xx);
  
  t0 = ts_xx.time(1);
  tt = ts_xx.time - ts_xx.time(1);
  xx = ts_xx.data;
  yy = ts_yy.data;
  
  scatter(hca,xx,yy,5,tt)
  
  hca.XLabel.String = sprintf('B_%s (nT)',comps(1));
  hca.YLabel.String = sprintf('B_%s (nT)',comps(2));
  
  hcb = colorbar('peer',hca);
  hb(ihb) = hcb;
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 30*[-1 1];
  hca.XLim = 30*[-1 1];
end
if 1 % (By,Bz)
  hca = h2(isub); isub = isub + 1;
  
  ts_xx = eval(sprintf('%sB%g_srvy.tlim(tint_scatter).y;',cs,ic));
  ts_yy = eval(sprintf('%sB%g_srvy.z;',cs,ic));
  ts_yy = ts_yy.resample(ts_xx);
  
  t0 = ts_xx.time(1);
  tt = ts_xx.time - ts_xx.time(1);
  xx = ts_xx.data;
  yy = ts_yy.data;
  
  scatter(hca,xx,yy,5,tt)
  
  hca.XLabel.String = sprintf('B_%s (nT)',comps(2));
  hca.YLabel.String = sprintf('B_%s (nT)',comps(3));
  
  hcb = colorbar('peer',hca);
  hb(ihb) = hcb;
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 30*[-1 1];
  hca.XLim = 30*[-1 1];
end
if 1 % (Bx,Ez)
  hca = h2(isub); isub = isub + 1;
  
  ts_xx = eval(sprintf('%sB%g_srvy.tlim(tint_scatter).x;',cs,ic));
  ts_yy = eval(sprintf('%sE%g_fast.z;',cs,ic));
  ts_yy = ts_yy.resample(ts_xx);
  
  t0 = ts_xx.time(1);
  tt = ts_xx.time - ts_xx.time(1);
  xx = ts_xx.data;
  yy = ts_yy.data;
  
  scatter(hca,xx,yy,5,tt)
  
  hca.XLabel.String = sprintf('B_%s (nT)',comps(1));
  hca.YLabel.String = sprintf('E_%s (mV/m)',comps(3));
  
  hcb = colorbar('peer',hca);
  hb(ihb) = hcb;
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 30*[-1 1];
  hca.XLim = 30*[-1 1];
end
if 1 % (Bx,Ez+vixB)
  hca = h2(isub); isub = isub + 1;
  
  ts_xx = eval(sprintf('%sB%g_srvy.tlim(tint_scatter).x;',cs,ic));
  ts_yy = eval(sprintf('%sE%g_fast.resample(%sVi%g) + %sVi%g.cross(%sB%g_srvy.resample(%sVi%g))*1e3*1e-9*1e3;',cs,ic,cs,ic,cs,ic,cs,ic,cs,ic));
  ts_yy = ts_yy.resample(ts_xx).z;
  
  t0 = ts_xx.time(1);
  tt = ts_xx.time - ts_xx.time(1);
  xx = ts_xx.data;
  yy = ts_yy.data;
  
  scatter(hca,xx,yy,5,tt)
  
  hca.XLabel.String = sprintf('B_%s (nT)',comps(1));
  hca.YLabel.String = sprintf('(E+v_ixB)_%s (mV/m)',comps(3));
  
  hcb = colorbar('peer',hca);
  hb(ihb) = hcb;
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 30*[-1 1];
  hca.XLim = 30*[-1 1];
end
if 1 % (By,Ez)
  hca = h2(isub); isub = isub + 1;
  
  ts_xx = eval(sprintf('%sB%g_srvy.tlim(tint_scatter).y;',cs,ic));
  ts_yy = eval(sprintf('%sE%g_fast.z;',cs,ic));
  ts_yy = ts_yy.resample(ts_xx);
  
  t0 = ts_xx.time(1);
  tt = ts_xx.time - ts_xx.time(1);
  xx = ts_xx.data;
  yy = ts_yy.data;
  
  scatter(hca,xx,yy,5,tt)
  
  hca.XLabel.String = sprintf('B_%s (nT)',comps(2));
  hca.YLabel.String = sprintf('E_%s (mV/m)',comps(3));
  
  hcb = colorbar('peer',hca);
  hb(ihb) = hcb;
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 30*[-1 1];
  hca.XLim = 30*[-1 1];
end
if 1 % (Ey,Ez)
  hca = h2(isub); isub = isub + 1;
  
  ts_xx = eval(sprintf('%sE%g_fast.tlim(tint_scatter).y;',cs,ic));
  ts_yy = eval(sprintf('%sE%g_fast.z;',cs,ic));
  ts_yy = ts_yy.resample(ts_xx);
  
  t0 = ts_xx.time(1);
  tt = ts_xx.time - ts_xx.time(1);
  xx = ts_xx.data;
  yy = ts_yy.data;
  
  scatter(hca,xx,yy,5,tt)
  
  hca.XLabel.String = sprintf('E_%s (mV/m)',comps(2));
  hca.YLabel.String = sprintf('E_%s (mV/m)',comps(3));
  
  hcb = colorbar('peer',hca);
  hb(ihb) = hcb;
  
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.YLim = 30*[-1 1];
  hca.XLim = 30*[-1 1];
end


% Add markings to time panels
hmark = irf_pl_mark(h1(1),tint_scatter);
drawnow

% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];
colormap(hca,cmap_time)


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

%% B time-shift plots for fronts/islands

%% Figure for FPI
ic = 1;

npanels = 9;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

isub = 0;
zoomy = [];
link_clim = [];
link_clim2 = [];

if 1 % B gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gsmB?.x,gsmB?.y,gsmB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 
if 1 % i DEF omni
  isub = isub + 1;
  link_clim = [link_clim isub];
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};  
  irf_legend(hca,{sprintf('mms%g_dis_dist_brst',ic)},[0.02 0.1],'interpreter','none') 
end
if 1 % i psd x, no bg removed
  isub = isub + 1;
  link_clim2 = [link_clim2 isub];
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  irf_legend(hca,{'f_i(v_y)'},[0.02 0.1],'interpreter','tex')
  hca.YLabel.String = {'v_{iy}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF omni, error
  isub = isub + 1;
  link_clim = [link_clim isub];
  hca = irf_panel('i DEF omni error');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDistErr?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};   
  irf_legend(hca,{sprintf('mms%g_dis_disterr_brst',ic)},[0.02 0.1],'interpreter','none')
  
end
if 0 % i DEF omni, onecount
  isub = isub + 1;  
  hca = irf_panel('i DEF omni onecount');  
  
  specrec.p_label = 'log_{10} counts'; 
  c_eval('specrec.p = sum(iPDist?_onecount.data,[3 4],''omitnan'');',ic)
  c_eval('specrec.t = iPDist?_onecount.time.epochUnix;',ic)
  c_eval('specrec.f = iPDist?_onecount.depend{1};',ic)
  c_eval('specrec.f_label = ''E_i (eV)'';',ic)
      
  c_eval('[hout,hcb] = irf_spectrogram(hca,specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};   
  irf_legend(hca,{sprintf('(mms%g_dis_dist_brst/mms%g_dis_disterr_brst)^2',ic,ic),'summed over polar and azimuthal angles'}',[0.02 0.02],'color',[0 0 0],'interpreter','none')
  
end
if 0 % i DEF omni, bg removed
  isub = isub + 1;
  link_clim = [link_clim isub];
  hca = irf_panel('i DEF omn bg removed');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_nobg.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};    
  irf_legend(hca,{sprintf('mms%g_dis_dist_brst < 1.01 mms%g_dis_disterr_brst removed',ic,ic)},[0.02 0.1],'interpreter','none') 
end
if 0 % i psd x, bg removed
  isub = isub + 1;
  link_clim2 = [link_clim2 isub];
  hca = irf_panel('iLine x bg removed');
  c_eval('if1D = if1Dx?_nobg;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  irf_legend(hca,{'f_i(v_x)'},[0.02 0.1],'interpreter','tex')
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i DEF omni, elim
  isub = isub + 1;
  link_clim = [link_clim isub];
  hca = irf_panel('i DEF omni elim');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.elim([50 1e6]).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};    
  irf_legend(hca,{sprintf('E > %g eV',50)},[0.02 0.1],'interpreter','none') 
end
if 1 % i psd x, elim 50
  isub = isub + 1;
  link_clim2 = [link_clim2 isub];
  hca = irf_panel('iLine y bg elim 50');
  c_eval('if1D = if1Dy?_050;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  irf_legend(hca,{sprintf('f_i(v_y), E > %g eV',50)},[0.02 0.1],'interpreter','tex')
  hca.YLabel.String = {'v_{iy}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i DEF omni, elim
  isub = isub + 1;
  link_clim = [link_clim isub];
  hca = irf_panel('i DEF omni elim 100');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.elim([100 1e6]).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};    
  irf_legend(hca,{sprintf('E > %g eV',100)},[0.02 0.1],'interpreter','none') 
end
if 1 % i psd x, elim 50
  isub = isub + 1;
  link_clim2 = [link_clim2 isub];
  hca = irf_panel('iLine y bg elim 100');
  c_eval('if1D = if1Dy?_100;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  irf_legend(hca,{sprintf('f_i(v_y), E > %g eV',100)},[0.02 0.1],'interpreter','tex')
  hca.YLabel.String = {'v_{iy}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i DEF omni, elim
  isub = isub + 1;
  link_clim = [link_clim isub];
  hca = irf_panel('i DEF omni elim 200');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.elim([200 1e6]).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i','(eV)'};    
  irf_legend(hca,{sprintf('E > %g eV',200)},[0.02 0.1],'interpreter','none') 
end
if 1 % i psd x, elim 50
  isub = isub + 1;
  link_clim2 = [link_clim2 isub];
  hca = irf_panel('iLine y bg elim 200');
  c_eval('if1D = if1Dy?_200;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  irf_legend(hca,{sprintf('f_i(v_y), E > %g eV',200)},[0.02 0.1],'interpreter','tex')
  hca.YLabel.String = {'v_{iy}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

hlinks = linkprop(h(link_clim),{'CLim','YLim'});
hlinks2 = linkprop(h(link_clim2),{'CLim'});
irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Energy partition at different times
tint_tseries = tint;
%n_dists = 5;
%tint(1) + linspace(0,tint(2)-tint(1),n_dists);
tint_distributions = irf_time('2017-07-25T22:06:30.00Z','utc>EpochTT'):30:irf_time('2017-07-25T22:10:00.00Z','utc>EpochTT');
n_dists = tint_distributions.length;

cmap = irf_colormap('waterfall');
color_distributions = interp1(1:size(cmap,1),cmap,linspace(1,size(cmap,1),n_dists));


npanels = 6; % TSeries panels
nrows = 3; % Side panels
ncols = 1;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.7,'vertical'); % horizontal
isub = 1;
zoomy = [];

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
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Te/i par perp Ti/Tref
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('1234'))
  refTi = 1;
  c_eval('irf_plot(hca,{Te?par.tlim(tint),Te?perp,Ti?par.tlim(tint)/refTi,Ti?perp.tlim(tint)/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_{i,||}/' num2str(refTi,'%.0f')],['T_{i,\perp}/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
  hca.YLim = [0 1e4];
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % ve_par
    hold(hca,'on')
    c_eval('vexb = gseVe?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
end
if 1 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('specrec = eis_omni?.specrec(''energy'');',ic)
  % add a nan energy level to get the yscale in 10^x format
  %specrec.p = [specrec.p, nan(size(specrec.t))];
  %specrec.f = [specrec.f, 10e4*ones(size(specrec.t))];
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  %hca.YLim(2) = specrec.f(1,end-1);
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',ic)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex'; 
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('1234'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end
if 0 % Ti par perp
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
if 0 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  c_eval('if1D = if1Dx?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z');
  c_eval('if1D = if1Dz?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % iPitch
  isub = isub + 1;
  hca = irf_panel('iPitchj');
  c_eval('iPitch = ePitch?.elim([5000 10000]);',ic)
  irf_spectrogram(hca,iPitch.specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h1(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h1(ii).FontSize = 12;
end

irf_plot_axis_align
%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h1,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h1(zoomy),'y')
for idist = 1:n_dists
  for ipanel = 1:npanels
    hmark(ipanel,idist) = irf_pl_mark(h1(ipanel),tint_distributions(idist),color_distributions(idist,:))
  end
end

% Energy spectograms for different times
isub = 1;
if 1 % ion f
  hca = h2(isub); isub = isub + 1;
  c_eval('dist = iPDist?.omni;',ic)
  idist_all = [];
  for id = 1:n_dists
    time = tint_distributions(id);    
    [~,idist] = min(abs(dist.time-time));
    idist_all(id) = idist;
  end
  dist_tmp = dist(idist_all);    
  toplot = dist_tmp;
  set(hca,'ColorOrder',color_distributions)
  loglog(hca,dist_tmp.depend{1}(1,:),toplot.data,'-');   
  set(hca,'ColorOrder',color_distributions)
  hca.XLabel.String = 'E_i (eV)';
  hca.YLabel.String = {'Ion phase space density','(s^3/cm^6)'};  
end
if 1 % ion deflux
  hca = h2(isub); isub = isub + 1;
  c_eval('dist = iPDist?.deflux.omni;',ic)
  idist_all = [];
  for id = 1:n_dists
    time = tint_distributions(id);    
    [~,idist] = min(abs(dist.time-time));
    idist_all(id) = idist;
  end
  dist_tmp = dist(idist_all);    
  toplot = dist_tmp;
  set(hca,'ColorOrder',color_distributions)
  loglog(hca,dist_tmp.depend{1}(1,:),toplot.data,'-');   
  set(hca,'ColorOrder',color_distributions)
  hca.XLabel.String = 'E_i (eV)';
  hca.YLabel.String = {'Ion differential energy flux','(keV/(cm^2 s sr keV))'};  
end

%% Long time overview
cmap = irf_colormap('waterfall');
doFilt = 1;
ffilt = 1;

npanels = 8;
h = irf_plot(npanels);
isub = 1;
zoomy = [];

if 1 % B gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?_srvy.x,gseB?_srvy.y,gseB?_srvy.z,gseB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 0 % B lmn fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?_srvy.x,lmnB?_srvy.y,lmnB?_srvy.z,lmnB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{LMN}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'  [L; M; N] = ',...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(1,1),R(1,2),R(1,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(2,1),R(2,2),R(2,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(3,1),R(3,2),R(3,3))}',[1.002 0.9],'fontsize',12,'color','k');  
end
if 0 % E gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  if doFilt
    c_eval('irf_plot(hca,{gseE?_fast.filt(0,ffilt,[],5).x,gseE?_fast.filt(0,ffilt,[],5).y,gseE?_fast.filt(0,ffilt,[],5).z},''comp'');',ic)
    irf_legend(hca,sprintf('f < %g Hz',ffilt),[0.02 0.07])
  else
    c_eval('irf_plot(hca,{gseE?_fast.x,gseE?_fast.y,gseE?_fast.z},''comp'');',ic)
  end
    
  hca.YLabel.String = {'E_{GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 1 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?_fast.x,gseVi?_fast.y,gseVi?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % V oplus gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V oplus gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{vOp?_srvy.x,vOp?_srvy.y,vOp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{O+}^{GSM}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % V hplus gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V hplus gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{vHp?_srvy.x,vHp?_srvy.y,vHp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{H+}^{GSM}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % VExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % E + VixB gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E + VixB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))    
  c_eval('irf_plot(hca,{gseEVixB?_fast.x,gseEVixB?_fast.y,gseEVixB?_fast.z},''comp'');',ic)    
  hca.YLabel.String = {'E_{GSE} + V_ixB','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % n log scale
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 0 % n log scale, no burst
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{fast}','n_{O+}','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 0 % nO/nHp
  hca = irf_panel('nOp/nHp log scale');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{nOp1_srvy/nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{fast}','n_{O+}','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n_{O^+}/n_{H^+}'};
  %hca.YScale = 'log';
  %hca.YLim(1) = 1e-5;
end

elim_feeps = [8e04 Inf]; % [8e04 Inf]
if 1 % FEEPS Pitch 0 allE
%  feeps_pa2.elim([9.7+04 Inf])
  isub = isub + 1;
  hca = irf_panel('feeps Pitch all E');  
  %elim = [0 32000];
  c_eval('irf_spectrogram(hca,feeps_pa?_srvy.elim(elim_feeps).specrec(''pitchangle''));',ic)
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 1 % i PEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?_srvy.elim(elim_feeps).specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval('irf_plot(hca,{betai?_srvy_hpca,beta?_srvy_fpi},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
end

if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_fast.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % hplus DEF omni
  isub = isub + 1;
  hca = irf_panel('hplus DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_Hp_fast.dpflux(-1).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*1*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');    
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end  
  hca.YLabel.String = {'E_{H+}','(eV)'};     
  hca.YLabel.Interpreter = 'tex';
end
if 1 % oplus DEF omni
  isub = isub + 1;
  hca = irf_panel('oplus DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_Op_fast.dpflux(-1).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*16*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');    
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end  
  hca.YLabel.String = {'E_{O+}','(eV)'};     
  hca.YLabel.Interpreter = 'tex';
end

if 0 % e PEF feeps omni
  isub = isub + 1;
  hca = irf_panel('e DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ele_omni?_srvy.elim(elim_feeps).specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_e^{FEEPS}','(eV)'};   
end

irf_zoom(h,'x',tint_fast)
irf_plot_axis_align

%% Calculate Lp from p and vy
%%%%%%%
%%%%%%% Coordinate transform in rec_stop_start.lmn_trans
%%%%%%%

ic = 1;
units = irf_units;
q = units.e;
% Lp = p/(q*n*B*vy_perp)
pp = Pi1perp.data*1e-9;               % Pa
vv = lmnVi1perp.y.data*1e3;           % m/s
nn = ne1.resample(ni1).data*1e6;      % m^-3
bb = lmnB1.x.resample(ni1).data*1e-9; % T 

Lp1 = irf.ts_scalar(gseVi1.time,pp./(q*vv.*nn.*bb))*1e-3; % km
Lp1.data(Lp1.data > 10000,:) = NaN;
Lp1.data(Lp1.data < 0,:) = NaN;

LpB1 = irf.ts_scalar(gseVi1.time,Lp1.data.*lmnB1.x.resample(ni1).data); % km*nT

% Plot
npanels = 10;
h = irf_plot(npanels);
doFilt = 0;
ffilt = 1;

if 1 % B lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?.x,lmnB?.y,lmnB?.z,lmnB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B^{LMN}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'  [L; M; N] = ',...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(1,1),R(1,2),R(1,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(2,1),R(2,2),R(2,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]  (GSE)',R(3,1),R(3,2),R(3,3))},[0.02 1.05],'fontsize',12,'color','k');  
end
if 1 % E lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  if doFilt
    c_eval('irf_plot(hca,{lmnE?.filt(0,ffilt,[],5).x,lmnE?.filt(0,ffilt,[],5).y,lmnE?.filt(0,ffilt,[],5).z},''comp'');',ic)    
    irf_legend(hca,sprintf('f < %g Hz',ffilt),[0.02 0.07])
  else
    c_eval('irf_plot(hca,{lmnE?.x,lmnE?.y,lmnE?.z},''comp'');',ic)    
  end
  hca.YLabel.String = {'E^{LMN}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end 
if 1 % Vi lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVi?.x,lmnVi?.y,lmnVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVi?perp.x,lmnVi?perp.y,lmnVi?perp.z,lmnVi?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_{\perp,L}','v_{\perp,M}','v_{\perp,N}','v_{||}'},[0.98 0.9],'fontsize',12);
end
if 1 % VExB lmn
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{lmnVExB?.x,lmnVExB?.y,lmnVExB?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}^{LMN}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end

if 1 % Pi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pi');
  set(hca,'ColorOrder',mms_colors('1234'))
  
  c_eval('irf_plot(hca,{gsePi?.trace/3,Pi?perp,Pi?par},''comp'');',ic)
  hca.YLabel.String = {'P_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'P_{i,||}','P_{i,\perp}','P_{i,||}'}',[1.01 0.90],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
  %hca.YLim = [0 1e4];
end
if 0 % n log scale
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('1234'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.8])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_fast.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
    irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % Lp
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{Lp?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L_p','(km)'};
  hca.YLim = [0 10000];
end
if 1 % Lp*Bx
  hca = irf_panel('Lp*Bx');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{LpB?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L_pB_x','(km nT)'};
  %hca.YLim = [0 10000];
end
if 1 % Lp*Bx/(B0/2)
  hca = irf_panel('Lp*Bx/(B0)');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{LpB?/(B0_*1e9)},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L = L_pB_x/(B_0)','(km)'};
  %hca.YLim = [0 10000];
end

irf_zoom(h,'x',tint)
irf_plot_axis_align

%% Scatter plots of 

%% Estimate beam with of initial field-aligned flow
tint_beam = irf.tint('2017-07-25T22:05:59.00Z/2017-07-25T22:06:00.00Z');

vint = [-Inf Inf];
%c_eval('if1Dx? = iPDist?.reduce(''1D'',[1 0 0],''vint'',vint);',ic)
%vdf = if1Dx1.tlim(tint(1)+0.2*[-1 1]);
T = 0.150;
c_eval('vdf = iPDist?_nobg.tlim(tint_beam(1)+0.5*T*[-1 1]).reduce(''1D'',[1 0 0],''vint'',vint);',ic)

[data,timeseries] = funFitVDF(vdf,'plot',1);

%% Figure: Part A, reconnection cessation 
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 10;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 1 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-3:1:2];
  hca.YLabel.Interpreter = 'tex';
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end

elim_feeps = [8e4 Inf];
if 1 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 1 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,feeps_ion_omni1234.elim(elim_feeps).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 1 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  eval(sprintf('irf_plot(hca,{%sJcurl.x,%sJcurl.y,%sJcurl.z},''comp'');',cs,cs,cs))
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
if 0 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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

if 0 % EIS Pitch
  isub = isub + 1;
  hca = irf_panel('eis Pitch');  
  irf_spectrogram(hca,eis_pa2.elim([14 4500]*1e3).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('xyza'))  
  hca.YLim = [0 5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);


hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [0 4]*0.99;

% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
else
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

drawnow
h(6).YLim(1) = h(7).YLim(2);
h(2).YLim(1) = 0.001;

drawnow
hmark = irf_pl_mark(h,time_mark,'k');
c_eval('hmark(?).LineStyle = '':'';',1:5)
c_eval('hmark(?).LineWidth = 1;',1:5)

%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);




%% Scatter plots of Vi, Ve, VExB

h = setup_subplots(1,2);
isub = 1;

ic = 1;
cs = 'lmn';
vlim_min = 5000*0.99;

%for comp = ['x','y','z']
for comp = ['y']
  for isp = ['i','e']
    if 1
      hca = h(isub); isub = isub + 1;
      %comp = 'x';
      %isp = 'i';
      yy = eval(sprintf('%sV%s%g.(comp)',cs,isp,ic));
      xx = eval(sprintf('%sVExB%g.(comp).resample(yy)',cs,ic));
      cc = eval(sprintf('%sV%s%g.time - %sV%s%g.time(1)',cs,isp,ic,cs,isp,ic));
      scatter(hca,xx.data,yy.data,1,cc)
      hcb = colorbar('peer',hca);
      hca.XGrid = 'on';
      hca.YGrid = 'on';
      hca.XLabel.String= 'v_{ExB} (km/s)';
      hca.YLabel.String= sprintf('v_{%s%s} (km/s)',isp,comp);
      colormap(hca,irf_colormap('waterfall'))
      %axis(hca,'equal')
      if 0 && max(abs(hca.YLim)) > vlim_min
%         hca.XLim = vlim_min*[-1 1];
%         hca.YLim = vlim_min*[-1 1];
      end
      if 1 % set limits as percentiles
        hca.XLim = prctile(xx.data,[0.5 99.5]);
        hca.YLim = prctile(yy.data,[0.5 99.5]);
      end
      axis(hca,'equal')
      axis(hca,'square')
    end
  end
end


for comp = ['x','y','z']
  for isp = ['i','e']
    if 0
      hca = h(isub); isub = isub + 1;
      %comp = 'x';
      %isp = 'e';
      yy = eval(sprintf('%sV%s%g.(comp)',cs,isp,ic));
      xx = eval(sprintf('%sVExB%g.(comp).resample(yy)',cs,ic));
      cc = eval(sprintf('%sV%s%g.time - %sV%s%g.time(1)',cs,isp,ic,cs,isp,ic));
      scatter(hca,xx.data,yy.data,1,cc)
      hcb = colorbar('peer',hca);
    end
  end
end

%% Scatter plots of Te/Ti/Bx/n

h = setup_subplots(2,2);
isub = 1;

ic = 1;
cs = 'lmn';
vlim_min = 5000*0.99;

if 1 % (Ti,Te)
  hca = h(isub); isub = isub + 1;
%   yy = eval(sprintf('%sV%s%g.(comp)',cs,isp,ic));
%   xx = eval(sprintf('%sVExB%g.(comp).resample(yy)',cs,ic));
%   cc = eval(sprintf('%sV%s%g.time - %sV%s%g.time(1)',cs,isp,ic,cs,isp,ic));
  c_eval('xx = gseTi?.trace/3;',ic)
  c_eval('yy = gseTe?.resample(gseTi?).trace/3;',ic)
  c_eval('cc = gseTi?.time - gseTi?.time(1);',ic);
  scatter(hca,xx.data,yy.data,7,cc)
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('pasteljet'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'T_i (eV)';
  hca.YLabel.String = 'T_e (eV)';
  hca.YLim(1) = 0;
  hca.XLim(1) = 0;
  hca.YLim(2) = 1e4;
end
if 1 % (Pi,Pe)
  hca = h(isub); isub = isub + 1;
%   yy = eval(sprintf('%sV%s%g.(comp)',cs,isp,ic));
%   xx = eval(sprintf('%sVExB%g.(comp).resample(yy)',cs,ic));
%   cc = eval(sprintf('%sV%s%g.time - %sV%s%g.time(1)',cs,isp,ic,cs,isp,ic));
  c_eval('xx = gsePi?.trace/3;',ic)
  c_eval('yy = gsePe?.resample(gsePi?).trace/3;',ic)
  c_eval('cc = gsePi?.time - gsePi?.time(1);',ic);
  scatter(hca,xx.data,yy.data,7,cc)
  hcb = colorbar('peer',hca);
  colormap(hca,pic_colors('thermal'))
  colormap(hca,pic_colors('candy4'))
  colormap(hca,irf_colormap('waterfall'))
  colormap(hca,'parula')
  colormap(hca,pic_colors('pasteljet'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XLabel.String = 'P_i (nPa)';
  hca.YLabel.String = 'P_e (nPa)';
  hca.YLim(1) = 0;
  hca.XLim(1) = 0;
  %hca.Color = [0 0 0];
end
if 1
  hca = h(isub); isub = isub + 1;
%   yy = eval(sprintf('%sV%s%g.(comp)',cs,isp,ic));
%   xx = eval(sprintf('%sVExB%g.(comp).resample(yy)',cs,ic));
%   cc = eval(sprintf('%sV%s%g.time - %sV%s%g.time(1)',cs,isp,ic,cs,isp,ic));
  c_eval('xx1 = gseB?.resample(gsePe?).x;',ic)  
  c_eval('yy1 = gsePe?.resample(gsePe?).trace/3;',ic)
  c_eval('xx2 = gseB?.resample(gsePi?).x;',ic)
  c_eval('yy2 = gsePi?.resample(gsePi?).trace/3;',ic)
  c_eval('cc = gsePi?.time - gsePi?.time(1);',ic);
  hs1 = scatter(hca,xx1.data,yy1.data,1,'k');
  hold(hca,'on')
  hs2 = scatter(hca,xx2.data,yy2.data,1,'r');
  hold(hca,'off')
  hcb = colorbar('peer',hca);
  %colormap(hca,pic_colors('pasteljet'))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %legend(hca,{'e','i'},'box','off')
  irf_legend(hca,{'e','i'},[0.98 0.98])
end


%% Long time overview, MVA
cmap = irf_colormap('waterfall');
cmap = pic_colors('pasteljet');
doFilt = 1;
ffilt = 1;
ic = 1;

npanels = 6;
h = irf_plot(npanels);
isub = 1;
zoomy = [];

if 1 % B gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?_srvy.x,gseB?_srvy.y,gseB?_srvy.z,gseB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 0 % B lmn fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{lmnB?_srvy.x,lmnB?_srvy.y,lmnB?_srvy.z,lmnB?_srvy.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{LMN}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'  [L; M; N] = ',...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(1,1),R(1,2),R(1,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(2,1),R(2,2),R(2,3)),...
                  sprintf('[%5.2f, %5.2f, %5.2f]',R(3,1),R(3,2),R(3,3))}',[1.002 0.9],'fontsize',12,'color','k');  
end
if 0 % E gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  if doFilt
    c_eval('irf_plot(hca,{gseE?_fast.filt(0,ffilt,[],5).x,gseE?_fast.filt(0,ffilt,[],5).y,gseE?_fast.filt(0,ffilt,[],5).z},''comp'');',ic)
    irf_legend(hca,sprintf('f < %g Hz',ffilt),[0.02 0.07])
  else
    c_eval('irf_plot(hca,{gseE?_fast.x,gseE?_fast.y,gseE?_fast.z},''comp'');',ic)
  end
    
  hca.YLabel.String = {'E_{GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % Vi gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?_fast.x,gseVi?_fast.y,gseVi?_fast.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i^{GSE}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % V oplus gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V oplus gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{vOp?_srvy.x,vOp?_srvy.y,vOp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{O+}^{GSM}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % V hplus gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V hplus gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{vHp?_srvy.x,vHp?_srvy.y,vHp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{H+}^{GSM}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % VExB gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?_srvy.x,gseVExB?_srvy.y,gseVExB?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % E + VixB gse fast
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E + VixB gse');
  set(hca,'ColorOrder',mms_colors('xyza'))    
  c_eval('irf_plot(hca,{gseEVixB?_fast.x,gseEVixB?_fast.y,gseEVixB?_fast.z},''comp'');',ic)    
  hca.YLabel.String = {'E_{GSE} + V_ixB','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end 
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy*16,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}*16','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % n log scale
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?,ne?_fast,nOp1_srvy*16,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}*16','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 0 % n log scale, no burst
  hca = irf_panel('n log scale');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{fast}','n_{O+}','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YScale = 'log';
  hca.YLim(1) = 1e-5;
end
if 1 % n lin scale, no burst
  hca = irf_panel('n lin scale o burst');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{ne?_fast,nOp1_srvy*16,nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  irf_legend(hca,{'n_e^{fast}','n_{O+} (m_{O+}/m_p)','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YLim(1) = 0;
end
if 0 % nO/nHp
  hca = irf_panel('nOp/nHp log scale');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{nOp1_srvy/nHp1_srvy},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{fast}','n_{O+}','n_{H+}'}',[1.02 0.98])
  hca.YLabel.String = {'n_{O^+}/n_{H^+}'};
  %hca.YScale = 'log';
  %hca.YLim(1) = 1e-5;
end
if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_fast.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XMinorGrid = 'off';
  hca.YMinorGrid = 'off';
  hca.Layer = 'top';
end
if 1 % hplus DEF omni
  isub = isub + 1;
  hca = irf_panel('hplus DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_Hp_fast.dpflux(-1).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*1*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');    
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end  
  hca.YLabel.String = {'E_{H+}','(eV)'};     
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XMinorGrid = 'off';
  hca.YMinorGrid = 'off';
  hca.Layer = 'top';
end
if 1 % oplus DEF omni
  isub = isub + 1;
  hca = irf_panel('oplus DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?_Op_fast.dpflux(-1).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*16*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');    
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end  
  hca.YLabel.String = {'E_{O+}','(eV)'};     
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XMinorGrid = 'off';
  hca.YMinorGrid = 'off';
  hca.Layer = 'top';
end
if 0 % oplus DEF omni
  isub = isub + 1;
  hca = irf_panel('oplus DEF omni EIS');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?_oplus.dpflux(-1).deflux.specrec,''log'');',2)  
  set(hca,'yscale','log');  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*16*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');    
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end  
  hca.YLabel.String = {'E_{O+}/{EIS}','(eV)'};     
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.XMinorGrid = 'off';
  hca.YMinorGrid = 'off';
  hca.Layer = 'top';
end

if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval('irf_plot(hca,{betai?_srvy_hpca,beta?_srvy_fpi},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
end

irf_zoom(h,'x',tint_fast)
irf_plot_axis_align

ipref = npanels-2;
if 1 % link spectrograms
  %%
  linkprop(h(ipref+(0:2)),{'CLim'})
  h(ipref).CLim = [3.5 5.5];
  h(ipref).YLim = h(ipref+1).YLim;
end
%% Figure: Current sheet structure, energetic plasma
ic = 1;

units = irf_units;
q = units.e;
% Lp = p/(q*n*B*vy_perp)
c_eval('pp = Pi?perp.data*1e-9;',ic)               % Pa
c_eval('vv = lmnVi?perp.y.data*1e3;',ic)           % m/s
c_eval('nn = ne?.resample(ni?).data*1e6;',ic)      % m^-3
c_eval('bb = lmnB?.x.resample(ni?).data*1e-9;',ic) % T 

c_eval('Lp? = irf.ts_scalar(gseVi?.time,pp./(q*vv.*nn.*bb))*1e-3;',ic) % km
c_eval('Lp?.data(Lp?.data > 10000,:) = NaN;',ic)
c_eval('Lp?.data(Lp?.data < 0,:) = NaN;',ic)

c_eval('LpB? = irf.ts_scalar(gseVi?.time,Lp?.data.*lmnB?.x.resample(ni?).data);',ic) % km*nT


time_mark = EpochTT('2017-07-25T22:09:08.00Z');

Etop_fpi = iPDist1.ancillary.energy(1,end)+iPDist1.ancillary.delta_energy_plus(1,end);

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
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
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 1 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
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
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
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

if 0 % e DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('e DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ele_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_e^{FEEPS}','(eV)'};   
end
if 0 % EIS Pitch 1 lowE
  isub = isub + 1;
  hca = irf_panel('eis Pitch 1');  
  elim = [0 32000];
  irf_spectrogram(hca,eis_pa2.elim(elim).specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 0 % EIS Pitch 2 highE
  isub = isub + 1;
  hca = irf_panel('eis Pitch 2');  
  elim = [32000 100000];
  irf_spectrogram(hca,eis_pa2.elim(elim).specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end

if 0 % EIS Pitch 0 allE oxygen
  isub = isub + 1;
  hca = irf_panel('eis Pitch ox');  
  %elim = [0 32000];
  irf_spectrogram(hca,eis_pa2_op.specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
elim_eis = [32000 Inf];
if 0 % EIS Pitch 0 allE
  isub = isub + 1;
  hca = irf_panel('eis Pitch all E');  
  
  irf_spectrogram(hca,eis_pa3.elim(elim_eis).specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.elim(elim_eis).specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
elim_feeps = [8e4 Inf];
if 0 % FEEPS Pitch 0 allE
%  feeps_pa2.elim([9.7+04 Inf])
  isub = isub + 1;
  hca = irf_panel('feeps Pitch all E');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_pa2.elim([8e04 Inf]).specrec('pitchangle'));  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.elim([8e04 Inf]).specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 0 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,feeps_ion_omni1234.elim(elim_feeps).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % FEEPS ele, single sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps ele pitch mms ');  
  elim_feeps = [0 Inf];
  c_eval('irf_spectrogram(hca,feeps_ele_omni?.elim(elim_feeps).specrec(''energy''));',ic)  
  hca.YLabel.String = {'E^{FEEPS}_e','eV'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])  
end
if 0 % FEEPS ele, single sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps ele pitch mms lines');  
  elim_feeps = [0 Inf];
  c_eval('irf_plot(hca,feeps_ele_omni?.elim(elim_feeps));',ic)  
  hca.YLabel.String = {'E^{FEEPS}_e','eV'};  
  hca.YLabel.Interpreter = 'tex';
  hca.YScale= 'log';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])  
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
if 0 % Ti par perp
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

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 1 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end
if 1 % Lp
  hca = irf_panel('Lp');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{Lp?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L_p','(km)'};
  hca.YLim = [0 10000];
end
if 1 % Lp*Bx
  hca = irf_panel('Lp*Bx');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{LpB?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L_pB_x','(km nT)'};
  %hca.YLim = [0 10000];
end
if 1 % Lp*Bx/(B0/2)
  hca = irf_panel('Lp*Bx/(B0)');
  set(hca,'ColorOrder',mms_colors('matlab'))
  %c_eval('irf_plot(hca,{ne?,ni?,ne?_fast,ni?_fast,nOp1_srvy,nHp1_srvy},''comp'')',ic)
  c_eval('irf_plot(hca,{LpB?/(B0_*1e9)},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))
  %irf_legend(hca,{'n_e^{brst}','n_i^{brst}','n_e^{fast}','n_i^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  %irf_legend(hca,{'n_e^{brst}','n_e^{fast}','n_{O+}','n_{H+}'},[0.98 0.98])
  hca.YLabel.String = {'L = L_pB_x/(B_0)','(km)'};
  %hca.YLim = [0 10000];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
h(1).Title.String = irf_ssub('MMS ?',ic);


if 0
drawnow
h(4).YLim(1) = h(5).YLim(2);
h(2).YLim(1) = 0.001;

drawnow
hmark = irf_pl_mark(h,time_mark,'k');
c_eval('hmark(?).LineStyle = '':'';',1:5)
c_eval('hmark(?).LineWidth = 1;',1:5)

hat=annotation('textarrow',[0.593 0.593],[.3 .32],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
end

%% Figure: Reconnection onset, density cavity
ic = 1;

npanels = 7;
h = irf_plot(npanels);  
iisub = 0;
cmap = colormap(pic_colors('candy4'));
timeline = tint(1):0.5:tint(2);
doResample = 0;

isub = 0;
zoomy = [];

cs = 'lmn';
cs_str = 'LMN';
if 1 % B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z,%sB?.abs},''comp'');',cs,cs,cs,cs),ic)    
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3),'|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, 4sc average
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_plot(hca,{PBav,Piav,Peav,PDiav,Piav + Peav.resample(Piav) + PBav.resample(Piav)},'comp')
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i}','P_{e}','P_{i,dyn,x}','P_{B}+P_i+P_e'},[0.98 0.9],'fontsize',12);  
end
if 0 % Pressures, single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3 + gsePe?.resample(gsePi?).trace/3,'...
                         'PDi?,',...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_{i+e}','P_{dyn,i}','P_i','P_e','P_{tot}'},[0.98 0.9],'fontsize',12);  
end
if 0 % Je single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Je');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sJe?.x,%sJe?.y,%sJe?.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'J_e','()'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 0 % Ji single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ji');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sJi?.x,%sJi?.y,%sJi?.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'J_i','()'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 0 % J single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sJ?.x,%sJ?.y,%sJ?.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  if doResample
    irf_plot(hca,{gseJcurl.x.resample(timeline),gseJcurl.y.resample(timeline),gseJcurl.z.resample(timeline)},'comp');  
  else
    irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');  
  end
  hca.YLabel.String = {'J^{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);    
end
if 0 % JxB single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JxB');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sJxB?.x,%sJxB?.y,%sJxB?.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'JxB','()'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 0 % JxB/ne single sc
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JxB/ne');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sJxBne?_mVm.x,%sJxBne?_mVm.y,%sJxBne?_mVm.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'JxB/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 0 % JxB, gradB, curvB, 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    scale = 1e0*1e0; units = 'T Am^{-2}';
    scale = 1e9*1e9; units = 'nT nAm^{-2}';
    
    hca = irf_panel(['JxB ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))   
    irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale},'comp');  
    %irf_plot(hca,{gseJxB.(comp)*scale,-gseDivPb.(comp)*scale,-gseCurvB.(comp)*scale,-gseDivPb.(comp)*scale-gseCurvB.(comp)*scale},'comp');  
    hca.YLabel.String = {['JxB_' comp],['(' units ')']};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  end
end
if 0 % JxB, 1 panels: x,y,z  
  isub = isub + 1;
  zoomy = [zoomy isub];      
  hca = irf_panel('JxB xyz');
  scale = 1e9*1e9; units = 'nT nAm^{-2}';
  set(hca,'ColorOrder',pic_colors('matlab'))   
  irf_plot(hca,{gseJxB.x*scale,gseJxB.y*scale,gseJxB.z*scale},'comp');      
  hca.YLabel.String = {['JxB'],['(' units ')']};
  set(hca,'ColorOrder',pic_colors('matlab'))
  %irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  irf_legend(hca,{'x','y','z'}',[1.01 0.9],'fontsize',12);    
end
if 0 % JxB curl, 1 panels: LMN 
  isub = isub + 1;
  zoomy = [zoomy isub];      
  hca = irf_panel('JxB lmn');
  scale = 1e9*1e9; units = 'nT nAm^{-2}';
  set(hca,'ColorOrder',pic_colors('xyza'))   
  eval(sprintf('irf_plot(hca,{%sJxB.x*scale,%sJxB.y*scale,%sJxB.z*scale},''comp'');',cs,cs,cs))
  hca.YLabel.String = {['JxB'],['(' units ')']};
  set(hca,'ColorOrder',pic_colors('xyza'))
  %irf_legend(hca,{'JxB','-\nabla B^2/2\mu_0','-\mu_0^{-1}B\cdot\nabla B'}',[1.01 0.9],'fontsize',12);  
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)}',[1.01 0.9],'fontsize',12);    
end
if 0 % electron mom eq. 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',mms_colors('xyza'))
    to_plot_gseVexBav = gseVexBav;
    to_plot_gseVexBav.data(abs(to_plot_gseVexBav.data)>100) = NaN;
    irf_plot(hca,{gseEav.(comp),-1*to_plot_gseVexBav.(comp),-1*gseGradPene.resample(gseVi1).(comp),-1*(gseEav.(comp).resample(to_plot_gseVexBav)+to_plot_gseVexBav.(comp))},'comp');  
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'E','-v_exB','-divPe/ne','-(E+v_exB)'}',[1.02 0.9],'fontsize',12);  
  end
end
if 1 % ohm's law, 3 panels: x,y,z
  for comp = ['x','y','z']
    isub = isub + 1;
    zoomy = [zoomy isub];
    %comp = 'x';
    hca = irf_panel(['ele eom ', comp]);
    set(hca,'ColorOrder',pic_colors('matlab'))
    if doResample
      to_plot_gseVixBav = gseVixBav.resample(timeline);
    else
      to_plot_gseVixBav = gseVixBav;
    end
    to_plot_gseVixBav.data(abs(to_plot_gseVixBav.data)>100) = NaN;
    if doResample
      irf_plot(hca,{...
        gseEav.(comp).resample(timeline),...
        to_plot_gseVixBav.(comp).resample(timeline),...
        -1*gseGradPene.(comp).resample(timeline),...
        gseJxBne_mVm.(comp).resample(timeline),...
        +1*(gseEav.(comp).resample(timeline)+to_plot_gseVixBav.(comp).resample(timeline))},'comp');  
    else
      irf_plot(hca,{gseEav.(comp),to_plot_gseVixBav.(comp),-1*gseGradPene.(comp).resample(gseVi1),gseJxBne_mVm.(comp),+1*(gseEav.(comp).resample(to_plot_gseVixBav)+to_plot_gseVixBav.(comp))},'comp');  
    end
    hca.YLabel.String = {['E_' comp],'(mV/m)'};
    set(hca,'ColorOrder',pic_colors('matlab'))
    irf_legend(hca,{'E','v_ixB','-divPe/ne','JxB/ne','(E+v_ixB)'}',[1.01 0.9],'fontsize',12);  
  end
end
if 0 % vExB gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gsmVExB?.x,gsmVExB?.y,gsmVExB?.z},''comp'');',ic)  
  c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVExB?.resample(gsmVi?).y,gsmVExB?.resample(gsmVi?).z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x,gsmVi?perp.x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,x}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.y gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).y,gsmVi?perp.y},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,y}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.z gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).z,gsmVi?perp.z},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{\perp,z}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
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
if 0 % E gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmE?.x,gsmE?.y,gsmE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Vi perp,par
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
if 0 % Ve gse
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve gsm
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x.tlim(tint),gsmVe?.y.tlim(tint),gsmVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve filt
  hca = irf_panel('Ve gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 3;
  c_eval('irf_plot(hca,{gsmVe?.filt(0,3,[],3).x.tlim(tint),gsmVe?.filt(0,3,[],3).y.tlim(tint),gsmVe?.filt(0,3,[],3).z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e (GSM)','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve  gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVe?.x,gsmVe?.y,gsmVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVe?perp.x,%sVe?perp.y,%sVe?perp.z},''comp'');',cs,cs,cs),ic) 
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVe?par.tlim(tint)},''comp'');',cs),ic)
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
if 0 % Ve perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x,gseVe?perp.y,gseVe?perp.z,gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {['v_{e} mms' num2str(ic)],'(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve perp par av
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{gseVeperpav.x,gseVeperpav.y,gseVeperpav.z,gseVeparav},'comp');
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp x','\perp y','\perp z','||'},[0.98 0.9],'fontsize',12);
end

if 0 % vExB.x , Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi x');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).x,%sVi?perp.x},''comp'');',cs,cs,cs),ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {sprintf('v_{\\perp,z}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x , Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi y');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).y,%sVi?perp.y},''comp'');',cs,cs,cs),ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {sprintf('v_{\\perp,y}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x , Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB Vi z');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).z,%sVi?perp.z},''comp'');',cs,cs,cs),ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {sprintf('v_{\\perp,z}^{%s}',cs),'(km/s)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'ExB','V_i'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ve',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVe?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('irf_plot(hca,{%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x gsm, Vi resample
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVExB?.resample(gsmVi?).x},''comp'');',ic)  
  %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end


if 1 % Te par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('1234'))
  refTi = 5;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
end

if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % ePDist pa
  isub = isub + 1;
  hca = irf_panel('e PA deflux lowe');  
  eint = [100 1000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 
  isub = isub + 1;
  hca = irf_panel('e PA deflux highe');  
  eint = [100 40000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

h(1).Title.String = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));

irf_zoom(h,'x',tint_cavity)
drawnow
irf_zoom(h(zoomy),'y')
irf_zoom(h,'y')
irf_plot_axis_align
%tmark_eis = tint(1):20:tint(2);

%% Plot 2D electron distributions
ic = 1;
%[h,h2] = initialize_combined_plot(8,3,2,0.5,'vertical');
h = irf_plot(8);

%tint_figure = tint_cavity + [-135 15];
tint_figure = tint;

times_dist = EpochTT(['2017-07-25T22:09:30.00Z';...
                      '2017-07-25T22:09:45.00Z';...
                      '2017-07-25T22:09:50.00Z';...
                      '2017-07-25T22:09:52.20Z';...
                      '2017-07-25T22:09:52.60Z';...
                      '2017-07-25T22:09:53.00Z']);

times_dist = EpochTT(['2017-07-25T22:09:30.00Z';...
                      '2017-07-25T22:09:45.00Z';...
                      '2017-07-25T22:09:50.00Z';...
                      '2017-07-25T22:09:52.10Z';...
                      '2017-07-25T22:09:52.80Z';...
                      '2017-07-25T22:09:53.00Z']);

cs = 'lmn';
cs_str = 'LMN';

cs = 'gse';
cs_str = 'xyz';


comps = ['x','y','z'];
comp_str = 'LMN';
comp_str = 'xyz';

comps = ['x'];
comps_str = 'x';
%comp_str = 'L';

if 0
cs = 'gse';
cs_str = 'xyz';

comps = ['x','y','z'];
comp_str = 'xyz';
end

zoomy = []; isub = 1;
if 1 % B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z,%sB?.abs},''comp'');',cs,cs,cs,cs),ic)    
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3),'|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % E
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3)},[0.98 0.9],'fontsize',12);  
end
if 1 % Epar
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Epar');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{gseE?par},''comp'');',cs,cs,cs),ic)    
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 1 % vExB.x/y/z , Vi resample, ve
  colors_xyz = [0.6 0.6 0.6; 0 0 0; mms_colors('2')];
  for icomp = 1:numel(comps) % vExB.x , Vi resample, ve  
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ve',comp]);
    set(hca,'ColorOrder',colors_xyz)
    c_eval(sprintf('ve = %sVe?perp.(comp);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVe?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}',comp_str(icomp)),'(km/s)'};
    set(hca,'ColorOrder',colors_xyz)
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % Vi.x/y/z
  for comp = comps % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['Vi',comp]);
    set(hca,'ColorOrder',mms_colors('2'))    
    c_eval(sprintf('irf_plot(hca,{gseVi?.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{(%s)}',comp_str),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'v_i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end

if 1 % vepar
  isub = isub + 1;
  hca = irf_panel('vepar');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{gseVe?par},''comp'')',ic)
  hold(hca,'on')
  c_eval('irf_patch(hca,{gseVe?par,0})',ic)
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('12'))  
  hca.YLabel.String = {'v_{e||}','(km/s)'};
  hca.YLabel.Interpreter = 'tex';
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.tlim(tint_figure).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 1 % ePDist pa 
  isub = isub + 1;
  hca = irf_panel('e PA deflux 2');  
  eint = [100 40000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint_figure).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f') ' eV'],[0.02 0.1],'color',0*[1 1 1],'fontsize',11)  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)  
  hca.CLim = [6.2 7];
  
  if 1 % plot trap angle
    %%
    hold(hca,'on')
    irf_plot(hca,trapping_angle(20,90,gseB1.tlim(EpochTT(['2017-07-25T22:09:35';'2017-07-25T22:09:52'])).abs),'k--')
    irf_plot(hca,trapping_angle(13,90,gseB1.tlim(EpochTT(['2017-07-25T22:09:52';'2017-07-25T22:09:59'])).abs),'k--')
    hold(hca,'off')
  end
  
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
end
if 0 % ePDist pa 
  isub = isub + 1;
  hca = irf_panel('e PA deflux');  
  eint = [5000 40000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint_figure).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f') ' eV'],[0.02 0.1],'color',0*[1 1 1],'fontsize',11)
  
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  
  hca.CLim = [6.2 7];
  
  if 1 % plot trap angle
    %%
    hold(hca,'on')
    irf_plot(hca,trapping_angle(20,90,gseB1.tlim(EpochTT(['2017-07-25T22:09:35';'2017-07-25T22:09:52'])).abs),'k--')
    irf_plot(hca,trapping_angle(13,90,gseB1.tlim(EpochTT(['2017-07-25T22:09:52';'2017-07-25T22:09:59'])).abs),'k--')
    hold(hca,'off')
  end
  
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
end


if 0 % Teperp, Tepar
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Te perp par');
  set(hca,'ColorOrder',mms_colors('1234'))
  refTi = 1;
  c_eval('irf_plot(hca,{Te?par.tlim(tint),Te?perp},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [10 400];  
  %irf_zoom(hca,'y')
  hca.YLim = [0 1e4];
end

if 1 % JdotE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JdotE');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{JdotE?.resample(gseJi?)},''comp'');',ic)    
  hca.YLabel.String = {'J\cdot E','(...)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
end
if 0 % fe(vpar)
  isub = isub + 1;
  hca = irf_panel('fe(vpar)');
  c_eval('fred = ef1D?_par;',ic)
  irf_spectrogram(hca,fred.specrec('velocity_1D'),'lin');  
  hca.YLim = fred.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{e||}','(10^3 km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end

drawnow
irf_zoom(h,'x',tint_figure)
drawnow
%irf_zoom(h(zoomy),'y')
irf_plot_axis_align(h)
h(1).Title.String = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));
drawnow

c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
irf_zoom(h,'x',tint_cavity)
irf_zoom(h,'y')

hca = irf_panel('e PA deflux 2');  hca.YLim = [0 180];
%% Distributions
if exist('hmark'); delete(hmark); end
times_dist = EpochTT(['2017-07-09T17:33:17.00Z';...
                      '2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:21.00Z';...
                      '2017-07-09T17:33:22.00Z']);
                    

times_dist = EpochTT(['2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:18.50Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:19.50Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:20.50Z']);

times_dist = EpochTT(['2017-07-09T17:33:19.00Z']) + 3*0.03*[-2:1:3];
                    
times_dist_exact = EpochTT([]);
for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  c_eval('v1 = mean(dmpaB?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  c_eval('v2 = mean(dslE?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  v2 = cross(v1,cross(v2,v1)); v2 = v2/norm(v2);
  v3 = cross(v1,v2);
  vlim = 70;
  if 0 % vE,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v2,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_E (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{E} (10^3 km/s)';
  end
  if 1 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
end

c_eval('h2(?).CLim = [-14.5 -9.5];',1:numel(h2))
hmark = irf_pl_mark(h,tocolumn(times_dist_exact_eu),'r');
c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))

%% Plot 2D ion distributions to show off cold ions
ic = 1;
[h,h2] = initialize_combined_plot(6,2,2,0.5,'vertical');
%h = irf_plot(10);

%tint_figure = tint_cavity + [-135 15];
tint_figure = tint_action;

times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.00Z';...
                      '2017-07-25T22:10:08.00Z';...
                      '2017-07-25T22:10:09.20Z']);

cs = 'lmn';
cs_str = 'LMN';
elim = [600 40000];

comps = ['x','y','z'];
comp_str = 'LMN';
comps = ['x'];
comp_str = 'x';

if 0
cs = 'gse';
cs_str = 'xyz';

comps = ['x','y','z'];
comp_str = 'xyz';
end

zoomy = []; isub = 1;
if 1 % B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z,%sB?.abs},''comp'');',cs,cs,cs,cs),ic)    
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3),'|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % vExB.x/y/z , Vi resample, ve
  for icomp = 1:numel(comps) % vExB.x , Vi resample, ve  
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ve',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVe?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}',comp_str(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % Vi.x/y/z
  for comp = comps % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['Vi',comp]);
    set(hca,'ColorOrder',mms_colors('2'))    
    c_eval(sprintf('irf_plot(hca,{gseVi?.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{(%s)}',comp_str),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'v_i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end

if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % elim
    %%
    hold(hca,'on')
    c_eval('timeline = iPDist?.time;',ic)    
    irf_plot(hca,irf.ts_scalar(timeline,700*ones([timeline.length,1])),'k')
    hold(hca,'off')
  end
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.XGrid = 'off';
  hca.YGrid = 'off'; 
end
if 1 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  c_eval('if1D = if1Dx?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z');
  c_eval('if1D = if1Dz?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end

colormap(pic_colors('candy4'))
drawnow
%irf_zoom(h,'x',tint_figure)
drawnow
irf_zoom(h(zoomy),'y')
irf_plot_axis_align(h)
%h(1).Title.String = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));
drawnow

%% Distributions
if exist('hmark'); delete(hmark); end
times_dist = EpochTT(['2017-07-09T17:33:17.00Z';...
                      '2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:21.00Z';...
                      '2017-07-09T17:33:22.00Z']);
                    

times_dist = EpochTT(['2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:18.50Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:19.50Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:20.50Z']);

times_dist = EpochTT(['2017-07-09T17:33:19.00Z']) + 3*0.03*[-2:1:3];




times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:09.00Z';...
                      '2017-07-25T22:10:15.00Z';...
                      %'2017-07-25T22:10:29.00Z';...
                      %'2017-07-25T22:10:49.00Z';...
                      %'2017-07-25T22:10:27.00Z';...
                      '2017-07-25T22:10:29.00Z';...
                    ]);
                    
                  
times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.00Z';...
                      '2017-07-25T22:10:08.00Z';...
                      '2017-07-25T22:10:09.20Z']);
times_dist_exact = EpochTT([]);

%elim = [700 Inf];

for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  v1 = [1 0 0];
  v2 = [0 1 0];
  v3 = [0 0 1];
  v1 = L;
  v2 = M;
  v3 = N;
  % dugcheck, slightly rotate system
%   rotang = 30;
%   v2n = v2*cosd(rotang) + v3*sind(rotang);
%   v3n = -v2*sind(rotang) + v3*cosd(rotang);
  vlim = 2300;
  vg = -vlim:100:vlim;
  if 1 % x,z
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+2*0.15*0.5*[-1 1]+0);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist(1).time.epochUnix;
    f2d = ddist.reduce('2D',v2,v3n);

    f2d.plot_plane(hca,'smooth',0);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = 'v_z (km/s)';
  end
  if 0 % x,y
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+0.15*0.5*[-1 1]);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vg',vg);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  colormap(hca,pic_colors('pasteljet'))
end
if 0
for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  c_eval('v1 = mean(dmpaB?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  c_eval('v2 = mean(dslE?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  v2 = cross(v1,cross(v2,v1)); v2 = v2/norm(v2);
  v3 = cross(v1,v2);
  vlim = 70;
  if 0 % vE,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v2,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_E (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{E} (10^3 km/s)';
  end
  if 1 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
end
end

hlinks = linkprop(h2,{'CLim','XLim','YLim'});

%c_eval('h2(?).CLim = [-14.5 -9.5];',1:numel(h2))
hmark = irf_pl_mark(h,tocolumn(times_dist_exact_eu),'r');
c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))

%% MMS SWT: Figure: Overview
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');
fontsize = 16;

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
  hca.YLim = [-300 1200]*0.99;
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1e-2 1e3];
end

elim_feeps = [8e4 Inf];
if 0 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  specrec = feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle');
  %specrec.f = f*1e-3;
  irf_spectrogram(hca,specrec);
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 1 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234'); 
  specrec = feeps_ion_omni1234.elim(elim_feeps).specrec('energy');
  specrec.f = specrec.f*1e-3;
  specrec.f_label =  {'E_i (keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(keV)'};   
end
if 1 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 1 % e DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % ve_par
    hold(hca,'on')
    c_eval('vexb = gseVe?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.me/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end


if 0 % i psd x,y,z, 3 panels
  for comp = ['L']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 5];
end
if 0 % Tepar/Teperp log
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  c_eval('irf_patch(hca,{tsdata,1});',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 2];
  hca.YLabel.Interpreter = 'tex';
end
if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  tsdata.data = log10(tsdata.data);
  c_eval('irf_patch(hca,{tsdata,0});',ic)  
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.301 0.301];
  hca.YLabel.Interpreter = 'tex';
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = 16;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = 16;',1:numel(h))

%hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-15 25];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];


hca1 = irf_panel('feeps omni mms 1234');  
hca2 = irf_panel('i DEF omni'); 
hca1.YLim(1) = hca2.YLim(2)*1e-3;


% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

%% MMS SWT: Figure: ongoing reconnection
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end
if 0 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 1 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1e-2 1e3];
end

elim_feeps = [8e4 Inf];
if 0 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  specrec = feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle');
  %specrec.f = f*1e-3;
  irf_spectrogram(hca,specrec);
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 1 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234'); 
  specrec = feeps_ion_omni1234.elim(elim_feeps).specrec('energy');
  specrec.f = specrec.f*1e-3;
  specrec.f_label =  {'E_i (keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(keV)'};   
end
if 0 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end


if 1 % i psd x,y,z, 3 panels
  for comp = ['L']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 5];
end
if 0 % Tepar/Teperp log
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  c_eval('irf_patch(hca,{tsdata,1});',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 2];
  hca.YLabel.Interpreter = 'tex';
end
if 1 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  tsdata.data = log10(tsdata.data);
  c_eval('irf_patch(hca,{tsdata,0});',ic)  
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.301 0.301];
  hca.YLabel.Interpreter = 'tex';
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 16;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = 16;',1:numel(h))

hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-15 25];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];


%hca1 = irf_panel('feeps omni mms 1234');  
%hca2 = irf_panel('i DEF omni'); 
%hca2.YLim(2) = hca1.YLim(1);
%h(end).YLim  = [-350 1150];

% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

%% MMS SWT: Figure: off-time current sheet structure, vi vExB
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

fontsize_leg = 14;

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
end 
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',fontsize_leg);
end
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize_leg);  
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-3:1:2];
  hca.YLabel.Interpreter = 'tex';
end
if 0 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E cs, resamp vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x.resample(gseVi?),%sE?.y.resample(gseVi?),%sE?.z.resample(gseVi?)},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
  irf_zoom(hca,'y')
end

elim_feeps = [8e4 Inf];
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 1 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval(sprintf('irf_plot(hca,{%sVe?perp.resample(%sVi?).(comp),%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',fontsize_leg);    
    %irf_legend(hca,{'v_e','v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  eval(sprintf('irf_plot(hca,{%sJcurl.x,%sJcurl.y,%sJcurl.z},''comp'');',cs,cs,cs))
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
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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

if 0 % EIS Pitch
  isub = isub + 1;
  hca = irf_panel('eis Pitch');  
  irf_spectrogram(hca,eis_pa2.elim([14 4500]*1e3).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('xyza'))  
  hca.YLim = [0 5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = 14;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);


hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [0 4]*0.99;

hca = irf_panel(['V ExB Vi ','x']); hca.YLim = [-100 250];
hca = irf_panel(['V ExB Vi ','y']); hca.YLim = [-250 699];
hca = irf_panel(['V ExB Vi ','z']); hca.YLim = [-300 300];
hca = irf_panel('E gsm'); hca.YLim = [-10 10];

c_eval('h(?).FontSize = 16;',1:numel(h))

% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

drawnow
%h(6).YLim(1) = h(7).YLim(2);
%h(2).YLim(1) = 0.001;

drawnow
%hmark = irf_pl_mark(h,time_mark,'k');
%c_eval('hmark(?).LineStyle = '':'';',1:5)
%c_eval('hmark(?).LineWidth = 1;',1:5)

%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
%hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);

%% MMS SWT: Figure: off-time current sheet structure, current
ic = 1;
ic = 1;

fontsize_leg = 14;

npanels = 4+3;
nrows = 1;
ncols = 1;
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
%h1 = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];
matlab_colors = pic_colors('matlab');
j_colors = [mms_colors('123'); matlab_colors(6,:)];
j_colors = [mms_colors('123'); matlab_colors(3,:)];
j_colors = [mms_colors('xyz'); 0 0 0];
%j_colors = matlab_colors;

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',fontsize_leg);
end 
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['x','y','z']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',j_colors)
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize_leg)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',j_colors)
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize_leg);  
  end
end
if 1 % Only  Jav ,curl
  for comp = ['x','y','z']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',j_colors([1 4],:))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize_leg)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',j_colors([1 4],:))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize_leg);  
  end
end
if 1 % Only  Jiav Jeav
  for comp = ['x','y','z']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s_ie',comp));
    set(hca,'ColorOrder',j_colors([2 3],:))
    if 0
      irf_plot(hca,{gseJiav.(comp),gseJeav.(comp)},'comp');
      %irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize_leg)
    elseif doResample    
      irf_plot(hca,{gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize_leg)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',j_colors([2 3],:))
    irf_legend(hca,{sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp)},[0.02 0.99],'fontsize',fontsize_leg);  
  end
end
if 0 % J, mms1-4,curl
  for comp = ['x','y','z']  
    species = 'i';    
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('123b'))
    %irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    eval(sprintf('irf_plot(hca,{gseJ%s1.(comp),gseJ%s2.(comp),gseJ%s3.(comp),gseJ%s4.(comp),gseJcurl.(comp)},''comp'');',species,species,species,species))
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s,%s}^{mms1}',species,comp),sprintf('J_{%s,%s}^{mms2}',species,comp),sprintf('J_{%s,%s}^{mms3}',species,comp),sprintf('J_{%s,%s}^{mms4}',species,comp),'curl'},[0.02 0.99],'fontsize',12);  
  end
end

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
irf_pl_mark(h1(1),tint_harris)

h1(2).YLim = [-20 20];
h1(3).YLim = [-14 19];
h1(4).YLim = [-12 4];
c_eval('h1(?).FontSize = 16;',1:numel(h1))
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 0.5;',1:numel(hl))





%
isub = 1;
if 0 % Jy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJav.tlim(tint_harris).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint_harris).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 1 % Jy,Jicurl, binning
  hca = h2(isub); isub = isub + 1;
  xx = gseJav.tlim(tint_harris).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint_harris).y.data;
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -20;
  jmax = -jmin;
  jstep = 0.5;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = '<J_{y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,(N))
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  hca.CLim = [0 35];
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))  
  plot(hca,j_edges,j_edges,'color',[0.5 0.5 0.5],'linestyle','--')
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = [jmin jmax];
  hca.YLim = [jmin jmax]; 
  axis(hca,'square')
  drawnow
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %axis(hca,'equal')
end


c_eval('h2(?).FontSize = 16;',1:numel(h2))

%% MMS SWT: Figure: onset 
% for feeps_ion_omni1234, see feeps_data.m
ic = 1;
time_mark = EpochTT('2017-07-25T22:09:08.00Z');
time_df = EpochTT('2017-07-25T22:10:06.50Z');

tint_figure = irf.tint('2017-07-25T22:07:00.00Z/2017-07-25T22:11:30.00Z');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

fontsize_leg = 14;

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-3:1:2];
  hca.YLabel.Interpreter = 'tex';
end
if 0% Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end

elim_feeps = [8e4 Inf];
if 1 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 1 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,feeps_ion_omni1234.elim(elim_feeps).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  eval(sprintf('irf_plot(hca,{%sJcurl.x,%sJcurl.y,%sJcurl.z},''comp'');',cs,cs,cs))
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
if 0 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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

if 0 % EIS Pitch
  isub = isub + 1;
  hca = irf_panel('eis Pitch');  
  irf_spectrogram(hca,eis_pa2.elim([14 4500]*1e3).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 1 % i psd x,y,z, 3 panels
  for comp = ['x']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('xyza'))  
  hca.YLim = [0 5];
end

if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{iL}','v_{iM}','v_{iN}'}',[1.02 0.9],'fontsize',12);
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);


hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca1 = irf_panel('feeps omni mms 1234');  
hca2 = irf_panel('i DEF omni'); 
hca2.YLim(2) = hca1.YLim(1);
h(end).YLim = [-350 1150];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [0 4]*0.99;


c_eval('h(?).FontSize = fontsize_leg;',1:numel(h))


drawnow
%h(6).YLim(1) = h(7).YLim(2);
%h(2).YLim(1) = 0.001;

drawnow
hmark = irf_pl_mark(h,time_mark,'k');
c_eval('hmark(?).LineStyle = '':'';',1:numel(hmark))
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))


hmark2 = irf_pl_mark(h,time_df,'k');
c_eval('hmark2(?).LineStyle = '':'';',1:numel(hmark2))
c_eval('hmark2(?).LineWidth = 1;',1:numel(hmark2))

%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
%hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
hat=annotation('textarrow',[0.593 0.593]+0.0094,[.56 .54]+0.08,'string',{'first energetic ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none');
hat=annotation('textarrow',[0.593 0.593]+0.104,[.4 .38]+0.01,'string',{'first cold ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none');
%

%% MMS SWT: detailed look at cold ions
% Plot 2D ion distributions to show off cold ions
ic = 1;
[h,h2] = initialize_combined_plot(7,2,1,0.5,'vertical');
c_eval('h(?).Position(1) = 0.10;',1:numel(h))
%h = irf_plot(10);

%tint_figure = tint_cavity + [-135 15];


tint_figure = tint_action;
tint_figure = irf.tint('2017-07-25T22:09:50.00Z/2017-07-25T22:10:20.00Z');

times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.00Z';...
                      '2017-07-25T22:10:08.00Z';...
                      '2017-07-25T22:10:09.20Z']);

fontsize_leg = 14;
                    
cs = 'gse';
cs_str = 'xyz';
elim = [600 40000];

comps = ['x','y','z'];
comp_str = 'LMN';
comps = ['x'];
comp_str = 'x';

if 0
cs = 'gse';
cs_str = 'xyz';

comps = ['x','y','z'];
comp_str = 'xyz';
end

zoomy = []; isub = 1;
if 1 % B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z,%sB?.abs},''comp'');',cs,cs,cs,cs),ic)    
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{cs_str(1),cs_str(2),cs_str(3),'|B|'}',[1.02 0.9],'fontsize',fontsize_leg);  
end 
if 1 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98],'fontsize',fontsize_leg)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
end

if 1 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.tlim(tint_figure).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % elim
    %%
    hold(hca,'on')
    c_eval('timeline = iPDist?.time;',ic)    
    irf_plot(hca,irf.ts_scalar(timeline,700*ones([timeline.length,1])),'k')
    hold(hca,'off')
  end
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YLabel.Interpreter = 'tex';  
  hca.XGrid = 'off';
  hca.YGrid = 'off'; 
end

if 1 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  c_eval('if1D = if1Dx?_700.tlim(tint_figure);',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iL}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?_700.tlim(tint_figure);',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iM}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z');
  c_eval('if1D = if1Dz?_700.tlim(tint_figure);',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iN}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end


hlinks = linkprop([irf_panel('iLine x'),irf_panel('iLine y'),irf_panel('iLine z')],{'CLim'});
colormap(pic_colors('candy4'))
drawnow
irf_zoom(h,'x',tint_figure)
drawnow
%irf_zoom(h(zoomy),'y')
irf_plot_axis_align(h)
%h(1).Title.String = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));
drawnow

% Distributions
if exist('hmark'); delete(hmark); end
if exist('hmark2'); delete(hmark); end
times_dist = EpochTT(['2017-07-09T17:33:17.00Z';...
                      '2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:21.00Z';...
                      '2017-07-09T17:33:22.00Z']);
                    

times_dist = EpochTT(['2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:18.50Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:19.50Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:20.50Z']);

times_dist = EpochTT(['2017-07-09T17:33:19.00Z']) + 3*0.03*[-2:1:3];




times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:09.00Z';...
                      '2017-07-25T22:10:15.00Z';...
                      %'2017-07-25T22:10:29.00Z';...
                      %'2017-07-25T22:10:49.00Z';...
                      %'2017-07-25T22:10:27.00Z';...
                      '2017-07-25T22:10:29.00Z';...
                    ]);
                    
                  
times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.00Z';...
                      '2017-07-25T22:10:08.00Z';...
                      '2017-07-25T22:10:09.20Z']);
                    
                    
times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.50Z']);
                    
times_dist = EpochTT(['2017-07-25T22:10:03.00Z';...
                      '2017-07-25T22:10:07.50Z']);
                    
times_dist_exact = EpochTT([]);
clear times_dist_exact_eu hmarkb hmark_tmp times_dist_exact_eu_edges

%elim = [700 Inf];

for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  v1 = [1 0 0];
  v2 = [0 1 0];
  v3 = [0 0 1];
  v1 = R(1,:);
  v2 = R(2,:);
  v3 = R(3,:);
  x_str = 'v_L (km/s)';
  y_str = 'v_M (km/s)';
  z_str = 'v_N (km/s)';
  % bugcheck, slightly rotate system
%   rotang = 30;
%   v2n = v2*cosd(rotang) + v3*sind(rotang);
%   v3n = -v2*sind(rotang) + v3*cosd(rotang);
%   v2 = v2n; v3 = v3n;
  vlim = 2300;
  vg = -vlim:100:vlim;
  if 1 % x,z
    
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+6*0.15*0.5*[-1 1]+0);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = mean(ddist.time.epochUnix);
    times_dist_exact_eu_edges(id,1) = ddist.time(1).epochUnix-0.03*0.5;
    times_dist_exact_eu_edges(id,2) = ddist.time(end).epochUnix+0.03*0.5;
    f2d = ddist.elim([00 Inf]).reduce('2D',v1,v3);

    f2d.plot_plane(hca,'smooth',0);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = x_str;
    hca.YLabel.String = z_str;
    
    hmark_tmp = irf_pl_mark(h,ddist.time,'k','facealpha',0.1);
    drawnow
    %c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))

  end
  if 0 % x,y
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+0.15*0.5*[-1 1]);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vg',vg);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  colormap(hca,pic_colors('pasteljet'))
end
if 0
for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  c_eval('v1 = mean(dmpaB?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  c_eval('v2 = mean(dslE?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  v2 = cross(v1,cross(v2,v1)); v2 = v2/norm(v2);
  v3 = cross(v1,v2);
  vlim = 70;
  if 0 % vE,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v2,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_E (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{E} (10^3 km/s)';
  end
  if 1 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
end
end

hlinks = linkprop(h2,{'CLim','XLim','YLim'});
h2(1).CLim = [-10 -7.];

c_eval('h(?).FontSize = 16;',1:numel(h))
c_eval('h2(?).FontSize = 16;',1:numel(h2))
c_eval('h2(?).Position(1) = 0.58;',1:numel(h2))

h(end).XTickLabelRotation = 0;
h2(1).Title.String = 'Hot pre-existing plasma sheet';
h2(2).Title.String = 'Post-DF cold and heated ions';
%c_eval('h2(?).CLim = [-14.5 -9.5];',1:numel(h2))
%hmark = irf_pl_mark(h,tocolumn(times_dist_exact_eu),'r');
%c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))
hmark2 = irf_pl_mark(h,tocolumn(times_dist_exact_eu_edges(:)),'r');
c_eval('hmark2(?).LineWidth = 0.5; hmark2(?).Color = [0 0 0]; hmark2(?).LineStyle = ''--'';',1:numel(hmark2))

%% MMS SWT: Harris current sheet (?)
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
units = irf_units;

% Parameters for Harris fit.
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L 
%    = -(B0*((Bx/B0)^2 - 1))/L
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;
mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);
% mf_L_B = solve(Bx,L);
% mf_L_J = solve(Jy,L);


cmap_time = irf_colormap('waterfall');
%cmap_time = colormap('parula');
cmap_time = pic_colors('candy4');
cmap_count = pic_colors('candy4');
cmap_count = flipdim(colormap('bone'),1);
cmap_count = colormap('hot');
cmap_count = flipdim(pic_colors('thermal'),1);

iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

% Parameters for current Harris current sheet
L = [1000 1500 2000 2500 4000]*1e3;
L = [1000 2000 3000 4000]*1e3;
colors_L = pic_colors('matlab');
B0 = 21e-9;
J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  
fontsize = 12;

% Prepare figure
npanels = 5;
nrows = 4;
ncols = 2;

h1 = irf_plot(2);
h1(1).Position = [0.15 0.75 0.7 0.2];
h1(2).Position = [0.15 0.55 0.7 0.2];
h2(1) = axes('position',[0.15 0.10 0.3 0.3]);
h2(2) = axes('position',[0.55 0.10 0.3 0.3]);

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z},'comp');
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end 

if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize);  
  end
end
if 1 % J, Jeav, Jiav, Jav
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp)},[0.02 0.99],'fontsize',fontsize);  
  end
end

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
hmark = irf_pl_mark(h1(1),tint_harris);

% Non-TSeries panels.
isub = 1;

for species = {'','i','e'}
  if 0 % Jy,Jcurl
    hca = h2(isub); isub = isub + 1;
    eval(sprintf('xx = gseJ%sav.tlim(tint).y.data;',species{1}))
    yy = gseJcurl.resample(gseJav).tlim(tint).y.data;
    irem = find(xx == 0);
    xx(irem) = [];
    yy(irem) = [];
    p = polyfit(xx,yy,1);
    J_edges1 = linspace(-10,10,41);
    J_edges2 = linspace(-10,10,39);
    [N edges mid loc] = histcn([xx(:) yy(:)],J_edges1,J_edges2);
    N(N==0) = NaN;
    pcolor(hca,mid{:},log10(N)')
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log_{10}(counts)';
    colormap(pic_colors('candy4'))
    %plot(hca,xx,yy,'.')
    %histcn_plot(hca,xx,yy)
    hca.XLabel.String = sprintf('<J_{%sy}^{FPI}>',species{1});
    hca.YLabel.String = 'J_y^{curl}';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
    plot(hca,[-100 100],p(2)+[-100 100]*p(1))
    hca.XLim = xlim;
    hca.YLim = ylim;
    hold(hca,'off')
    %hca.XLim = 40*[-1 1];
    %hca.YLim = 40*[-1 1]; 
    hca.XGrid = 'on';
    hca.YGrid = 'on'; 
  end
end
if 0 % Jey,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJiav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJiav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{i,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 0 % Jiy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJeav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{e,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 1 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = 'J_{y}^{curl}';
  dt = ts_yy.time(2)-ts_yy.time(1);
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N*dt)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  colormap(hca,cmap_count)
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'time (s)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  clear hl
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1,'color',colors_L(ileg,:));
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 1 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.05:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); 
  %ts_yy = gseJcurl;
  %y_str = '<J_{y}^{curl}>';
  y_str = 'J_{y}^{curl}';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,4,tt,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  colormap(hca,cmap_time)
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1,'color',colors_L(ileg,:));
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  hca.Box = 'on';
end
if 0 % Bx,JxBz, binned
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  %ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  %irem = find(xx == 0);
  %xx(irem) = [];
  %yy(irem) = [];
  %p = polyfit(xx,yy,1);
  JxB_edges = linspace(-2,12,70);
  B_edges = linspace(0,24,71);
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end
if 0 % Bx,JxBz, scatter
  hca = h2(isub); isub = isub + 1;
  
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data*1e18;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  cc = ts_yy.tlim(tint_harris).time - ts_yy.tlim(tint_harris).time(1);
  scatter(hca,xx.^2,yy,1,cc)
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time since tt';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  if 0
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  end
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end

h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];
colormap(hca,cmap_time)


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

h1(1).FontSize = fontsize;
h1(2).FontSize = fontsize;
h2(1).FontSize = fontsize;
h2(2).FontSize = fontsize;

hlegs.Position(1) = 0.3;
hlegs.FontSize = 16;

fontsize_labels = 16;


irf_legend(h1(1),'(a)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h1(2),'(b)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h2(1),'(c)',[-0.1 0.98],'fontsize',fontsize_labels);
irf_legend(h2(2),'(d)',[-0.1 0.98],'fontsize',fontsize_labels);

h1(1).FontSize = 16;
h1(2).FontSize = 16;
h2(1).FontSize = 16;
h2(2).FontSize = 16;

%% MMS SWT: Harris current sheet (?), extra poanel with P
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
units = irf_units;

% Parameters for Harris fit.
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L 
%    = -(B0*((Bx/B0)^2 - 1))/L
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy
mf_L = @(Bx,B0,Jy) -(B0.*((Bx./B0).^2 - 1))./Jy;
mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);
% mf_L_B = solve(Bx,L);
% mf_L_J = solve(Jy,L);


cmap_time = irf_colormap('waterfall');
%cmap_time = colormap('parula');
cmap_time = pic_colors('candy4');
cmap_count = pic_colors('candy4');
cmap_count = flipdim(colormap('bone'),1);
cmap_count = colormap('hot');
cmap_count = flipdim(pic_colors('thermal'),1);

iisub = 0;
cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];

% Parameters for current Harris current sheet
L = [1000 1500 2000 2500 4000]*1e3;
L = [1000 2000 3000 4000]*1e3;
colors_L = pic_colors('matlab');
B0 = 21e-9;
J_edges = linspace(-2,12,70);
B_edges = linspace(0,22,71);
  
fontsize = 16;

% Prepare figure
npanels = 5;
nrows = 4;
ncols = 2;

h1 = irf_plot(3);
h1(1).Position = [0.15 0.75+0.05 0.7 0.17];
h1(2).Position = [0.15 0.58+0.05 0.7 0.17];
h1(3).Position = [0.15 0.41+0.05 0.7 0.17];
h2(1) = axes('position',[0.15 0.10 0.3 0.25]);
h2(2) = axes('position',[0.55 0.10 0.3 0.25]);

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{lmnBav.x,lmnBav.y,lmnBav.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize);
end 

if 1 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize);  
  hca(1).YLim = [0 0.499];
end
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize);  
  end
end
if 1 % J, Jeav, Jiav, Jav
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp)}',[1.02 0.90],'fontsize',fontsize);  
  end
  hca.YLim = 0.99*[-10 20];
end

irf_zoom(h1,'x',tint)
%irf_zoom(h1(zoomy),'y')
irf_plot_axis_align
hmark = irf_pl_mark(h1(1),tint_harris);

% Non-TSeries panels.
isub = 1;

for species = {'','i','e'}
  if 0 % Jy,Jcurl
    hca = h2(isub); isub = isub + 1;
    eval(sprintf('xx = gseJ%sav.tlim(tint).y.data;',species{1}))
    yy = gseJcurl.resample(gseJav).tlim(tint).y.data;
    irem = find(xx == 0);
    xx(irem) = [];
    yy(irem) = [];
    p = polyfit(xx,yy,1);
    J_edges1 = linspace(-10,10,41);
    J_edges2 = linspace(-10,10,39);
    [N edges mid loc] = histcn([xx(:) yy(:)],J_edges1,J_edges2);
    N(N==0) = NaN;
    pcolor(hca,mid{:},log10(N)')
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log_{10}(counts)';
    colormap(pic_colors('candy4'))
    %plot(hca,xx,yy,'.')
    %histcn_plot(hca,xx,yy)
    hca.XLabel.String = sprintf('<J_{%sy}^{FPI}>',species{1});
    hca.YLabel.String = 'J_y^{curl}';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
    plot(hca,[-100 100],p(2)+[-100 100]*p(1))
    hca.XLim = xlim;
    hca.YLim = ylim;
    hold(hca,'off')
    %hca.XLim = 40*[-1 1];
    %hca.YLim = 40*[-1 1]; 
    hca.XGrid = 'on';
    hca.YGrid = 'on'; 
  end
end
if 0 % Jey,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJiav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJiav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{i,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 0 % Jiy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJeav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{e,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 1 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = 'J_{y}^{curl}';
  dt = ts_yy.time(2)-ts_yy.time(1);
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},(N*dt)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  colormap(hca,cmap_count)
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'time (s)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  clear hl
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1,'color',colors_L(ileg,:));
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
end
if 1 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.05:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); 
  %ts_yy = gseJcurl;
  %y_str = '<J_{y}^{curl}>';
  y_str = 'J_{y}^{curl}';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,6,tt,'o')  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  colormap(hca,cmap_time)
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1,'color',colors_L(ileg,:));
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  hca.Box = 'on';
end
if 0 % Bx,JxBz, binned
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  %ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  %irem = find(xx == 0);
  %xx(irem) = [];
  %yy(irem) = [];
  %p = polyfit(xx,yy,1);
  JxB_edges = linspace(-2,12,70);
  B_edges = linspace(0,24,71);
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end
if 0 % Bx,JxBz, scatter
  hca = h2(isub); isub = isub + 1;
  
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data*1e18;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  cc = ts_yy.tlim(tint_harris).time - ts_yy.tlim(tint_harris).time(1);
  scatter(hca,xx.^2,yy,1,cc)
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time since tt';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  if 0
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  end
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end

h2(2).YTick = h2(1).YTick;
h2(2).XTick = h2(1).XTick;


% Clone colorbar for time and place on top of panel 1.
hca = h1(1);
h_pos = hca.Position;
hb = colorbar('peer',hca,'location','northoutside');
hb.YTick = [];
colormap(hca,cmap_time)


xlim = hca.XLim;
xmark = [min(hmark.XData) max(hmark.XData)];
x1_rel = xmark(1)/diff(xlim);
x2_rel = xmark(2)/diff(xlim);

hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
drawnow
hb.Position(2) = hca.Position(2) + hca.Position(4);

h1(1).FontSize = fontsize;
h1(2).FontSize = fontsize;
h2(1).FontSize = fontsize;
h2(2).FontSize = fontsize;

hlegs.Position(1) = 0.3;
hlegs.FontSize = 16;

fontsize_labels = 16;


irf_legend(h1(1),'(a)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h1(2),'(b)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h1(3),'(c)',[-0.05 0.98],'fontsize',fontsize_labels);
irf_legend(h2(1),'(d)',[-0.1 0.98],'fontsize',fontsize_labels);
irf_legend(h2(2),'(e)',[-0.1 0.98],'fontsize',fontsize_labels);

h1(1).FontSize = 16;
h1(2).FontSize = 16;
h1(3).FontSize = 16;
h2(1).FontSize = 16;
h2(2).FontSize = 16;

%% EGU 2023: Figure: Overview
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

dt_resample = 1;
timeline = tint(1):dt_resample:tint(2);

npanels = 5;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');
fontsize = 16;
fontsize_leg = 15;

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x ','B_y ','B_z '},[0.05 0.9],'fontsize',12);
end 
if 1 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?.resample(timeline)},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.05 0.9])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  c_eval('irf_plot(hca,{gseVHp?_srvy.x,gseVHp?_srvy.y,gseVHp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{ix} ','v_{iy} ','v_{iz} '},[0.05 0.9],'fontsize',fontsize_leg);
  hca.YLim = [-300 1200]*0.99;
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1e-2 1e3];
end

if 1 % i psd x,y,z, 3 panels
  for comp = ['x']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end

elim_feeps = [8e4 Inf];

if 0 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  specrec = feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle');
  %specrec.f = f*1e-3;
  irf_spectrogram(hca,specrec);
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 0 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234'); 
  specrec = feeps_ion_omni1234.elim(elim_feeps).specrec('energy');
  specrec.f = specrec.f*1e-3;
  specrec.f_label =  {'E_i (keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(keV)'};   
end
if 1 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
end
if 0 % e DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % ve_par
    hold(hca,'on')
    c_eval('vexb = gseVe?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.me/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end



if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 5];
end
if 0 % Tepar/Teperp log
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  c_eval('irf_patch(hca,{tsdata,1});',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 2];
  hca.YLabel.Interpreter = 'tex';
end
if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp patch');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_patch(hca,{Te?par/Te?perp,1});',ic) 
  c_eval('tsdata = Te?par.resample(gseVi?)/Te?perp.resample(gseVi?);',ic)
  tsdata.data = log10(tsdata.data);
  c_eval('irf_patch(hca,{tsdata,0});',ic)  
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.301 0.301];
  hca.YLabel.Interpreter = 'tex';
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = 16;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = 16;',1:numel(h))

%hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-15 25];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];


%hca1 = irf_panel('feeps omni mms 1234');  
%hca2 = irf_panel('i DEF omni'); 
%hca1.YLim(1) = hca2.YLim(2)*1e-3;


% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end
%%
c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

c_eval('h(?).LineWidth = 1;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('ht(?).FontWeight = ''bold'';',1:numel(ht))
c_eval('h(?).FontWeight = ''bold'';',1:numel(h))

%% EGU 2023, current sheet pressure
% Since the total pressure is changing, it is likely that B0 is also
% changing. So recalculate a varying B0 based on the total pressure.
% B^2/2mu0 + Pi + Pe = Ptot = B0^2/2*mu0 (RHS is asymptotical lobe field)
%  B0 = (2*mu0*Ptot)^0.5
tsPtot = irf.ts_scalar(gsePi1.time,gsePe1.resample(gsePi1).trace.data/3+PB1.resample(gsePi1).data+gsePi1.trace.data/3);
tsB0 = irf.ts_scalar(tsPtot.time,sqrt(tsPtot.data*1e-9*2*units.mu0)*1e9); tsB0.name = 'B_0'; tsB0.units = 'nT'; % nT

cmap = colormap(pic_colors('candy4'));
dt_resample = 0.5;
timeline = tint(1):dt_resample:tint(2);
doResample = 1;
isub = 0;
zoomy = [];
matlab_colors = pic_colors('matlab');
j_colors = [mms_colors('123'); matlab_colors(6,:)];
j_colors = [mms_colors('123'); matlab_colors(3,:)];
j_colors = [mms_colors('xyz'); 0 0 0];

plot_colors = [0.8 0.8 0.8; j_colors([4 3],:)];
leg_colors = [0.6 0.6 0.6; j_colors([4 3],:)];
  
fontsize = 16;
fontsize_leg = 14;

h = irf_plot(4);

if 1 % B gse + tsB0 from total pressure
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{lmnBav.x,lmnBav.y,lmnBav.z,tsB0},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','B_0 = 2\mu_0 P_{tot}^{1/2}'},[0.02 0.09],'fontsize',fontsize);
end 
if 0 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{lmnBav.x,lmnBav.y,lmnBav.z},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'}',[0.02 0.02],'fontsize',fontsize);
end 

if 1 % Pressures, PB, Pi, Pe
  colors = mms_colors('matlab');
  colors(4,:) = [0 0 0];
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',colors)  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',colors)
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'},[0.02 0.95],'fontsize',fontsize);  
  hca(1).YLim = [0 0.499];
end
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k','fontsize',fontsize)
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',fontsize);  
  end
end

if 1 % Only  Jav ,curl
  for comp = comps      
    isub = isub + 1;    
    hca = irf_panel(sprintf('J%s',comp));    
    set(hca,'ColorOrder',plot_colors)
    irf_plot(hca,{gseJav.(comp),gseJav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');        
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};    
    set(hca,'ColorOrder',leg_colors)    
    irf_legend(hca,{sprintf('<J_{%s}^{FPI}> (%g ms)',comp,30),sprintf('<J_{%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('J_{%s}^{curl} (%g s)',comp,dt_resample)}',[0.01 0.99],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('J_{%s}^{curl} (%g s)',comp,dt_resample)},[0.02 0.99],'fontsize',fontsize_leg);      
    %hl = findobj(hca,'type','line');    
    %c_eval('hl(?).LineWidth = 0.5;',1:2)
  end  
end
if 1 % Only  Jiav Jeav
  for comp = comps
    isub = isub + 1;    
    hca = irf_panel(sprintf('J%s_ie',comp));            
    set(hca,'ColorOrder',plot_colors)
    irf_plot(hca,{gseJeav.(comp),gseJeav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),},'comp');    
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',leg_colors)     
    
    irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,30),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.99],'fontsize',fontsize_leg);          
    set(hca,'ColorOrder',leg_colors(3,:))     
    irf_legend(hca,{sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.07],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.05],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)},[0.02 0.99],'fontsize',fontsize_leg);  
    
    %hl = findobj(hca,'type','line');    
    %c_eval('hl(?).LineWidth = 0.5;',1:2)    
  end
end


irf_zoom(h,'x',tint)
%irf_zoom(h1(zoomy),'y')
irf_plot_axis_align


h(1).YLim = [-18 35]*0.99;
h(3).YLim = [-15 30]*0.99;
h(4).YLim = [-20 20]*0.99;

%c_eval('h1(?).Position(1) = 0.08;',1:numel(h1))
c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

c_eval('h(?).LineWidth = 1;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('ht(?).FontWeight = ''bold'';',1:numel(ht))
c_eval('h(?).FontWeight = ''bold'';',1:numel(h))


%% Timeline of Harris sheet thickness, taking into account potentially varying B0
tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');
tint_L = tint_harris + [0 15];
units = irf_units;

% Parameters for Harris fit. Want an expression Bx(Jy)
syms z z0 L B0
Bx = B0*tanh((z-z0)/L);
Jy = gradient(Bx,z)/units.mu0; % -(B0*(tanh((z - z0)/L)^2 - 1))/L
JxBz = -Jy*Bx;

% Jy = -(B0*(tanh((z - z0)/L)^2 - 1))/L/mu0 
%    = -(B0*((Bx/B0)^2 - 1))/L/mu0
%     
% L  = -(B0*((Bx/B0)^2 - 1))/Jy/mu0
mf_L = @(Bx,B0,Jy) (B0./Jy/units.mu0).*(1-(Bx./B0).^2);

mf_Bx = matlabFunction(Bx);
mf_Jy = matlabFunction(Jy);
mf_JxBz = matlabFunction(JxBz);

mf_L_B = solve(Bx,L);
mf_L_J = solve(Jy,L);

% Calculate L directly from J and Bx
B0 = 22e-9;
dt_resample = 0.5;
dt_L = 1;
timeline = tint(1):dt_resample:tint(2);

dataL = mf_L(gseBav.x.data*1e-9,B0,Jcurl.resample(gseBav).y.data)*1e-3; % m -> km
dataL(abs(dataL)>prctile(dataL,95)) = NaN;
dataL(Jcurl.resample(gseBav).y.data<0) = NaN;
tsL = irf.ts_scalar(gseBav.time,dataL);
tsL = tsL.tlim(tint_L);

% Data is quite scattered, so make specrec of data with binning
%L_edges = linspace(0,max(dataL),20);
L_edges = linspace(0,10000*0.99,20);
t_edges = (tint(1)+-dt_L*0.5):dt_L:(tint(2)+dt_L*0.5);
[N,edges,mid,loc] = histcn([tsL.data,tsL.time-t_edges(1)], L_edges, t_edges-t_edges(1));
specrecL.p_label = 'counts';
specrecL.f_label = 'L (km)';
specrecL.p = N';
specrecL.t = irf_time(t_edges(1) + mid{2},'EpochTT>epoch');
specrecL.f = mid{1};
Lpeak = nan(numel(mid{2}),1);
for it = 1:numel(mid{2})
  [PKS,LOCS] = findpeaks(N(:,it),'NPeaks',1,'SortStr','descend','minpeakprominence',5);
  if not(isempty(LOCS))
    Lpeak(it) = mid{1}(LOCS);
  end
end
tsLpeak = irf.ts_scalar(t_edges(1) + mid{2},Lpeak);

% Since the total pressure is changing, it is likely that B0 is also
% changing. So recalculate a varying B0 based on the total pressure.
% B^2/2mu0 + Pi + Pe = Ptot = B0^2/2*mu0 (RHS is asymptotical lobe field)
%  B0 = (2*mu0*Ptot)^0.5
tsPtot = irf.ts_scalar(gsePi1.time,gsePe1.resample(gsePi1).trace.data/3+PB1.resample(gsePi1).data+gsePi1.trace.data/3);
tsB0 = irf.ts_scalar(tsPtot.time,sqrt(tsPtot.data*1e-9*2*units.mu0)*1e9); tsB0.name = 'B_0'; tsB0.units = 'nT'; % nT

dataL_varB0 = mf_L(gseBav.x.data*1e-9,tsB0.resample(gseBav).data*1e-9,Jcurl.resample(gseBav).y.data)*1e-3; % m -> km
dataL_varB0(abs(dataL_varB0)>prctile(dataL_varB0,95)) = NaN;
dataL_varB0(Jcurl.resample(gseBav).y.data<0) = NaN;
tsL_varB0 = irf.ts_scalar(gseBav.time,dataL_varB0);
tsL_varB0 = tsL_varB0.tlim(tint_L);
% Data is quite scattered, so make specrec of data with binning
%L_edges = linspace(0,max(dataL_varB0),20);
L_edges = linspace(0,10000*0.99,20);
t_edges = (tint(1)+-dt_L*0.5):dt_L:(tint(2)+dt_L*0.5);
[N_varB0,edges,mid,loc] = histcn([tsL_varB0.data,tsL_varB0.time-t_edges(1)], L_edges, t_edges-t_edges(1));
specrecL_varB0.p_label = 'counts';
specrecL_varB0.f_label = 'L (km)';
specrecL_varB0.p = N_varB0';
specrecL_varB0.p(specrecL_varB0.p==0) = NaN;
specrecL_varB0.t = irf_time(t_edges(1) + mid{2},'EpochTT>epoch');
specrecL_varB0.f = mid{1};

Lpeak_varB0 = nan(numel(mid{2}),1);
for it = 1:numel(mid{2})
  [PKS,LOCS] = findpeaks(N_varB0(:,it),'NPeaks',1,'SortStr','descend','minpeakprominence',5);
  if not(isempty(LOCS))
    Lpeak_varB0(it) = mid{1}(LOCS);
  end
end
tsLpeak_varB0 = irf.ts_scalar(t_edges(1) + mid{2},Lpeak_varB0);
%tsL = tsL.resample(timeline);

% Setup figure
fontsize = 18;
fontsize_leg = 14;
npanels = 2;
nrows = 1;
ncols = 1;

gca = gcf;
doResize = 0;
if not(isempty(gca))
  fig_position = get(gca, 'Position');
  doResize = 1;
end
[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.5,'vertical');
c_eval('h1(?).Position(1) = 0.10;',1:numel(h1))

if doResize
  fig = h1.Parent;
  set(fig,'position',fig_position)
end
iisub = 0;
cmap = colormap(pic_colors('candy4'));
doResample = 1;
isub = 0;
zoomy = [];

L = [1000 2000 3000 4000]*1e3;

J_edges = linspace(-2,12,50);
B_edges = linspace(0,24,71);
B2_edges = linspace(0,24.^2,51);
  

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B lmn');
  set(hca,'ColorOrder',mms_colors('xyzab'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  %irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z,tsB0.resample(tsB0.time(1):dt_resample:tsB0.time(end))},'comp');
  irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z,tsB0.resample(tsB0.time(1):dt_resample:tsB0.time(end)),gseJcurl.(comp).resample(timeline)},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyzab'))
  irf_legend(hca,{'B_x ','B_y ','B_z ','B_0 = 2\mu_0P_{tot}^{1/2}        ','  J_y^{curl}'},[0.02 0.07],'fontsize',fontsize_leg);
end 
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_plot(hca,{gseJcurl.(comp)},'comp');  
    
    hca.YLabel.String = {sprintf('J_%s^{curl}',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    %irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',12);  
  end
end
if 0 % J, Jeav, Jiav, Jav ,curl
  for comp = ['y']      
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    if doResample    
      irf_plot(hca,{gseJav.(comp).resample(timeline),gseJiav.(comp).resample(timeline),gseJeav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');
      irf_legend(hca,{sprintf('resampled to %g s (for visibility)',dt_resample)},[0.02 0.05],'color','k')
    else
      irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    end
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s}^{FPI}',comp),sprintf('J_{i%s}^{FPI}',comp),sprintf('J_{e%s}^{FPI}',comp),sprintf('J_{%s}^{curl}',comp)},[0.02 0.99],'fontsize',12);  
  end
end
if 0 % J, mms1-4,curl
  for comp = ['x','y','z']  
    species = 'i';    
    isub = isub + 1;
    %zoomy = [zoomy isub];
    hca = irf_panel(sprintf('J%s',comp));
    set(hca,'ColorOrder',mms_colors('1234a'))
    %irf_plot(hca,{gseJ1.(comp),gseJ2.(comp),gseJ3.(comp),gseJ4.(comp),gseJcurl.(comp)},'comp');  
    eval(sprintf('irf_plot(hca,{gseJ%s1.(comp),gseJ%s2.(comp),gseJ%s3.(comp),gseJ%s4.(comp),gseJcurl.(comp)},''comp'');',species,species,species,species))
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234a'))
    irf_legend(hca,{sprintf('J_{%s,%s}^{mms1}',species,comp),sprintf('J_{%s,%s}^{mms2}',species,comp),sprintf('J_{%s,%s}^{mms3}',species,comp),sprintf('J_{%s,%s}^{mms4}',species,comp),'curl'},[0.02 0.99],'fontsize',12);  
  end
end

if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize);  
  hca(1).YLim = [0 0.499];
end
if 0 % B0
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B0');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{tsB0},'comp');
  hca.YLabel.String = {'B_0','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'const B_0'},[0.98 0.9],'fontsize',12);
end 
if 0 % L
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('L');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{tsL},'comp');
  hca.YLabel.String = {'L','(km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'const B_0'},[0.98 0.9],'fontsize',12);
end 
if 0 % L specrec
  isub = isub + 1;
  hca = irf_panel('L specrec');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,specrecL,''log'');',2)  
  set(hca,'yscale','lin');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]); 
  hcb.YLabel.String = 'log_{10} counts';
  colormap(hca,cmap)   
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  irf_legend(hca,{sprintf('B_0 = %.0f nT',B0*1e9)},[0.98 0.98],'color','k','fontsize',fontsize_leg)
  hold(hca,'on')
  irf_plot(hca,tsLpeak,'k')
  hold(hca,'off')
  hca.YLabel.String = {'L','(km)'}; 
end
if 1 % L specrec
  isub = isub + 1;
  hca = irf_panel('L specrec var B0');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,specrecL_varB0,''log'');',2)  
  set(hca,'yscale','lin');
  hcb.YLabel.String = 'log_{10} counts';
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);    
  colormap(hca,cmap) 
  hca.YLabel.Interpreter = 'tex';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  irf_legend(hca,{'B_0 = 2\mu_0P_{tot}^{1/2}'},[0.98 0.98],'color','k','fontsize',fontsize_leg)
  %colormap(hca,irf_colormap('waterfall'))
  colormap(hca,pic_colors('candy4'))
  hold(hca,'on')
  irf_plot(hca,tsLpeak_varB0,'k')
  hold(hca,'off')
  hca.YLabel.String = {'L','(km)'};
end
if 0 % compare L with const and varying B0
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('L const and varying');
  set(hca,'ColorOrder',mms_colors('21'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{tsLpeak,tsLpeak_varB0},'comp');
  hca.YLabel.String = {'L','(km)'};
  set(hca,'ColorOrder',mms_colors('21'))
  irf_legend(hca,{sprintf('B_0 = %.0f nT',B0*1e9),'B_0 = 2\mu_0P_{tot}^{1/2}'}',[0.98 0.98],'fontsize',fontsize_leg);
  %irf_legend(hca,{sprintf('resampled to %g s (for noise reduction)',dt_L)},[0.02 0.98],'color','k')
end 

irf_zoom(h1,'x',tint)
irf_zoom(h1(zoomy),'y')
drawnow
h1(end).YLabel.Position(1) = -0.10;
drawnow
irf_plot_axis_align
%hmark = irf_pl_mark(h1(1),tint_harris);
c_eval('h1(?).FontSize = fontsize;',1:numel(h1))
c_eval('h1(?).YLabel.Position(1) = -0.12;',1:numel(h1))

% Non-TSeries panels.
isub = 1;
if 0 % Bx,Jy
  hca = h2(isub); isub = isub + 1;
  plot(hca,abs(gseBav.tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'.')  
end
clear hl
for species = {'','i','e'}
  if 0 % Jy,Jcurl
    hca = h2(isub); isub = isub + 1;
    eval(sprintf('xx = gseJ%sav.tlim(tint).y.data;',species{1}))
    yy = gseJcurl.resample(gseJav).tlim(tint).y.data;
    irem = find(xx == 0);
    xx(irem) = [];
    yy(irem) = [];
    p = polyfit(xx,yy,1);
    J_edges1 = linspace(-10,10,41);
    J_edges2 = linspace(-10,10,39);
    [N edges mid loc] = histcn([xx(:) yy(:)],J_edges1,J_edges2);
    N(N==0) = NaN;
    pcolor(hca,mid{:},log10(N)')
    shading(hca,'flat')
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = 'log_{10}(counts)';
    colormap(pic_colors('candy4'))
    %plot(hca,xx,yy,'.')
    %histcn_plot(hca,xx,yy)
    hca.XLabel.String = sprintf('<J_{%sy}^{FPI}>',species{1});
    hca.YLabel.String = 'J_y^{curl}';
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hold(hca,'on')
    xlim = hca.XLim;
    ylim = hca.YLim;
    plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
    plot(hca,[-100 100],p(2)+[-100 100]*p(1))
    hca.XLim = xlim;
    hca.YLim = ylim;
    hold(hca,'off')
    %hca.XLim = 40*[-1 1];
    %hca.YLim = 40*[-1 1]; 
    hca.XGrid = 'on';
    hca.YGrid = 'on'; 
  end
end
if 0 % Jey,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJiav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJiav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{i,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 0 % Jiy,Jicurl
  hca = h2(isub); isub = isub + 1;
  xx = gseJeav.tlim(tint).y.data;
  yy = gseJcurl.resample(gseJeav).tlim(tint).y.data;
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = '<J_{e,y}^{FPI}>';
  hca.YLabel.String = 'J_y^{curl}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  hca.XLim = 40*[-1 1];
  hca.YLim = 40*[-1 1];  
end
if 0 % Bx,Jiy, binning with B^2
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_L).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_L).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:).^2 yy(:)],B2_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1},mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'counts';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';  
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
  %colormap(hca,irf_colormap('waterfall'))
  colormap(pic_colors('candy4'))
end
if 1 % Bx,Jiy, binning with B^2
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = 'J_{y}^{curl}';
  yy = ts_yy.tlim(tint_L).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_L).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:).^2 yy(:)],B2_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1},mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  %hcb.YLabel.String = 'counts';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  clim = hca.CLim;
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  % mf_L(B,B0,Jy)
  J_edges_pos = J_edges(J_edges>0);
  B_grid = linspace(B2_edges(1),B2_edges(end),1000);
  J_grid = linspace(0,J_edges_pos(end),1001);
  %B_grid = logspace(log10(B2_edges(1)),log10(B2_edges(end)),500);
  %J_grid = logspace(log10(0.0000001),log10(J_edges_pos(end)),501);
  
  if 1 % B0 as defined above
  [BB2,JJ] = meshgrid(B_grid,J_grid);
  LL = mf_L(sqrt(BB2),B0*1e9,JJ);
  [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:1000:6000],'k','linewidth',1.5);
  % [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:500:6000],'k',"ShowText",true,'LabelFormat','L = %g km');
  clabel(cc,hh,'LabelSpacing',400,'Color','k','FontWeight','bold','fontsize',fontsize);
  end
  B02 = 25;
  [BB2,JJ] = meshgrid(B_grid,J_grid);
  LL = mf_L(sqrt(BB2),B02,JJ);
  [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:1000:6000],'r','linewidth',1.5);
  % [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:500:6000],'k',"ShowText",true,'LabelFormat','L = %g km');
  clabel(cc,hh,'LabelSpacing',600,'Color','r','FontWeight','bold','fontsize',fontsize);
  
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{sprintf('B_0 = %g nT',B0*1e9),sprintf('B_0 = %g nT',B02)}',[0.98 0.98],'fontsize',fontsize_leg)
  
%  hca_pos = hca.Position;
%  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
%  hlegs.Title.String = 'L = ';  
%  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.CLim = clim;
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
  colormap(hca,irf_colormap('waterfall'))
  %colormap(pic_colors('candy4'))
  hca.Layer = 'top';
end

if 0 % Bx,Jiy
  hca = h2(isub); isub = isub + 1;
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  %hcb.YLabel.String = 'log_{10}(counts)';
  hcb.YLabel.String = 'counts';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
   
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  
  hca_pos = hca.Position;
  hlegs = legend(hl,legs,'fontsize',fontsize,'location','northoutside','box','off','orientation','horizontal');
  hlegs.Title.String = 'L = ';  
  hca.Position = hca_pos;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hlegs.Title.String = 'L = ';
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);  
  colormap(hca,irf_colormap('waterfall'))
end
if 0 % Bx,Jiy, time
  hca = h2(isub); isub = isub + 1;
  timeline = tint_harris(1):0.02:tint_harris(2);
  ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  tt = ts_yy.resample(ts_yy).tlim(tint_harris).time - ts_yy.resample(ts_yy).tlim(tint_harris).time(1);
  start_time = irf_time(ts_yy.resample(ts_yy).tlim(tint_harris).time(1),'EpochTT>utc_HH:MM:SS');
  scatter(hca,xx.^2,yy,1,tt)  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = sprintf('seconds since %s',start_time);
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  
  % Add L lines
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = B0
    for L_ = L  
      ileg = ileg + 1;
      %legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      legs{ileg} = sprintf('L = %g km',L_*1e-3);
      legs{ileg} = sprintf('%g km',L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_Jy(B0_,L_,z_,z0_)*1e9*(1+0*1/5),'k:');
    end
  end
  %hca_pos = hca.Position;
  %hlegs = legend(hl,legs,'fontsize',10,'location','northoutside','box','off','orientation','horizontal');
  %hlegs.Title.String = 'L = ';  
  %hca.Position = hca_pos;
  %hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  %hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  %hlegs.Box = 'on';
  hold(hca,'off')
  
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
  
  hca.XLim = B_edges([1 end]).^2;
  hca.YLim = J_edges([1 end]);
  colormap(hca,irf_colormap('waterfall'))
  
end
if 0 % Bx,JxBz, binned
  hca = h2(isub); isub = isub + 1;
  %plot(hca,abs(gseBav.resample(gseJi1).tlim(tint_harris).x.data),gseJi1.tlim(tint_harris).y.data,'.')
  %plot(hca,abs(gseBav.resample(gseJcurl).tlim(tint_harris).x.data),gseJcurl.tlim(tint_harris).y.data,'k.')
  %ts_yy = gseJiav; y_str = '<J_{iy}^{FPI}>';
  %ts_yy = gseJcurl.resample(timeline); y_str = '<J_{y}^{curl}>';
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  
  %irem = find(xx == 0);
  %xx(irem) = [];
  %yy(irem) = [];
  %p = polyfit(xx,yy,1);
  JxB_edges = linspace(-2,12,70);
  B_edges = linspace(0,24,71);
  [N edges mid loc] = histcn([xx(:) yy(:)],B_edges,J_edges);
  N(N==0) = NaN;
  pcolor(hca,mid{1}.^2,mid{2},log10(N)')
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'log_{10}(counts)';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end
if 0 % Bx,JxBz, scatter
  hca = h2(isub); isub = isub + 1;
  
  ts_yy = gseJxB; y_str = '-J_{y}xB_x';
  %ts_yy = gseJeav; y_str = '<J_{ey}^{FPI}>';
  %ts_yy = gseJav.resample(timeline); y_str = '<J_{y}^{FPI}>';
  yy = ts_yy.tlim(tint_harris).y.data*1e18;  
  xx = gseBav.resample(ts_yy).tlim(tint_harris).x.data;  
  cc = ts_yy.tlim(tint_harris).time - ts_yy.tlim(tint_harris).time(1);
  scatter(hca,xx.^2,yy,1,cc)
  shading(hca,'flat')    
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time since tt';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
    
  if 0
  hold(hca,'on')
  % mf_Bx = @(B0,L,z,z0)B0.*tanh((z-z0)./L)
  % mf_Jy = @(B0,L,z,z0)-(B0.*(tanh((z-z0)./L).^2-1.0))./L
  legs = {};
  ileg = 0;
  z0_ = 00e3;
  for B0_ = [21]*1e-9
    for L_ = [1500 2000 2500 4000]*1e3  
      ileg = ileg + 1;
      legs{ileg} = sprintf('B_0 = %g nT, L = %g km',B0_*1e9,L_*1e-3);
      z_ = L_*linspace(0.0,3,20);      
      hl(ileg) = plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'linewidth',1);
      plot(hca,(mf_Bx(B0_,L_,z_,z0_)*1e9).^2,mf_JxB(B0_,L_,z_,z0_)*1e9*(1+1/5),'k:');
    end
  end
  hlegs = legend(hl,legs,'fontsize',10,'location','northeast','box','off');
  %hca_pos = hca.Position;
  %hlegs.Location = 'northoutside';
  %hca.Position = hca_pos;
  hlegs.Position(2) = hca.Position(2) + hca.Position(4);
  hlegs.Position = [0.8521    0.3418    0.1180    0.1044];
  hlegs.Box = 'on';
  hold(hca,'off')
  end
  hca.XLabel.String = 'B_x^2 (nT^2)';
  hca.YLabel.String = sprintf('%s (nA/m^2)',y_str);
end

%h2(2).YTick = h2(1).YTick;
%h2(2).XTick = h2(1).XTick;

h2(1).Position = [0.58 0.25 0.2 0.5];
c_eval('h2(?).FontSize = fontsize;',1:numel(h2))
c_eval('h2(?).Position(3) = 0.3; h2(?).Position(4) = 0.6;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

leg_str = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
for ih = 1:numel(h1)
  irf_legend(h1(ih),leg_str{1},[0.01 0.98],'color','k','fontsize',fontsize_leg);
  leg_str = leg_str(2:end);  
end
for ih = 1:numel(h2)
  irf_legend(h2(ih),leg_str{1},[0.02 0.98],'color','k','fontsize',fontsize_leg);
  leg_str = leg_str(2:end);
end

if 1 % add time marking on top of panel a to show were the data comes from
    %%
    drawnow
    % Clone colorbar for time and place on top of panel 1.
    hca = h1(1);
    h_pos = hca.Position;
    hb = colorbar('peer',hca,'location','northoutside');
    hb.YTick = [];
    colormap(hb,irf_colormap('waterfall'))
    colormap(hb,[0 0 0]+0.8)
    
    hmark = irf_pl_mark(hca,tint_L.epochUnix','k');


    xlim = hca.XLim;
    xmark = [min(hmark.XData) max(hmark.XData)];
    x1_rel = xmark(1)/diff(xlim);
    x2_rel = xmark(2)/diff(xlim);

    hb.Position(1) = hca.Position(1) + hca.Position(3)*x1_rel;
    hb.Position(3) = hca.Position(3)*(x2_rel-x1_rel);
    drawnow
    hb.Position(2) = hca.Position(2) + hca.Position(4);
    delete(hmark)
end

%
h = findobj(gcf,'type','axes'); h = h(end:-1:1);
  c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

c_eval('h(?).LineWidth = 1;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('ht(?).FontWeight = ''bold'';',1:numel(ht))
c_eval('h(?).FontWeight = ''bold'';',1:numel(h))

%% EGU 2023: Figure: onset 
% for feeps_ion_omni1234, see feeps_data.m
ic = 1;
time_mark = EpochTT('2017-07-25T22:09:08.00Z');
time_df = EpochTT('2017-07-25T22:10:06.50Z');

tint_figure = irf.tint('2017-07-25T22:07:00.00Z/2017-07-25T22:11:30.00Z');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

fontsize_leg = 14;

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x  ','B_y  ','B_z'},[0.02 0.12],'fontsize',12);
end 
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',12);  
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-3:1:2];
  hca.YLabel.Interpreter = 'tex';
end
if 0% Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,L}','v_{i\perp,M}','v_{i\perp,N}','v_{i||}'}',[1.02 0.9],'fontsize',12);
end

elim_feeps = [8e4 Inf];
if 1 % FEEPS Pitch all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta^{FEEPS}_i','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  hca.YTick = 0:60:180;
end
if 1 % i DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps omni mms 1234');  
  [hout,hcb] = irf_spectrogram(hca,feeps_ion_omni1234.elim(elim_feeps).specrec('energy'),'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 1 % i DEF omni, FPI
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i^{FPI}','(eV)'};   
end
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','M','N'];
  comps = ['x','y','z'];
  for icomp = 1:3 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  eval(sprintf('irf_plot(hca,{%sJcurl.x,%sJcurl.y,%sJcurl.z},''comp'');',cs,cs,cs))
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
if 0 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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

if 0 % EIS Pitch
  isub = isub + 1;
  hca = irf_panel('eis Pitch');  
  irf_spectrogram(hca,eis_pa2.elim([14 4500]*1e3).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 1 % i psd x,y,z, 3 panels
  for comp = ['x']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('xyza'))  
  hca.YLim = [0 5];
end

if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{ix}  ','v_{iy}  ','v_{iz}'},[0.05 0.9],'fontsize',12);
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);


hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca1 = irf_panel('feeps omni mms 1234');  
hca2 = irf_panel('i DEF omni'); 
hca1.YLim(1) = hca2.YLim(2);
h(end).YLim = [-350 1150];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [0 4]*0.99;


c_eval('h(?).FontSize = fontsize_leg;',1:numel(h))


drawnow
%h(6).YLim(1) = h(7).YLim(2);
%h(2).YLim(1) = 0.001;

drawnow
hmark = irf_pl_mark(h,time_mark,'k');
c_eval('hmark(?).LineStyle = '':'';',1:numel(hmark))
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))


hmark2 = irf_pl_mark(h,time_df,'k');
c_eval('hmark2(?).LineStyle = '':'';',1:numel(hmark2))
c_eval('hmark2(?).LineWidth = 1;',1:numel(hmark2))

%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
%hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
hat=annotation('textarrow',[0 0]+0.4820,[.56 .54]+0.08,'string',{'first energetic ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none','fontweight','bold');
%hat=annotation('textarrow',[0.593 0.593]+0.104,[.4 .38]+0.01,'string',{'first cold ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none');
hat=annotation('textarrow',[0 0]+0.628,[.4 .38]+0.01,'string',{'first cold ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none','fontweight','bold');
%

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
  c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

c_eval('h(?).LineWidth = 1;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('ht(?).FontWeight = ''bold'';',1:numel(ht))
c_eval('h(?).FontWeight = ''bold'';',1:numel(h))

c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))

%% MMS SWT: detailed look at cold ions
% Plot 2D ion distributions to show off cold ions
ic = 1;
[h,h2] = initialize_combined_plot(6,2,1,0.55,'vertical');
c_eval('h(?).Position(1) = 0.10;',1:numel(h))
%h = irf_plot(10);

%tint_figure = tint_cavity + [-135 15];


tint_figure = tint_action;
tint_figure = irf.tint('2017-07-25T22:09:50.00Z/2017-07-25T22:10:20.00Z');

times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.00Z';...
                      '2017-07-25T22:10:08.00Z';...
                      '2017-07-25T22:10:09.20Z']);

fontsize_leg = 14;
                    
cs = 'gse';
cs_str = 'xyz';
elim = [600 40000];

comps = ['x','y','z'];
comp_str = 'LMN';
comps = ['x'];
comp_str = 'x';

if 0
cs = 'gse';
cs_str = 'xyz';

comps = ['x','y','z'];
comp_str = 'xyz';
end

zoomy = []; isub = 1;
if 1 % B
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z,%sB?.abs},''comp'');',cs,cs,cs,cs),ic)    
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x ','B_y ','B_z ','|B|'},[0.98 0.12],'fontsize',fontsize_leg);  
end 
if 1 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98],'fontsize',fontsize_leg)
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 1 % Vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{ix} ','v_{iy} ','v_{iz}'},[0.02 0.98],'fontsize',fontsize_leg);
end

if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.tlim(tint_figure).omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 1 % elim
    %%
    hold(hca,'on')
    c_eval('timeline = iPDist?.time;',ic)    
    irf_plot(hca,irf.ts_scalar(timeline,700*ones([timeline.length,1])),'k')
    hold(hca,'off')
  end
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YLabel.Interpreter = 'tex';  
  hca.XGrid = 'off';
  hca.YGrid = 'off'; 
end

if 1 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  c_eval('if1D = if1Dx?_700.tlim(tint_figure);',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);  
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.x,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{ix}','(km/s)'}; 
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?_700.tlim(tint_figure);',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 1 % i psd z
  isub = isub + 1;
  hca = irf_panel('iLine z');
  c_eval('if1D = if1Dz?_700.tlim(tint_figure);',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
end


hlinks = linkprop([irf_panel('iLine x'),irf_panel('iLine y'),irf_panel('iLine z')],{'CLim'});
colormap(pic_colors('candy4'))
drawnow
irf_zoom(h,'x',tint_figure)
drawnow
%irf_zoom(h(zoomy),'y')
irf_plot_axis_align(h)
%h(1).Title.String = sprintf('L = [%5.2f,%5.2f,%5.2f], M = [%5.2f,%5.2f,%5.2f], N = [%5.2f,%5.2f,%5.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));
drawnow

% Distributions
if exist('hmark'); delete(hmark); end
if exist('hmark2'); delete(hmark); end
times_dist = EpochTT(['2017-07-09T17:33:17.00Z';...
                      '2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:21.00Z';...
                      '2017-07-09T17:33:22.00Z']);
                    

times_dist = EpochTT(['2017-07-09T17:33:18.00Z';...
                      '2017-07-09T17:33:18.50Z';...
                      '2017-07-09T17:33:19.00Z';...
                      '2017-07-09T17:33:19.50Z';...
                      '2017-07-09T17:33:20.00Z';...
                      '2017-07-09T17:33:20.50Z']);

times_dist = EpochTT(['2017-07-09T17:33:19.00Z']) + 3*0.03*[-2:1:3];




times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:09.00Z';...
                      '2017-07-25T22:10:15.00Z';...
                      %'2017-07-25T22:10:29.00Z';...
                      %'2017-07-25T22:10:49.00Z';...
                      %'2017-07-25T22:10:27.00Z';...
                      '2017-07-25T22:10:29.00Z';...
                    ]);
                    
                  
times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.00Z';...
                      '2017-07-25T22:10:08.00Z';...
                      '2017-07-25T22:10:09.20Z']);
                    
                    
times_dist = EpochTT(['2017-07-25T22:10:00.00Z';...
                      '2017-07-25T22:10:07.50Z']);
                    
times_dist = EpochTT(['2017-07-25T22:10:03.00Z';...
                      '2017-07-25T22:10:07.50Z']);
                    
times_dist_exact = EpochTT([]);
clear times_dist_exact_eu hmarkb hmark_tmp times_dist_exact_eu_edges

%elim = [700 Inf];

for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  v1 = [1 0 0];
  v2 = [0 1 0];
  v3 = [0 0 1];
  v1 = R(1,:);
  v2 = R(2,:);
  v3 = R(3,:);
  x_str = 'v_L (km/s)';
  y_str = 'v_M (km/s)';
  z_str = 'v_N (km/s)';
  % bugcheck, slightly rotate system
%   rotang = 30;
%   v2n = v2*cosd(rotang) + v3*sind(rotang);
%   v3n = -v2*sind(rotang) + v3*cosd(rotang);
%   v2 = v2n; v3 = v3n;
  vlim = 2300;
  vg = -vlim:100:vlim;
  if 1 % x,z
    
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+6*0.15*0.5*[-1 1]+0);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = mean(ddist.time.epochUnix);
    times_dist_exact_eu_edges(id,1) = ddist.time(1).epochUnix-0.03*0.5;
    times_dist_exact_eu_edges(id,2) = ddist.time(end).epochUnix+0.03*0.5;
    f2d = ddist.elim([00 Inf]).reduce('2D',v1,v3);

    f2d.plot_plane(hca,'smooth',0);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = x_str;
    hca.YLabel.String = z_str;
    
    hmark_tmp = irf_pl_mark(h,ddist.time,'k','facealpha',0.1);
    drawnow
    %c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))

  end
  if 0 % x,y
    c_eval('ddist = iPDist?.elim(elim).tlim(times_dist(id)+0.15*0.5*[-1 1]);',ic)
%    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vg',vg);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_x (km/s)';
    hca.YLabel.String = 'v_y (km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  colormap(hca,pic_colors('pasteljet'))
end
if 0
for id = 1:times_dist.length
  hca = h2(id);
  dt = 0.05;
  c_eval('v1 = mean(dmpaB?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  c_eval('v2 = mean(dslE?.tlim(times_dist(id)+dt*[-1 1]).norm.data,1);',ic)
  v2 = cross(v1,cross(v2,v1)); v2 = v2/norm(v2);
  v3 = cross(v1,v2);
  vlim = 70;
  if 0 % vE,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v2,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_E (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
  if 0 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v2,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{E} (10^3 km/s)';
  end
  if 1 % vB,vExB
    c_eval('ddist = ePDist?.tlim(times_dist(id)+0.03*0.5*[-1 1]);',ic)
    ddist = ddist(1).elim([100 Inf]);
    times_dist_exact_eu(id) = ddist.time.epochUnix;
    f2d = ddist.reduce('2D',v1,v3,'vint',20e3*[-1 1]);

    f2d.plot_plane(hca);
    axis(hca,'square')
    hca.XLim = 0.99*vlim*[-1 1];
    hca.YLim = 0.99*vlim*[-1 1];
    hca.Layer = 'top';
    hca.XLabel.String = 'v_B (10^3 km/s)';
    hca.YLabel.String = 'v_{ExB} (10^3 km/s)';
  end
end
end

hlinks = linkprop(h2,{'CLim','XLim','YLim'});
h2(1).CLim = [-10 -7.];

c_eval('h(?).FontSize = 16;',1:numel(h))
c_eval('h2(?).FontSize = 16;',1:numel(h2))
c_eval('h2(?).Position(1) = 0.58;',1:numel(h2))

h(end).XTickLabelRotation = 0;
h2(1).Title.String = 'Hot pre-existing plasma sheet';
h2(2).Title.String = 'Post-DF cold and heated ions';
%c_eval('h2(?).CLim = [-14.5 -9.5];',1:numel(h2))
%hmark = irf_pl_mark(h,tocolumn(times_dist_exact_eu),'r');
%c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))
hmark2 = irf_pl_mark(h,tocolumn(times_dist_exact_eu_edges(:)),'r');
c_eval('hmark2(?).LineWidth = 0.5; hmark2(?).Color = [0 0 0]; hmark2(?).LineStyle = ''--'';',1:numel(hmark2))

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
  c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

c_eval('h(?).LineWidth = 1;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('ht(?).FontWeight = ''bold'';',1:numel(ht))
c_eval('h(?).FontWeight = ''bold'';',1:numel(h))

%% MMS SWT: Figure: off-time current sheet structure, vi vExB
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 4;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

fontsize_leg = 14;

isub = 0;
zoomy = [];
cs = 'lmn';
if 1 % B cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval(sprintf('irf_plot(hca,{%sB?.x,%sB?.y,%sB?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {sprintf('B_{%s}',cs),'(nT)'};
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'B_x ','B_y ','B_z'},[0.02 0.12],'fontsize',fontsize_leg);
end 
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i,\perp,||}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i\perp,x}  ','v_{i\perp,y}  ','v_{i\perp,z}  ','v_{i||}'},[0.02 0.12],'fontsize',fontsize_leg);
end
if 0 % ne
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e'},[0.98 0.98])
  hca.YLabel.String = {'n','(cm^{-3})'};
end
if 0 % Pressures, PB, Pi, Pe
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize_leg);  
end
if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,1/0.3},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-3:1:2];
  hca.YLabel.Interpreter = 'tex';
end
if 0 % E cs
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x,%sE?.y,%sE?.z},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E cs, resamp vi
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sE?.x.resample(gseVi?),%sE?.y.resample(gseVi?),%sE?.z.resample(gseVi?)},''comp'');',cs,cs,cs),ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',fontsize_leg);
  irf_legend(hca,{'E_x ','E_y ','E_z'},[0.02 0.12],'fontsize',fontsize_leg);
  irf_zoom(hca,'y')
end

elim_feeps = [8e4 Inf];
if 0 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    ve.data(abs(ve.data)>5000) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve,%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'v_e','v_i','ExB'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
  end
end
if 1 % vExB.x/y/z , Vi resample, 3 panels
  strcomp = ['L','y','N'];
  comps = ['x','y','z'];
  for icomp = 2 % vExB.x , Vi resample, ve
    comp = comps(icomp);
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',mms_colors('21'))
    c_eval(sprintf('irf_plot(hca,{%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs),ic)  
    %c_eval(sprintf('irf_plot(hca,{%sVe?perp.resample(%sVi?).(comp),%sVExB?.resample(%sVi?).(comp),%sVi?perp.(comp)},''comp'');',cs,cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{i\\perp,(%s)}^{%s}',comp,cs),'(km/s)'};
    hca.YLabel.String = {sprintf('v_{i\\perp,%s}',strcomp(icomp)),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('21'))
    %irf_legend(hca,{'v_{ExB}','v_i'},[0.98 0.9],'fontsize',12);
    irf_legend(hca,{'v_{ExB}','v_i'}',[0.07 0.9],'fontsize',fontsize_leg);    
    %irf_legend(hca,{'v_e','v_{ExB}','v_i'}',[1.02 0.9],'fontsize',12);    
  end
end

if 0 % E+vxB.x/y/z , Vi resample, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve
    isub = isub + 1;    
    hca = irf_panel(['E + VxB ',comp]);
    set(hca,'ColorOrder',mms_colors('123'))
    c_eval(sprintf('evexb = %sE?perp.(comp).resample(%sVi?)+gseVexB?.(comp).resample(%sVi?);',cs,cs,cs),ic)
    c_eval(sprintf('irf_plot(hca,{evexb,%sE?perp.(comp).resample(%sVi?)+gseVixB?.(comp)},''comp'');',cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('E+v_{(%s)}\times B',comp,cs),'(km/s)'};
    set(hca,'ColorOrder',mms_colors('123'))
    irf_legend(hca,{'e','i'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'x'},[0.98 0.9],'fontsize',12);
    hca.YLim = [-5 5];
  end
end
if 0 % J curl
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyza'))
  eval(sprintf('irf_plot(hca,{%sJcurl.x,%sJcurl.y,%sJcurl.z},''comp'');',cs,cs,cs))
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
if 0 % E gsm lowpass filt
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E gsm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  ffilt = 1;
  c_eval('irf_plot(hca,{gsmE?.filt(0,1,[],3).x,gsmE?.filt(0,1,[],3).y,gsmE?.filt(0,1,[],3).z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{sprintf('f < %.0f',ffilt)},[0.1 0.9],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'y')
end
if 0 % Vi gsm
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gsmVi?.x,gsmVi?.y,gsmVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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

if 0 % EIS Pitch
  isub = isub + 1;
  hca = irf_panel('eis Pitch');  
  irf_spectrogram(hca,eis_pa2.elim([14 4500]*1e3).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,eis_omni?.specrec(''energy''),''log'');',2)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  if exist('Etop_fpi','var')
    hold(hca,'on')
    c_eval('tint_tmp = eis_omni?.time;',2)
    irf_plot(hca,irf.ts_scalar(tint_tmp([1 end]),Etop_fpi*[1 1]),'k')
    hold(hca,'off')
  end
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{EIS}','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF feeps omni
  isub = isub + 1;
  hca = irf_panel('i DEF feeps omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,feeps_ion_omni?.specrec(''energy''),''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_i^{FEEPS}','(eV)'};   
end
if 0 % i DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_i','(eV)'};   
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
if 0 % Ti par perp
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

if 0 % i psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?;',comp),ic)
    irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
    hca.YLim = if1D.depend{1}(1,[1 end]);  
    if 1 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{i%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end
if 0 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % scPot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end

if 0 % i psd y
  isub = isub + 1;
  hca = irf_panel('iLine y');
  c_eval('if1D = if1Dy?;',ic)
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 0 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.y,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iy}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
  hca.CLim = [-5 -1];
end

if 0 % Tepar/Teperp
  isub = isub + 1;  
  hca = irf_panel('Tepar/Teperp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par/Te?perp},''comp'');',ic)  
  hca.YLabel.String = 'T_{e||}/T_{e\perp}';
  set(hca,'ColorOrder',mms_colors('xyza'))  
  hca.YLim = [0 5];
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = 14;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);


hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [0 4]*0.99;

%hca = irf_panel(['V ExB Vi ','x']); hca.YLim = [-100 250];
hca = irf_panel(['V ExB Vi ','y']); hca.YLim = [-250 699];
%hca = irf_panel(['V ExB Vi ','z']); hca.YLim = [-300 300];
hca = irf_panel('E gsm'); hca.YLim = [-10 10]*0.99;

c_eval('h(?).FontSize = 16;',1:numel(h))

% text: 'Moderate beta...'
if 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]-dy);
elseif 0
  annotation('textarrow',[0.3 0.3],[.7 .72],'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]);
  annotation('textarrow',[0.63 0.63],[.7 .72],'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63],[.68 .64]);  
elseif 0
  %%
  dy = 0.08;
  dx = 0.07; 
  annotation('textarrow',[0.3 0.3],[.7 .72]-dy,'string',{'moderate \beta with v_{i||}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.68 .64]-dy);
  annotation('textarrow',[0.63 0.63]-dx,[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.63 0.63]-dx,[.68 .64]-dy);
end

drawnow
%h(6).YLim(1) = h(7).YLim(2);
%h(2).YLim(1) = 0.001;

drawnow
%hmark = irf_pl_mark(h,time_mark,'k');
%c_eval('hmark(?).LineStyle = '':'';',1:5)
%c_eval('hmark(?).LineWidth = 1;',1:5)

%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
%hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);


h = findobj(gcf,'type','axes'); h = h(end:-1:1);
  c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

c_eval('h(?).LineWidth = 1;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

ht = findall(gcf,'type','text');
c_eval('ht(?).FontSize = fontsize;',1:numel(ht))
c_eval('ht(?).FontWeight = ''bold'';',1:numel(ht))
c_eval('h(?).FontWeight = ''bold'';',1:numel(h))
