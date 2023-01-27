tint_harris = irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:09:30.00Z');

%% Figure: Overview 1
ic = 1;

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));

fontsize_leg = 13;
fontsize = 14;
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
  %zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
  hca.YLim = [-17 24];
end 
if 1 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  colors_plot = [0.7 0.7 0.7; mms_colors('1')];
  colors_leg = [0.6 0.6 0.6; mms_colors('1')];
  set(hca,'ColorOrder',colors_plot)
  c_eval('irf_plot(hca,{ne?,nHp?_srvy},''comp'')',ic)
  set(hca,'ColorOrder',colors_leg)
  irf_legend(hca,{'FPI: n_e','HPCA: n_{H^+}'}',[1.02 0.9],'fontsize',fontsize_leg)
  hca.YLabel.String = {'n','(cm^{-3})'};
  hca.YLim = [0 0.32*0.99];
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
 % zoomy = [zoomy isub];
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',fontsize_leg);
  irf_legend(hca,{'v_x','v_y','v_z'}',[1.02 0.9],'fontsize',fontsize_leg);
  irf_legend(hca,{'FPI'},[0.04 0.9],'color','k','fontsize',fontsize_leg);
  hca.YLim = 0.99*[-600 1650];
end
if 1 % Vi gse hpca
  isub = isub + 1;
  %zoomy = [zoomy isub];
  hca = irf_panel('Vi gse hpca');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVHp?_srvy.x,gseVHp?_srvy.y,gseVHp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'v_x','v_y','v_z'}',[1.02 0.9],'fontsize',fontsize_leg);
  %irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',fontsize_leg);
  irf_legend(hca,{'HPCA'},[0.04 0.9],'color','k','fontsize',fontsize_leg);
  ytick = h(isub-1).YTick;
  hca.YTick = -3000:(ytick(2)-ytick(1)):3000;
  hca.YLim = 0.99*[-600 1650];
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
  
  if 0 % elim
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
  c_eval('if1D = if1Dx?_700;',ic)
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
  c_eval('if1D = if1Dy?_700;',ic)
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
  c_eval('if1D = if1Dz?_700;',ic)
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
if 1 % e DEF omni
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

if 1 % Pressures, PB, Pi, Pe
  isub = isub + 1;  
  hca = irf_panel('Pressure');  
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval(['irf_plot(hca,{PB?,'...
                         'gsePi?.trace/3,'...
                         'gsePe?.trace/3,'...
                         'gsePe?.resample(gsePi?).trace/3+PB?.resample(gsePi?)+gsePi?.trace/3},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('matlab'))
  irf_legend(hca,{'P_B','P_i','P_e','P_{tot}'}',[1.02 0.9],'fontsize',fontsize_leg);  
  hca.YLim(1) = 0;
  hca.YLim(2) = 0.5;
end
if 0 % Te par perp Ti/Tref
  isub = isub + 1;
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 1;
  c_eval('irf_plot(hca,{facTe?.resample(facTi?).trace/3,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_e',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);  
  hca.YLim = [0 10000];  
  %irf_zoom(hca,'y')
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);


h(end).XTickLabelRotation = 0;

%% Figure: Overview 1&2 together, full page
ic = 1;

time_mark = EpochTT('2017-07-25T22:09:08.00Z');
time_df = EpochTT('2017-07-25T22:10:06.00Z');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
feeps_ion_omni1234 = feeps_ion_omni1;
%feeps_ion_omni1234.data = (feeps_ion_omni1.data + feeps_ion_omni2.resample(feeps_ion_omni1).data + feeps_ion_omni3.resample(feeps_ion_omni1).data + feeps_ion_omni4.resample(feeps_ion_omni1).data)/4;
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

fontsize = 14;
fontsize_leg = 13;

npanels = 10;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'gse';
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
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
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
if 1 % Vi xyz
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi');
  colors_plot = [mms_colors('xyza'); mms_colors('x').^0.5];
  set(hca,'ColorOrder',colors_plot)
  c_eval(sprintf('irf_plot(hca,{%sVi?.x,%sVi?.y,%sVi?.z},''comp'');',cs,cs,cs),ic)  
  c_eval('irf_plot(hca,{gseVHp?_srvy.x,gseVHp?_srvy.y,gseVHp?_srvy.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',colors_plot)
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{ix}','v_{iy}','v_{iz}'}',[1.02 0.9],'fontsize',fontsize_leg);

  irf_legend(hca,{'FPI'},[.17 0.77],'fontsize',fontsize_leg,'color','k');
  irf_legend(hca,{'HPCA'},[.03 0.8],'fontsize',fontsize_leg,'color','k');
end
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  colors_plot = [mms_colors('xyza'); mms_colors('x').^0.5];
  set(hca,'ColorOrder',colors_plot)
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  c_eval('irf_plot(hca,{gseVHp?perp_srvy.x,gseVHp?perp_srvy.y,gseVHp?perp_srvy.z,gseVHp?par_srvy},''comp'');',ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',colors_plot)
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i{\perp}x}','v_{i{\perp}y}','v_{i{\perp}z}','v_{i||}'}',[1.02 0.99],'fontsize',fontsize_leg);
  irf_legend(hca,{'FPI'},[.17 0.77],'fontsize',fontsize_leg,'color','k');
  irf_legend(hca,{'HPCA'},[.03 0.8],'fontsize',fontsize_leg,'color','k');
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

if 1 % i psd x,y,z, 3 panels
  for comp = ['x']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?_700;',comp),ic)
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
    %hca.XGrid = 'on';
    %hca.YGrid = 'on';
    %hca.Layer = 'top';
    %irf_legend(hca,sprintf('E > %.0f eV',if1D.ancillary.energy(1,1)-if1D.ancillary.delta_energy_minus(1,1)),[0.02 0.08],'color','k','fontsize',fontsize_leg)
  end
end

if 1 % beta
  hca = irf_panel('beta');
  betalow = 1/0.3;
  betalow = 3;
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,betalow},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  %hp.FaceColor = [0.7 0.5 0.0];
  hp.EdgeColor = hp.FaceColor*0;
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1.01e-2 1e3];
  irf_legend(hca,['\beta < ' num2str(betalow,'%.2f')],[0.02 0.62],'color',hp.FaceColor,'fontsize',fontsize_leg)
end

elim_feeps = [8e4 Inf];
if 1 % FEEPS Pitch all sc
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
  irf_legend(hca,'four-spacecraft average',[0.04 0.95],'fontsize',fontsize_leg)
  hca.Color = [1 1 1]*0.9;
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
  irf_legend(hca,'four-spacecraft average',[0.04 0.95],'fontsize',fontsize_leg)
  hca.Color = [1 1 1]*0.9;
end
if 1% i DEF omni, FPI
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



if 0 % e DEF pitch
  isub = isub + 1;
  hca = irf_panel('i DEF pitch');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPitch?.elim([200 Inf]).deflux.specrec,''log'');',ic)  
  %set(hca,'yscale','log');
  set(hca,'ytick',0:30:180);  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};   
  colormap(hca,cmap) 
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
  tsdata = tsdata.tlim(irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:11:00.00Z'));
  if 1
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN; % doesn't work well with patch
    tsdata.data = log10(tsdata.data);
    c_eval('hp = irf_patch(hca,{tsdata,0});',ic)  
    hp.FaceColor = [0 0 0];
    hp.EdgeColor = [0 0 0];
  else
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN;
    tsdata.data = log10(tsdata.data);
    irf_plot(hca,tsdata,'k')
  end
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.201 0.201];
  hca.YLabel.Interpreter = 'tex';
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize_leg;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = fontsize;',1:numel(h))

linkprop([irf_panel('feeps omni mms 1234'),irf_panel('feeps pitch mms 1234')],{'CLim'})

h(end).XTickLabelRotation = 0;

%hca = irf_panel('feeps omni mms 1234'); hca.CLim = [-1 100];
%hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-16 25];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];

h(6).YLim(1) = h(7).YLim(2)*1e-3;
%hca1 = irf_panel('feeps omni mms 1234');  
%hca2 = irf_panel('i DEF omni'); 
%hca2.YLim(2) = hca1.YLim(1);
%h(end).YLim  = [-350 1150];

drawnow
hmark = irf_pl_mark(h,time_mark,'k');
c_eval('hmark(?).LineStyle = '':'';',1:numel(hmark))
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))


hmark2 = irf_pl_mark(h,time_df,'k');
c_eval('hmark2(?).LineStyle = '':'';',1:numel(hmark2))
c_eval('hmark2(?).LineWidth = 1;',1:numel(hmark2))

% text: 'Moderate beta...'
if 1
    %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.67 .66]-dy,'string',{'moderate \beta','with v_{ix}'},'fontsize',fontsize_leg,'horizontalalignment','center');
  %annotation('arrow',[0.3 0.3],[.69 .66]-dy);
  annotation('textarrow',[0.55 0.55],[.67 .66]-dy,'string',{'moderate \beta','without v_{ix}'},'fontsize',fontsize_leg,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  %annotation('arrow',[0.57 0.57],[.69 .66]-dy);
  annotation('textarrow',[0.56 0.6],[.415 .433],'string',{'field-aligned','energetic ions','appearing'},'fontsize',fontsize_leg,'horizontalalignment','right');
elseif 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.71 .73]-dy,'string',{'moderate \beta with v_{ix}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.69 .66]-dy);
  annotation('textarrow',[0.57 0.57],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.57 0.57],[.69 .66]-dy);
elseif 0
  %%
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


h(end).XTickLabelRotation = 0;


%% MMS SWT: Figure: ongoing reconnection
ic = 1;

time_mark = EpochTT('2017-07-25T22:09:08.00Z');
time_df = EpochTT('2017-07-25T22:10:06.00Z');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
feeps_ion_omni1234 = feeps_ion_omni1;
%feeps_ion_omni1234.data = (feeps_ion_omni1.data + feeps_ion_omni2.resample(feeps_ion_omni1).data + feeps_ion_omni3.resample(feeps_ion_omni1).data + feeps_ion_omni4.resample(feeps_ion_omni1).data)/4;
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

fontsize = 14;
fontsize_leg = 13;

npanels = 8;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'gse';
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
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
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
  colors_plot = [mms_colors('xyza'); mms_colors('x').^0.5];
  set(hca,'ColorOrder',colors_plot)
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par,gseVHp?perp_srvy.x},''comp'');',cs,cs,cs,cs),ic)  
  c_eval('irf_plot(hca,{gseVHp?perp_srvy.x,gseVHp?perp_srvy.y,gseVHp?perp_srvy.z,gseVHp?par_srvy},''comp'');',ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',colors_plot)
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i{\perp}x}','v_{i{\perp}y}','v_{i{\perp}z}','v_{i||}'}',[1.02 0.9],'fontsize',fontsize_leg);

  irf_legend(hca,{'FPI'},[.17 0.77],'fontsize',fontsize_leg,'color','k');
  irf_legend(hca,{'HPCA'},[.03 0.8],'fontsize',fontsize_leg,'color','k');
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

if 1 % i psd x,y,z, 3 panels
  for comp = ['x']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?_700;',comp),ic)
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
    %hca.XGrid = 'on';
    %hca.YGrid = 'on';
    %hca.Layer = 'top';
    %irf_legend(hca,sprintf('E > %.0f eV',if1D.ancillary.energy(1,1)-if1D.ancillary.delta_energy_minus(1,1)),[0.02 0.08],'color','k','fontsize',fontsize_leg)
  end
end

if 1 % beta
  hca = irf_panel('beta');
  betalow = 1/0.3;
  betalow = 3;
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,betalow},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  %hp.FaceColor = [0.7 0.5 0.0];
  hp.EdgeColor = hp.FaceColor*0;
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1.01e-2 1e3];
  irf_legend(hca,['\beta < ' num2str(betalow,'%.2f')],[0.02 0.62],'color',hp.FaceColor,'fontsize',fontsize_leg)
end

elim_feeps = [8e4 Inf];
if 1 % FEEPS Pitch all sc
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
  irf_legend(hca,'four-spacecraft average',[0.04 0.95],'fontsize',fontsize_leg)
  hca.Color = [1 1 1]*0.9;
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
  irf_legend(hca,'four-spacecraft average',[0.04 0.95],'fontsize',fontsize_leg)
  hca.Color = [1 1 1]*0.9;
end
if 1% i DEF omni, FPI
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



if 0 % e DEF pitch
  isub = isub + 1;
  hca = irf_panel('i DEF pitch');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPitch?.elim([200 Inf]).deflux.specrec,''log'');',ic)  
  %set(hca,'yscale','log');
  set(hca,'ytick',0:30:180);  
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};   
  colormap(hca,cmap) 
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
  tsdata = tsdata.tlim(irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:11:00.00Z'));
  if 1
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN; % doesn't work well with patch
    tsdata.data = log10(tsdata.data);
    c_eval('hp = irf_patch(hca,{tsdata,0});',ic)  
    hp.FaceColor = [0 0 0];
    hp.EdgeColor = [0 0 0];
  else
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN;
    tsdata.data = log10(tsdata.data);
    irf_plot(hca,tsdata,'k')
  end
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.201 0.201];
  hca.YLabel.Interpreter = 'tex';
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize_leg;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = fontsize;',1:numel(h))

linkprop([irf_panel('feeps omni mms 1234'),irf_panel('feeps pitch mms 1234')],{'CLim'})

h(end).XTickLabelRotation = 0;

%hca = irf_panel('feeps omni mms 1234'); hca.CLim = [-1 100];
%hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-16 25];
%hca = irf_panel('V ExB Vi x'); hca.YLim = [-200 1000]*0.99;
%hca = irf_panel('V ExB Vi y'); hca.YLim = [-500 800]*0.99;
%hca = irf_panel('V ExB Vi z'); hca.YLim = [-500 500]*0.99;
%hca = irf_panel('Pressure');  hca.YLim = [0 0.5];
%hca = irf_panel('beta'); hca.YLim = [1e-2 1e3];
%hca = irf_panel('Tepar/Teperp'); hca.YLim = [0.5 1.5];

h(6).YLim(1) = h(7).YLim(2)*1e-3;
%hca1 = irf_panel('feeps omni mms 1234');  
%hca2 = irf_panel('i DEF omni'); 
%hca2.YLim(2) = hca1.YLim(1);
%h(end).YLim  = [-350 1150];

drawnow
hmark = irf_pl_mark(h,time_mark,'k');
c_eval('hmark(?).LineStyle = '':'';',1:numel(hmark))
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))


hmark2 = irf_pl_mark(h,time_df,'k');
c_eval('hmark2(?).LineStyle = '':'';',1:numel(hmark2))
c_eval('hmark2(?).LineWidth = 1;',1:numel(hmark2))

% text: 'Moderate beta...'
if 1
    %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.67 .66]-dy,'string',{'moderate \beta','with v_{ix}'},'fontsize',fontsize_leg,'horizontalalignment','center');
  %annotation('arrow',[0.3 0.3],[.69 .66]-dy);
  annotation('textarrow',[0.55 0.55],[.67 .66]-dy,'string',{'moderate \beta','without v_{ix}'},'fontsize',fontsize_leg,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  %annotation('arrow',[0.57 0.57],[.69 .66]-dy);
  annotation('textarrow',[0.56 0.6],[.415 .433],'string',{'field-aligned','energetic ions','appearing'},'fontsize',fontsize_leg,'horizontalalignment','right');
elseif 0 % pressure panel included
  %%
  dy = 0.08;
  annotation('textarrow',[0.3 0.3],[.71 .73]-dy,'string',{'moderate \beta with v_{ix}'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.3 0.3],[.69 .66]-dy);
  annotation('textarrow',[0.57 0.57],[.7 .72]-dy,'string',{'moderate \beta without v_{i||}'},'fontsize',12,'horizontalalignment','center');
  %annotation('textarrow',[0.63 0.63],[.7 .72]-dy,'string',{'moderate \beta, v_{i||} \sim 0'},'fontsize',12,'horizontalalignment','center');
  annotation('arrow',[0.57 0.57],[.69 .66]-dy);
elseif 0
  %%
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

%% MMS SWT: Figure: off-time current sheet structure, current
ic = 1;
ic = 1;

comps = ['x','y','z'];
comps = ['y'];

fontsize_leg = 14;
fontsize = 14;

npanels = 1+numel(comps)*2;
nrows = 2;
ncols = 1;
%[h1,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
h1 = irf_plot(npanels);
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

plot_colors = [0.8 0.8 0.8; j_colors([4 3],:)];
leg_colors = [0.6 0.6 0.6; j_colors([4 3],:)];
%j_colors = matlab_colors;

if 1 % B gse
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.02 0.9],'fontsize',fontsize_leg);
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
  for comp = comps      
    isub = isub + 1;    
    hca = irf_panel(sprintf('J%s',comp));    
    set(hca,'ColorOrder',plot_colors)
    irf_plot(hca,{gseJav.(comp),gseJav.(comp).resample(timeline),gseJcurl.(comp).resample(timeline)},'comp');        
    hca.YLabel.String = {sprintf('J_%s',comp),'(nA/m^2)'};    
    set(hca,'ColorOrder',leg_colors)    
    irf_legend(hca,{sprintf('<J_{%s}^{FPI}> (%g ms)',comp,dt_resample*1e3),sprintf('<J_{%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('J_{%s}^{curl} (%g s)',comp,dt_resample)}',[0.01 0.99],'fontsize',fontsize_leg);  
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
    
    irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,dt_resample*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.99],'fontsize',fontsize_leg);          
    set(hca,'ColorOrder',leg_colors(3,:))     
    irf_legend(hca,{sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.07],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)}',[0.01 0.05],'fontsize',fontsize_leg);  
    %irf_legend(hca,{sprintf('<J_{e%s}^{FPI}> (%g ms)',comp,cadence*1e3),sprintf('<J_{e%s}^{FPI}> (%g s)',comp,dt_resample),sprintf('<J_{i%s}^{FPI}> (%g s)',comp,dt_resample)},[0.02 0.99],'fontsize',fontsize_leg);  
    
    %hl = findobj(hca,'type','line');    
    %c_eval('hl(?).LineWidth = 0.5;',1:2)    
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
%irf_pl_mark(h1(1),tint_harris)

h1(1).YLim = [-18 25]*0.99;
h1(2).YLim = [-15 30]*0.99;
h1(3).YLim = [-20 20]*0.99;

%c_eval('h1(?).Position(1) = 0.08;',1:numel(h1))
c_eval('h1(?).Layer = ''top'';',1:numel(h1))
c_eval('h1(?).FontSize = fontsize;',1:numel(h1))
hl = findobj(gcf,'type','line');
%c_eval('hl(?).LineWidth = 0.5;',1:numel(hl))

%annotation('textarrow',[0.33 0.38],[.18 .23],'string',{'<J_{ey}^{FPI}> changes','direction'}','fontsize',fontsize_leg,'horizontalalignment','center');
%annotation('textarrow',[0.36 0.38],[.18 .23],'string',{'<J_{ey}^{FPI}> changes direction','<J_{ey}^{FPI}> remains positive'}','fontsize',fontsize_leg,'horizontalalignment','center');

if 0 % binnin plots of current
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
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  xx = gseJav.tlim(tint_harris).y.data;
  yy = gseJcurl.resample(timeline).tlim(tint_harris).y.data;
  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -17;
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
  pcolor(hca,j_center,j_center,(N)')
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
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
end

if 1 % Jey,Jiy, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  xx = gseJeav.tlim(tint_harris).y.data;
  yy = gseJiav.resample(timeline).tlim(tint_harris).y.data;  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -17;
  jmax = -jmin;
  jstep = 0.5;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = '<J_{ey}^{FPI}>';
  hca.YLabel.String = '<J_{iy}^{FPI}>';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  hca.CLim = [0 35];
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))  
  plot(hca,j_edges,-j_edges,'color',[0.5 0.5 0.5],'linestyle','--')
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
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
  irf_legend(hca,{'data from';'yellow-shaded';'interval'},[0.98 0.08],'fontsize',fontsize_leg,'color',[0 0 0])
end


c_eval('h2(?).FontSize = fontsize;',1:numel(h2))
c_eval('h2(?).Box = ''on'';',1:numel(h2))
end

%% MMS SWT: Figure: off-time current sheet structure, vi vExB
ic = 1;

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 6;
%npanels = 1+numel(comps)*2;
nrows = 2;
ncols = 1;
%[h,h2] = initialize_combined_plot(npanels,nrows,ncols,0.6,'vertical');
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

j_colors = [mms_colors('xyz'); 0 0 0];

plot_colors = [0.8 0.8 0.8; j_colors([4 3],:)];
leg_colors = [0.6 0.6 0.6; j_colors([4 3],:)];

fontsize = 14;
fontsize_leg = 13;

isub = 0;
zoomy = [];
cs = 'gse';
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
  irf_legend(hca,{'x','y','z'}',[1.02 0.9],'fontsize',fontsize_leg);
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
  irf_legend(hca,{'v_{i{\perp}x}','v_{i{\perp}y}','v_{i{\perp}z}','v_{i||}'}',[1.02 0.98],'fontsize',fontsize_leg);
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
  irf_legend(hca,{'x','y','z'}',[1.02 0.9],'fontsize',fontsize_leg);
  irf_zoom(hca,'y')
end

elim_feeps = [8e4 Inf];
if 1 % vExB.x/y/z , Vi resample, ve, 3 panels
  for comp = ['x','y','z'] % vExB.x , Vi resample, ve  
    isub = isub + 1;
    zoomy = [zoomy isub];
    hca = irf_panel(['V ExB Vi ',comp]);
    set(hca,'ColorOrder',plot_colors)
    c_eval(sprintf('ve = %sVe?perp.(comp).resample(%sVi?);',cs,cs),ic);
    c_eval(sprintf('ne = ne?.resample(%sVi?);',cs),ic);    
    %ve.data(abs(ve.data)>5000) = NaN;
    ve.data(abs(ne.data)<0.04,:) = NaN;
    c_eval(sprintf('irf_plot(hca,{ve.resample(%sVi?),%sVi?perp.(comp),%sVExB?.resample(%sVi?).(comp)},''comp'');',cs,cs,cs,cs),ic)  
    %c_eval('irf_plot(hca,{gsmVExB?.x},''comp'');',ic)  
    hca.YLabel.String = {sprintf('v_{\\perp%s}',comp),'(km/s)'};
    set(hca,'ColorOrder',leg_colors)
    irf_legend(hca,{'v_e','v_i','ExB'}',[1.02 0.9],'fontsize',fontsize_leg);    
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
  h(ii).FontSize = fontsize;
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

%
hca = irf_panel(['V ExB Vi ','x']); hca.YLim = [-100 250]; hca.YLim = [-150 450]*0.99;
hca = irf_panel(['V ExB Vi ','y']); hca.YLim = [-300 700]*0.99;
hca = irf_panel(['V ExB Vi ','z']); hca.YLim = [-300 400]*0.99;
%hca = irf_panel('E gsm'); hca.YLim = [-15 15]*0.99;

c_eval('h(?).FontSize = fontsize;',1:numel(h))

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

h(end).XTickLabelRotation = 0;
%hat=annotation('textarrow',[0.593 0.593],[.2 .22],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
%hat=annotation('textarrow',[0.593 0.593],[.245 .26],'string',{'first energetic ions appearing'}','fontsize',12,'horizontalalignment','center','TextBackgroundColor',[1 1 1]);
if 0 % binning plots of current
isub = 1;
if 0 % ve, vExB scatter
  hca = h2(isub); isub = isub + 1;
  c_eval('xx = gseVe?.tlim(tint_harris).y.data;',ic)
  c_eval('yy = gseVExB?.resample(gseVe?).tlim(tint_harris).y.data;',ic)
  p = polyfit(xx,yy,1);
  plot(hca,xx,yy,'.')
  hca.XLabel.String = 'v_{ey}^{FPI}';
  hca.YLabel.String = 'v_{ExB}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  %plot(hca,[-100 100],[-100 100],'color',[0.5 0.5 0.5])
  %plot(hca,[-100 100],p(2)+[-100 100]*p(1))
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'off')
  %hca.XLim = 40*[-1 1];
  %hca.YLim = 40*[-1 1];  
end
if 1 % ve, vExB, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  timeline = timeline(1):0.5:timeline(end);
  cadence = timeline(2)-timeline(1);
  c_eval('xx = gseVe?.resample(timeline).tlim(tint_harris).y.data;',ic)
  c_eval('yy = gseVExB?.resample(timeline).tlim(tint_harris).y.data;',ic)
  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -1000;
  jmax = -jmin;
  jstep = 20;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = 'v_{ey}^{FPI}';
  hca.YLabel.String = 'v_{ExB}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,log10(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  %hca.CLim = [0 35];
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
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
  
end
if 1 % vi, vExB, binning
  hca = h2(isub); isub = isub + 1;
  timeline = gseJeav.time;
  cadence = timeline(2)-timeline(1);
  c_eval('xx = gseVi?.tlim(tint_harris).y.data;',ic)
  c_eval('yy = gseVExB?.resample(gseVi?).tlim(tint_harris).y.data;',ic)
  
  p = polyfit(xx,yy,1);
  %plot(hca,xx,yy,'.')
  jmin = -1000;
  jmax = -jmin;
  jstep = 20;
  j_edges = jmin:jstep:jmax;
  j_center = (jmin+jstep*0.5):jstep:(jmax-jstep*0.5);
  [N,edges,mid,loc] = histcn([xx yy], j_edges, j_edges);
  
  hca.XLabel.String = 'v_{iy}^{FPI}';
  hca.YLabel.String = 'v_{ExB}';
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hold(hca,'on')
  xlim = hca.XLim;
  ylim = hca.YLim;
  pcolor(hca,j_center,j_center,log10(N)')
  shading(hca,'flat')
  hcb = colorbar(hca,'peer',hca);
  hcb.YLabel.String = 'counts';
  %hca.CLim = [0 35];
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
  irf_legend(hca,sprintf('%g ms',cadence*1e3),[0.02 0.98],'fontsize',fontsize_leg)
end


c_eval('h2(?).FontSize = fontsize;',1:numel(h2))
c_eval('h2(?).Box = ''on'';',1:numel(h2))
colormap(irf_colormap('waterfall'))
end

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
fontsize = 13;
fontsize_leg = 12;
npanels = 4;
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
  set(hca,'ColorOrder',mms_colors('xyza'))  
  %c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  irf_plot(hca,{gseBav.x,gseBav.y,gseBav.z,tsB0.resample(tsB0.time(1):dt_resample:tsB0.time(end))},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'B_x','B_y','B_z','B_0 = 2\mu_0P_{tot}^{1/2}'}',[1.02 0.98],'fontsize',fontsize_leg);
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
if 1 % L specrec
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
if 1 % compare L with const and varying B0
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
  ts_yy = gseJcurl; y_str = '<J_{y}^{curl}>';
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
  [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:1000:6000],'k');
  % [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:500:6000],'k',"ShowText",true,'LabelFormat','L = %g km');
  clabel(cc,hh,'LabelSpacing',400,'Color','k','FontWeight','normal','fontsize',fontsize_leg);
  end
  B02 = 25;
  [BB2,JJ] = meshgrid(B_grid,J_grid);
  LL = mf_L(sqrt(BB2),B02,JJ);
  [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:1000:6000],'r');
  % [cc,hh] = contour(hca,BB2,JJ,LL*1e-3,[1000:500:6000],'k',"ShowText",true,'LabelFormat','L = %g km');
  clabel(cc,hh,'LabelSpacing',600,'Color','r','FontWeight','normal','fontsize',fontsize_leg);
  
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

fontsize = 14;
fontsize_leg = 13;
                    
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
  irf_legend(hca,{'B_x','B_y','B_z','|B|'}',[1.02 1.0],'fontsize',fontsize_leg);  
  %irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.98 0.12],'fontsize',fontsize_leg);  
end 
if 1 % ne
  isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?},''comp'')',ic)
  set(hca,'ColorOrder',mms_colors('12'))
  %irf_legend(hca,{'n_e'},[1.02 0.98],'fontsize',fontsize_leg)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
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
  irf_legend(hca,{'v_{ix}','v_{iy}','v_{iz}'}',[1.02 0.99],'fontsize',fontsize_leg);
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
    hleg = irf_legend(hca,{'Lower energy limit','for reduced distributions'}',[0.20 0.6],'color','k','fontsize',12);
    c_eval('hleg(?).HorizontalAlignment = ''center'';',1:numel(hleg))
    c_eval('hleg(?).VerticalAlignment = ''top'';',1:numel(hleg))
    
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
  x_str = 'v_x (km/s)';
  y_str = 'v_y (km/s)';
  z_str = 'v_z (km/s)';
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

c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h2(?).FontSize = fontsize;',1:numel(h2))
hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
c_eval('hb(?).FontSize = fontsize_leg;',1:numel(hb))
c_eval('h2(?).Position(1) = 0.56;',1:numel(h2))

h(1).YLim = [-17 24];

h(end).XTickLabelRotation = 0;
h2(1).Title.String = 'I. Hot pre-existing plasma sheet';
h2(2).Title.String = 'II. Post-DF cold and heated ions';
h2(1).Title.FontWeight = "normal";
h2(2).Title.FontWeight = "normal";
irf_legend(h(1),'I',[0.43 1.02],'color','k','fontsize',fontsize)
irf_legend(h(1),'II',[0.59 1.02],'color','k','fontsize',fontsize)

%c_eval('h2(?).CLim = [-14.5 -9.5];',1:numel(h2))
%hmark = irf_pl_mark(h,tocolumn(times_dist_exact_eu),'r');
%c_eval('hmark(?).LineWidth = 0.5; hmark(?).Color = [0 0 0]; hmark(?).LineStyle = ''-'';',1:numel(hmark))
hmark2 = irf_pl_mark(h,tocolumn(times_dist_exact_eu_edges(:)),'r');
c_eval('hmark2(?).LineWidth = 0.5; hmark2(?).Color = [0 0 0]; hmark2(?).LineStyle = ''--'';',1:numel(hmark2))

%
leg_str = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)'};
for ih = 1:numel(h)
  irf_legend(h(ih),leg_str{1},[0.01 0.98],'color','k','fontsize',fontsize_leg );
  leg_str = leg_str(2:end);
end
for ih = 1:numel(h2)
  irf_legend(h2(ih),leg_str{1},[0.02 0.98],'color','k','fontsize',fontsize_leg );
  leg_str = leg_str(2:end);
end

%% The two DFs or DF/island
ic = 1;

time_df1 = irf_time('2017-07-25T22:10:06.00Z','utc>EpochTT');
time_df2 = irf_time('2017-07-25T22:10:28.00Z','utc>EpochTT');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
feeps_ion_omni1234 = feeps_ion_omni1;
%feeps_ion_omni1234.data = (feeps_ion_omni1.data + feeps_ion_omni2.resample(feeps_ion_omni1).data + feeps_ion_omni3.resample(feeps_ion_omni1).data + feeps_ion_omni4.resample(feeps_ion_omni1).data)/4;
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

fontsize = 14;
fontsize_leg = 13;

npanels = 9;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

isub = 0;
zoomy = [];
cs = 'gse';
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
  irf_legend(hca,{'B_x','B_y','B_z'}',[1.02 0.9],'fontsize',fontsize_leg);
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
if 1 % Vi perp,par
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Vi perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval(sprintf('irf_plot(hca,{%sVi?perp.x,%sVi?perp.y,%sVi?perp.z,%sVi?par},''comp'');',cs,cs,cs,cs),ic)  
  hca.YLabel.String = {'v_{i}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'v_{\perp,x}','v_{\perp,y}','v_{\perp,z}','v_{||}'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'v_{i{\perp}x}','v_{i{\perp}y}','v_{i{\perp}z}','v_{i||}'}',[1.02 0.9],'fontsize',fontsize_leg);
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
if 0 % beta
  hca = irf_panel('beta');
  betalow = 1/0.3;
  betalow = 3;
  set(hca,'ColorOrder',mms_colors('1'))  
  c_eval('irf_plot(hca,{beta?},''comp'')',ic)
  hold(hca,'on')   
  hp = irf_patch(hca,{beta1,betalow},'smaller');
  hp.EdgeColor = 'none';
  hp.FaceColor = colors(1,:);
  hold(hca,'off')
  
  hca.YScale = 'log';
  %irf_legend(hca,{'hpca','fpi'},[0.98 0.98])
  hca.YLabel.String = '\beta';
  hca.YTick = 10.^[-5:1:5];
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1.01e-2 1e3];
  irf_legend(hca,['\beta < ' num2str(betalow,'%.2f')],[0.02 0.62],'color',hp.FaceColor,'fontsize',fontsize_leg)
end
if 1 % JotE
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('JdotE');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{JdotE?.resample(gseVi?)},''comp'');',ic)  
  hca.YLabel.String = {'J\cdot E','(nW/m^3)'}; % 
  set(hca,'ColorOrder',mms_colors('xyza'))
end

elim_feeps = [8e4 Inf];
if 1 % FEEPS Pitch all sc
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
  specrec = feeps_ion_omni1.elim(elim_feeps).specrec('energy');
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
  for comp = ['x']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('if1D = if1D%s?_700;',comp),ic)
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
    %hca.XGrid = 'on';
    %hca.YGrid = 'on';
    %hca.Layer = 'top';
    %irf_legend(hca,sprintf('E > %.0f eV',if1D.ancillary.energy(1,1)-if1D.ancillary.delta_energy_minus(1,1)),[0.02 0.08],'color','k','fontsize',fontsize_leg)
  end
end

if 1 % e DEF feeps omni all sc
  % for data see feeps_data.m
  isub = isub + 1;
  hca = irf_panel('feeps e omni mms 1234'); 
  specrec = feeps_ele_omni1.elim(elim_feeps).specrec('energy');
  specrec.f = specrec.f*1e-3;
  specrec.f_label =  {'E_e (keV)'};
  [hout,hcb] = irf_spectrogram(hca,specrec,'log');
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  hca.YLabel.String = {'E_e^{FEEPS}','(keV)'};   
end
if 1 % e DEF omni
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
  tsdata = tsdata.tlim(irf.tint('2017-07-25T22:05:30.00Z/2017-07-25T22:11:00.00Z'));
  if 1
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN; % doesn't work well with patch
    tsdata.data = log10(tsdata.data);
    c_eval('hp = irf_patch(hca,{tsdata,0});',ic)  
    hp.FaceColor = [0 0 0];
    hp.EdgeColor = [0 0 0];
  else
    tsdata.data(ne1.resample(tsdata).data<0.01) = NaN;
    tsdata.data = log10(tsdata.data);
    irf_plot(hca,tsdata,'k')
  end
  hca.YLabel.String = 'log_{10}(T_{e||}/T_{e\perp})';
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [-0.201 0.201];
  hca.YLabel.Interpreter = 'tex';
end
%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize_leg;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint)
%irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
tmark_eis = tint(1):20:tint(2);
%irf_pl_mark(h(3),tmark_eis.EpochUnix')
%h(1).Title.String = irf_ssub('MMS ?',ic);

c_eval('h(?).FontSize = fontsize;',1:numel(h))


h(end).XTickLabelRotation = 0;

irf_zoom(h,'x',tint_action)
%irf_zoom(h,'y')
c_eval('hmark_df!(?) = irf_pl_mark(h(?),time_df!.epochUnix,''k''); hmark_df!(?).LineStyle = ''--'';',1:numel(h),1:2)


%hca = irf_panel('feeps omni mms 1234'); hca.CLim = [-1 100];
hca = irf_panel('Vi perp'); hca.YLim = [-750 1150];
hca = irf_panel('B gsm'); hca.YLim = 0.99*[-16 25];
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
  %%
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


%% MMS SWT: Figure: cessation and onset 
% for feeps_ion_omni1234, see feeps_data.m
ic = 1;
time_mark = EpochTT('2017-07-25T22:09:08.00Z');
time_df = EpochTT('2017-07-25T22:10:06.50Z');

tint_figure = irf.tint('2017-07-25T22:07:00.00Z/2017-07-25T22:11:30.00Z');

%c_eval('Etop_fpi = iPDist?.ancillary.energy(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)
c_eval('Etop_fpi = iPDist?.depend{1}(1,end)+iPDist?.ancillary.delta_energy_plus(1,end);',ic)

npanels = 8;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
colors = mms_colors('matlab');

fontsize_leg = 14;

isub = 0;
zoomy = [];
cs = 'gse';
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
%  feeps_pa2.elim([9.7+04 Inf])
  isub = isub + 1;
  hca = irf_panel('feeps pitch mms 1234');  
  %elim = [0 32000];
  irf_spectrogram(hca,feeps_ion_pa1234.elim(elim_feeps).specrec('pitchangle'));
  hca.YLim = [0 180];
  hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('%.0f < E_i < %.0f keV',elim(1)*1e-3,elim(2)*1e-3),[0.05 0.98])
  irf_legend(hca,'four-spacecraft average',[0.05 0.98])
end
if 1 % i DEF feeps omni all sc
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

elim_eis = [00 4500]*1e3;
if 0 % EIS Pitch all
  isub = isub + 1;
  hca = irf_panel('eis Pitch all');  
  irf_spectrogram(hca,eis_pa1234.elim(elim_eis).specrec('pitchangle'));  
  hca.YLim = [0 180];
  %hca.YLabel.String = {'\theta','(deg)'};  
  hca.YLabel.Interpreter = 'tex';
end
if 0 % i DEF EIS omni all
  isub = isub + 1;
  hca = irf_panel('i DEF eis omni all');  
  [hout,hcb] = irf_spectrogram(hca,eis_omni1234.elim(elim_eis).specrec('energy'),'log');
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
    c_eval(sprintf('if1D = if1D%s?_700;',comp),ic)
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
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',fontsize_leg)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

%irf_zoom(h(1:iisub),'x',fastTint)
%irf_zoom(h,'x',tint_figure)
irf_zoom(h,'x',tint)
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

%h1 = irf_panel('feeps omni mms 1234');  
%h2 = irf_panel('i DEF omni');
%h1.YLim(1) = h2.YLim(2);

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
%hat=annotation('textarrow',[0.593 0.593]+0.0094,[.56 .54]+0.08,'string',{'first energetic ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none');
%hat=annotation('textarrow',[0.593 0.593]+0.104,[.4 .38]+0.01,'string',{'first cold ions appearing'}','fontsize',fontsize_leg,'horizontalalignment','center','TextBackgroundColor','none');
%

