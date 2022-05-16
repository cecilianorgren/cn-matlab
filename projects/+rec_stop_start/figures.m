tint_b = irf.tint('2017-07-25T22:09:00.82Z/2017-07-25T22:11:31.82Z');

tint_c = irf.tint('2017-07-25T22:09:00.82Z/2017-07-25T22:11:00.00Z');
tint_d = irf.tint('2017-07-25T22:09:30.82Z/2017-07-25T22:11:00.00Z');

tint_action = irf.tint('2017-07-25T22:09:40.00Z/2017-07-25T22:10:50.00Z');
tint_cavity = irf.tint('2017-07-25T22:09:46.00Z/2017-07-25T22:09:56.00Z');

%% Figure: Overview 1
ic = 2;

npanels = 11;
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
if 1 % vExB gsm, Vi resample
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
if 1 % E gsm
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
if 1 % Te/i par perp Ti/Tref
  isub = isub + 1;
  zoomy = [zoomy isub];
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
if 1 % iPitch
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
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Figure: V ExB
ic = 2;

Etop_fpi = iPDist2.ancillary.energy(1,end)+iPDist2.ancillary.delta_energy_plus(1,end);

npanels = 6;
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
if 1 % vExB gse
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

%% Figure: momentum equation / ohms law
ic = 2;

npanels = 6;
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
if 0 % JxB, 3 panels: x,y,z
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

%% Figure: V ExB
ic = 1;

Etop_fpi = iPDist2.ancillary.energy(1,end)+iPDist2.ancillary.delta_energy_plus(1,end);

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

if 1 % i psd x,y,z, 3 panels
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
ic = 2;

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
  c_eval('irf_plot(hca,{gseB?_fast.x,gseB?_fast.y,gseB?_fast.z,gseB?_fast.abs},''comp'');',ic)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);  
end 
if 1 % B gse fast
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

npanels = 7;
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
if 1 % V oplus gse
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
if 1 % n log scale, no burst
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
if 0 % hplus DEF omni
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
if 0 % oplus DEF omni
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


irf_zoom(h,'x',tint_fast)
irf_plot_axis_align
