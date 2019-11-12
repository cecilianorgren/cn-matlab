mms_20170706_005403.load_data
mms_20170706_005403.prepare_data_single_sc

%% Goddard talk
ic = 1;
npanels = 3;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

if 1 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,ePDist1.deflux.omni.specrec,'log');
  hold(hca,'on')
  lineScpot = irf_plot(hca,scPot1,'k');
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_e','(eV)'};   
  hca.CLim = [4 8];
end
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
if 0 % Ve
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x*1e-3,gseVe?.y*1e-3,gseVe?.z*1e-3},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % Te par perp
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
if 1 % Te par perp
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{(facTe?.xx+facTe?.yy+facTe?.zz)/3},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  %irf_legend(hca,{'T_{e}','T_{e,\perp}'},[0.98 0.9],'fontsize',12);  
end
if 0 % Te par perp
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


irf_zoom(h,'x',tint)
irf_zoom(h(zoomy),'y')
irf_plot_axis_align
colormap('jet');
