%% Separatrix streaming
% Make reduced distribution
%tintZoom = irf.tint('2017-07-18T01:40:30.00Z',20);
tint = irf.tint('2017-07-09T04:41:43.00Z/2017-07-09T04:43:32.00Z');
%tintZoom = irf.tint('2017-07-06T08:18:00.00Z',13);
strTintZoom = [irf_time(tintZoom(1),'epochtt>utc_yyyymmdd_HHMMSS') '_' irf_time(tintZoom(2),'epochtt>utc_HHMMSS')];

eint = [000 40000];
vint = [-Inf Inf];

eDist = ePDist1.tlim(tintZoom).elim(eint);
iDist = iPDist1.tlim(tintZoom).elim(eint);
ve = gseVe1.tlim(eDist.time).resample(eDist);
vi = gseVi1.tlim(iDist.time).resample(iDist);
scpot_margin = 1.5; % keep in mind that this also affects the velocity at lower energies
scpot_lim = scPot1.resample(eDist)*scpot_margin;
eLine = dmpaB1.resample(eDist).norm;
iLine = dmpaB1.resample(iDist).norm;

ePara = dmpaB1.resample(eDist).norm;
ePerp1 = ePara.cross(irf.ts_vec_xyz(ePara.time,repmat([1 0 0],ePara.length,1))).norm;
ePerp2 = ePara.cross(ePerp1).norm;

iPara = dmpaB1.resample(iDist).norm;
iPerp1 = iPara.cross(irf.ts_vec_xyz(iPara.time,repmat([1 0 0],iPara.length,1))).norm;
iPerp2 = iPara.cross(ePerp1).norm;

tic; if1Dx = iDist.reduce('1D',[1 0 0],'vint',vint); toc % reduced distribution along B
tic; if1Dy = iDist.reduce('1D',[0 1 0],'vint',vint); toc % reduced distribution along B
tic; if1Dz = iDist.reduce('1D',[0 0 1],'vint',vint); toc % reduced distribution along B

%tic; ef1D = eDist.reduce('1D',eLine,'vint',vint,'scpot',scpot_lim); toc % reduced distribution along B
tic; if1D = iDist.reduce('1D',iLine,'vint',vint); toc % reduced distribution along B
lineVe = ve.dot(eLine); % projection of Vi on B
lineVi = vi.dot(iLine); % projection of Vi on B


%% Plot
ic = 1;
npanels = 9;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];

if 0 % B absh
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

  if 1 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?.abs.resample(iPDist?);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'k')
    hold(hca,'off')
  end
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
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
if 1 % i psd vpar
  isub = isub + 1;
  hca = irf_panel('iLine');
  irf_spectrogram(hca,if1D.specrec('velocity_1D'));
  hold(hca,'on')
  irf_plot(hca,{lineVi},'comp')
  %irf_plot(hca,gseVi1)
  hold(hca,'off')
  hca.YLim = if1D.depend{1}(1,[1 end]);
  hca.YLabel.String = 'v_{i,||} (km/s)'; 
end
if 1 % i psd x
  isub = isub + 1;
  hca = irf_panel('iLine x');
  irf_spectrogram(hca,if1Dx.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  
  if 1 % % Vi
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
  irf_spectrogram(hca,if1Dy.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 1 % % Vi
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
  irf_spectrogram(hca,if1Dz.specrec('velocity_1D'));  
  hca.YLim = if1D.depend{1}(1,[1 end]);
  if 1 % % Vi
    hold(hca,'on')    
    irf_plot(hca,gseVi1.z,'k')
    hold(hca,'off')
  end
  hca.YLabel.String = {'v_{iz}','(km/s)'};  
  hca.YLabel.Interpreter = 'tex';
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
if 0 % eDEF omni
  isub = isub + 1;
  hca = irf_panel('eDEF');
  [hout,hcb] = irf_spectrogram(hca,eDist.deflux.omni.specrec,'log');
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
if 1 % E 
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E ');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
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
%colormap('jet');
colormap(pic_colors('candy4'))

%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];

%% Plot fred, electrons
ic = 1;
npanels = 7;
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
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{gseVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(km/s)'};  
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
if 0 % ne
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
colormap(cn.cmap('blue_white'));
colormap(irf_panel('fe reduced * v'),cn.cmap('blue_red'))
hca = irf_panel('phase velocity');
hca.CLim = [-5 -2];

hca = irf_panel('fe reduced * v'); hca.CLim = min(abs(h(4).CLim))*[-1 1];
%colormap(irf_panel('fe reduced * v^2'),cn.cmap('white_blue'))

%h=irf_plot({gseB1,gseVi1,iPDist1.deflux.omni.specrec('energy'),f1D.specrec('velocity_1D')}); h(3).YScale = 'log'; %h(4).YLim = [-1000 1000];

%% Plot, including f proj and v phi and vtrap
ic = 1;
npanels = 6;
h = irf_plot(npanels); 
isub = 0;
zoomy = [];
tint_zoom = tintZoom;
%tint_zoom = irf.tint('2017-07-06T01:38:00.00Z/2017-07-18T01:39:00.00Z');

vmin = tsVphpar-tsVtrap;
vmax = tsVphpar+tsVtrap;
c_eval('vmin? = tsVphpar-tsVtrap?;',1:4)
c_eval('vmax? = tsVphpar+tsVtrap?;',1:4)
  
  
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
  hca = irf_panel('eLine');
  %irf_plot(hca,ef1D.specrec('velocity_1D'));
  irf_spectrogram(hca,ef1D.specrec('velocity_1D','10^3 km/s'));
  %hold(hca,'on')
  %irf_plot(hca,{lineVe},'comp')
  %set(hca,'ColorOrder',mms_colors('122'))
  %irf_plot(hca,{tsVphpar,vmin,vmax},'comp');
  %irf_plot(hca,gseVi1)
  %hold(hca,'off')
  %hca.YLim = ef1D.depend{1}(1,[1 end]);
  hca.YLabel.String = {'v_{e,||}','(10^3 km/s)'}; 
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
if 0 % ne
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
  hh = irf_plot(hca,{tsPhi});
  c_eval('hh(?).Color = mms_colors(''?'');',1:4)
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
  c_eval('vmin? = tsVphpar-tsVtrap?;',1:4)
  c_eval('vmax? = tsVphpar+tsVtrap?;',1:4)
  
  %set(hca,'ColorOrder',mms_colors('111223344'))
  set(hca,'ColorOrder',mms_colors('111111111'))
  vscale = 1e-3;
  irf_plot(hca,{tsVphpar*vscale,vmin1*vscale,vmax1*vscale,vmin2*vscale,vmax2*vscale,vmin3*vscale,vmax3*vscale,vmin4*vscale,vmax4*vscale},'comp');
  
  %irf_patch(hca,{vmin,vmax})
  %hca.YLim = sort(real([max([vmax1.data; vmax2.data; vmax3.data; vmax4.data]) min([vmin1.data; vmin2.data; vmin3.data; vmin4.data])]));
  hold(hca,'off')
  hca.YLabel.String = {'v_{||}','(10^3 km/s)'};  
  set(hca,'ColorOrder',mms_colors('111'))
  
  hleg = irf_legend(hca,{'- v_{ph}+v_{trap}';'- v_{ph}         ';'- v_{ph}-v_{trap}'},[0.95 0.99],'fontsize',12);  
  
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

%% Plot 1D plot
%if exist('fig') && ~isempty(fig) , close(fig); end
units = irf_units;
iesw = 11;
fig = figure(22);
fig.Position = [729   765   442   260];
vph = esw_data{9}(iesw);
c_eval('phi(?) = esw_data{9+?}(iesw);',1:4)
vtrap = sqrt(2*units.e*phi/units.me)*1e-3; % km/s
vtrap = max(vtrap);
time = EpochTT(esw_data{5}{iesw});

tindPitch = find(abs(ePitch1.time-time) == min(abs(ePitch1.time-time)));
tindF1D =   find(abs(ef1D.time   -time) == min(abs(ef1D.time   -time)));
epitch = ePitch1(tindPitch);
ef1d = ef1D(tindF1D);

ylim = [0 1.4e-3];

if 0
  xlim = [1e1 1e4];
  ylim = [1e-33 1e-25];
  figure(10)
  hca = subplot(2,1,1);
  plot(hca,epitch.depend{1},squeeze(epitch.data))
  hca.XScale = 'log';
  hca.YScale = 'log';
  hca.XLim = xlim;
  hca.YLim = ylim;
  hold(hca,'on')

  hold(hca,'off')
end
hca = subplot(1,1,1);
hline = plot(hca,ef1d.depend{1}*1e-3,squeeze(ef1d.data),'color',mms_colors('1'));
hline.LineWidth = 2;
hca.XLabel.String = 'v_{||} (10^3 km/s)';
hca.YLim = ylim;
%hca.YScale = 'log';
hold(hca,'on')
hvph = plot(hca,-vph*1e-3*[1 1],hca.YLim,'color',mms_colors('4')); 
hvph.LineWidth = 2;
hpatch = patch(hca,[-vph-vtrap -vph-vtrap -vph+vtrap -vph+vtrap]*1e-3,[hca.YLim hca.YLim([2 1])],'c');
hpatch.FaceAlpha = 0.2;
hpatch.EdgeColor = hpatch.FaceColor;
hpatch.EdgeAlpha = hpatch.FaceAlpha;
hold(hca,'off')
hca.XLim = [-30 30];
hca.YLabel.String = ['f_{e,red} (' ef1d.units ')'];

hca.Title.String = time.utc;

%hca.XScale = 'log';
%hca.YScale = 'log';
%hca.XLim = xlim;

%irf_legend()
set(hca,'ColorOrder',[0 0 0;mms_colors('4');[0 1 1]])
irf_legend(hca,{'f_e';'vph';'v_{trap}'},[0.95 0.95],'fontsize',14)
hca.FontSize = 14;
cn.print(sprintf('fred_e_vtrap_1time_%s',time.utc),'path',eventPath)

