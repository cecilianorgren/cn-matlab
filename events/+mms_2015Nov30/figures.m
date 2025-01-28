%% Figure 1 v2: Overview: Ions larges scale, zoom in on electrons 
tint = irf.tint('2015-11-30T00:21:44.00Z/2015-11-30T00:27:00.00Z');
npanels = 8;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x,gseB?.y,gseB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('xy'))
  c_eval('irf_plot(hca,{ni?,ne?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('xy'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  c_eval('irf_plot(hca,{dbcsVi?.x,dbcsVi?.y,dbcsVi?.z},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % Ve
  iisub = iisub+1;
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVi?.x,gseVi?.y,gseVi?.z},''comp'');',ic)
  c_eval('irf_plot(hca,{dbcsVe?.x,dbcsVe?.y,dbcsVe?.z},''comp'');',ic)
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % E
  iisub = iisub+1;
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x,gseE?.y,gseE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end

if 1 % B scm  
  hca = irf_panel('scm B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?scm.x,gseB?scm.y,gseB?scm.z},''comp'');',ic)
  hca.YLabel.String = {'B_{SCM}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1
  hca = irf_panel('pfft E');
  c_eval('irf_spectrogram(hca,pfftE?);',ic)
  hca.YScale ='log';
  hca.YLabel.String = 'f (Hz)';
end
if 1
  hca = irf_panel('pfft B');
  c_eval('irf_spectrogram(hca,pfftB?);',ic)
  hca.YScale='log';
  hca.YLabel.String = 'f (Hz)';
end

irf_zoom(h,'x',tint)
irf_zoom(h(1:5),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);


%%
if 1 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni error');  
  c_eval('irf_spectrogram(hca,iPDistError?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni counts');  
  c_eval('irf_spectrogram(hca,iPDistCounts?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  %hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 1 % iPDist pa 32
  iisub = iisub+1;
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,18).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap)
end

hca = irf_panel('delete for space');

if 0 % Ve n>1cc
  c_eval('Veplot = gseVe?;',ic)
  c_eval('neplot = ne?;',ic)
  tmpdata = Veplot.data;
  indLowNe = find(neplot.data<1);
  tmpdata(indLowNe,:) = repmat([NaN NaN NaN],numel(indLowNe),1);
  Veplot.data =tmpdata;
  c_eval('Veplot = gseVe?.clone(gseVe?.time,tmpdata);',ic)
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{Veplot.x.tlim(tint),Veplot.y.tlim(tint),Veplot.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'V_e_x','V_e_y','V_e_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'n>1 cm^{-3}'},[0.08 0.95],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [-1200 700];  
end
if 1 % Ve all n
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'V_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_e_x','V_e_y','V_e_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.15],'fontsize',12);
  hca.YLim = [-900 400];  
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap) 
end
if 0 % e DEF omni 32
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % ePDist omni 32
  hca = irf_panel('e DEF omni');
  c_eval('irf_spectrogram(hca,ePDist?.omni(''e'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux');  
  eint = [0 200];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [200 1000];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all e');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(eTint).pitchangles(dmpaB?,20).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

if 0 % i DEF omni 64
  hca = irf_panel('i DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iDEFomni64_?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hout,irf_colormap('space'))  
end
if 0 % i DEF omni 32
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(tint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % J curl  
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'4 sc average'},[0.98 0.9],'fontsize',12,'color','k');
end

%load('caa/cmap.mat');
for ii = [4 5 8 9]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',iTint)
irf_zoom(h(iisub+2:npanels),'x',eTint)
%irf_zoom(h([1:3 6 10]),'y')
hca = irf_panel('Te'); irf_zoom(hca,'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
hca = irf_panel('delete for space'); delete(hca) 

%% Figure 1: 1 sc overview 
tint = irf.tint('2015-11-30T00:21:44.00Z/2015-11-30T00:26:44.00Z');
tintZoom = irf.tint('2015-11-30T00:24:20.00Z/2015-11-30T00:24:35.00Z');
npanels = 11;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 12];
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.e64.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 1 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.deflux.tlim(tint).pitchangles(dmpaB?,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
  %hca.YTick
end
if 0 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 0 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 200];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 64
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [200 1000];  
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h,'x',tint)
irf_zoom(h([1:3 6:7 11]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Figure with rh panels for particle distributions
%% Plot detailed electron distribution for selected times
% Initialize figure
[h1,h2] = initialize_combined_plot(8,2,2,0.4,'vertical');
cmap = 'jet';
ic = 1;
tintZoom = irf.tint('2015-11-30T00:24:26.00Z/2015-11-30T00:24:28.00Z');

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1200 700];  
end
if 0% Ve perp
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  ylabel(hca,{'v_{e,\perp}','(km/s)'},'interpreter','tex');
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
end
if 0 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 0 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
end
if 1 % eDist omni 64
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom).pitchangles(dmpaB?,pas);',ic)
  hca = irf_panel('e DEF 1');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom).pitchangles(dmpaB?,pas);',ic)
  hca = irf_panel('e DEF 2');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom).pitchangles(dmpaB?,pas);',ic)
  hca = irf_panel('e DEF 3');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.1],'fontsize',12,'color',[0 0 0]);
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:3 8]),'y')

for ip = 4:7
  h1(ip).YLim = [10 1000];
  h1(ip).XGrid = 'off';
  h1(ip).YGrid = 'off';
end

irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

%% Plot single time particle distributions, 1 sc, 4 projections,
time = irf_time('2015-11-30T00:24:27.13Z','utc>epochtt')+0*0.03;
if exist('hmark'); delete(hmark); end
hmark = irf_pl_mark(h1,time.epochUnix','green');
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)

% Projection coordinate system
c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
c_eval('hatExB = cross(hatE,hatB);',ic)

y = hatE;
z = hatB;
x = cross(y,z)/norm(cross(y,z));

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 6];  
  

timeUTC = time.utc;      
isub = 1;

% Perpendicular plane    
hca = h2(isub); isub = isub + 1; 
xyz = [x;y;z]; vlabels = {'v_{ExB}','v_E','v_B'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
%mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';
hca.Title.String = timeUTC(1:23);

% B plane 1
hca = h2(isub); isub = isub + 1; 
xyz = [x;z;-y]; vlabels = {'v_{ExB}','v_B','-v_{E}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';

% B plane 2
hca = h2(isub); isub = isub + 1;
xyz = [y;z;x]; vlabels = {'v_E','v_B','v_{ExB}'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
hca.Title.String = '';

% Perp plane, but with markings
vectors = {hatExB,'ExB'; hatE,'E';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
hca = h2(isub); isub = isub + 1; 
xyz = [x;y;z]; vlabels = {'v_{ExB}','v_E','v_B'};
mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
hca.Title.String = '';

for ii = 1:4
  colormap(h2(ii),strCMap)
end