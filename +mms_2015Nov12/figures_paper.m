pathPaper =  '/Users/Cecilia/Dropbox (IRFU)/MMS_Nov12/';

%% Overview, fast + burst data
fastTint = irf.tint('2015-11-12T06:10:00.00Z/2015-11-12T07:30:00.00Z');
brstTint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
npanels = 7;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 0 % B
  iisub = iisub+1;
  hca = irf_panel('B fast');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?srvy.x,gseB?srvy.y,gseB?srvy.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 0 % J curl  
  iisub = iisub+1;
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurlsrvy.x.tlim(fastTint),gseJcurlsrvy.y.tlim(fastTint),gseJcurlsrvy.z.tlim(fastTint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  %hca.YLim = [-1200 1500];
end
if 0 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?fast.x,gseVi?fast.y,gseVi?fast.z},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-300 220];
end
if 0 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?fast},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni counts');  
  c_eval('irf_spectrogram(hca,iPDist?fast.deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end


if 1 % B
  hca = irf_panel('B brst');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 1 % Vi  
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 12];
end
if 1 % J  
  hca = irf_panel('J zoom');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Ve  
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.e64.omni.deflux.specrec,''log'');',ic)  
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
if 0 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux');  
  eint = [0 200];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(brstTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [10 1000];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(brstTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f') ' eV'],[0.99 0.90],'color',0*[1 1 1])
  %hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YLabel.String = {'Pitchangle','e^- (\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(brstTint).pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all e');  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(brstTint).pitchangles(dmpaB?,20).deflux.specrec(''pa''),''log'');',ic)
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
  c_eval('irf_plot(hca,{facPe?.xx.tlim(brstTint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(brstTint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'4 sc average'},[0.98 0.9],'fontsize',12,'color','k');
end

%load('caa/cmap.mat');
%for ii = [4 5 8 9]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  %colormap(h(ii),'jet')
%end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

if iisub, 
  irf_zoom(h(1:iisub),'x',fastTint); 
  irf_zoom(h(iisub+2:npanels),'x',brstTint)
else
  irf_zoom(h(1:npanels),'x',brstTint)
end
%irf_zoom(h([1:3 6 10]),'y')

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
h(iisub).XLabel.String = '';
irf_pl_mark(h(1:iisub),brstTint, 'green')


%%
hca = irf_panel('Te'); irf_zoom(hca,'y')
hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
hca = irf_panel('J curl'); hca.YLim = [-500 600];
hca = irf_panel('Vi'); hca.YLim = [-400 230];
hca = irf_panel('Ve'); hca.YLim = [-900 400];
hca = irf_panel('Te'); hca.YLim = [30 110];
irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

%% Overview, only fast
fastTint = irf.tint('2015-11-12T06:10:00.00Z/2015-11-12T07:30:00.00Z');
brstTint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
npanels = 4;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B fast');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?srvy.x,gseB?srvy.y,gseB?srvy.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % J curl  
  iisub = iisub+1;
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseJcurlsrvy.x.tlim(fastTint),gseJcurlsrvy.y.tlim(fastTint),gseJcurlsrvy.z.tlim(fastTint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  %hca.YLim = [-1200 1500];
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?fast.x,gseVi?fast.y,gseVi?fast.z},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-300 220];
end
if 0 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?fast},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni counts');  
  c_eval('irf_spectrogram(hca,iPDist?fast.deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % ePDist omni 64
  iisub = iisub+1;
  hca = irf_panel('e DEF omni counts');  
  c_eval('irf_spectrogram(hca,ePDist?fast.deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end


%load('caa/cmap.mat');
for ii = [4]
  colormap(h(ii),'jet')
end
legends = {'a)','b) curlometer','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.95],'color',[0 0 0])
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',fastTint)
hca = irf_panel('B fast'); hca.YLim = [-40 40];
hca = irf_panel('J curl'); hca.YLim = [-450 600];
hca = irf_panel('Vi'); hca.YLim = [-450 300];
%irf_zoom(h([1:3 6 10]),'y')
%hca = irf_panel('i DEF omni counts');  hca.YLim = []
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
h(iisub).XLabel.String = '';
irf_pl_mark(h(1:iisub),brstTint, 'green')

%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
%% Overview, only burst data
fastTint = irf.tint('2015-11-12T06:10:00.00Z/2015-11-12T07:30:00.00Z');
brstTint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
tintBCS = irf.tint('2015-11-12T07:19:19.80Z/2015-11-12T07:19:22.70Z');
npanels = 9;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');

if 1 % B
  hca = irf_panel('B brst');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
  if 1 % plot ion scale
    %%
    T0 = tintBCS(1)+0.5*(tintBCS.stop-tintBCS.start);
    vBCS = 70; % km/s
    rploc = 250; % km
    Lploc = 90; % km
    rplocT = rploc/vBCS;
    LplocT = Lploc/vBCS;
    tintRP = T0 + 0.5*rplocT*[-1 1];   
    tintLP = T0 + 0.5*LplocT*[-1 1];   
    tsRP = irf.ts_scalar(tintRP,27*[1 1]');
    tsLP = irf.ts_scalar(tintLP,20*[1 1]');
    hold(hca,'on') 
    colors = mms_colors('matlab');
    hh1 = irf_plot(hca,tsRP,'linewidth',2,'color',[0 0 0]);
    irf_legend(hca,'\rho_p',[0.45 0.98])
    hh2 = irf_plot(hca,tsLP,'linewidth',2,'color',[0 0 0]);
    irf_legend(hca,'L_p',[0.41 0.80])
    hold(hca,'off')
  end
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 12];
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(brstTint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 1 % Ve x B
  hca = irf_panel('E + Ve x B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseEVexB?_new.x.tlim(tint),gseEVexB?_new.y.tlim(tint),gseEVexB?_new.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E + v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % J  
  hca = irf_panel('J zoom');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(brstTint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % e DEF omni 64
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
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [10 400];  
  try 
    c_eval('irf_spectrogram(hca,ePitch?.elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(brstTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',brstTint)
%irf_zoom(h([1:3 6 10]),'y')
%hca = irf_panel('Te'); irf_zoom(hca,'y')
%hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
hca = irf_panel('n'); hca.YLim = [0 11];
hca = irf_panel('J zoom'); hca.YLim = [-500 800];
hca = irf_panel('Ve'); hca.YLim = [-900 500];
hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
hca = irf_panel('Ve x B'); hca.YLim = [-5 5];
%hca = irf_panel('Ti'); hca.YLim = [300 700];
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
c_eval('hmark(?) = irf_pl_mark(h(?),tintBCS, ''yellow'');',1:7)
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
for ii = 1:npanels;
  h(ii).FontSize = 12;
end
%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))
%% Overview, only burst data, started 12-12-2017
fastTint = irf.tint('2015-11-12T06:10:00.00Z/2015-11-12T07:30:00.00Z');
brstTint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
tintBCS = irf.tint('2015-11-12T07:19:19.80Z/2015-11-12T07:19:22.70Z');
timeMP = irf_time('2015-11-12T07:19:17.80Z','utc>epochtt');
npanels = 13;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
cmap = colormap('jet');

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};  
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
  if 1 % plot ion scales
    %%
    T0 = tintBCS(1)+0.5*(tintBCS.stop-tintBCS.start);
    vBCS = 70; % km/s
    rploc = 250; % km
    Lploc = 90; % km
    rplocT = rploc/vBCS;
    LplocT = Lploc/vBCS;
    tintRP = T0 + 0.5*rplocT*[-1 1];   
    tintLP = T0 + 0.5*LplocT*[-1 1];   
    tsRP = irf.ts_scalar(tintRP,27*[1 1]');
    tsLP = irf.ts_scalar(tintLP,20*[1 1]');
    hold(hca,'on') 
    colors = mms_colors('matlab');
    hh1 = irf_plot(hca,tsRP,'linewidth',2,'color',[0 0 0]);
    irf_legend(hca,'\rho_p',[0.45 0.98],'fontsize',14)
    hh2 = irf_plot(hca,tsLP,'linewidth',2,'color',[0 0 0]);
    irf_legend(hca,'L_p',[0.41 0.80],'fontsize',14)
    hold(hca,'off')
  end
  hca.YLabel.String = {'B','(nT)'};  
end
if 1 % J
  hca = irf_panel('J');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?.x.tlim(tint),gseJ?.y.tlim(tint),gseJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-490 900];  
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-910 780];  
  hca.YTick = [-1200:400:1200];
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(brstTint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [21 105];
  %hca.YTick
  %irf_zoom(hca,'y')
end
if 1 % P, 1sc, given above
  hca = irf_panel(sprintf('P sc %g',ic));
  Pishift = 0.4;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',ic)
  hca.YLabel.String = {sprintf('P',ic),'(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.98 0.9],'fontsize',12);
  hca.YLim = [0.0 0.3];  
end
if 1 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{beta?},''comp'');',ic)
  hca.YLabel.String = {'\beta'};
  hca.YScale = 'log';
  hca.YLim = [0.09 200];
  hca.YMinorTick = 'on';
  %hca.YMinorTick = [0.1:0.1:10 20:10:100 200];
  hca.YTick = [0.1 1 10 100];
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 9.9];
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
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('plotE = gseE?perp - 0*gseE?.filt(1,0,[],3);',ic)
  
  irf_plot(hca,{plotE.x.tlim(tint),plotE.y.tlim(tint),plotE.z.tlim(tint)},'comp');
  %c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-5 5];  
end
if 1 % J x B
  hca = irf_panel('E + Ve x B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJxB?.x.tlim(tint),gseJxB?.y.tlim(tint),gseJxB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J x B','()'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 1 % diff P
  hca = irf_panel('gradPi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  scale = 1e9*1e-3;
  irf_plot(hca,{gseGradPe.x.tlim(tint)*scale,gseGradPe.y.tlim(tint)*scale,gseGradPe.z.tlim(tint)*scale},'comp');
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'div(Pe)','()'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-3 3]*scale;  
end
if 0 % Ve x B
  hca = irf_panel('E + Ve x B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseEVexB?_new.x.tlim(tint),gseEVexB?_new.y.tlim(tint),gseEVexB?_new.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E + v_e x B','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % Pe par perp
  hca = irf_panel('Pe');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{facPe?.xx.tlim(brstTint),(facPe?.yy+facPe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'P_e','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'P_{||}','P_{\perp}'},[0.98 0.9],'fontsize',12);
end
if 1 % i DEF omni 64
  hca = irf_panel('i DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]); 
  hca.YLabel.String = {'E_i','(eV)'};   
  colormap(hca,cmap) 
end
if 1 % e DEF omni 64
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
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [10 1000];  
  try 
    c_eval('irf_spectrogram(hca,ePitch?.elim(eint).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(brstTint).e64.pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f') ' eV'],[0.99 0.90],'color',0*[1 1 1],'fontsize',14)
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h(1:npanels),'x',brstTint)
%irf_zoom(h([1:3 6 10]),'y')
if 0
%hca = irf_panel('Te'); irf_zoom(hca,'y')
%hca = irf_panel('e PA e64 deflux lowe');  hca.YLim = [0 180];
hca = irf_panel('n'); hca.YLim = [0 11];
hca = irf_panel('J zoom'); hca.YLim = [-500 800];
hca = irf_panel('Ve'); hca.YLim = [-900 500];
hca = irf_panel('Te'); hca.YLim = [20 110]; hca.YTick = [20:20:120];
hca = irf_panel('gradPe'); hca.YLim = [-2.2 2.2];
hca = irf_panel('B brst'); hca.YLabel.String = {'B','(nT)'}
hca = irf_panel('Ve x B'); hca.YLim = [-5 5];
end
%hca = irf_panel('Ti'); hca.YLim = [300 700];
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
colors = mms_colors('matlab');
c_eval('hmark_mp = irf_pl_mark(h(?),timeMP,0*colors(1,:));',1:npanels)
irf_legend(h(1),'MP',[0.25 1.05],'fontsize',14);
irf_legend(h(1),'BCS',[0.35 1.05],'fontsize',14);
hleg = irf_legend(h(1),'unperturbed MSH',[0.97 1.05],'fontsize',14);
hleg.HorizontalAlignment = 'center';

if 1 % irf_pl_hmark
  %%
imark = 1;
clear hmark;

c_eval('hmark(imark) = irf_pl_mark(h(?),tintBCS, ''yellow''); imark = imark + 1;',1:(npanels-3))
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h(ii).FontSize = 14;
end

hca = irf_panel('e DEF omni 64'); hca.YGrid = 'off'; hca.XGrid = 'off';
%irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

%% Constrained and unconstrained minimum variance magnetic field (average), 2 panels
% tint_bcs
c_eval('mvaB?_unc = gseB?*mva_mean''; mvaB?_unc.name = ''B LMN (unconstr.)'';')
c_eval('mvaB?_Bn0 = gseB?*mva_mean_td''; mvaB?_Bn0.name = ''B LMN (<BN>=0)'';')

mvaAvB_unc = (mvaB1_unc+mvaB2_unc.resample(mvaB1_unc.time)+mvaB3_unc.resample(mvaB1_unc.time)+mvaB4_unc.resample(mvaB1_unc.time))/4; mvaAvB_unc.name = 'avg B LMN  (unconstr.)';
mvaAvB_Bn0 = (mvaB1_Bn0+mvaB2_Bn0.resample(mvaB1_Bn0.time)+mvaB3_Bn0.resample(mvaB1_Bn0.time)+mvaB4_Bn0.resample(mvaB1_Bn0.time))/4; mvaAvB_Bn0.name = 'avg B LMN (<BN>=0)';

h = irf_plot(3);

hca = irf_panel('B LMN (unconstr.)');
set(hca,'ColorOrder',mms_colors('xyz'))  
c_eval('irf_plot(hca,{mvaAvB_unc.x.tlim(tint_bcs),mvaAvB_unc.y.tlim(tint_bcs),mvaAvB_unc.z.tlim(tint_bcs)},''comp'');',ic)
hca.YLabel.String = {'B (nT)'};
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.02 0.2],'fontsize',12);
 
hca = irf_panel('B LMN (<BN>=0)');
set(hca,'ColorOrder',mms_colors('xyz'))  
c_eval('irf_plot(hca,{mvaAvB_Bn0.x.tlim(tint_bcs),mvaAvB_Bn0.y.tlim(tint_bcs),mvaAvB_Bn0.z.tlim(tint_bcs)},''comp'');',ic)
hca.YLabel.String = {'B (nT)'};
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L_{B_N = 0}','M_{B_N = 0}','N_{B_N = 0}'},[0.02 0.2],'fontsize',12);
 
hca = irf_panel('remove'); hca.Visible = 'off';

for ih = 1:2
  h(ih).YLim = [-13 13];
  h(ih).YTick = [-10:5:10];
end
irf_zoom(h(1:2),'x',tint_bcs)
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);

%% Timing: B, ve, only the shifted field 
npanels = 5;
h = irf_plot(npanels);
tintZoom = irf.tint('2015-11-12T07:19:19.80Z/2015-11-12T07:19:22.70Z');
tintZoom = irf.tint('2015-11-12T07:19:20.50Z/2015-11-12T07:19:22.00Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.6;
hcf = gcf; hcf.Position = figurePostition;

dt0 = [0 0 0 0];
dt1 = [0.00     0.05     -0.10     -0.08]; v1 = 70*[-0.86,0.37,0.35]; % good for second peak and BL
dt2 = [0.00     0.05     -0.05     -0.03]; v2 = 120*[-0.82 0.49 0.30];% first part

dt = dt2; vv = v2; vvn = vv/norm(vv);
pshift = 0;

hca = irf_panel('invisible'); hca.Visible = 'off'; pshift = pshift + 1;

if 0 % BL 
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BL dt
  hca = irf_panel('BL dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('12341'))  
  irf_legend(hca,{['\Delta t = [ ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4)),'] s'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['v = ' num2str(norm(vv),'%.0f') '\times[' num2str(vvn(1),'%.2f') ', ' num2str(vvn(2),'%.2f') ', ' num2str(vvn(3),'%.2f') '] km/s']},[0.98 0.7],'fontsize',12);
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.02 0.3],'fontsize',12);
end
if 0 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BM dt
  hca = irf_panel('BM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{M}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  %irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 0 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % VeL
  hca = irf_panel('VeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.x.tlim(tint),mvaVe2.x.tlim(tint),mvaVe3.x.tlim(tint),mvaVe4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % VeL dt
  hca = irf_panel('VeL dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.x.tlim(tint),mvaVe2.x.tlim(tint),mvaVe3.x.tlim(tint),mvaVe4.x.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'V_{e,L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('1234'))  
  %irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 0 % VeM
  hca = irf_panel('VeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.y.tlim(tint),mvaVe2.y.tlim(tint),mvaVe3.y.tlim(tint),mvaVe4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,M}','(nT)'};
end
if 1 % VeM dt
  hca = irf_panel('VeM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.y.tlim(tint),mvaVe2.y.tlim(tint),mvaVe3.y.tlim(tint),mvaVe4.y.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'V_{e,M}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  %irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 0 % VeN
  hca = irf_panel('VeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.z.tlim(tint),mvaVe2.z.tlim(tint),mvaVe3.z.tlim(tint),mvaVe4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,N}','(nT)'};
end
if 0 % ne dt
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % VeL
  hca = irf_panel('VeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).abs,mvaVe2perp.tlim(tint).abs,mvaVe3perp.tlim(tint).abs,mvaVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 0 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end


irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

if norm(vv)>100
  htl = add_length_on_top(h(2),norm(vv),0.2);
else
  htl = add_length_on_top(h(2),norm(vv),0.1);
end
%htl.XTickLabel{end} = ['    ' htl.XTickLabel{end} ' km']
htl.XLabel.String = 'km';
htl.FontSize = 14;


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
legends = {'e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
for ii = 1:(npanels-pshift)
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14)
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Timing: B, ve
npanels = 8;
h = irf_plot(npanels);
tintZoom = irf.tint('2015-11-12T07:19:19.80Z/2015-11-12T07:19:22.70Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;


dt1 = [0.00     0.05     -0.10     -0.08]; % good for second peak
dt2 = [0.00     0.05     -0.05     -0.03]; % first part

dt = dt1;

if 1 % BL 
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BL dt
  hca = irf_panel('BL dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BM dt
  hca = irf_panel('BM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{M}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 0 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % VeL
  hca = irf_panel('VeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.x.tlim(tint),mvaVe2.x.tlim(tint),mvaVe3.x.tlim(tint),mvaVe4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % VeL dt
  hca = irf_panel('VeL dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.x.tlim(tint),mvaVe2.x.tlim(tint),mvaVe3.x.tlim(tint),mvaVe4.x.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'V_{e,L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 1 % VeM dt
  hca = irf_panel('VeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.y.tlim(tint),mvaVe2.y.tlim(tint),mvaVe3.y.tlim(tint),mvaVe4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,M}','(nT)'};
end
if 1 % VeM
  hca = irf_panel('VeM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.y.tlim(tint),mvaVe2.y.tlim(tint),mvaVe3.y.tlim(tint),mvaVe4.y.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'V_{e,M}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4))},[0.98 0.9],'fontsize',12);
end
if 0 % VeN
  hca = irf_panel('VeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.z.tlim(tint),mvaVe2.z.tlim(tint),mvaVe3.z.tlim(tint),mvaVe4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,N}','(nT)'};
end

if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % VeL
  hca = irf_panel('VeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).abs,mvaVe2perp.tlim(tint).abs,mvaVe3perp.tlim(tint).abs,mvaVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 0 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end


irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% 2 timings: BLM, veLM
[out_Bn0,l_Bn0,v_Bn0] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z')),'<Bn>=0');
L_Bn0 = v_Bn0(1,:); M_Bn0 = -v_Bn0(2,:); N_Bn0 = -v_Bn0(3,:);
coordLabels_Bn0 = {'L','M','N'};
lmn_Bn0 = [L_Bn0;M_Bn0;N_Bn0];

[out_Bn20,l_Bn20,v_Bn20] = irf_minvar(gseB1.tlim(irf.tint('2015-11-12T07:19:20.116Z/2015-11-12T07:19:22.136Z')),'td');
L_Bn20 = v_Bn20(1,:); M_Bn20 = -v_Bn20(2,:); N_Bn20 = -v_Bn20(3,:);
coordLabels_Bn20 = {'L','M','N'};
lmn_Bn20 = [L_Bn20;M_Bn20;N_Bn20];

npanels = 8;
pshift = 3;
h = irf_plot(npanels + pshift);
tintZoom = irf.tint('2015-11-12T07:19:19.80Z/2015-11-12T07:19:22.70Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;


%[0.00      0.05     -0.13     -0.11]

dt1 = [0.00     0.05     -0.10     -0.08]; % good for second peak
%dt1 = [0.00      0.05     -0.13     -0.11];
dt2 = [0.00     0.03     -0.05     -0.03]; % first part
dt2 = dt2*0;

v1 = 71*[-0.13 0.24 -0.96]; n1 = v1/norm(v1);
v2 = 120*[-0.24 0.16 -0.96]; n2 = v2/norm(v2);
n_Bn0 = N_Bn0*lmn';
n_Bn20 = N_Bn20*lmn';
n_ = lmn(3,:)*lmn';


hca = irf_panel('delete1'); hca.Visible = 'off';
hca = irf_panel('delete2'); hca.Visible = 'off';
hca = irf_panel('delete3'); hca.Visible = 'off';

if 1 % BL dt
  hca = irf_panel('BL dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp','dt',dt2);
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.02 0.6],'fontsize',12);
  irf_legend(hca,{['dt = ' num2str(dt2(1))],num2str(dt2(2)),num2str(dt2(3)),num2str(dt2(4))},[0.98 0.9],'fontsize',12);
end
if 1 % BL 
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp','dt',dt1);
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{['dt = ' num2str(dt1(1))],num2str(dt1(2)),num2str(dt1(3)),num2str(dt1(4))},[0.98 0.9],'fontsize',12);
end

if 1 % BM dt
  hca = irf_panel('BM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt2);
  hca.YLabel.String = {'B_{M}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt2(1))],num2str(dt2(2)),num2str(dt2(3)),num2str(dt2(4))},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt1);
  hca.YLabel.String = {'B_{M}','(nT)'};
  irf_legend(hca,{['dt = ' num2str(dt1(1))],num2str(dt1(2)),num2str(dt1(3)),num2str(dt1(4))},[0.98 0.9],'fontsize',12);
end

if 0 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp','dt',dt1);
  hca.YLabel.String = {'B_{N}','(nT)'};
  irf_legend(hca,{['dt = ' num2str(dt1(1))],num2str(dt1(2)),num2str(dt1(3)),num2str(dt1(4))},[0.98 0.9],'fontsize',12);
end

if 1 % VeL dt
  hca = irf_panel('VeL dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.x.tlim(tint),mvaVe2.x.tlim(tint),mvaVe3.x.tlim(tint),mvaVe4.x.tlim(tint)},'comp','dt',dt2);
  hca.YLabel.String = {'V_{e,L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt2(1))],num2str(dt2(2)),num2str(dt2(3)),num2str(dt2(4))},[0.98 0.9],'fontsize',12);
end
if 1 % VeL
  hca = irf_panel('VeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.x.tlim(tint),mvaVe2.x.tlim(tint),mvaVe3.x.tlim(tint),mvaVe4.x.tlim(tint)},'comp','dt',dt1);
  hca.YLabel.String = {'V_{e,L}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{['dt = ' num2str(dt1(1))],num2str(dt1(2)),num2str(dt1(3)),num2str(dt1(4))},[0.98 0.9],'fontsize',12);
end

if 1 % VeM
  hca = irf_panel('VeM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.y.tlim(tint),mvaVe2.y.tlim(tint),mvaVe3.y.tlim(tint),mvaVe4.y.tlim(tint)},'comp','dt',dt2);
  hca.YLabel.String = {'V_{e,M}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  irf_legend(hca,{['dt = ' num2str(dt2(1))],num2str(dt2(2)),num2str(dt2(3)),num2str(dt2(4))},[0.98 0.9],'fontsize',12);
end
if 1 % VeM dt
  hca = irf_panel('VeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.y.tlim(tint),mvaVe2.y.tlim(tint),mvaVe3.y.tlim(tint),mvaVe4.y.tlim(tint)},'comp','dt',dt1);
  hca.YLabel.String = {'V_{e,M}','(nT)'};
  irf_legend(hca,{['dt = ' num2str(dt1(1))],num2str(dt1(2)),num2str(dt1(3)),num2str(dt1(4))},[0.98 0.9],'fontsize',12);
end

if 0 % VeN
  hca = irf_panel('VeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.z.tlim(tint),mvaVe2.z.tlim(tint),mvaVe3.z.tlim(tint),mvaVe4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'V_{e,N}','(nT)'};
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 2; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% 2 timings: Plot sc positions
if exist('h2'); delete(h2); end
nrows = 4;
ncols = 3;
h2(1) = subplot(nrows,ncols,1); %h2(1).Position(2) = h2(1).Position(2)+0.02;
h2(2) = subplot(nrows,ncols,2); %h2(2).Position(2) = h2(2).Position(2)+0.02;
h2(3) = subplot(nrows,ncols,3); %h2(2).Position(2) = h2(2).Position(2)+0.02;

mms_marker={{'ks','markersize',10},{'rd','markersize',10},...
	{'go','markersize',10,'color',[0 0.6 0]},{'bv','markersize',10}};
mms_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
	{'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
mms_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
	{'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};
sc_list = 1:4;

x = {mvaRR1,mvaRR2,mvaRR3,mvaRR4};

hca = h2(1);
hold(hca,'on');
for ic=1:4
  % put Cluster markers  
  plot(hca,x{ic}(3),x{ic}(1),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.0*[-10 10];
hca.YLim = 1.0*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'L (km)';

hca = h2(2);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.0*[-10 10];
hca.YLim = 1.0*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'M (km)';
%hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
%hleg.Box = 'off';
%hleg.Position(1) = hleg.Position(1)+0.2;

hca = h2(3);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(2),x{ic}(1),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.0*[-10 10];
hca.YLim = 1.0*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'M (km)';
hca.YLabel.String = 'L (km)';
%hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
%hleg.Box = 'off';
%hleg.Position(1) = hleg.Position(1)+0.2;


% plot planes
x = 50*[-1 1];
y = 50*[-1 1];
z = 50*[-1 1];

funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);

funL = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funM = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funN = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);

if exist('hmshN1'); delete(hmshN1); end
if exist('hmshN2'); delete(hmshN2); end
if exist('hmspN1'); delete(hmspN1); end
if exist('hmspN2'); delete(hmspN2); end

hca = h2(1); subplot(nrows,ncols,1);
hold(hca,'on');
quiver(hca,0,0,n1(3),n1(1),4,'tag','quiver','color',[0 0 0])
quiver(hca,0,0,n2(3),n2(1),4,'tag','quiver','color',[0 0 0])
quiver(hca,0,0,n_Bn0(3),n_Bn0(1),4,'tag','quiver','color',[1 0 0])
quiver(hca,0,0,n_Bn20(3),n_Bn20(1),4,'tag','quiver','color',[0 1 0])


plotN = hca.XLim;
plotL1 = funX(z,[0 0],n1);
plotL2 = funX(z,[0 0],n2);
%hmshN1 = plot(hca,z,plotL1,'k-.');
%hmspN1 = plot(hca,z,plotL2,'k-');
hmshN1 = plot(hca,z+0,funX([0 0],z,n1),'k-.');
hmspN1 = plot(hca,z-0,funX([0 0],z,n2),'k-');


%legend([hmshN1 hmspN1],{'v1','v2'},'location','eastoutside')
ht = text(9.5,9,'MSH'); ht.Rotation = 55; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;
ht = text(-11.2,8,'MSP'); ht.Rotation = -80; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;

hold(hca,'off');

y = [-10 10];
hca = h2(2);
hold(hca,'on');
quiver(hca,0,0,n1(3),n1(2),4,'tag','quiver','color',[0 0 0])
quiver(hca,0,0,n2(3),n2(2),4,'tag','quiver','color',[0 0 0])
quiver(hca,0,0,n_Bn0(3),n_Bn0(2),4,'tag','quiver','color',[1 0 0])
quiver(hca,0,0,n_Bn20(3),n_Bn20(2),4,'tag','quiver','color',[0 1 0])

hmshN1 = plot(hca,z+0,funY([0 0],z,n1),'k-.');
hmspN2 = plot(hca,z-0,funY([0 0],z,n2),'k-');
hold(hca,'off');


hca = h2(3);
hold(hca,'on');
quiver(hca,0,0,n1(2),n1(1),4,'tag','quiver','color',[0 0 0])
quiver(hca,0,0,n2(2),n2(1),4,'tag','quiver','color',[0 0 0])
quiver(hca,0,0,n_Bn0(2),n_Bn0(1),4,'tag','quiver','color',[1 0 0])
quiver(hca,0,0,n_Bn20(2),n_Bn20(1),4,'tag','quiver','color',[0 1 0])

hmshN1 = plot(hca,z+0,funX(z,[0 0],n1),'k-.');
hmspN2 = plot(hca,z-0,funX(z,[0 0],n2),'k-');
hold(hca,'off');


%h2(1).Position(2) = h2(1).Position(2)+0.02;
%h2(2).Position(2) = h2(2).Position(2)+0.02;
irf_legend(h2(1),'a)',[-0.15 1.1],'color',[0 0 0])
irf_legend(h2(2),'b)',[-0.15 1.1],'color',[0 0 0])
irf_legend(h2(3),'c)',[-0.15 1.1],'color',[0 0 0])
%
for ii = 1:3;
  h2(ii).YLabel.FontSize = 12;
  h2(ii).XLabel.FontSize = 12;
  h2(ii).FontSize = 12;
  h2(ii).Position(2) = h2(ii).Position(2)+0.03;
end

%% Figure: eletron pa plots, mva fields
npanels = 10;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)
tref = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % V ExB
  hca = irf_panel('VExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x.tlim(tint),mvaVExB?.y.tlim(tint),mvaVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = 0.5*[-1100 1100];  
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
  %hca.YTick
end
if 1 % eDist omni 64 par
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
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
  if 1 % plot perpendicular energy baed on magnetic moment conservation
    tref4 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength4 = [0 -0.30];
    tref3 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength3 = [0 -0.35];
    tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
    tref1 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength1 = [0 -0.40];
    c_eval('tref = tref?;',ic)
    c_eval('tintB = tref+tlength?;',ic)
    tintB = tintZoom;
    c_eval('muref = 0.5*units.me*mvaVe?perp.abs2.resample(tref).data*10^6/units.eV/mvaB?.abs.resample(tref).data;',ic)
    c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*muref);')
    hold(hca,'on')
    c_eval('lineWperp = irf_plot(hca,Wperp,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 apar
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all low E');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all mid E');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [100 300];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all mid E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [20 500];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % potential in normal direction
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaIntE?N.tlim(tint)*80},''comp'');',ic)
  hca.YLabel.String = {'\phi','(mV/m s km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % length scales
  hca = irf_panel('length scales');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{re?,Le?},''comp'');',ic)
  hca.YLabel.String = {'length scales','(km)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'r_e','L_e'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:4 10:end]),'y')
for  ii = 5:7; h(ii).CLim = [5 8.5]; end
for  ii = 5:7; h(ii).YLim = [10 1000]; end
irf_plot_axis_align

hca = irf_panel('B'); hca.YLim = [-14 14];
hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.2 8.2];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.2 8.2];

h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

hca = irf_panel('gradPe'); hca.YLim = [-2.5 2.5];
hca = irf_panel('E'); hca.YLim = [-3 3];
irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Figure: Forces, pressures, and ohms law
 
ic = 1;
npanels = 8;
h = irf_plot(npanels);
pshift=0;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 0 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 0 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.abs.tlim(tint),mvaB2.abs.tlim(tint),mvaB3.abs.tlim(tint),mvaB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % PB
  hca = irf_panel('Pressure, B');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{PB1.tlim(tint),PB2.tlim(tint),PB3.tlim(tint),PB4.tlim(tint)},'comp');
  hca.YLabel.String = {'P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Pe
  hca = irf_panel('Pressure, e');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaPe1.trace.tlim(tint)/3,mvaPe2.trace.tlim(tint)/3,mvaPe3.trace.tlim(tint)/3,mvaPe4.trace.tlim(tint)/3},'comp');
  hca.YLabel.String = {'P_e','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Pe perp
  hca = irf_panel('Pressure, e perp');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{0.5*(mvaPe1.yy+mvaPe1.yy),...
                0.5*(mvaPe2.yy+mvaPe2.yy),...    
                0.5*(mvaPe3.yy+mvaPe3.yy),...    
                0.5*(mvaPe4.yy+mvaPe4.yy)},'comp');
  hca.YLabel.String = {'P_{e,\perp}','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Pi
  hca = irf_panel('Pressure, i');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaPi1.trace.tlim(tint)/3,mvaPi2.trace.tlim(tint)/3,mvaPi3.trace.tlim(tint)/3,mvaPi4.trace.tlim(tint)/3},'comp');
  hca.YLabel.String = {'P_i','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Pressures
  hca = irf_panel('Pressure');
  set(hca,'ColorOrder',mms_colors('1234b'))  
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'mvaPi?.trace/3,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 1 % Pressures, e i B
  hca = irf_panel('Pressure, rel');
  set(hca,'ColorOrder',mms_colors('1234bx'))  
  c_eval('refPi = mean(mvaPi?.trace.tlim(tintZoom).data)/3;',ic)
  c_eval('refPi = mvaPi?.trace/3; refPi = refPi.tlim(tintZoom).data(1);',ic)
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-refPi,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234bx'))
  irf_legend(hca,{'P_e',['P_i-' num2str(refPi,'%.2f') ' nPa'],'P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 1 % Pressures, e i B
  icc = ic;
  ic = 3;
  hca = irf_panel('Pressure, rel, mms_');
  set(hca,'ColorOrder',mms_colors('1234bx'))  
  c_eval('refPi = mean(mvaPi?.trace.tlim(tintZoom).data)/3;',ic)
  c_eval('refPi = mvaPi?.trace/3; refPi = refPi.tlim(tintZoom).data(1);',ic)
  refPi = 0.42;
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-refPi,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234bx'))
  irf_legend(hca,{'P_e',['P_i-' num2str(refPi,'%.2f') ' nPa'],'P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
  ic = icc;
end
if 0 % Pressures, e i B, perps incl
  hca = irf_panel('Pressure, rel');
  set(hca,'ColorOrder',mms_colors('1234bx'))  
  c_eval('refPi = mean(mvaPi?.trace.tlim(tintZoom).data)/3;',ic)
  c_eval('refPi = mvaPi?.trace/3; refPi = refPi.tlim(tintZoom).data(1);',ic)
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),0.5*(facPi?.yy+facPi?.zz)-refPi,'...
                         'mvaPi?.trace/3-refPi,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234bx'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_{i\perp}',['P_i-' num2str(refPi,'%.2f') ' nPa'],'P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
  ic = icc;
end
if 0 % Pressures
  hca = irf_panel('Pressure, rel');
  set(hca,'ColorOrder',mms_colors('1234b'))  
  c_eval('refPi = mean(mvaPi?.trace.tlim(tintZoom).data)/3;',ic)
  c_eval('refPi = mvaPi?.trace/3; refPi = refPi.tlim(tintZoom).data(1);',ic)
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'mvaPi?.trace/3-refPi,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}',['P_i-' num2str(refPi,'%.2f') ' nPa'],'P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 0 % Pressures
  hca = irf_panel('Pressure, 4sc');
  set(hca,'ColorOrder',mms_colors('1234'))  
  if 1
  irf_plot(hca,{mvaPe1.trace/3+PB1.resample(mvaPe1)+mvaPi1.resample(mvaPe1).trace/3,...    
                mvaPe2.trace/3+PB2.resample(mvaPe2)+mvaPi1.resample(mvaPe2).trace/3,...    
                mvaPe3.trace/3+PB3.resample(mvaPe3)+mvaPi1.resample(mvaPe3).trace/3,...    
                mvaPe4.trace/3+PB4.resample(mvaPe4)+mvaPi1.resample(mvaPe4).trace/3},'comp');
  else
    irf_plot(hca,{0.5*(mvaPe1.yy+mvaPe1.yy)+PB1.resample(mvaPe1),...    
                  0.5*(mvaPe2.yy+mvaPe2.yy)+PB2.resample(mvaPe2),...    
                  0.5*(mvaPe3.yy+mvaPe3.yy)+PB3.resample(mvaPe3),...    
                  0.5*(mvaPe4.yy+mvaPe4.yy)+PB4.resample(mvaPe4)},'comp');
  end
  hca.YLabel.String = {'P_e+P_i+P_B','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).abs,mvaVe2perp.tlim(tint).abs,mvaVe3perp.tlim(tint).abs,mvaVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 0 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % Ohm's law x
  hca = irf_panel('Ohm L');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x,mvaAvE.x+mvaOhmVexB.resample(mvaAvE).x+mvaOhmGradPe.resample(mvaAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_L','(mV/m)'};
end
if 0 % Ohm's law y
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y,-mvaOhmVixB.y,mvaOhmJxB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law z
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z,-mvaOhmVixB.z,mvaOhmJxB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law residual 3 components
  hca = irf_panel('Ohm residual');
  set(hca,'ColorOrder',mms_colors('xyz'))
  residual = mvaAvE+mvaOhmVexB.resample(mvaAvE)+mvaOhmGradPe.resample(mvaAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  hca.YLabel.String = {'E+v_exB+\nabla P_e','(mV/m)'};
  ylabel({'E+v_exB+\nabla P_e','(mV/m)'},'interpreter','tex');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end


irf_zoom(h,'x',tintZoom)
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

%% Figure: Crescent distributions
% Initialize figure
h1 = irf_plot(12);
cmap = 'jet';
ic = 3;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [20 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  hca.CLim = [7.3 8.3];
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all mid E');
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
  %hca.CLim = [7.3 8.3];
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1]),'y')
for  ii = 6:8; h1(ii).CLim = [5 8.5]; end
for  ii = 6:8; h1(ii).YLim = [10 1000]; end
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

h1 = h1(1:3);

ipanel = 0;
pshift = 0; legshift = 0;
for ii = 1:3
  ipanel = ipanel + 1
  irf_legend(h1(ipanel),legends{ipanel},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h1(ipanel).FontSize = 14;  
  h1(ipanel).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

%% plot crescents
if exist('hmark'); delete(hmark); end
c_eval('times = ePDist?.time;',ic)
tind = 902:1:(960+9*2);
tind_mms1 = [903 905 907 909 911 913];
tind_mms3 = [900 901 902 907 908 909];
tind_mms3 = [900 907 901 908 902 909]; % sorted into columns
legends = {'a)','b)','c)','d)','g)','e)','h)','f)','i)'};
c_eval('tind = tind_mms?;',ic)

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  
nrows = 3;
ncols = 2;
rowshift = 1;

for ii = 1:nrows*ncols
  h2(ii) = subplot(rowshift+nrows,ncols,rowshift*ncols+ii);
end

isub = 1;
for ii=1:nrows*ncols  
  it = tind(ii);
  time = times(it);
  timeUTC = time.utc; 
  hmark(ii) = irf_pl_mark(h1(1),time.epochUnix','green');
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
  x = cross(y,z)/norm(cross(y,z)); % ExB
  y = cross(z,x); % Bx(ExB)
  perp1 = x;
  perp2 = y;
  par = z;

  %x = L;
  %y = M;
  %z = N;

  if 1
    %vectors = {hatExB,'ExB'; hatE,'E';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
    %vectors = {[1 0 0],'x';[0 1 0],'y';[0 0 1],'z';L,'L';M,'M';N,'N'};
    vectors = {L,'L';M,'M';N,'N'};
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_E','v_B'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);            
    hca.Title.String = '';
    hca.Title.String = timeUTC(12:23);
    hca.FontWeight = 'normal';
    hca.FontSize = 12;
    colormap(hca,strCMap);
  end
  
  if 0
    hca = h2(isub); isub = isub + 1; 
    set(hca,'ColorOrder',colors)
    indPA = [1 7 8 12];
    c_eval('plot(hca,ePitch?(it).depend{1},squeeze(ePitch?(it).convertto(''s^3/km^6'').data(1,:,indPA)))',ic)
    hca.YScale = 'log'; hca.XScale = 'log';
    hca.YLabel.String = 'f_e (s^3/km^6)'; hca.XLabel.String = 'E (eV)';
    hca.YLim = [1e-1 1e5]; hca.XLim = [1e1 1e3];
    hca.YTick = [1e-1 1e0 1e1 1e2 1e3 1e4 1e5];
    %legend(hca,num2str(ePitch3.depend{2}(indPA)'),'location','southwest')
    hca.Title.String = timeUTC(12:23);
  end
  
end
ipanel = 0;
for ii = 1:(nrows*ncols)
  originalPosition{ii} = h2(ii).Position;  
  h2(ii).XGrid = 'on';
  h2(ii).YGrid = 'on';
  ipanel = ipanel + 1;
  irf_legend(h2(ipanel),legends{ipanel+3},[0.01 0.95],'color',[0 0 0],'fontsize',14);  
end

irf_zoom(h1,'x',tintZoom)
h1(3).XLabel.String = '';


hcf = gcf;
hCB = findall(hcf,'type','ColorBar'); 
delete(hCB(2:end)); hCB = hCB(1);
hCBx = hCB.Position(1);
hCB.Position(1) = hCBx+0.01;
xyfactor = 1.2;
dy = 0.055;
leftdx = -0.04;
rightdx = 0.05;

ip = 1;
h2(ip).Position(1) = [originalPosition{ip}(1)-leftdx];
h2(ip).Position(2) = [originalPosition{ip}(2)-dy];
h2(ip).Position(3) = [originalPosition{ip}(3)-00];
h2(ip).Position([3 4]) = [originalPosition{ip}([3 4])*xyfactor];
h2(ip).XLabel = []; h2(ip).XTick = [];

ip = 2;
h2(ip).Position(1) = [originalPosition{ip}(1)-rightdx];
h2(ip).Position(2) = [originalPosition{ip}(2)-dy];
h2(ip).Position(3) = [originalPosition{ip}(3)+00];
h2(ip).Position([3 4]) = [originalPosition{ip}([3 4])*xyfactor];
h2(ip).XLabel = []; h2(ip).XTick = [];
h2(ip).YLabel = []; h2(ip).YTick = [];

ip = 3;
h2(ip).Position(1) = [originalPosition{ip}(1)-leftdx];
h2(ip).Position(2) = [originalPosition{ip}(2)-dy];
h2(ip).Position(3) = [originalPosition{ip}(3)+00];
h2(ip).Position([3 4]) = [originalPosition{ip}([3 4])*xyfactor];
h2(ip).XLabel = []; h2(ip).XTick = [];

ip = 4;
h2(ip).Position(1) = [originalPosition{ip}(1)-rightdx];
h2(ip).Position(2) = [originalPosition{ip}(2)-dy];
h2(ip).Position(3) = [originalPosition{ip}(3)+00];
h2(ip).Position([3 4]) = [originalPosition{ip}([3 4])*xyfactor];
h2(ip).XLabel = []; h2(ip).XTick = [];
h2(ip).YLabel = []; h2(ip).YTick = [];

ip = 5;
h2(ip).Position(1) = [originalPosition{ip}(1)-leftdx];
h2(ip).Position(2) = [originalPosition{ip}(2)-dy];
h2(ip).Position(3) = [originalPosition{ip}(3)+00];
h2(ip).Position([3 4]) = [originalPosition{ip}([3 4])*xyfactor];

ip = 6;
h2(ip).Position(1) = [originalPosition{ip}(1)-rightdx];
h2(ip).Position(2) = [originalPosition{ip}(2)-dy];
h2(ip).Position(3) = [originalPosition{ip}(3)+00];
h2(ip).Position([3 4]) = [originalPosition{ip}([3 4])*xyfactor];
h2(ip).YLabel = []; h2(ip).YTick = [];

%% Figure 2: 4 sc, MVA
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:22.50Z');
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.50Z');
npanels = 8;
pshift = 4;
h = irf_plot(npanels+pshift);


scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

c_eval('hca = irf_panel(''delete1?''); hca.Visible = ''off'';',1:pshift)

if 1 % BX
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');      
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp N
  hca = irf_panel('Ve perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).z,mvaVe2perp.tlim(tint).z,mvaVe3perp.tlim(tint).z,mvaVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,N}','(km/s)'},'interpreter','tex');
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{mvaEVexB1.abs,mvaEVexB2.abs,mvaEVexB3.abs,mvaEVexB4.abs},'comp'); 
  hca.YLabel.String = {'|E+v_{e}xB|','(mV/m)'};
  %hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  %hca.YLim = [-10 10];
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P',['(10^{-3} ' mvaGradPe.units ')']};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end


irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 3; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];
hl = add_length_on_top(h(1+pshift),70,0.1);
hl.XLabel.String = 'km';
%% 2 timings: Plot sc positions
if exist('h2'); delete(h2); end
nrows = 4;
ncols = 3;
h2(1) = subplot(nrows,ncols,1); %h2(1).Position(2) = h2(1).Position(2)+0.02;
h2(2) = subplot(nrows,ncols,2); %h2(2).Position(2) = h2(2).Position(2)+0.02;
h2(3) = subplot(nrows,ncols,3); %h2(2).Position(2) = h2(2).Position(2)+0.02;

mms_marker={{'ks','markersize',10},{'rd','markersize',10},...
	{'go','markersize',10,'color',[0 0.6 0]},{'bv','markersize',10}};
mms_marker_small={{'ks','markersize',8},{'rd','markersize',8},...
	{'go','markersize',8,'color',[0 0.6 0]},{'bv','markersize',8}};
mms_marker_shaded={{'ks','color',[0.3 0.3 0.3]},...
	{'rd','color',[1 0.3 0.3]},{'go','color',[.3 1 .3]},{'bv','color',[.3 .3 1]}};
sc_list = 1:4;

x = {mvaRR1,mvaRR2,mvaRR3,mvaRR4};

hca = h2(1);
hold(hca,'on');
for ic=1:4
  % put Cluster markers  
  plot(hca,x{ic}(3),x{ic}(1),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.0*[-10 10];
hca.YLim = 1.0*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'L (km)';

hca = h2(2);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.0*[-10 10];
hca.YLim = 1.0*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'M (km)';
%hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
%hleg.Box = 'off';
%hleg.Position(1) = hleg.Position(1)+0.2;

hca = h2(3);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(2),x{ic}(1),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.0*[-10 10];
hca.YLim = 1.0*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'M (km)';
hca.YLabel.String = 'L (km)';
%hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
%hleg.Box = 'off';
%hleg.Position(1) = hleg.Position(1)+0.2;


% plot planes
x = 50*[-1 1];
y = 50*[-1 1];
z = 50*[-1 1];

funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);

funL = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funM = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funN = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);

if exist('hmshN1'); delete(hmshN1); end
if exist('hmshN2'); delete(hmshN2); end
if exist('hmspN1'); delete(hmspN1); end
if exist('hmspN2'); delete(hmspN2); end

hca = h2(1); subplot(nrows,ncols,1);
hold(hca,'on');
quiver(hca,0,0,n1(3),n1(1),4,'tag','quiver','color',[0 0 0])
%quiver(hca,0,0,n2(3),n2(1),4,'tag','quiver','color',[0 0 0])
%quiver(hca,0,0,n_Bn0(3),n_Bn0(1),4,'tag','quiver','color',[1 0 0])
%quiver(hca,0,0,n_Bn20(3),n_Bn20(1),4,'tag','quiver','color',[0 1 0])


plotN = hca.XLim;
plotL1 = funX(z,[0 0],n1);
%plotL2 = funX(z,[0 0],n2);
%hmshN1 = plot(hca,z,plotL1,'k-.');
%hmspN1 = plot(hca,z,plotL2,'k-');
hmshN1 = plot(hca,z+0,funX([0 0],z,n1),'k-.');
%hmspN1 = plot(hca,z-0,funX([0 0],z,n2),'k-');


%legend([hmshN1 hmspN1],{'v1','v2'},'location','eastoutside')
%ht = text(9.5,9,'<- MSH'); ht.Rotation = 0; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;
%ht = text(-11.2,8,'MSP ->'); ht.Rotation = 0; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;
irf_legend(h2(1),'<- MSH',[0.01 0.99],'fontsize',11)
irf_legend(h2(1),'MSP ->',[0.99 0.99],'fontsize',11)

hold(hca,'off');

y = [-10 10];
hca = h2(2);
hold(hca,'on');
quiver(hca,0,0,n1(3),n1(2),4,'tag','quiver','color',[0 0 0])
%quiver(hca,0,0,n2(3),n2(2),4,'tag','quiver','color',[0 0 0])
%quiver(hca,0,0,n_Bn0(3),n_Bn0(2),4,'tag','quiver','color',[1 0 0])
%quiver(hca,0,0,n_Bn20(3),n_Bn20(2),4,'tag','quiver','color',[0 1 0])

hmshN1 = plot(hca,z+0,funY([0 0],z,n1),'k-.');
%hmspN2 = plot(hca,z-0,funY([0 0],z,n2),'k-');
hold(hca,'off');


hca = h2(3);
hold(hca,'on');
quiver(hca,0,0,n1(2),n1(1),4,'tag','quiver','color',[0 0 0])
%quiver(hca,0,0,n2(2),n2(1),4,'tag','quiver','color',[0 0 0])
%quiver(hca,0,0,n_Bn0(2),n_Bn0(1),4,'tag','quiver','color',[1 0 0])
%quiver(hca,0,0,n_Bn20(2),n_Bn20(1),4,'tag','quiver','color',[0 1 0])

hmshN1 = plot(hca,z+0,funX(z,[0 0],n1),'k-.');
%hmspN2 = plot(hca,z-0,funX(z,[0 0],n2),'k-');
hold(hca,'off');


%h2(1).Position(2) = h2(1).Position(2)+0.02;
%h2(2).Position(2) = h2(2).Position(2)+0.02;
irf_legend(h2(1),'a)',[-0.15 1.1],'color',[0 0 0])
irf_legend(h2(2),'b)',[-0.15 1.1],'color',[0 0 0])
irf_legend(h2(3),'c)',[-0.15 1.1],'color',[0 0 0])
%
for ii = 1:3;
  h2(ii).YLabel.FontSize = 12;
  h2(ii).XLabel.FontSize = 12;
  h2(ii).FontSize = 12;
  h2(ii).Position(2) = h2(ii).Position(2)+0.03;
end
%% Figure 2: 4 sc: eletron pa
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:22.50Z');
npanels = 8;
h = irf_plot(npanels);

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = colormap('jet');

if 1 % BX
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end

elim = [10 400];
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E mms1');
  %elim = [20 400];
  try
    irf_spectrogram(hca,ePitch1.tlim(tint).elim(elim).deflux.specrec('pa'),'log');
  catch
    irf_spectrogram(hca,ePDist1.tlim(tint).pitchangles(dmpaB1,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E mms2');
  %elim = [20 400];
  try
    irf_spectrogram(hca,ePitch2.tlim(tint).elim(elim).deflux.specrec('pa'),'log');
  catch
    irf_spectrogram(hca,ePDist2.tlim(tint).pitchangles(dmpaB2,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E mms3');
  %elim = [20 400];
  try
    irf_spectrogram(hca,ePitch3.tlim(tint).elim(elim).deflux.specrec('pa'),'log');
  catch
    irf_spectrogram(hca,ePDist3.tlim(tint).pitchangles(dmpaB3,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E mms4');
  %elim = [20 400];
  try
    irf_spectrogram(hca,ePitch4.tlim(tint).elim(elim).deflux.specrec('pa'),'log');
  catch
    irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:4]),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

for ii = 4:8, h(ii).CLim = [7.3 8.2]; end
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];
 
%% Figure: Electron adiabaticity
npanels = 9;
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';
Eref = [40 400]; colorEref = [0 0 0];
Efactor = 2;
tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)
tref = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)

% prepare energy levels for adiabaticity
tref4 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength4 = [0 -0.30];
tref3 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength3 = [0 -0.35];
tref2 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength2 = [-0.50 0];
tref1 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength1 = [0 -0.40];
c_eval('tref = tref?;',ic)
c_eval('tintB = tref+tlength?;',ic)
tintB = tintZoom;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp)',ic)
Wrefperp = 100;
Wrefpar = 100;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wpar = irf.ts_scalar(Wperp.time,Wreftot-Wperp.data);',ic)
Efactorpar = 1;
Efactorperp = 1;

zoomE = [];
zoomPA = [];
zoomY = [];
isub = 0;
if 1 % B
  hca = irf_panel('B'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 1 % magnetic moment
  hca = irf_panel('mu'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{mag_mom?},''comp'');',ic)
  c_eval('mag_mom_units = mag_mom?.units;',ic)
  hca.YLabel.String = {'\mu',['(' mag_mom_units ')']};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 1 % length scales
  hca = irf_panel('length scales');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('123'))
  if 1 % estimate magnetic length scale
    c_eval('LB? = abs(mvaB?.x.resample(mvaJ1)*1e9/(mvaJ?.y*1e9*units.mu0))*1e-3; Lm?.units = ''km'';',ic)
    c_eval('irf_plot(hca,{re?,Le?,LB?},''comp'');',ic)
    irf_legend(hca,{'r_e','L_e','L_B'},[0.95 0.9],'fontsize',12);
  else
    c_eval('irf_plot(hca,{re?,Le?},''comp'');',ic)
    irf_legend(hca,{'r_e','L_e'},[0.95 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'length','(km)'};
  set(hca,'ColorOrder',mms_colors('123'))
  
  %hca.YLim = [0 12];
end
if 0 % Te par perp
  hca = irf_panel('Te');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
  %hca.YTick
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation   
    %%
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 2'); isub = isub + 1;  zoomE = [zoomE isub];    
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    %c_eval('lineWperp = irf_plot(hca,Wperp,''k'');',ic)
    
    c_eval('lineWperp = irf_plot(hca,Wperp*Efactorperp,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorperp)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
      
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
  if 1 % plot reference energy level
    hold(hca,'on')    
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactor)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 1 % ePDist pa low energies
  %%
  hca = irf_panel('e PA e32 deflux all low E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));',ic)
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));',ic)
  elim = [10 300];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
    
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa mid energies
  hca = irf_panel('e PA e32 deflux all mid E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [100 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all mid E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [20 500];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist Par/Apar
  hca = irf_panel('e Anis e32 par apar');  isub = isub + 1;  zoomE = [zoomE isub];
  try
    c_eval('plotAnis = eAnis?_app;',ic)
    %plotAnis.data(plotAnis.data<1e-5) = NaN;
    plotAnis.data(plotAnis.data>1e4) = NaN;
    c_eval('[hc,hcb] = irf_spectrogram(hca,plotAnis.tlim(tint).specrec(''e''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YScale = 'log';
  hcb.YLabel.String = 'f_{||}/f_{\perp}'
  hca.YLabel.String = {'E_e','(\circ)'};   
  hca.YLim = [10 1000];
  %cc = colormap('jet');  
  colormap(hca,cn.cmap('bluered3'))
  %colormap(hca,'jet')
  hca.CLim = [-4 4];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E new
  hca = irf_panel('E perp par new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp_new.x.tlim(tint),mvaE?perp_new.y.tlim(tint),mvaE?perp_new.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp_new.x.tlim(tint),mvaE?perp_new.y.tlim(tint),mvaE?perp_new.z.tlim(tint),mvaE?par_new.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % E par
  hca = irf_panel('E par');
  try
    set(hca,'ColorOrder',mms_colors('b12'))
    c_eval('irf_plot(hca,{E?par.x,E?par.y,gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('b12'))
    irf_legend(hca,{'error','E_{||}','E_{||}'},[0.98 0.9],'fontsize',12);
  catch
    set(hca,'ColorOrder',mms_colors('2'))
    c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'E_{||}'},[0.98 0.9],'fontsize',12);
  end
end
if 1 % E  
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % E par
  hca = irf_panel('E par');
  try
    set(hca,'ColorOrder',mms_colors('b12'))
    c_eval('irf_plot(hca,{E?par.x,E?par.y,gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('b12'))
    irf_legend(hca,{'error','E_{||}','E_{||}'},[0.98 0.9],'fontsize',12);
  catch
    set(hca,'ColorOrder',mms_colors('2'))
    c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
    hca.YLabel.String = {'E_{||}','(mV/m)'};
    set(hca,'ColorOrder',mms_colors('2'))
    irf_legend(hca,{'E_{||}'},[0.98 0.9],'fontsize',12);
  end
end
if 0 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % potential in normal direction
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaIntE?N.tlim(tint)*80},''comp'');',ic)
  hca.YLabel.String = {'\phi','(mV/m s km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:3 8:end]),'y')
for  ii = zoomE; h(ii).YLim = [10 1000]; end
for  ii = zoomE(1:3); h(ii).CLim = [5.5 8.4]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end
irf_plot_axis_align
%hca = irf_panel('E par'); irf_zoom(hca,'y');
hca = irf_panel('length scales'); hca.YLim = [0 8]; hca.YTick = [2:2:8];
hca = irf_panel('B'); hca.YLim = [-15 18];
hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.4 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.5 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.2 8.2];
hca = irf_panel('E perp par');  hca.YLim = [-4 4];
hca = irf_panel('E par');  hca.YLim = [-3 3];
%hca = irf_panel('e Anis e32 par apar'); hca.YLim = [10 1000]; hca.CLim = [-4 4];

h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];
 
%% Figure: Electron adiabaticity, compare more pitchangle bins.
npanels = 10;
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 3;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';
Eref = [40 400]; colorEref = [0 0 0];
Efactor = 2;
tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)
tref = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)

% prepare energy levels for adiabaticity
tref4 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength4 = [0 -0.30];
tref3 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength3 = [0 -0.35];
tref2 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength2 = [-0.50 0];
tref1 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength1 = [0 -0.40];
c_eval('tref = tref?;',ic)
c_eval('tintB = tref+tlength?;',ic)
tintB = tintZoom;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp)',ic)
Wrefperp = 100;
Wrefpar = 100;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wpar = irf.ts_scalar(Wperp.time,Wreftot-Wperp.data);',ic)
Efactorpar = 1;
Efactorperp = 1;

zoomE = [];
zoomPA = [];
zoomY = [];
isub = 0;
if 1 % B
  hca = irf_panel('B'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 0 % magnetic moment
  hca = irf_panel('mu'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{mag_mom?},''comp'');',ic)
  c_eval('mag_mom_units = mag_mom?.units;',ic)
  hca.YLabel.String = {'\mu',['(' mag_mom_units ')']};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 0 % length scales
  hca = irf_panel('length scales');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('123'))
  if 1 % estimate magnetic length scale
    c_eval('LB? = abs(mvaB?.x.resample(mvaJ1)*1e9/(mvaJ?.y*1e9*units.mu0))*1e-3; Lm?.units = ''km'';',ic)
    c_eval('irf_plot(hca,{re?,Le?,LB?},''comp'');',ic)
    irf_legend(hca,{'r_e','L_e','L_B'},[0.95 0.9],'fontsize',12);
  else
    c_eval('irf_plot(hca,{re?,Le?},''comp'');',ic)
    irf_legend(hca,{'r_e','L_e'},[0.95 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'length','(km)'};
  set(hca,'ColorOrder',mms_colors('123'))
  
  %hca.YLim = [0 12];
end
if 0 % Te par perp
  hca = irf_panel('Te');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
  %hca.YTick
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation   
    %%
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 30-60'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [30 60];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation   
    %%
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 60-90'); isub = isub + 1;  zoomE = [zoomE isub];    
  pas = [60 90];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    %c_eval('lineWperp = irf_plot(hca,Wperp,''k'');',ic)
    
    c_eval('lineWperp = irf_plot(hca,Wperp*Efactorperp,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorperp)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 2'); isub = isub + 1;  zoomE = [zoomE isub];    
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    %c_eval('lineWperp = irf_plot(hca,Wperp,''k'');',ic)
    
    c_eval('lineWperp = irf_plot(hca,Wperp*Efactorperp,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorperp)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 90-120'); isub = isub + 1;  zoomE = [zoomE isub];    
  pas = [90 120];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    %c_eval('lineWperp = irf_plot(hca,Wperp,''k'');',ic)
    
    c_eval('lineWperp = irf_plot(hca,Wperp*Efactorperp,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorperp)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 120-150'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [120 150];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
      
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
  if 1 % plot reference energy level
    hold(hca,'on')    
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactor)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
      
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
  if 1 % plot reference energy level
    hold(hca,'on')    
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactor)},[0.77 0.2])
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 1 % ePDist pa low energies
  %%
  hca = irf_panel('e PA e32 deflux all low E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));',ic)
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));',ic)
  elim = [10 300];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
    
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa mid energies
  hca = irf_panel('e PA e32 deflux all mid E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [100 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all mid E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [20 500];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist Par/Apar
  %%
  hca = irf_panel('e Anis e32 par apar');  isub = isub + 1;  zoomE = [zoomE isub];
  try
    c_eval('plotAnis = eAnis?_app;',ic)
    %plotAnis.data(plotAnis.data<1e-5) = NaN;
    plotAnis.data(plotAnis.data>1e4) = NaN;
    c_eval('[hc,hcb] = irf_spectrogram(hca,plotAnis.tlim(tint).specrec(''e''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YScale = 'log';
  hcb.YLabel.String = 'f_{||}/f_{\perp}'
  %hca.YLabel.String = {'E_e','(\circ)'};   
  hca.YLim = [10 1000];
  %cc = colormap('jet');  
  colormap(hca,cn.cmap('bluered3'))
  %colormap(hca,'jet')
  hca.CLim = [-4 4];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E new
  hca = irf_panel('E perp par new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp_new.x.tlim(tint),mvaE?perp_new.y.tlim(tint),mvaE?perp_new.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp_new.x.tlim(tint),mvaE?perp_new.y.tlim(tint),mvaE?perp_new.z.tlim(tint),mvaE?par_new.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % E  
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:3 8:end]),'y')
for  ii = zoomE; h(ii).YLim = [10 1000]; end
for  ii = zoomE(1:7); h(ii).CLim = [5.5 8.4]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end
irf_plot_axis_align
%hca = irf_panel('E par'); irf_zoom(hca,'y');
hca = irf_panel('length scales'); hca.YLim = [0 8]; hca.YTick = [2:2:8];
hca = irf_panel('B'); hca.YLim = [-15 18];
hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.4 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.5 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.2 8.2];
hca = irf_panel('E perp par');  hca.YLim = [-4 4];
hca = irf_panel('E par');  hca.YLim = [-3 3];
%hca = irf_panel('e Anis e32 par apar'); hca.YLim = [10 1000]; hca.CLim = [-4 4];

h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Combined overview and adiabaticity plot
npanels = 11;
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';
Eref = [40 400]; colorEref = [0 0 0];

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)
tref = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)

% prepare energy levels for adiabaticity
tref4 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength4 = [0 -0.30];
tref3 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength3 = [0 -0.35];
tref2 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength2 = [-0.50 0];
tref1 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength1 = [0 -0.40];
c_eval('tref = tref?;',ic)
c_eval('tintB = tref+tlength?;',ic)
tintB = tintZoom;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp)',ic)
Wrefperp = 100;
Wrefpar = 100;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wpar = irf.ts_scalar(Wperp.time,Wreftot-Wperp.data);',ic)
Efactorpar = 1;
Efactorperp = 1;

zoomE = [];
zoomPA = [];
zoomY = [];
isub = 0;
if 1 % B
  hca = irf_panel('B'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 1 % n
  hca = irf_panel('n'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 1 % Ve
  hca = irf_panel('Ve'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Ve
  hca = irf_panel('Ve perp par'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Te par perp
  hca = irf_panel('Te');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
  %hca.YTick
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 0 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation   
    %%
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    irf_plot(hca,irf.ts_scalar(tref,Wrefpar),'*')
    irf_legend(hca,{'E_{\mu}'},[0.70 0.2])
    hold(hca,'off')    
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 2'); isub = isub + 1;  zoomE = [zoomE isub];    
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    c_eval('lineWperp = irf_plot(hca,Wperp*Efactorperp,''k'');',ic)  
    irf_plot(hca,irf.ts_scalar(tref,Wrefpar),'*')
    irf_legend(hca,{'E_{\mu}'},[0.70 0.2])
    hold(hca,'off')    
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)      
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
  if 0 % plot reference energy level
    hold(hca,'on')    
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    irf_plot(hca,irf.ts_scalar(tref,Wrefpar),'*')
    irf_legend(hca,{'E_{\mu}'},[0.70 0.2])
    hold(hca,'off')        
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 1 % ePDist pa low energies
  %%
  hca = irf_panel('e PA e32 deflux all low E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
    
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa mid energies
  hca = irf_panel('e PA e32 deflux all mid E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [100 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all mid E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [20 500];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist Par/Apar
  hca = irf_panel('e Anis e32 par apar');  isub = isub + 1;  zoomE = [zoomE isub];
  try
    c_eval('plotAnis = eAnis?_app;',ic)
    %plotAnis.data(plotAnis.data<1e-5) = NaN;
    plotAnis.data(plotAnis.data>1e4) = NaN;
    c_eval('[hc,hcb] = irf_spectrogram(hca,plotAnis.tlim(tint).specrec(''e''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YScale = 'log';
  hcb.YLabel.String = 'f_{||}/f_{\perp}'
  hca.YLabel.String = {'E_e','(\circ)'};   
  hca.YLim = [10 1000];
  %cc = colormap('jet');  
  colormap(hca,cn.cmap('bluered3'))
  %colormap(hca,'jet')
  hca.CLim = [-4 4];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % magnetic moment
  hca = irf_panel('mu'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{mag_mom?},''comp'');',ic)
  c_eval('mag_mom_units = mag_mom?.units;',ic)
  hca.YLabel.String = {'\mu',['(' mag_mom_units ')']};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 1 % length scales
  hca = irf_panel('length scales');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('123'))
  if 1 % estimate magnetic length scale
    c_eval('LB? = abs(mvaB?.x.resample(mvaJ1)*1e9/(mvaJ?.y*1e9*units.mu0))*1e-3; Lm?.units = ''km'';',ic)
    c_eval('irf_plot(hca,{re?,Le?,LB?},''comp'');',ic)
    irf_legend(hca,{'r_e','L_e','L_B'},[0.95 0.9],'fontsize',12);
  else
    c_eval('irf_plot(hca,{re?,Le?},''comp'');',ic)
    irf_legend(hca,{'r_e','L_e'},[0.95 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'length','(km)'};
  set(hca,'ColorOrder',mms_colors('123'))
  
  %hca.YLim = [0 12];
end
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % potential in normal direction
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaIntE?N.tlim(tint)*80},''comp'');',ic)
  hca.YLabel.String = {'\phi','(mV/m s km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:end]),'y')
for  ii = zoomE; h(ii).YLim = [10 1000]; end
for  ii = zoomE(1:3); h(ii).CLim = [6 8.2]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end
irf_plot_axis_align
%hca = irf_panel('E par'); irf_zoom(hca,'y');
hca = irf_panel('length scales'); hca.YLim = [0 8]; hca.YTick = [2:2:8];
hca = irf_panel('B'); hca.YLim = [-15 18];
hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.4 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.5 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.2 8.2];
%hca = irf_panel('E perp par');  hca.YLim = [-4 4];
%hca = irf_panel('E par');  hca.YLim = [-3 3];
%hca = irf_panel('e Anis e32 par apar'); hca.YLim = [10 1000]; hca.CLim = [-4 4];

h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

hkm = add_length_on_top(h(1),70,0.1);
labels = hkm.XTickLabel;
labelss = {'','','','','','',labels{1:6},'','','','','','',''}
hkm.XTickLabel = labelss;

h(1).Title.String = '';
hkm.XLabel.String = 'km                      ';
hkm.FontSize = 12;

set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'paperpositionmode','auto');
    set(gcf,'color','white');
    

irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Combined overview and adiabaticity plot, Matsuyama
npanels = 8;
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';
Eref = [40 400]; colorEref = [0 0 0];

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)
tref = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)

% prepare energy levels for adiabaticity
tref4 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength4 = [0 -0.30];
tref3 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength3 = [0 -0.35];
tref2 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength2 = [-0.50 0];
tref1 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength1 = [0 -0.40];
c_eval('tref = tref?;',ic)
c_eval('tintB = tref+tlength?;',ic)
tintB = tintZoom;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp)',ic)
Wrefperp = 100;
Wrefpar = 100;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wpar = irf.ts_scalar(Wperp.time,Wreftot-Wperp.data);',ic)
Efactorpar = 1;
Efactorperp = 1;

zoomE = [];
zoomPA = [];
zoomY = [];
isub = 0;
if 1 % B
  hca = irf_panel('B'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.01 0.2],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 1 % Te par perp
  hca = irf_panel('Te');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [20 120];
  %hca.YTick
end
if 1 % magnetic moment
  hca = irf_panel('mu'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{mag_mom?},''comp'');',ic)
  c_eval('mag_mom_units = mag_mom?.units;',ic)
  hca.YLabel.String = {'\mu',['(' mag_mom_units ')']};
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 0 % Curvature of magnetic field
  hca = irf_panel('curv B'); isub = isub + 1;  
  set(hca,'colororder',mms_colors('xyzb'))
  irf_plot(hca,{mvaCurvB.x,mvaCurvB.y,mvaCurvB.z},'comp');  
  set(hca,'colororder',mms_colors('xyzb'))
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)  
  hca.YLabel.String = {'curv B','(1/km)'};
end
if 0 % n
  hca = irf_panel('n'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 0 % Ve
  hca = irf_panel('Ve'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Ve
  hca = irf_panel('Ve perp par'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 0 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation   
    %%
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    irf_plot(hca,irf.ts_scalar(tref,Wrefpar),'*')
    irf_legend(hca,{'E_{\mu}'},[0.70 0.2])
    hold(hca,'off')    
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 2'); isub = isub + 1;  zoomE = [zoomE isub];    
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
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
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    c_eval('lineWperp = irf_plot(hca,Wperp*Efactorperp,''k'');',ic)  
    irf_plot(hca,irf.ts_scalar(tref,Wrefpar),'*')
    irf_legend(hca,{'E_{\mu}'},[0.70 0.2])
    hold(hca,'off')    
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)      
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
  if 0 % plot reference energy level
    hold(hca,'on')    
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    c_eval('lineWperp = irf_plot(hca,Wpar*Efactorpar,''k'');',ic)  
    irf_plot(hca,irf.ts_scalar(tref,Wrefpar),'*')
    irf_legend(hca,{'E_{\mu}'},[0.70 0.2])
    hold(hca,'off')        
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 1 % ePDist pa low energies
  %%
  hca = irf_panel('e PA e32 deflux all low E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';
    
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa mid energies
  hca = irf_panel('e PA e32 deflux all mid E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [100 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all mid E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [20 500];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist Par/Apar
  hca = irf_panel('e Anis e32 par apar');  isub = isub + 1;  zoomE = [zoomE isub];
  try
    c_eval('plotAnis = eAnis?_app;',ic)
    %plotAnis.data(plotAnis.data<1e-5) = NaN;
    plotAnis.data(plotAnis.data>1e4) = NaN;
    c_eval('[hc,hcb] = irf_spectrogram(hca,plotAnis.tlim(tint).specrec(''e''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YScale = 'log';
  hcb.YLabel.String = 'f_{||}/f_{\perp}'
  hca.YLabel.String = {'E_e','(\circ)'};   
  hca.YLim = [10 1000];
  %cc = colormap('jet');  
  colormap(hca,cn.cmap('bluered3'))
  %colormap(hca,'jet')
  hca.CLim = [-4 4];
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  try
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),E?par.y.tlim(tint)},''comp'');',ic)
  catch
    c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par.tlim(tint)},''comp'');',ic)
  end
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % potential in normal direction
  hca = irf_panel('phi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaIntE?N.tlim(tint)*80},''comp'');',ic)
  hca.YLabel.String = {'\phi','(mV/m s km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % length scales
  hca = irf_panel('length scales');  isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  if 1 % estimate magnetic length scale
    c_eval('LB? = abs(mvaB?.x.resample(mvaJ1)*1e9/(mvaJ?.y*1e9*units.mu0))*1e-3; Lm?.units = ''km'';',ic)
    c_eval('irf_plot(hca,{re?,Le?,LB?,curvBradius},''comp'');',ic)
    set(hca,'ColorOrder',mms_colors('1234'))
    %irf_legend(hca,{'r_e','L_e','L_B','R_c'},[0.95 0.9],'fontsize',12);
    irf_legend(hca,{'r_e','L_e','L_B'},[0.95 0.9],'fontsize',12);
  else
    c_eval('irf_plot(hca,{re?,Le?},''comp'');',ic)
    set(hca,'ColorOrder',mms_colors('1234'))
    irf_legend(hca,{'r_e','L_e'},[0.95 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'length','(km)'};
  set(hca,'ColorOrder',mms_colors('123'))
  
  hca.YLim = [0 12];
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:end]),'y')
for  ii = zoomE; h(ii).YLim = [10 1000]; end
for  ii = zoomE(1:3); h(ii).CLim = [6 8.2]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end
irf_plot_axis_align
%hca = irf_panel('E par'); irf_zoom(hca,'y');
hca = irf_panel('length scales'); hca.YLim = [0 8]; hca.YTick = [2:2:8];
hca = irf_panel('B'); hca.YLim = [-15 18];
hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.4 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.5 8.1];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.2 8.2];
%hca = irf_panel('E perp par');  hca.YLim = [-4 4];
%hca = irf_panel('E par');  hca.YLim = [-3 3];
%hca = irf_panel('e Anis e32 par apar'); hca.YLim = [10 1000]; hca.CLim = [-4 4];

h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

hkm = add_length_on_top(h(1),70,0.1);
labels = hkm.XTickLabel;
labelss = {'','','','','','',labels{1:6},'','','','','','',''}
hkm.XTickLabel = labelss;

h(1).Title.String = '';
hkm.XLabel.String = 'km                      ';
hkm.FontSize = 12;

set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'paperpositionmode','auto');
    set(gcf,'color','white');
    

irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Figure: the miniature current sheet
npanels = 4;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
tintZoomTiny = irf.tint('2015-11-12T07:19:22.00Z/2015-11-12T07:19:23.00Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 3;
pshift = 0;
%scrsz = get(groot,'ScreenSize');
%figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
%hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

v_tinycs = 38*[0.89,0.21,-0.41];
n_tinycs = v_tinycs/norm(v_tinycs);
dt = [0.00      0.21      0.28      0.33];
c_eval('dt = dt - dt(?);',ic)

hca = irf_panel('remove'); hca.Visible = 'off'; h = h(2:4);
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp','dt',dt);      
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
  irf_legend(hca,{['dt = ' num2str(dt(1),'%.2f')],num2str(dt(2),'%.2f'),num2str(dt(3),'%.2f'),num2str(dt(4),'%.2f')},[0.1 0.9],'fontsize',14);
end

if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.9],'fontsize',14,'color',[0 0 0]);
  irf_legend(hca,{irf_ssub('MMS?',ic)},[0.08 0.9],'fontsize',14,'color',[0 0 0]);
end
if 1 % E 4sc dt
  hca = irf_panel('E N');
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_plot(hca,{mvaE1.z.tlim(tint),mvaE2.z.tlim(tint),mvaE3.z.tlim(tint),mvaE4.z.tlim(tint)},'comp','dt',dt);
  irf_plot(hca,{gseE1.dot(n_tinycs).tlim(tint),gseE2.dot(n_tinycs).tlim(tint),gseE3.dot(n_tinycs).tlim(tint),gseE4.dot(n_tinycs).tlim(tint)},'comp','dt',dt);      
  ylabel(hca,{'E_{N,loc}','(mv/m)'},'interpreter','tex');
  irf_legend(hca,{['dt = ' num2str(dt(1),'%.2f')],num2str(dt(2),'%.2f'),num2str(dt(3),'%.2f'),num2str(dt(4),'%.2f')},[0.1 0.9],'fontsize',14);
end

irf_zoom(h,'x',tintZoomTiny)

irf_zoom(h([1 3]),'y')

irf_plot_axis_align

hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.2 8.2]; hca.YLim = [0 180];
hca = irf_panel('Ve par'); hca.YLim = [-400 900];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels-1
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end
irf_plot_axis_align
hl = add_length_on_top(h(1),38,0.5);
hl.XLabel.String = 'km';
hl.FontSize = 14;

%% 2 figure: one with TSeries markings and one with pitchangle plots
tind = 2+[890 895 900 905 910 915 920];
tind = +5+[890 895 900 905 910 915];
tind = +5+[890 900 910];
%tind = 0+[915 920 925];
ic = 1;
c_eval('ePitch = ePitch?;',ic)

%% Par apar and perp in separate subplots
figure(42)
nRows = 1;
nCols = 3;
c_eval('h(?) = subplot(nRows,nCols,?);',1:nRows*nCols);

colorOrder = mms_colors('blues');
colorOrder = colormap('jet');
colorOrder = cn.cmap('discrete_jet');
posPleg = [0.98 0.98];

isub = 1;
if 1
  hca = h(isub); isub = isub + 1;
  paind = 1;
  set(hca,'ColorOrder',colorOrder)
  hp = plot(hca,ePitch.depend{1}(tind,:)',squeeze(ePitch.data(tind,:,paind))');
  %c_eval('hp(?).Color = colorOrder(?,:);',1:numel(hp))
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,'pa = 0{\circ}',posPleg)
end
if 1
  hca = h(isub); isub = isub + 1;
  paind = 8;
  set(hca,'ColorOrder',colorOrder)
  hp=plot(hca,ePitch.depend{1}(tind,:)',squeeze(ePitch.data(tind,:,paind))');
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,'pa = 90{\circ}',posPleg)
end
if 1
  hca = h(isub); isub = isub + 1;
  paind = size(ePitch.data(tind,:,:),3);
  set(hca,'ColorOrder',colorOrder)
  hp=plot(hca,ePitch.depend{1}(tind,:)',squeeze(ePitch.data(tind,:,paind))');
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,'pa = 180{\circ}',posPleg)
end

c_eval('lineColors(?,:) = hp(?).Color;',1:numel(hp));

c_eval('h(?).YTick = 10.^[0:5];',1:3)
c_eval('h(?).YLim = [1e0 1e5];',1:3)
c_eval('h(?).XLim = [10 1000];',1:3)

%% Par apar and perp in subplots by time
figure(43)
nRows = 1;
nCols = 3;
nPanels = nRows*nCols;
c_eval('h(?) = subplot(nRows,nCols,?);',1:nPanels);

posPleg = [0.98 0.98];

isub = 1;
for it = 1:numel(tind);
  hca = h(isub); isub = isub + 1;
  paind = [1 8 15];
  hp = plot(hca,ePitch.depend{1}(tind(it),:)',squeeze(ePitch.data(tind(it),:,paind))');
  %c_eval('hp(?).Color = colorOrder(?,:);',1:numel(hp))
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,ePitch(tind(it)).time.utc('HH:MM:SS:mmm'),posPleg,'fontsize',14)
end
hleg = legend(h(2),{'0\circ','90\circ','180\circ'},'location','southwest'); hleg.Box = 'off';
c_eval('h(?).YTick = 10.^[0:5];',1:nPanels)
c_eval('h(?).YLim = [1e0 1e5];',1:nPanels)
c_eval('h(?).XLim = [10 1000];',1:nPanels)
c_eval('h(?).FontSize = 14;',1:nPanels)

%% TSeries
figure(41)
npanels = 2;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; 
figurePostition(3)=figurePostition(3)*0.5; 
figurePostition(4)=figurePostition(4)*0.5;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % eDist omni 64
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all low E');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?_new.x.tlim(tint),mvaE?_new.y.tlim(tint),mvaE?_new.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1]),'y')
irf_plot_axis_align

% Mark time
hmark = irf_pl_mark(h,ePitch.time(tind).epochUnix);
c_eval('hmark(!,?).Color = lineColors(?,:);',1:numel(hp),1:4);

h(2).XGrid = 'off';
h(2).YGrid = 'off';

%% Plot par apar and perp in different subplots for a few times
tind = 2+[890 895 900 905 910 915 920];
tind = +5+[890 895 900 905 910 915];
tind = +4+[890 895 900 905 910];
%tind = 0+[915 920 925];
ic = 4;
c_eval('ePitch = ePitch?;',ic)

npanels = 4;
[h,h2] = initialize_combined_plot(npanels,3,1,0.6,'vertical');


tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; 
%figurePostition(3)=figurePostition(3)*0.5; 
%figurePostition(4)=figurePostition(4)*0.5;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % E new
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?_new.x.tlim(tint),mvaE?_new.y.tlim(tint),mvaE?_new.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E new + par
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?_new.x.tlim(tint),mvaE?_new.y.tlim(tint),mvaE?_new.z.tlim(tint),mvaE?par_new.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','||'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % eDist omni 64
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all low E');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end

h(1).Title.String = irf_ssub('MMS ?',ic);
irf_zoom(h,'x',tintZoom)
irf_zoom(h([1 2 3]),'y')
irf_plot_axis_align

% Mark time
lineColors = mms_colors('matlab');
c_eval('hmark = irf_pl_mark(h,ePitch.time(tind(?)).epochUnix,lineColors(?,:));',1:numel(tind))

c_eval('h(?).XGrid = ''off'';',1:npanels)
h(4).XGrid = 'off';
h(4).YGrid = 'off';


colorOrder = mms_colors('blues');
colorOrder = colormap('jet');
colorOrder = cn.cmap('discrete_jet');
posPleg = [0.98 0.98];

% PLOT PITCHANGLE DISTRIBUTIONS
isub = 1;
if 1
  hca = h2(isub); isub = isub + 1;
  paind = 1;
  set(hca,'ColorOrder',colorOrder)
  hp = plot(hca,ePitch.depend{1}(tind,:)',squeeze(ePitch.data(tind,:,paind))');
  %c_eval('hp(?).Color = colorOrder(?,:);',1:numel(hp))
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,'pa = 0{\circ}',posPleg)
end
if 1
  hca = h2(isub); isub = isub + 1;
  paind = 8;
  set(hca,'ColorOrder',colorOrder)
  hp=plot(hca,ePitch.depend{1}(tind,:)',squeeze(ePitch.data(tind,:,paind))');
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,'pa = 90{\circ}',posPleg)
end
if 1
  hca = h2(isub); isub = isub + 1;
  paind = size(ePitch.data(tind,:,:),3);
  set(hca,'ColorOrder',colorOrder)
  hp=plot(hca,ePitch.depend{1}(tind,:)',squeeze(ePitch.data(tind,:,paind))');
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  irf_legend(hca,'pa = 180{\circ}',posPleg)
end

c_eval('lineColors(?,:) = hp(?).Color;',1:numel(hp));

c_eval('h2(?).YTick = 10.^[0:5];',1:3)
c_eval('h2(?).YLim = [1e1 2e4];',1:3)
c_eval('h2(?).XLim = [10 1000];',1:3)
