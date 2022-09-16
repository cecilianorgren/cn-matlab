%% Overview, MSP > MSH > MSP, MVA ions and electrons
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere
npanels = 7;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.3; figurePostition(4)=figurePostition(4)*0.5;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 0 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'|B|','(nT)'};  
end
if 0 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curl'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'moments'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1000 1000];
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % P, 1sc, given above
  hca = irf_panel(sprintf('P sc %g',ic));
  Pishift = 0.4;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',ic)
  %%
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1'))
  irf_legend(hca,{'P_e'},[0.05 0.9],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('2'))
  irf_legend(hca,{sprintf('P_i - %.2f nPa',Pishift)},[0.05 0.6],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('3'))
  irf_legend(hca,{'P_B'},[0.05 0.2],'fontsize',12);
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.01 0.2],'fontsize',14);
  hca.YLim = [0.0 0.3];  
end
if 0 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.e64.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 500];  
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.8 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 80];
  %hca.YTick
end
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)
  hca.YLim = [10 1000];
  hca.CLim = [5 8.3];
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 100];
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
if 0 % ePDist pa 32
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [100 1000];  
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
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  c_eval('gseE?_lowf = gseE?.resample(gseVe?);',ic)
  hca = irf_panel('E lowf');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?_lowf.x.tlim(tint),gseE?_lowf.y.tlim(tint),gseE?_lowf.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % B sum spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bsum');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Bsum;
  specrec.f_label='';
  specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.0,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.0,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.0,'color','w')
  irf_plot(hca,fpe,'linewidth',1.0,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-8 -4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

%% Overview, MSP > MSH > MSP, MVA, mostly ions
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere
npanels = 7;
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
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 1 % B abs
  hca = irf_panel('abs B');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  hca.YLim = [0 10];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'moments'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1000 1000];
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
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
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)
  hca.YLim = [10 1000];
  hca.CLim = [5 8.3];
end
irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

%% Overview, MSH > MSP, MVA electrons
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere
npanels = 10;
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
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 1 % B
  hca = irf_panel('Curv B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaCurvB.x.tlim(tint),mvaCurvB.y.tlim(tint),mvaCurvB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'Curv B','(1/km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-0.05 0.02];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'moments'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1000 1000];
end
if 0 % B abs
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(nT)'};  
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 500];  
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(nT)'};    
  hca.YLim = [-500 500];  
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-300 300];  
end
if 1 % Te par perp
  hca = irf_panel('Te per par');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.6 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 80];
  %hca.YTick
end
if 1 % Te par perp, Ti/10
  hca = irf_panel('T');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,gseTi?.trace/3/10},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}','T_i/10'},[0.6 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 80];
  %hca.YTick
end
if 0 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curl'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 1 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)
  hca.YLim = [10 1000];
  hca.CLim = [5 8.3];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 120];
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
if 1 % ePDist pa 32
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [120 1000];  
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
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF par');
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF par');
  pas = [75 115];
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF par');
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  c_eval('gseE?_lowf = gseE?.resample(gseVe?);',ic)
  hca = irf_panel('E lowf');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?_lowf.x.tlim(tint),gseE?_lowf.y.tlim(tint),gseE?_lowf.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

tint = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:55.00Z'); irf_zoom(h,'x',tint)

if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
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
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  c_eval('gseE?_lowf = gseE?.resample(gseVe?);',ic)
  hca = irf_panel('E lowf');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?_lowf.x.tlim(tint),gseE?_lowf.y.tlim(tint),gseE?_lowf.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

%tint = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:55.00Z'); irf_zoom(h,'x',tint)

%% Overview, MSH > MSP, MVA electrons
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosheath-magnetosphere
npanels = 9;
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
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 0 % B curv
  hca = irf_panel('Curv B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaCurvB.x.tlim(tint),mvaCurvB.y.tlim(tint),mvaCurvB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'Curv B','(1/km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-0.05 0.02];
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'moments'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1000 1000];
end
if 0 % B abs
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(nT)'};  
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 500];  
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(nT)'};    
  hca.YLim = [-500 500];  
end
if 1 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-300 300];  
end
if 1 % E+VexB perp
  hca = irf_panel('E VexB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaEVexB?.x.tlim(tint),mvaEVexB?.y.tlim(tint),mvaEVexB?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E+v_{e}xB','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-300 300];  
end
comps_str = ['L','M','N'];
comps = ['x','y','z'];
for icomp = 1:numel(comps)
  comp = comps(icomp);
  comp_str = comps_str(icomp);
  if 1 % ve/iperp +VExB perp
    hca = irf_panel(['E VexB ' comp]);
    set(hca,'ColorOrder',mms_colors('xyza'))
    %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
    c_eval('irf_plot(hca,{mvaVe?perp.(comp).tlim(tint),mvaVi?perp.(comp).tlim(tint),mvaVExB?.(comp).tlim(tint)},''comp'');',ic)  
    hca.YLabel.String = {['v_{\perp' comp_str '}'],'(km/s)'};
    set(hca,'ColorOrder',mms_colors('xyza'))
    irf_legend(hca,{'v_e','v_i','v_{ExB}'},[0.98 0.9],'fontsize',12);  
    %hca.YLim = [-300 300];  
  end
end
if 0 % Te par perp, Ti/10
  hca = irf_panel('T');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,gseTi?.trace/3/10},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}','T_i/10'},[0.6 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 80];
  %hca.YTick
end
if 1 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curl'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)
  hca.YLim = [10 1000];
  hca.CLim = [5 8.3];
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 120];
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
if 0 % ePDist pa 32
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [120 1000];  
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
if 0 % eDist omni 64 par
  hca = irf_panel('e DEF par');
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 1000];
end
if 0 % eDist omni 64 perp
  hca = irf_panel('e DEF perp');
  pas = [75 115];
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 1000];
end
if 0 % eDist omni 64 apar
  hca = irf_panel('e DEF apar');
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 1000];
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
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

%tint = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:55.00Z'); irf_zoom(h,'x',tint)


%% Electron adiabaticity, prepare data
ic = 1;
units = irf_units;
%c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,15);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('Te?perp = (facTe?.yy + facTe?.zz)/2;',ic);
c_eval('Te?par = facTe?.xx;',ic);
c_eval('vte?perp = irf.ts_scalar(gseB?.time,sqrt(units.eV*Te?perp.resample(gseB?).data*2/units.me)*1e-3); vte?perp.units =''km/s'';',ic);
c_eval('vte?par = irf.ts_scalar(gseB?.time,sqrt(units.eV*Te?par.resample(gseB?).data*2/units.me)*1e-3); vte?par.units =''km/s'';',ic);
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
c_eval('mag_mom? = irf.ts_scalar(gseB?.time,0.5*units.me*vte?perp.data.^2*10^6./(gseB?.abs.data*1e-9)*1e9);  mag_mom?.units = ''nAm^2''; mag_mom?.name = ''magnetic moment'';',ic)


[out,l,v] = irf_minvar(gseB1.tlim(irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z')));
L = v(1,:); M = v(2,:); N = v(3,:);
coordLabels = {'L','M','N'};
lmn = [L;M;N];

disp(sprintf('L = [%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L,M,N))
% Rotate data
c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',ic)
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
%c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)

%% Electron adiabaticity
npanels = 7;
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintZoom = irf.tint('2015-10-16T10:33:43.00Z/2015-10-16T10:33:50.00Z');
tint = tintZoom;

tintB = tintZoom;

h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'parula';

Etotref = 100;
tref = irf_time('2015-10-16T10:33:45.00Z','utc>epochtt');
tref = irf_time('2015-10-16T10:33:44.50Z','utc>epochtt');
tref = irf_time('2015-10-16T10:33:45.00Z','utc>epochtt');
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = mag_mom?.resample(tref).data;',ic)
murefs = [(1:3:20)]*1e5*1e-6;
c_eval('tsMuref = irf.ts_scalar(mag_mom?.time,repmat(murefs,mag_mom?.length,1)); tsMuref.units = ''nAm^2'';',ic)
Eref_J = Bref*1e-9*murefs*1e-9; % J % perpendicular reference 
Eref_eV = Eref_J/units.e; % eV
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data/Bref*Eref_eV);',ic)
c_eval('A = Te?par.resample(tref)/Te?perp.resample(tref); A = A.data;',ic) % Epar/Eperp
c_eval('Wpar = irf.ts_scalar(mvaB?.tlim(tintB).time,(A+1-mvaB?.tlim(tintB).abs.data/Bref)*Eref_eV);',ic)

% prepare energy levels for adiabaticity
if 0
%c_eval('tintB = tref+tlength?;',ic)
tintB = tint;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp)',ic)
Wrefperp = [50 100 200 500 1000]; nWperp = numel(Wrefperp);
Wrefpar = 40;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
%c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wperp! = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp(!)/Bref);',ic,1:nWperp)
c_eval('Wpar = irf.ts_scalar(Wperp?.time,Wreftot-Wperp?.data);',1:nWperp)
Efactorpar = 1;
Efactorperp = Wrefperp/Wrefperp(1);
end

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
  set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
  c_eval('irf_plot(hca,mag_mom?);',ic)
  hold(hca,'on')
  set(hca,'ColorOrder',[mms_colors('matlab')])
  irf_plot(hca,tsMuref);
  hold(hca,'off')
  c_eval('mag_mom_units = mag_mom?.units;',ic)
  ylabel(hca,{'magnetic moment',['(' mag_mom_units ')']});
  %irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 1 % temperature
  hca = irf_panel('Te'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{Te?par,Te?perp},''comp'');',ic)  
  hca.YLabel.String = {'T_{e}',['(eV)']};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'||','\perp'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 1 % thermal velocity
  hca = irf_panel('vte'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{vte?par,vte?perp},''comp'');',ic)  
  hca.YLabel.String = {'v_{te}',['(km/s)']};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'||','\perp'},[0.8 0.9],'fontsize',12);
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
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 0 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
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
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wperp);
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
      
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)      
    hold(hca,'off')
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  if 0 % plot parallel energy based on magnetic moment conservation    
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
if 0 % ePDist pa low energies
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
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,16).elim(elim).deflux.specrec(''pa''),''log'');',ic)
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
if 0 % E new
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
irf_zoom(h,'y')
for  ii = zoomE; h(ii).YLim = [10 1000]; end
for  ii = zoomE(1:3); h(ii).CLim = [5.5 8.7]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end
irf_plot_axis_align


h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  irf_pl_mark(h(ii+pshift),tref,'k')
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Electron adiabaticity, Teperp
npanels = 7;
tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintZoom = irf.tint('2015-10-16T10:33:41.00Z/2015-10-16T10:33:50.00Z');

tintB = tintZoom;

h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

Etotref = 100;
tref = irf_time('2015-10-16T10:33:45.00Z','utc>epochtt');
tref = irf_time('2015-10-16T10:33:44.50Z','utc>epochtt');
tref = irf_time('2015-10-16T10:33:45.00Z','utc>epochtt');
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = mag_mom?.resample(tref).data;',ic)
murefs = (1:3:20)*1e5*1e-6;
c_eval('tsMuref = irf.ts_scalar(mag_mom?.time,repmat(murefs,mag_mom?.length,1)); tsMuref.units = ''nAm^2'';',ic)
Eref_J = Bref*1e-9*murefs*1e-9; % J % perpendicular reference 
Eref_eV = Eref_J/units.e; % eV
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data/Bref*Eref_eV);',ic)
c_eval('A = Te?par.resample(tref)/Te?perp.resample(tref); A = A.data;',ic) % Epar/Eperp
c_eval('Wpar = irf.ts_scalar(mvaB?.tlim(tintB).time,(A+1-mvaB?.tlim(tintB).abs.data/Bref)*Eref_eV);',ic)

% prepare energy levels for adiabaticity
if 0
%c_eval('tintB = tref+tlength?;',ic)
tintB = tint;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp)',ic)
Wrefperp = [50 100 200 500 1000]; nWperp = numel(Wrefperp);
Wrefpar = 40;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
%c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wperp! = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp(!)/Bref);',ic,1:nWperp)
c_eval('Wpar = irf.ts_scalar(Wperp?.time,Wreftot-Wperp?.data);',1:nWperp)
Efactorpar = 1;
Efactorperp = Wrefperp/Wrefperp(1);
end

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
  set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
  c_eval('irf_plot(hca,mag_mom?);',ic)
  hold(hca,'on')
  set(hca,'ColorOrder',[mms_colors('matlab')])
  irf_plot(hca,tsMuref);
  hold(hca,'off')
  c_eval('mag_mom_units = mag_mom?.units;',ic)
  hca.YLabel.String = {'\mu',['(' mag_mom_units ')']};
  set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
  %irf_legend(hca,{'n_e','n_i'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 1 % temperature
  hca = irf_panel('Te'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{Te?perp,Te?par},''comp'');',ic)  
  hca.YLabel.String = {'T_{e}',['(eV)']};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'\perp','||'},[0.8 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end
if 1 % thermal velocity
  hca = irf_panel('vte'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{vte?perp,vte?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{te}',['(km/s)']};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'\perp','||'},[0.8 0.9],'fontsize',12);
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
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 0 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
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
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 1 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wperp);
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
      
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)      
    hold(hca,'off')
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  if 0 % plot parallel energy based on magnetic moment conservation    
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
if 0 % ePDist pa low energies
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
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,16).elim(elim).deflux.specrec(''pa''),''log'');',ic)
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
if 0 % E new
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
irf_zoom(h,'y')
for  ii = zoomE; h(ii).YLim = [10 1000]; end
for  ii = zoomE(1:3); h(ii).CLim = [5.5 8.7]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end
irf_plot_axis_align


h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  irf_pl_mark(h(ii+pshift),tref,'k')
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% electrons and electron distributions, MSH > MSP, 

tintZoom = irf.tint('2015-10-16T10:33:43.00Z/2015-10-16T10:33:52.00Z');
tint = tintZoom;

npanels = 5;

cmap = 'jet';

% Initialize figure
ic = 1;
[h,h2] = initialize_combined_plot(npanels,3,1,0.6,'vertical');

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.3; figurePostition(4)=figurePostition(4)*0.5;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 0 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 1 % B and B abs
  hca = irf_panel('B + Babs');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'L','M','N','||'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 29];
end
if 0 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 0 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'|B|','(nT)'};  
end
if 0 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'curl'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 0 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'moments'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1000 1000];
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % P, 1sc, given above
  hca = irf_panel(sprintf('P sc %g',ic));
  Pishift = 0.4;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',ic)
  %%
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1'))
  irf_legend(hca,{'P_e'},[0.05 0.9],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('2'))
  irf_legend(hca,{sprintf('P_i - %.2f nPa',Pishift)},[0.05 0.6],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('3'))
  irf_legend(hca,{'P_B'},[0.05 0.2],'fontsize',12);
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.01 0.2],'fontsize',14);
  hca.YLim = [0.0 0.3];  
end
if 0 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.e64.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 500];  
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx,(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.9 0.95],'fontsize',14);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [15 55];
  
  %hca.YTick
end
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)
  hca.YLim = [10 1000];
  hca.CLim = [5 8.3];
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 100];
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
if 0 % ePDist pa 32
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [100 1000];  
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
endif 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 0 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
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
if 0 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [0 30];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 0 % plot reference energy level
    hold(hca,'on')      
    c_eval('lineEref = irf_plot(hca,irf.ts_scalar(mvaVe?.time,repmat(Eref,mvaVe?.length,1)),''k'');',ic)  
    for iii = 1: numel(lineEref) lineEref(iii).Color = colorEref; lineEref(iii).LineWidth = 1.5; end
    hold(hca,'off')    
  end
  if 0 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
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
  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  if 0 % plot perpendicular energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wperp);
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 1200];
  hca.CLim = [6 9];
end
if 0 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); isub = isub + 1;  zoomE = [zoomE isub];
  pas = [150 180];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
      
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.tlim(tintZoom).deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    set(hca,'ColorOrder',[0 0 0;mms_colors('matlab')])
    c_eval('lineScpot = irf_plot(hca,scPot?);',ic)      
    hold(hca,'off')
  end
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('matlab'))
    irf_plot(hca,Wpar);
    hold(hca,'off')
    irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  if 0 % plot parallel energy based on magnetic moment conservation    
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
if 1 % ePDist pa mid energies
  hca = irf_panel('e PA e32 deflux all low E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  
  elim = [100 800];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,16).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  
  hca.XGrid = 'off';
  hca.YGrid = 'off';
    
  hca.YLabel.String = {'\theta_{PA}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
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
if 1 % B sum spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bsum');
  specrec=struct('t',polarization.t);
  specrec.f=frequency;
  specrec.p=Bsum;
  specrec.f_label='';
  specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.0,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.0,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.0,'color','w')
  irf_plot(hca,fpe,'linewidth',1.0,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-9 -2])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Plot single time particle distributions, 1 sc, 4 projections,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('ePitch = ePitch?.convertto(''s^3/km^6'');',ic)

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [-1 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 786;%:1:920;953; 911;
tind = 1160;
tind = 1240;
tind = 1451:1500;%1350:1450;
for it = 1670:1730;%%:2:1381;%2260:2300%:1710%:2370;1420:1530;2340:2370;tind;589;%:650;550:650;%1240;tind;
  time = times(it);
  %it
  if exist('hmark'); delete(hmark); end
  c_eval('hmark(?) = irf_pl_mark(h(?),time.epochUnix'',''black'');',1:npanels)
  
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);  
  

  timeUTC = time.utc;      
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  % Perpendicular plane    
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{\perp,1}','v_{\perp,2}','v_{||}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;par;-perp2]; vlabels = {'v_{\perp,1}','v_{||}','-v_{perp2}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  if 0 % B plane 2
    hca = h2(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{perp2}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';
  end
  % Pitchangle distribution
  hca = h2(isub); isub = isub + 1;
  fE = ePitch.depend{1}(it,:);
  fv = sqrt(fE*units.eV*2/units.me)/1000;
  fref = cn.maxwellian(fv,40,15,0,'e',1);
  plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])),'linewidth',1.5);
  hold(hca,'on')
  plot(hca,fE,fref*1e18,'k--')
  hold(hca,'off')
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  hca.XLim = [10 1000];
  legend(hca,{'0','90','180'},'location','southwest') %eastoutside
  hca.YTick = 10.^[-4:5];
  hca.YLim = [1e-1 1e5];
  axis(hca,'square')

  % ExB plane, with markings
  if 0
    vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    hca.Title.String = '';
  end
  for ii = 1:3
    colormap(h2(ii),strCMap)
  end
  c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',3:npanels)
   
  cn.print(['e-dist_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',[eventPath 'tseries-dist/'])
end

%% Few panels to show basic structure of fluxropes
tintZoom = irf.tint('2015-10-16T10:33:43.00Z/2015-10-16T10:33:52.00Z');
tint = tintZoom;
npanels = 3;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.2; figurePostition(4)=figurePostition(4)*0.4;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 1 % AV B and B abs
  hca = irf_panel('B + Babs');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_plot(hca,{mvaAvB.x,mvaAvB.y,mvaAvB.z,mvaAvB.abs},'comp');
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',14);
  hca.YLim = [-25 29];
end
if 0 % B and B abs
  hca = irf_panel('B + Babs');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z,mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'L','M','N','||'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 29];
end
if 0 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 0 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'|B|','(nT)'};  
end
if 1 % J curl 
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint)},'comp');
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',14);
  irf_legend(hca,{'curlometer'},[0.05 0.9],'color',[0 0 0],'fontsize',14);
  hca.YScale = 'lin';
  hca.YLim = [-999 599];
end
if 0 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?perp.x,mvaJ?perp.y,mvaJ?perp.z,mvaJ?par},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'moments'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1000 1000];
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 0 % P, 1sc, given above
  hca = irf_panel(sprintf('P sc %g',ic));
  Pishift = 0.4;
  set(hca,'ColorOrder',mms_colors('1234'))
   %irf_plot(hca,{avPe,avPi-Pishift,avPB,avPi/10},'comp');
  c_eval('irf_plot(hca,{mvaPe?.trace/3,mvaPi?.trace/3-Pishift,PB?},''comp'');',ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1'))
  irf_legend(hca,{'P_e'},[0.05 0.9],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('2'))
  irf_legend(hca,{sprintf('P_i - %.2f nPa',Pishift)},[0.05 0.6],'fontsize',12);
  set(hca,'ColorOrder',mms_colors('3'))
  irf_legend(hca,{'P_B'},[0.05 0.2],'fontsize',12);
  %irf_legend(hca,{'P_e',sprintf('P_i-%.2f',Pishift),'P_B','P_i/10'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'P_e',sprintf('P_i - %.2f nPa',Pishift),'P_B'},[0.01 0.2],'fontsize',14);
  hca.YLim = [0.0 1];  
end
if 0 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.e64.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-500 500];  
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,||}','(nT)'};    
  hca.YLim = [-500 500];  
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z},''comp'');',ic)  
  hca.YLabel.String = {'v_{e,\perp}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-300 300];  
end
if 0 % Ve perp par
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?perp.x,mvaVe?perp.y,mvaVe?perp.z,mvaVe?par},''comp'');',ic)  
  hca.YLabel.String = {'v_{e}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\perp,L','\perp,M','\perp,N','||'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-300 300];  
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx,(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.95 0.1],'fontsize',14);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [15 55];
  %hca.YTick
end
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.e64.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 1 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)
  hca.YLim = [10 1000];
  hca.CLim = [5 8.3];
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 100];
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
if 0 % ePDist pa 32
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [100 1000];  
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
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % E
  c_eval('gseE?_lowf = gseE?.resample(gseVe?);',ic)
  hca = irf_panel('E lowf');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?_lowf.x.tlim(tint),gseE?_lowf.y.tlim(tint),gseE?_lowf.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % B sum spectrogram
  isspec(isub) = 1; isub = isub + 1;
  hca=irf_panel('Bsum');
  specrec=struct('t',time);
  specrec.f=frequency;
  specrec.p=Bsum;
  specrec.f_label='';
  specrec.p_label={'log_{10}B^{2}','nT^2 Hz^{-1}'};
  irf_spectrogram(hca,specrec,'log','donotfitcolorbarlabel');
  hold(hca,'on');
  irf_plot(hca,ecfreq,'linewidth',1.0,'color','w')
  irf_plot(hca,ecfreq05,'linewidth',1.0,'color','w')
  irf_plot(hca,ecfreq01,'linewidth',1.0,'color','w')
  irf_plot(hca,fpe,'linewidth',1.0,'color','w')
  hold(hca,'off');
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3]);
  caxis(hca,[-8 -4])
  ylabel(hca,'f (Hz)','fontsize',12);
  colormap(hca,'jet');
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14)
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
end
