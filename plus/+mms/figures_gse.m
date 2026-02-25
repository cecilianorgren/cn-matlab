%% Figure: Electron adiabaticity
npanels = 9;
tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.10Z');
h = irf_plot(npanels);
ic = 1;

cmap = 'jet';

zoomE = [];
zoomPA = [];
zoomY = [];
isub = 0;
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
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,{re?,Le?,rp?,Lp?},''comp'');',ic)
  irf_legend(hca,{'r_e','L_e','r_p','L_p'},[0.95 0.9],'fontsize',12);  
  hca.YLabel.String = {'length','(km)'};
  set(hca,'ColorOrder',mms_colors('1234'))  
  hca.YScale = 'log';
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  ylabel(hca,{'E_e','(eV)'},'interpreter','tex'); 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64
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
end
if 1 % ePDist pa low energies
  %%
  hca = irf_panel('e PA e32 deflux all low E');  isub = isub + 1;  zoomPA = [zoomPA isub];
  elim = [80 30000];  
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)    
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
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
if 0 % E par
  hca = irf_panel('E par');
  set(hca,'ColorOrder',mms_colors('2'))
  c_eval('irf_plot(hca,{gseE?par},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('2'))
  irf_legend(hca,{'E_{||}'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ve
  isub = isub + 1;
  zoomY = [zoomY isub];
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x,gseVe?.y,gseVe?.z},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
