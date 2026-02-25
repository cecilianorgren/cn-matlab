%% LONGER OVERVIEW INTERVALS
h = irf_plot(4);
if 1 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x.tlim(tint),gseB2.x.tlim(tint),gseB3.x.tlim(tint),gseB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.y.tlim(tint),gseB2.y.tlim(tint),gseB3.y.tlim(tint),gseB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.z.tlim(tint),gseB2.z.tlim(tint),gseB3.z.tlim(tint),gseB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');   
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'X','Y','Z'},[0.98 0.2],'fontsize',12);  
end


%% SHORTER BURST INTERVALS
ic = 1;
tint = irf.tint('2015-10-16T10:25:40.00Z/2015-10-16T10:27:40.00Z');
TSTART = tic; mms_survey.load_data_overview; toc(TSTART);

%%
units = irf_units;
c_eval('gseJe?fast = -units.e*ne?fast*gseVe?fast*1e3*1e6*1e9; gseJe?fast.units = ''nA/m^2''; gseJe?fast.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi?fast = units.e*ne?fast*gseVi?fast.resample(ne?fast.time)*1e3*1e6*1e9; gseJi?fast.units = ''nA/m^2''; gseJi?fast.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ?fast = (gseJe?fast+gseJi?fast);',ic);

%% Make plot, fast data
tint = irf.tint('2015-11-12T06:05:00.00Z/2015-11-12T07:55:00.00Z');
%tint = irf.tint('2015-11-12T06:55:00.00Z/2015-11-12T07:32:00.00Z');
npanels = 7;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;
toZoom = [];
isub = 1;
% Plot
if 1 % B
  toZoom = [toZoom isub]; isub = isub + 1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{gseB?srvy.x.tlim(tint),gseB?srvy.y.tlim(tint),gseB?srvy.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % n
  toZoom = [toZoom isub]; isub = isub + 1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ni?fast,ne?fast},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'}; 
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_i','n_e'},[0.98 0.9],'fontsize',12);
end
if 0 % J
  toZoom = [toZoom isub]; isub = isub + 1;
  hca = irf_panel('J');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJ?fast.x.tlim(tint),gseJ?fast.y.tlim(tint),gseJ?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-300 220];
end
if 1 % Vi
  toZoom = [toZoom isub]; isub = isub + 1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?fast.x.tlim(tint),gseVi?fast.y.tlim(tint),gseVi?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % iPDist omni 64
  toZoom = [toZoom]; isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?fast.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 1 % Ve
  toZoom = [toZoom isub]; isub = isub + 1;
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?fast.x.tlim(tint),gseVe?fast.y.tlim(tint),gseVe?fast.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-300 220];
end
if 1 % ePDist omni 64
  toZoom = [toZoom]; isub = isub + 1;
  hca = irf_panel('e DEF omni');  
  c_eval('irf_spectrogram(hca,ePDist?fast.tlim(tint).deflux.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_e','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % iPDist pa 32
  toZoom = [toZoom]; isub = isub + 1;
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?fast.deflux.tlim(tint).pitchangles(dmpaB?srvy,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 1 % Ti par perp
  toZoom = [toZoom isub]; isub = isub + 1;
  hca = irf_panel('Ti');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTi?fast.xx.tlim(tint),(facTi?fast.yy+facTi?fast.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_i','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin'; %hca.YTick = [10:10:100 200:100:1000];  
end
if 0 % Walen ratio
  toZoom = [toZoom]; isub = isub + 1;
  hca = irf_panel('Walen angle');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{walAngle?i.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'\theta_W'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [0 180];  
end
if 0 % Walen angle
  toZoom = [toZoom]; isub = isub + 1;
  hca = irf_panel('Walen ratio');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{walR?i.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'R_W'};
  set(hca,'ColorOrder',mms_colors('1'))
  hca.YLim = [0 5];
end
if 0 % alfa
  hca = irf_panel('alpha');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{Alfa?.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'\alpha'};
  set(hca,'ColorOrder',mms_colors('1'))  
end


irf_zoom(h,'x',tint)
%irf_zoom(h([1:4 6 8]),'y')
irf_zoom(h(toZoom),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
hmark = irf_pl_mark(h,tintZoom.epochUnix');

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Figure 1: 1 sc overview 
%tint = irf.tint('2015-10-16T10:34:20.00Z/2015-10-16T10:35:30.00Z');
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

%% Figure 2: 4 sc: eletron pa
npanels = 8;
h = irf_plot(npanels);

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x.tlim(tint),gseB2.x.tlim(tint),gseB3.x.tlim(tint),gseB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.y.tlim(tint),gseB2.y.tlim(tint),gseB3.y.tlim(tint),gseB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.z.tlim(tint),gseB2.z.tlim(tint),gseB3.z.tlim(tint),gseB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E mms1');
  elim = [0 400];
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
  elim = [0 400];
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
  elim = [0 400];
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
  elim = [0 400];
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

irf_zoom(h,'x',tint)
irf_zoom(h([1:4]),'y')
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

%% Figure 2: 4 sc 
npanels = 11;
h = irf_plot(npanels);

pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x.tlim(tint),gseB2.x.tlim(tint),gseB3.x.tlim(tint),gseB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.y.tlim(tint),gseB2.y.tlim(tint),gseB3.y.tlim(tint),gseB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 1 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.z.tlim(tint),gseB2.z.tlim(tint),gseB3.z.tlim(tint),gseB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
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
  %lines = irf_plot(hca,{gseJ1.y,gseJ2.y,gseJ3.y,gseJ4.y,gseAvJ.y,gseJcurl.y},'comp');   
  lines = irf_plot(hca,{gseJ1.y,gseJ2.y,gseJ3.y,gseJ4.y},'comp');   
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
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp X
  hca = irf_panel('Ve perp X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).x,gseVe2perp.tlim(tint).x,gseVe3perp.tlim(tint).x,gseVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,x}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp Y
  hca = irf_panel('Ve perp Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).y,gseVe2perp.tlim(tint).y,gseVe3perp.tlim(tint).y,gseVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,y}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp Z
  hca = irf_panel('Ve perp Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).z,gseVe2perp.tlim(tint).z,gseVe3perp.tlim(tint).z,gseVe4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,z}','(km/s)'},'interpreter','tex');
end
if 1 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{gseEVexB1.abs,gseEVexB2.abs,gseEVexB3.abs,gseEVexB4.abs},'comp'); 
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
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P',['(10^{-3} ' gseGradPe.units ')']};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end


irf_zoom(h,'x',tint)
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

%% Figure 3: 4 sc, Ohm's law 
npanels = 9;
h = irf_plot(npanels);
pshift=0;

scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 0 % BX
  hca = irf_panel('BX');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.x.tlim(tint),gseB2.x.tlim(tint),gseB3.x.tlim(tint),gseB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{x}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % BY
  hca = irf_panel('BY');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.y.tlim(tint),gseB2.y.tlim(tint),gseB3.y.tlim(tint),gseB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{y}','(nT)'};
end
if 0 % BZ
  hca = irf_panel('BZ');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.z.tlim(tint),gseB2.z.tlim(tint),gseB3.z.tlim(tint),gseB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{z}','(nT)'};
end
if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.abs.tlim(tint),gseB2.abs.tlim(tint),gseB3.abs.tlim(tint),gseB4.abs.tlim(tint)},'comp');
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
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).abs,gseVe2perp.tlim(tint).abs,gseVe3perp.tlim(tint).abs,gseVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{gseJcurl.x,gseJcurl.y,gseJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ohm's law x
  hca = irf_panel('Ohm X');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{gseAvE.x,gseOhmGradPe.x,-gseOhmVexB.x,-gseOhmVixB.x,gseOhmJxB.x,gseAvE.x+gseOhmVexB.resample(gseAvE).x+gseOhmGradPe.resample(gseAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{gseAvE.x,gseOhmGradPe.x,-gseOhmVexB.x,-gseOhmVixB.x,gseOhmJxB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_x','(mV/m)'};
end
if 1 % Ohm's law y
  hca = irf_panel('Ohm Y');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseAvE.y,gseOhmGradPe.y,-gseOhmVexB.y,-gseOhmVixB.y,gseOhmJxB.y},'comp');
  hca.YLabel.String = {'E_y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law z
  hca = irf_panel('Ohm Z');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseAvE.z,gseOhmGradPe.z,-gseOhmVexB.z,-gseOhmVixB.z,gseOhmJxB.z},'comp');
  hca.YLabel.String = {'E_z','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law residual 3 components
  hca = irf_panel('Ohm residual');
  set(hca,'ColorOrder',mms_colors('xyz'))
  residual = gseAvE+gseOhmVexB.resample(gseAvE)+gseOhmGradPe.resample(gseAvE);
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

%% Plot EDR signatures
npanels = 10;
h = irf_plot(npanels);

if 1 % B abs
  hca = irf_panel('B abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseB1.abs.tlim(tint),gseB2.abs.tlim(tint),gseB3.abs.tlim(tint),gseB4.abs.tlim(tint)},'comp');
  hca.YLabel.String = {'|B|','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1perp.tlim(tint).abs,gseVe2perp.tlim(tint).abs,gseVe3perp.tlim(tint).abs,gseVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 1 % agyrotropy
  hca = irf_panel('agyrotropy Aphi');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{agyro1.tlim(tint).abs,agyro2.tlim(tint).abs,agyro3.tlim(tint).abs,agyro4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'A\phi';
end
if 1 % agyrotropy
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'Q';
end
if 1 % agyrotropy
  hca = irf_panel('Dng');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Dng1.tlim(tint).abs,Dng2.tlim(tint).abs,Dng3.tlim(tint).abs,Dng4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'D_{ng}';
end

%% Plot detailed electron distribution for selected times
% Initialize figure
[h1,h2] = initialize_combined_plot(8,2,2,0.4,'vertical');
cmap = 'jet';
ic = 3;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

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
if 1 % Ve perp
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
if 0 % ePDist pa 64
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [200 400];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).e64.pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
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
if 1 % Q
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'Q';
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:4 7:8]),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

%% Plot single time particle distributions, 1 sc, 4 projections,
time = irf_time('2015-10-16T10:34:32.83Z','utc>epochtt');
time = irf_time('2015-10-16T11:28:36.40Z','utc>epochtt');
time = irf_time('2015-09-15T11:20:31.74Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:20.60Z','utc>epochtt');
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
projclim = [0 5];  
  

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

%% Plot single time particle distributions, 1 sc, 4 projections,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('gradPeslow = gseGradPe.resample(gseVe?);',ic)

% Projection coordinate system

%y = hatE;
%z = hatB;
%x = cross(y,z)/norm(cross(y,z));

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 886:1:920;953; 911;
%tind = 920;
for it = tind
  time = times(it);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  
  c_eval('hatGradPe = gradPeslow.resample(time).data/gradPeslow.resample(time).abs.data;',ic)
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);
  
  c_eval('hatV = dbcsVe?(it).data/dbcsVe?(it).abs.data;',ic)

  

  timeUTC = time.utc;      
  isub = 1;

  % Perpendicular plane    
  hca = h2(isub); isub = isub + 1; 
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';hatGradPe,'\nabla P_e'};
  xyz = [x;y;z]; vlabels = {'x_{GSE}','y_{GSE}','z_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [x;z;-y]; vlabels = {'x_{GSE}','z_{GSE}','-y_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [y;z;x]; vlabels = {'y_{GSE}','z_{GSE}','x_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % ExB plane, with markings
  vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z';hatV,'v_e'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = '';

  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','hard_path',eventPath)
end

%% Waves
tintZoom = irf.tint('2015-11-12T07:19:18.00Z/2015-11-12T07:19:24.00Z');

npanels = 4;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
if 1
  hca = irf_panel('E spectrogram');
  c_eval('tE = irf_time(wavE?.t,''epoch>epochtt'');',ic)
  iE = tE.tlim(tintZoom);
  c_eval('plotWavE = wavE?; plotWavE.t = plotWavE.t(iE); plotWavE.p{1} = plotWavE.p{1}(iE,:);',ic)
  c_eval('irf_spectrogram(hca,plotWavE,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end
if 1
  hca = irf_panel('B spectrogram');
  c_eval('tB = irf_time(wavB?.t,''epoch>epochtt'');',ic)
  iB = tB.tlim(tintZoom);
  c_eval('plotWavB = wavB?; plotWavB.t = plotWavB.t(iE); plotWavB.p{1} = plotWavB.p{1}(iB,:);',ic)
  c_eval('irf_spectrogram(hca,plotWavB,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.99],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

%% Plot parallel and perpendicular ve and E
% Initialize figure
%[h1,h2] = initialize_combined_plot(8,2,2,0.4,'vertical');
h1 = irf_plot(8);
cmap = 'jet';
ic = 4;

tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

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
if 1 % Ve perp
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  ylabel(hca,{'v_{e}','(km/s)'},'interpreter','tex');
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
if 0 % ePDist pa 64
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [200 400];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).e64.pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
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
if 1 % E par
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
  %irf_zoom(hca,'y')
end
if 1 % EdotJ par
  hca = irf_panel('EdotJ par');  
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{EdotJ?par},''comp'');',ic)
  hca.YLabel.String = {'E\cdotJ_{||}','(nW/m^3)'};
  set(hca,'ColorOrder',mms_colors('1'))    
  %irf_zoom(hca,'y')
end
if 0 % Q
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'Q';
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:4 6:7]),'y')
hca = irf_panel('E par'); hca.YLim = [-4 4];
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

%% Quiver plot
h2 = subplot(1,1,1);
isub = 1;      
tintQuivers = irf.tint('2015-11-12T07:19:20.59Z/2015-11-12T07:19:21.61Z');
if 1 % GSE, 1 velocity
  hca = h2(isub); isub = isub +1;
  hold(hca,'on')
  %irf_pl_mark(h,tintQuivers.epochUnix')
  times = gseVe1.tlim(tintQuivers).time;
  gseVmatch = 68*[-0.88 0.33 0.35]; % GSE  
  % add motion tangential to the boundary
  gseV = gseVmatch-cross([0 1 0],gseVmatch/norm(gseVmatch))*100;
  c_eval('posR? = repmat(gseRR?,times.length,1)-(times-times(1))*gseV;')
  c_eval('posV?perp = gseVe?perp.resample(times).data;')
  c_eval('posV? = gseVe?.resample(times).data;')
  c_eval('posB? = gseB?.resample(times).data;')
  %c_eval('plot_quivers(hca,[posV?perp(:,1) posV?perp(:,2) posV?perp(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posV?(:,1) posV?(:,2) posV?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posB?(:,1) posB?(:,2) posB?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''b''))')
  hold(hca,'off')

  hca.XLabel.String = 'x_{GSE}';
  hca.YLabel.String = 'y_{GSE}';
  hca.ZLabel.String = 'z_{GSE}';
  %hca.YDir = 'normal';
  %hca.XDir = 'reverse';
  axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
end
if 0 % 3D
   
  times = mvaVe1.tlim(tintQuivers).time;
  
  gseVouteredge = 65.3*[0.37  0.32 -0.87];  
  gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
  gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
  gseMSP = 32*[-0.88 -0.42 0.24];
  
  lmnVmsh = gseVinneredge*[L' M' N'];
  lmnVmsp = gseMSP*[L' M' N'];
  
  vel_selection = 5;
  clear timesVUTC gseVdata
  gseV
  switch vel_selection
    case 1 % do velocities manually    
      tanV = irf_norm(cross(gseV,M));
      gseV = 55*[-0.90 -0.28 -0.33]-30*tanV;

      gseV = repmat(gseV,times.length,1);
      gseV(1:27,:) = repmat(-50*L,27,1);   
      gseV(28:31,:) = repmat(-50*N,4,1);   
      gseV(44:60,:) = repmat(gseVouteredge,17,1);  
    case 2 % define velocities at certain times, and then interpolate to other times      
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-25*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 3 % more 'vertical'
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 4
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*50;-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 5
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*dot(-cross(irf_norm(gseMSP),M),gseVinneredge)...
                                                                                 -cross(irf_norm(gseMSP),M)*33*0.75;      
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge*0.75; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
    case 6
      tanVinneredge = irf_norm(cross(gseVinneredge,M));    
      tanVouteredge = irf_norm(cross(gseVouteredge,M));    
      iv = 0;
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:26.00Z'; gseVdata(iv,:) = gseMSP;%-55*L; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:27.00Z'; gseVdata(iv,:) = gseMSP/3;%-55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:29.50Z'; gseVdata(iv,:) = -cross(irf_norm(gseMSP),M)*dot(-cross(irf_norm(gseMSP),M),gseVinneredge); 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.04Z'; gseVdata(iv,:) = -55*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.05Z'; gseVdata(iv,:) = -55*N; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.00Z'; gseVdata(iv,:) = gseVoutflow-65*L; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.08Z'; gseVdata(iv,:) = gseVinneredge-40*tanVinneredge; 
      iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.30Z'; gseVdata(iv,:) = gseVinneredge;-40*tanVinneredge; 
      %iv = iv + 1; timesVUTC(iv,:) = '2015-10-16T10:33:30.80Z'; gseVdata(iv,:) = gseVouteredge;+40*tanVouteredge;
      gseV = irf.ts_vec_xyz(timesVUTC,gseVdata);
      gseV = gseV.resample(times);
  end
  
  dt = times(2)-times(1);
  
  lmnV = gseV.data*[L' M' N'];
  
  
  c_eval('posR? = repmat(mvaRR?,times.length,1)-dt*cumsum(lmnV,1);')
  c_eval('posV? = mvaVe?.resample(times).data*0.6;')
  c_eval('posJ? = mvaJ?.resample(times).data;')
  c_eval('posB? = mvaB?.resample(times).data;')
  c_eval('posE? = mvaE?.resample(times).data;')
  c_eval('posRe? = mvaEVexB?.resample(times).data;')
  
  hca = h2(isub); isub = isub +1;
  sclist = 1:4;
  hold(hca,'on')
  %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
  c_eval('plot_quivers(hca,[posV?(:,3) -posV?(:,2) posV?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
  c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
  %c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
  %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
  hold(hca,'off')

  hca.XLabel.String = 'N (km)';
  hca.YLabel.String = 'M (km)';
  hca.ZLabel.String = 'L (km)';
  hca.ZDir = 'normal';
  hca.YDir = 'normal';
  hca.XDir = 'reverse';
  %axis(hca,'square')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';
  view(hca,[0 -1 0])
  %hca.YLim = [-20 20];
  %axis(hca,'equal')
  axis(hca,'equal')
  
  hca.Box = 'on';
  hca.ZLim = [-10 75];
  hca.XLim = [-8 27];
  tintUTCstart = tintQuivers(1).utc;
  tintUTCstop = tintQuivers(2).utc;
  hca.Title.String = ['V_e: ' tintUTCstop(12:22) ' - ' tintUTCstop(12:22)];
  hca.Title.String = ['V_e'];
  hleg=legend(hca,'MMS 1','MMS 2','MMS 3','MMS 4','location','northeast');
  fontsize = 17;
  hca.XLabel.FontSize = fontsize;
  hca.YLabel.FontSize = fontsize;
  hca.ZLabel.FontSize = fontsize;
  hca.Title.FontSize = fontsize;
  hca.FontSize = fontsize;
  hleg.FontSize = 10;
  if 0 % plot electric field arrows
    hca = h2(isub); isub = isub +1;
    hold(hca,'on')
    %c_eval('plot_quivers(hca,[posV?(:,3) posV?(:,1)],[posR?(:,3) posR?(:,1)],mms_colors(''?''))')
    c_eval('plot_quivers(hca,[posE?(:,3) -posE?(:,2) posE?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''?''))',sclist)
    c_eval('color? = mms_colors(''?''); color?=(color? + [1 1 1]*2)/3;',sclist)
    c_eval('plot_quivers(hca,[posJ?(:,3) -posJ?(:,2) posJ?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',sclist)
    %c_eval('plot_quivers(hca,[posB?(:,3) -posB?(:,2) posB?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],mms_colors(''b''))',1:4)
    %c_eval('plot_quivers(hca,[posRe?(:,3) -posRe?(:,2) posRe?(:,1)],[posR?(:,3) -posR?(:,2) posR?(:,1)],color?)',sclist)
    hold(hca,'off')

    hca.XLabel.String = 'N';
    hca.YLabel.String = 'M';
    hca.ZLabel.String = 'L';
    hca.ZDir = 'normal';
    hca.YDir = 'normal';
    hca.XDir = 'reverse';
    %axis(hca,'square')
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.ZGrid = 'on';
    view(hca,[0 -1 0])
    %hca.YLim = [-20 20];
    %axis(hca,'equal')
    %axis(hca,'equal')
    %hca.XLim = 20*[-1 1];
    hca.Title.String = 'E, J (faint)';
  end
  if 0 % Interpolate BM in LN-plane to get Hall field color surface
    posL = double([posR1(:,1); posR2(:,1); posR3(:,1); posR4(:,1)]);
    posM = double([posR1(:,2); posR2(:,2); posR3(:,2); posR4(:,2)]);
    posN = double([posR1(:,3); posR2(:,3); posR3(:,3); posR4(:,3)]);
    posB = double([posB1(:,2); posB2(:,2); posB3(:,2); posB4(:,2)]);
    dN = 2; dL = 2;
    [NN,LL] = meshgrid(min(posN):dN:max(posN),min(posL):dL:max(posL));
    fBM = griddata(posN,posL,posB,NN,LL);
    hold(hca,'on');
    mesh(hca,NN,NN*0+abs(min(posM)),LL,fBM);  
    hmcb = colorbar('peer',hca); 
    %plot3(hca,posN,posM,posL,posB,'o');
    hold(hca,'off');
  end
  
  
  if 1 % plot magnetosheath and magnetosphere boundary planes    
    gseVouteredge = 65.3*[0.37  0.32 -0.87];  
    gseVinneredge = 55*[-0.90 -0.28 -0.33]; % GSE
    gseVoutflow = 14.3*[-0.93 -0.13  0.35]; % GSE  
    gseMSP = 32*[-0.88 -0.42 0.24];

    lmnVmsh = gseVinneredge*[L' M' N'];
    lmnVmsp = gseMSP*[L' M' N'];
     
    mspN = irf_norm(lmnVmsp);
    mshN = irf_norm(lmnVmsh);

    x = 90*[-1 1];
    y = 90*[-1 1];
    z = 90*[-1 1];
    
    funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
    funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
    funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);


    if exist('hmshN1'); delete(hmshN1); end
    if exist('hmshN2'); delete(hmshN2); end
    if exist('hmspN1'); delete(hmspN1); end
    if exist('hmspN2'); delete(hmspN2); end
    if exist('ht1'); delete(ht1); end
    if exist('ht2'); delete(ht2); end
    if exist('ht3'); delete(ht3); end
    if exist('ht4'); delete(ht4); end
   
    hold(hca,'on')
    if 1
      hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+70,'k-.');
      hmspN1 = plot3(hca,funZ(x,y,mspN)-2.8,x*0,x,'k-');
    else
      hmshN1 = plot3(hca,funZ(x,y,mshN),x*0,x+72,'k-.');
      hmspN1 = plot3(hca,funZ(x,y,mspN)+20,x*0,x,'k-');
      hca.XLim = [0 50];
      hca.ZLim = [-10 40];
    end
    hold(hca,'off')

    ht1 = text(21.5,00,43,'MSH'); ht1.HorizontalAlignment = 'center'; ht1.FontSize = 13; ht1.Rotation = 55; 
    ht2 = text(2,0,50,'MSP'); ht2.HorizontalAlignment = 'center'; ht2.FontSize = 13; ht2.Rotation = -80;
    
    ht3 = text(16,0,-8,tintUTCstart(12:22)); ht3.HorizontalAlignment = 'center'; ht3.FontSize = 13;
    ht4 = text(16,0,72,tintUTCstop(12:22)); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
    
    hold(hca,'on') 
    quiver3(hca,17,0,4,0,0,5,2,'k')
    hold(hca,'off')
    ht4 = text(17,0,0,'time'); ht4.HorizontalAlignment = 'center'; ht4.FontSize = 13;
     
  end
end

%% Figure: eletron pa plots, mva fields
npanels = 10;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 4;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
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
if 1 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Te par perp
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
if 1 % eDist omni 64
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64
  pas = [65 105];
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.1],'fontsize',12,'color',[0 0 0]);
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
if 0 % ePDist pa 32
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
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 1 1]);
end
if 0 % ePDist pa 64
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
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E mms4');
  elim = [0 200];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all medium E mms4');
  elim = [200 1000];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:5]),'y')
for  ii = 6:8; h(ii).CLim = [5 8.5]; end
for  ii = 6:8; h(ii).YLim = [10 1000]; end
irf_plot_axis_align

h(1).Title.String = ['MMS' num2str(ic)];

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

%% Plot detailed electron distribution for selected times, MVA
% Initialize figure
[h1,h2] = initialize_combined_plot(9,2,2,0.4,'vertical');
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
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 1 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1200 700];  
end
if 1 % Ve perp par
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_{e}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  ylabel(hca,{'v_{e,\perp}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp
  hca = irf_panel('Ve perp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_{e,\perp}','(km/s)'};
  ylabel(hca,{'v_{e,\perp}','(km/s)'},'interpreter','tex');
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
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
  %c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint),mvaJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-800 1100];
end
if 0 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{mvaJcurl.x.tlim(tint),mvaJcurl.y.tlim(tint),mvaJcurl.z.tlim(tint),mvaJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
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
  elim = [20 200];
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
if 1 % ePDist pa 32
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
if 0 % ePDist pa 64
  %%
  hca = irf_panel('e PA e32 deflux high E');  
  elim = [200 400];
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).e64.pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
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
if 1 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')
end
if 0 % Q
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
  hca.YLabel.String = 'Q';
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.98 0.9],'fontsize',12);
end

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:4 8:9]),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

%% Plot single time particle distributions, 1 sc, 4 projections, LMN
time = irf_time('2015-10-16T10:34:32.83Z','utc>epochtt');
time = irf_time('2015-10-16T11:28:36.40Z','utc>epochtt');
time = irf_time('2015-09-15T11:20:31.74Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt');

c_eval('it = find(abs(ePDist?.time-time)==min(abs(ePDist?.time-time)));',ic);
c_eval('time = ePDist?(it).time;',ic)
timeUTC = time.utc;

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
perp1 = x;
perp2 = y;
par = z;

x = L;
y = M;
z = N;

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

isub = 1;

if 0 % Perpendicular plane    
  hca = h2(isub); isub = isub + 1; 
  xyz = [x;y;z]; vlabels = {'L','M','N'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);
end
if 0 % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [x;z;-y]; vlabels = {'L','N','-M'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
end
if 0 % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [y;z;x]; vlabels = {'M','N','L'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
end

if 1 % Perp plane, but with markings
  vectors = {hatExB,'ExB'; hatE,'E';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_E','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';
end
if 1
  hca = h2(isub); isub = isub + 1; 
  %c_eval('ploty = ePitch?(it).data(1,:,1)
  indPA = [1 8 15];
  c_eval('plot(hca,ePitch?(it).depend{1},squeeze(ePitch?(it).convertto(''s^3/km^6'').data(1,:,indPA)))',ic)
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = 'f_e (s^3/km^6)'; hca.XLabel.String = 'E (eV)';
  hca.YLim = [1e-1 1e5]; hca.XLim = [1e1 1e3];
  legend(hca,num2str(ePitch3.depend{2}(indPA)'),'location','southwest')
end

for ii = 1:3
  colormap(h2(ii),strCMap)
end

%% Plot single time particle distributions, 1 sc, 4 projections,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('gradPeslow = gseGradPe.resample(gseVe?);',ic)

% Projection coordinate system

%y = hatE;
%z = hatB;
%x = cross(y,z)/norm(cross(y,z));

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 886:1:920;953; 911;
%tind = 920;
for it = tind
  time = times(it);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  
  c_eval('hatGradPe = gradPeslow.resample(time).data/gradPeslow.resample(time).abs.data;',ic)
  c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
  c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);
  
  c_eval('hatV = dbcsVe?(it).data/dbcsVe?(it).abs.data;',ic)

  

  timeUTC = time.utc;      
  isub = 1;

  % Perpendicular plane    
  hca = h2(isub); isub = isub + 1; 
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';hatGradPe,'\nabla P_e'};
  xyz = [x;y;z]; vlabels = {'x_{GSE}','y_{GSE}','z_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [x;z;-y]; vlabels = {'x_{GSE}','z_{GSE}','-y_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [y;z;x]; vlabels = {'y_{GSE}','z_{GSE}','x_{GSE}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % ExB plane, with markings
  vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z';hatV,'v_e'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = '';

  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','hard_path',eventPath)
end

%%
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
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [20 200];
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
if 1 % ePDist pa 32
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

if exist('hmark'); delete(hmark); end
c_eval('times = ePDist?.time;',ic)
tind = 902:1:(960+9*2);

nrows = 3;
ncols = 3;
rowshift = 1;

for ii = 1:nrows*ncols
  h2(ii) = subplot(rowshift+nrows,ncols,rowshift*ncols+ii);
end
colors = ...
   [     0    0.4470    0.7410;...
    0.8500    0.3250    0.0980;...
    0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;...
    0.4660    0.6740    0.1880;...
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840];
  
isub = 1;
for ii=1:nrows*ncols  
  it = tind(ii);
  time = times(it);
  timeUTC = time.utc; 
  hmark = irf_pl_mark(h1(1),time.epochUnix','green');
  if 1
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
legend(h2(1),num2str(ePitch3.depend{2}(indPA)'),'location','southwest')


