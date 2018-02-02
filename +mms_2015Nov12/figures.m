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

cmap = colormap('jet');

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

%% Figure 3: 4 sc, Ohm's law, GSE 
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

%% Figure 3: 4 sc, Ohm's law, MVA 
ic = 1;
npanels = 10;
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
if 1 % Pressures
  hca = irf_panel('Pressure');
  set(hca,'ColorOrder',mms_colors('1234b'))  
  c_eval(['irf_plot(hca,{mvaPe?.trace/3,facPe?.xx,0.5*(facPe?.yy+facPe?.zz),'...
                         'mvaPi?.trace/3,PB?},''comp'');'],ic)
  hca.YLabel.String = {'P','(nPa)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'P_e','P_{e||}','P_{e\perp}','P_i','P_B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{['mms' num2str(ic)]},[0.05 0.9],'fontsize',12);
end
if 1 % Pressures
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
if 1 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp abs
  hca = irf_panel('Ve perp abs');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).abs,mvaVe2perp.tlim(tint).abs,mvaVe3perp.tlim(tint).abs,mvaVe4perp.tlim(tint).abs},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 1 % J curl
  hca = irf_panel('J curl');
  set(hca,'ColorOrder',mms_colors('xyz'))
  lines = irf_plot(hca,{mvaJcurl.x,mvaJcurl.y,mvaJcurl.z},'comp');   
  hca.YLabel.String = {'J_{curl}','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 1 % Ohm's law x
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
if 1 % Ohm's law y
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y,-mvaOhmVixB.y,mvaOhmJxB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law z
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
[h1,h2] = initialize_combined_plot(9,2,2,0.4,'vertical');
cmap = 'jet';
ic = 1;
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
if 1% Ve perp par
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','v_{e,||}'},[0.98 0.9],'fontsize',12);  
  hca.YLabel.String = {'v_{e}','(km/s)'};
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
if 1 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseGradPe.x*1e3,gseGradPe.y*1e3,gseGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
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

%irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:4 6:7]),'y')
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

%% Plot single time particle distributions, 1 sc, 4 projections,
time = irf_time('2015-10-16T10:34:32.83Z','utc>epochtt');
time = irf_time('2015-10-16T11:28:36.40Z','utc>epochtt');
time = irf_time('2015-09-15T11:20:31.74Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:20.60Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:17.70Z','utc>epochtt');
%time = irf_time('2015-11-12T07:19:18.85Z','utc>epochtt');
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

%x = [1 0 0];
%y = [0 1 0];
%z = [0 0 1];

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
ic = 1;

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
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)

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
if 1 % ePDist pa 32 low E
  hca = irf_panel('e PA e32 deflux all low E mms4');
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
  hca.CLim = [7.5 8.3];
end
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all medium E mms4');
  elim = [400 1000];
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
if 1 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h([1:4 9]),'y')
for  ii = 5:7; h(ii).CLim = [5 8.5]; end
for  ii = 5:7; h(ii).YLim = [10 1000]; end
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

hca = irf_panel('B'); hca.YLim = [-15 15];

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 13;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Plot detailed electron distribution for selected times, MVA
% Initialize figure
[h1,h2] = initialize_combined_plot(9,2,2,0.4,'vertical');
cmap = 'jet';
ic = 2;
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
time = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:33.10Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:21.45Z','utc>epochtt');

c_eval('it = find(abs(ePDist?.time-time)==min(abs(ePDist?.time-time)));',ic);
c_eval('time = ePDist?(it).time;',ic)
timeUTC = time.utc;

if exist('hmark'); delete(hmark); end
hmark = irf_pl_mark(h1,time.epochUnix','green');
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('dbcsVe?slow = dslE?.resample(gseVe?);',ic)

c_eval('hatVe = dbcsVe?.resample(time); hatVe = hatVe.data/norm(hatVe.data);',ic)

% Projection coordinate system
c_eval('hatE = dslE?slow.resample(time).data/dslE?slow.resample(time).abs.data;',ic)
c_eval('hatB = dmpaB?slow.resample(time).data/dmpaB?slow.resample(time).abs.data;',ic)
c_eval('hatExB = cross(hatE,hatB);',ic)

y = hatExB;
z = hatB;
x = cross(y,z)/norm(cross(y,z));
perp1 = hatExB;
perp2 = cross(hatExB,hatB)/norm(cross(hatExB,hatB));
par = hatB;

%x = L;
%y = M;
%z = N;

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 4.5];  
  

isub = 1;

if 1 % perp1 perp2
  hca = h2(isub); isub = isub + 1; 
  vectors = {hatVe,'Ve';hatB,'B'};
  xyz = [perp1;perp2;par]; vlabels = {'\perp 1','\perp 2','||'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);
end
if 1 % par perp1
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp2;par;perp1]; vlabels = {'\perp2','||','\perp1'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
end
if 1 % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [perp1;par;-perp2]; vlabels = {'\perp1','||','-\perp2'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
end

if 0 % Perp plane, but with markings
  vectors = {hatExB,'ExB'; hatE,'E';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_E','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';
end
if 1
  hca = h2(isub); isub = isub + 1; 
  %c_eval('ploty = ePitch?(it).data(1,:,1)
  indPA = [1 9 10 18];
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
ic = 4;
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

%% Figure: eletron pa plots, mva fields
npanels = 11;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 3;
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
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
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
  hca = irf_panel('e PA e32 deflux all low E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 200];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
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
if 1 % potential in normal direction
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
irf_zoom(h([1:4 9:end]),'y')
for  ii = 5:7; h(ii).CLim = [5 8.5]; end
for  ii = 5:7; h(ii).YLim = [10 1000]; end
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

%% Figure: eletron pa plots, mva fields + projections
npanels = 9;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
tintZoom = tint;
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
[h1,h2] = initialize_combined_plot(npanels,2,2,0.4,'vertical');
ic = 4;
pshift = 0;

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
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
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
  hca = irf_panel('e PA e32 deflux all low E mms4');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 200];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
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

irf_zoom(h1,'x',tintZoom)
irf_zoom(h1([1:4 9:end]),'y')
%hca = irf_panel('B'); hca.YLim = [-15 15];
for  ii = 5:7; h1(ii).CLim = [5 8.5]; end
for  ii = 5:7; h1(ii).YLim = [10 1000]; end
irf_plot_axis_align

h1(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h1(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h1(ii+pshift).FontSize = 12;  
  h1(ii+pshift).YLabel.FontSize = 11;
end
for ii = 5:8; h1(ii).XGrid = 'off'; h1(ii).YGrid = 'off'; end
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Plot single time particle distributions, 1 sc, 4 projections, LMN
time = irf_time('2015-10-16T10:34:32.83Z','utc>epochtt');
time = irf_time('2015-10-16T11:28:36.40Z','utc>epochtt');
time = irf_time('2015-09-15T11:20:31.74Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt');
time = irf_time('2015-11-12T07:19:33.10Z','utc>epochtt');
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

x = L;
y = M;
z = N;

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 1208;886:1:920;953; 911;906;911;
%tind = 886;
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

  % NL   
  hca = h2(isub); isub = isub + 1; 
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';hatGradPe,'\nabla P_e'};
  xyz = [N;L;M]; vlabels = {'N','L','M'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % ML
  hca = h2(isub); isub = isub + 1; 
  xyz = [M;L;-N]; vlabels = {'M','L','-N'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % NM
  hca = h2(isub); isub = isub + 1;
  xyz = [N;M;-L]; vlabels = {'N','M','-L'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % ExB plane, with markings
  vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';L,'L';M,'M';N,'N';hatV,'v_e'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = '';

  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  %cn.print(['lmn_eproj_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','hard_path',[eventPath 'lmn_eproj/mms' num2str(ic) '/'])
end

%% Plot single time particle distributions, 1 sc, skymaps,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
%c_eval('dist = ePDist?;',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('gradPeslow = gseGradPe.resample(gseVe?);',ic)

% Projection coordinate system

%y = hatE;
%z = hatB;
%x = cross(y,z)/norm(cross(y,z));

x = L;
y = M;
z = N;

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 886:1:920;953; 911;906;911;
%tind = 886;
tind = 905;
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

  energies = [40 70 120];

  timeUTC = time.utc;      
  isub = 1;
  
  
  % NL   
  hca = h2(isub); isub = isub + 1; 
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';hatGradPe,'\nabla P_e'};
  xyz = [N;L;M]; vlabels = {'N','L','M'};
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';};
  mms.plot_skymap(hca,dist,'tint',times(it),'flat','energy',energies(1),'vectors',vectors)
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % ML
  hca = h2(isub); isub = isub + 1; 
  xyz = [M;L;-N]; vlabels = {'M','L','-N'};
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';};
  mms.plot_skymap(hca,dist,'tint',times(it),'flat','energy',energies(2),'vectors',vectors)
  hca.Title.String = '';

  % NM
  hca = h2(isub); isub = isub + 1;
  xyz = [N;M;-L]; vlabels = {'N','M','-L'};
  vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B';hatV,'v_e';};
  mms.plot_skymap(hca,dist,'tint',times(it),'flat','energy',energies(3),'vectors',vectors)
   hca.Title.String = '';

  % ExB plane, with markings
  vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';L,'L';M,'M';N,'N';hatV,'v_e'};
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = '';

  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  for ii = 1:3
    h2(ii).YLim = [0 180];
    h2(ii).XLim = [0 360];
  end
  
  %cn.print(['lmn_skymap_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','hard_path',[eventPath 'lmn_flatskymap/mms' num2str(ic) '/'])
end

%% Figure 1: 1 sc overview electric field
%tint = irf.tint('2015-10-16T10:34:20.00Z/2015-10-16T10:35:30.00Z');
npanels = 11;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
if 1 % B gse
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 0 % B mva
  hca = irf_panel('B mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 12];
end
if 0 % Vi mva
  hca = irf_panel('Vi mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint),mvaVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
end
if 0 % Vi gse
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
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
if 0 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.deflux.tlim(tint).pitchangles(dmpaB?,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 0 % Ve mva
  hca = irf_panel('Ve mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','M'},[0.98 0.9],'fontsize',12);    
end
if 0 % Ve perp par mva
  hca = irf_panel('Ve perp par mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L,\perp','M,\perp','N,\perp','||'},[0.98 0.9],'fontsize',12);    
end
if 0 % Ve gse
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 1 % Ve perp par gse
  hca = irf_panel('Ve perp par gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x,\perp','y,\perp','z,\perp','||'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % V ExB LMN
  hca = irf_panel('V ExB mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x.tlim(tint),mvaVExB?.y.tlim(tint),mvaVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);    
end
if 1 % V ExB GSE
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.tlim(tint),gseVExB?.y.tlim(tint),gseVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % Te par perp
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
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
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
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseVe1par.tlim(tint),gseVe2par.tlim(tint),gseVe3par.tlim(tint),gseVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % E mva
  hca = irf_panel('E mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
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
if 0 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.x.tlim(tint),mvaE2.x.tlim(tint),mvaE3.x.tlim(tint),mvaE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{L}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 0 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.y.tlim(tint),mvaE2.y.tlim(tint),mvaE3.y.tlim(tint),mvaE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{M}','(mV/m)'};
end
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.z.tlim(tint),mvaE2.z.tlim(tint),mvaE3.z.tlim(tint),mvaE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{N}','(mV/m)'};
end
if 1 % E DSL X
  hca = irf_panel('E DSL X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.x.tlim(tint),dslE2.x.tlim(tint),dslE3.x.tlim(tint),dslE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{X,DSL}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % E DSL Y
  hca = irf_panel('E DSL Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{gseE1.y.tlim(tint),dslE2.y.tlim(tint),dslE3.y.tlim(tint),dslE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Y,DSL}','(mV/m)'};
end
if 1 % E DSL Z
  hca = irf_panel('E DSL Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.z.tlim(tint),dslE2.z.tlim(tint),dslE3.z.tlim(tint),dslE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Z,DSL}','(mV/m)'};
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % Ohm's law L, only electron terms
  hca = irf_panel('Ohm L');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x,mvaAvE.x+mvaOhmVexB.resample(mvaAvE).x+mvaOhmGradPe.resample(mvaAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_L','(mV/m)'};
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law LMN, only electron terms residul
  hca = irf_panel('Ohm LMN residual');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  residual = mvaAvE + mvaOhmVexB.resample(mvaAvE) + mvaOhmGradPe.resample(mvaAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E+vexB+divP/ne','(mV/m)'};
end
if 1 % Ohm's law L, only electron terms
  hca = irf_panel('Ohm X');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  irf_plot(hca,{gseAvE.x,gseOhmGradPe.x,-gseOhmVexB.x},'comp');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E_Z','(mV/m)'};
end
if 1 % Ohm's law Y
  hca = irf_panel('Ohm Y');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseAvE.y,gseOhmGradPe.y,-gseOhmVexB.y},'comp');
  hca.YLabel.String = {'E_Y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law Z
  hca = irf_panel('Ohm Z');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseAvE.z,gseOhmGradPe.z,-gseOhmVexB.z},'comp');
  hca.YLabel.String = {'E_Z','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 1 % Ohm's law LMN, only electron terms residul
  hca = irf_panel('Ohm XYZ residual');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  residual = gseAvE + gseOhmVexB.resample(gseAvE) + gseOhmGradPe.resample(gseAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'X','Y','Z'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E+vexB+divP/ne','(mV/m)'};
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

%% Figure 1: 1 sc overview electric field, spacecraft coord. sys.
%tint = irf.tint('2015-10-16T10:34:20.00Z/2015-10-16T10:35:30.00Z');
npanels = 9;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B_{DMPA}','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
end
if 1 % spacecraft potential
  hca = irf_panel('sc pot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{scPot1*(-1),scPot2*(-1),scPot3*(-1),scPot4*(-1)},'comp');
  hca.YLabel.String = {'-V_{SC}','(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))    
end
if 0 % Vi gse
  hca = irf_panel('Vi gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'V_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % Ve perp par mva
  hca = irf_panel('Ve perp par mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint),mvaVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L,\perp','M,\perp','N,\perp','||'},[0.98 0.9],'fontsize',12);    
end
if 0 % Ve gse
  hca = irf_panel('Ve gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Ve perp par gse
  hca = irf_panel('Ve perp par gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x,\perp','y,\perp','z,\perp','||'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % V ExB
  hca = irf_panel('V ExB mva');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVExB?.x.tlim(tint),mvaVExB?.y.tlim(tint),mvaVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);    
end
if 0 % V ExB
  hca = irf_panel('V ExB');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVExB?.x.tlim(tint),gseVExB?.y.tlim(tint),gseVExB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{ExB}','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
end
if 0 % Te par perp
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
  %irf_zoom(hca,'y')
end
if 1 % E DSL X
  hca = irf_panel('E DSL X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.x.tlim(tint),dslE2.x.tlim(tint),dslE3.x.tlim(tint),dslE4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{X,DSL}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % E DSL Y
  hca = irf_panel('E DSL Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.y.tlim(tint),dslE2.y.tlim(tint),dslE3.y.tlim(tint),dslE4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Y,DSL}','(mV/m)'};
end
if 1 % E DSL Z
  hca = irf_panel('E DSL Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{dslE1.z.tlim(tint),dslE2.z.tlim(tint),dslE3.z.tlim(tint),dslE4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'E_{Z,DSL}','(mV/m)'};
end
if 1 % ve x B GSE X
  hca = irf_panel('Ve x B GSE X');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVexB1.x.tlim(tint),-gseVexB2.x.tlim(tint),-gseVexB3.x.tlim(tint),-gseVexB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e x B_{X,GSE}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % Ve x B GSE Y
  hca = irf_panel('Ve x B GSE Y');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVexB1.y.tlim(tint),-gseVexB2.y.tlim(tint),-gseVexB3.y.tlim(tint),-gseVexB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e x B_{Y,GSE}','(mV/m)'};
end
if 1 % Ve x B GSE Z
  hca = irf_panel('Ve x B GSE Z');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{-gseVexB1.z.tlim(tint),-gseVexB2.z.tlim(tint),-gseVexB3.z.tlim(tint),-gseVexB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'-v_e x B_{Z,GSE}','(mV/m)'};
end
if 1 % gradPe GSE
  hca = irf_panel('gradPe/ne gse mv/m');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseOhmGradPe.x,gseOhmGradPe.y,gseOhmGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e/ne','(mv/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);    
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
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
if 0 % Ohm's law L, only electron terms
  hca = irf_panel('Ohm L');
  set(hca,'ColorOrder',mms_colors('1234ba'))
  if 0
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x,-mvaOhmVixB.x,mvaOhmJxB.x,mvaAvE.x+mvaOhmVexB.resample(mvaAvE).x+mvaOhmGradPe.resample(mvaAvE).x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B','-v_i\times B','J\times B/ne','residual'},[0.98 0.9],'fontsize',12);
  else
    irf_plot(hca,{mvaAvE.x,mvaOhmGradPe.x,-mvaOhmVexB.x},'comp');
    set(hca,'ColorOrder',mms_colors('1234ba'))
    irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  end
  hca.YLabel.String = {'E_L','(mV/m)'};
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.y,mvaOhmGradPe.y,-mvaOhmVexB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
end
if 0 % Ohm's law LMN, only electron terms residul
  hca = irf_panel('Ohm LMN residual');
  set(hca,'ColorOrder',mms_colors('xyz'))  
  residual = mvaAvE + mvaOhmVexB.resample(mvaAvE) + mvaOhmGradPe.resample(mvaAvE);
  irf_plot(hca,{residual.x,residual.y,residual.z},'comp');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLabel.String = {'E+vexB+divP/ne','(mV/m)'};
end

c_eval('h(?).YLim = [-4 4];',3:9)
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

%% Figure 2: 4 sc overview electric field potentials, spacecraft coord. sys.
tint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:35.00Z');
npanels = 10;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
ic = 1;
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B (DMPA)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,irf_ssub('MMS ?',ic),[0.02,0.1])
end
if 1 % spacecraft potential
  hca = irf_panel('sc pot');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{scPot1*(-1),scPot2*(-1),scPot3*(-1),scPot4*(-1)},'comp');  
  hca.YLabel.String = {irf_ssub('-V_{SC}',ic),'(V)'};
  set(hca,'ColorOrder',mms_colors('1234'))    
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.1],'fontsize',12);
end
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end
ic = 2;
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end
ic = 3;
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end
ic = 4;
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? DSL',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  hca.YLabel.String = {irf_ssub('DCV_{?}',ic),'(V)'};  
  set(hca,'ColorOrder',mms_colors('1234'))    
end


irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Figure 2: 1 sc overview electric field potentials, spacecraft coord. sys.
%tint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:35.00Z');
npanels = 5;
cmap = 'jet';

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Plot
ic = 4;
if 1 % B dmpa
  hca = irf_panel('B dmpa');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{dmpaB?.x.tlim(tint),dmpaB?.y.tlim(tint),dmpaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B (DMPA)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % spacecraft potential
  hca = irf_panel('sc pot');
  set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,scPot1*(-1));',ic)  
  ylabel(hca,{irf_ssub('-V_{SC}',ic),'(V)'},'interpreter','tex')
  %hca.YLabel.String = {irf_ssub('-V_{SC}',ic),'(V)'};    
end
if 1 % probe to spacecraft potentials
  hca = irf_panel(irf_ssub('dcv ?',ic));
  %set(hca,'ColorOrder',mms_colors('1234'))
  c_eval('irf_plot(hca,dcv?);',ic)
  ylabel(hca,{'DCV','(V)'},'interpreter','tex')
  irf_legend(hca,{'1','2','3','4','5','6'},[0.98 0.9],'fontsize',12);  
  %set(hca,'ColorOrder',mms_colors('1234'))    
end
if 1 % E1 dsl 
  hca = irf_panel(irf_ssub('E?',ic));
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{dslE?.x.tlim(tint),dslE?.y.tlim(tint),dslE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {irf_ssub('E_? (DSL)',ic),'(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
end
if 1 % difference between axial probes
  hca = irf_panel(irf_ssub('V6-V5',ic));
  baselength = 14.5;
  c_eval('axdiff = irf.ts_scalar(dcv?.time,(dcv?.data(:,6)-dcv?.data(:,5))/baselength*1e3);',ic)
  set(hca,'ColorOrder',mms_colors('z'))
  c_eval('irf_plot(hca,axdiff);',ic)
  ylabel(hca,{sprintf('(V_6 - V_5)/%.1f',baselength),'(mV/m)'},'interpreter','tex')
  %hca.YLabel.String = {sprintf('(V_6 - V_5)/%.1f',baselength),'(mV/m)'};    
end

h(1).Title.String = irf_ssub('MMS ?',ic);

irf_zoom(h,'x',tint)
%irf_zoom(h([1:5]),'y')
irf_plot_axis_align


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots
pshift = 0;
for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 12;  
  h(ii+pshift).YLabel.FontSize = 11;
end

%% Figure: Waves
tintSpec = '2015-11-12T07:19:05.000Z/2015-11-12T07:19:45.000Z';
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
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-25 35];
end
if 0 % n
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?},''comp'');',ic)
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 12];
end
if 0 % iPDist omni 64
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.tlim(tint).deflux.e64.omni.specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap) 
end
if 0 % iPDist pa 32
  hca = irf_panel('i PA deflux');  
  c_eval('irf_spectrogram(hca,iPDist?.deflux.tlim(tint).pitchangles(dmpaB?,18).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];
  colormap(hca,cmap) 
end
if 0 % Ve
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  hca.YLim = [-1100 1100];  
end
if 0 % Te par perp
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
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
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
if 1 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?_new.x.tlim(tint),gseE?_new.y.tlim(tint),gseE?_new.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E spectrogram par
  hca = irf_panel('E spectrogram par');
  c_eval('tE = irf_time(wavE?fac.t,''epoch>epochtt'');',ic)
  %iE = tE.tlim(tintZoom);
  iE = 1:10:tE.length;
  c_eval('plotWavE = wavE?fac; plotWavE.t = plotWavE.t(iE); plotWavE.p = {plotWavE.p{1}(iE,:)};',ic)
  c_eval('irf_spectrogram(hca,plotWavE,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.8],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end
if 1 % E spectrogram perp
  hca = irf_panel('E spectrogram perp');
  c_eval('tE = irf_time(wavE?fac.t,''epoch>epochtt'');',ic)
  %iE = tE.tlim(tintZoom);
  iE = 1:10:tE.length;
  c_eval('plotWavE = wavE?fac; plotWavE.t = plotWavE.t(iE); plotWavE.p = {0.5*(plotWavE.p{2}(iE,:)+plotWavE.p{3}(iE,:))};',ic)
  c_eval('irf_spectrogram(hca,plotWavE,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.8],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end
if 1 % B scm
  hca = irf_panel('B scm');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?scm.x.tlim(tint),gseB?scm.y.tlim(tint),gseB?scm.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B (SCM)','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % B spectrogram par
  hca = irf_panel('b spectrogram par');
  c_eval('tB = irf_time(wavB?fac.t,''epoch>epochtt'');',ic)
  %iE = tE.tlim(tintZoom);
  iB = 1:10:tB.length;
  c_eval('plotWavB = wavB?fac; plotWavB.t = plotWavB.t(iB); plotWavB.p = {plotWavB.p{1}(iB,:)};',ic)
  c_eval('irf_spectrogram(hca,plotWavB,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.8],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
end
if 1 % B spectrogram perp
  hca = irf_panel('B spectrogram perp');
  c_eval('tB = irf_time(wavB?fac.t,''epoch>epochtt'');',ic)
  %iE = tE.tlim(tintZoom);
  iB = 1:10:tB.length;
  c_eval('plotWavB = wavB?fac; plotWavB.t = plotWavB.t(iB); plotWavB.p = {0.5*(plotWavB.p{2}(iB,:)+plotWavB.p{3}(iB,:))};',ic)
  c_eval('irf_spectrogram(hca,plotWavB,''log'',''donotfitcolorbarlabel'');',ic)
  hca.YScale = 'log';
  hold(hca,'on')
  c_eval('hflh=irf_plot(hca,flh?,''black'');',ic)
  c_eval('hfce=irf_plot(hca,fce?,''cyan'');',ic)
  c_eval('hfpp=irf_plot(hca,fpp?,''green'');',ic)
  irf_legend(hca,'f_{lh}',[0.02 0.2],'color','black');
  irf_legend(hca,'f_{ce}',[0.02 0.99],'color','cyan');
  irf_legend(hca,'f_{pp}',[0.10 0.8],'color','green');
  hold(hca,'off')
  hca.YLabel.String = 'f [Hz]';
  hca.CLim = [-9 0]; 
  hca.YTick = [1e0 1e1 1e2 1e3 1e4];
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

%% Plot overview figure with focus on electrons, including single time electron distributions
npanels = 8;
tintZoom = tint;irf.tint('2016-12-25T16:00:53.00Z/2016-12-25T16:00:58.00Z');
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
tintZoom = irf.tint('2015-11-12T07:19:19.20Z/2015-11-12T07:19:23.10Z');
cmap = 'jet';
[h1,h2] = initialize_combined_plot(npanels,2,2,0.4,'vertical');
ic = 1;
iisub = 0;
cmap = colormap('jet');

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
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{ne?,ni?},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))  
  irf_legend(hca,{'n_e','n_i'},[0.98 0.8],'fontsize',12);
end
if 0 % J  
  hca = irf_panel('J fpi');
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
end
if 1 % Ve
  hca = irf_panel('Ve perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?perp.x.tlim(tint),gseVe?perp.y.tlim(tint),gseVe?perp.z.tlim(tint),gseVe?par.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x_{\perp}','y_{\perp}','z_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
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
if 1 % e DEF omni 64
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hold(hca,'on')
  hca.XGrid = 'off'; hca.YGrid = 'off';
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hca.YLabel.String = {'E_e','(eV)'};   
  colormap(hca,cmap) 
end
if 1 % ePDist pa 64
  hca = irf_panel('e PA e64 deflux lowe');  
  eint = [10 400];  
  c_eval('irf_spectrogram(hca,ePDist?.tlim(tint).pitchangles(dmpaB?,20).elim(eint).deflux.specrec(''pa''),''log'');',ic)
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  %irf_legend(hca,['E = [' num2str(eint(1),'%.0f') ' ' num2str(eint(2),'%.0f') ']'],[0.95 0.90],'color',0*[1 1 1])
  irf_legend(hca,[num2str(eint(1),'%.0f') '<E<' num2str(eint(2),'%.0f')],[0.99 0.90],'color',0*[1 1 1])
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('123'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.98 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % E
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?.x.tlim(tint),gseE?.y.tlim(tint),gseE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % E new
  hca = irf_panel('E');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseE?_new.x.tlim(tint),gseE?_new.y.tlim(tint),gseE?_new.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 1 % OhmgradPe
  hca = irf_panel('Ohm gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{gseOhmGradPe.x,gseOhmGradPe.y,gseOhmGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))  
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);      
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:npanels]  
  irf_legend(h1(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  nInd = nInd + 1;
end

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h1(1:npanels),'x',tintZoom + [-2 0])
irf_zoom(h1([1:3 6:npanels]),'y')
hca = irf_panel('E'); hca.YLim = [-4.9 4.9];
hca = irf_panel('e DEF omni 64'); hca.YLim = [10 2000];
hca = irf_panel('Ohm gradPe'); hca.YLim = [-2.4 2.2];

irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);
if 0
hmark = irf_pl_mark(h1(1:6),tintBCS, 'yellow');
for ii = 1:numel(hmark);
  hmark(ii).FaceAlpha = 0.5;
end
end
for ii = 1:npanels;
  h1(ii).FontSize = 12;
end

%% Plot single time particle distributions, 1 sc, 4 projections,
  
c_eval('dist = ePDist?.convertto(''s^3/km^6'');',ic)
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?_new.resample(gseVe?);',ic)
c_eval('ePitch = ePitch?;',ic)

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  

c_eval('times = ePDist?.time;',ic)
tind = 786;%:1:920;953; 911;
tind = 1160 ;
tind = 1240
tind = 781:800;881:930;
%tind = 880;
for it = tind;589%:650;550:650;%1240;tind;
  time = times(it);
  if exist('hmark'); delete(hmark); end
  hmark = irf_pl_mark(h1,time.epochUnix','green');
  
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
  xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  %mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);        
  hca.Title.String = '';
  hca.Title.String = timeUTC(1:23);

  % B plane 1
  hca = h2(isub); isub = isub + 1; 
  xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{||}','-v_{perp2}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % B plane 2
  hca = h2(isub); isub = isub + 1;
  xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{ExB}'};
  mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
  hca.Title.String = '';

  % Pitchangle distribution
  hca = h2(isub); isub = isub + 1;
  plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
  hca.YScale = 'log'; hca.XScale = 'log';
  hca.YLabel.String = ['f_e (' ePitch.units ')'];
  hca.XLabel.String = 'E (eV)';
  hca.XLim = [10 1000];
  legend(hca,{'0','90','180'})
  hca.YTick = 10.^[0:5];
  hca.YLim = [1e0 1e5];

  % ExB plane, with markings
  if 0
    vectors = {hatExB,'ExB'; hatE,'E';hatB,'B';[1 0 0],'x';[0 1 0],'y';[0 0 1],'z'};
    hca = h2(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{E,\perp}','v_B'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
    hca.Title.String = '';
  end
  for ii = 1:4
    colormap(h2(ii),strCMap)
  end
  cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',[eventPath '/time_dist/'])
end

%% Figure: reconnection signatures, mva fields
npanels = 12;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize')
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)

if 1 % B no abs
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 0 % B with abs
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
  c_eval('irf_plot(hca,{ne?,ni?*0.96},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 1 % J
  hca = irf_panel('J');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);     
end
if 0 % Ve
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
  hca.YLabel.String = {'v_e','(km/s)'};
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 0 % ePDist pa 32 low E
  hca = irf_panel('e PA e32 deflux all low E mms4');
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
  hca.CLim = [7.3 8.2];
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all medium E mms4');
  elim = [400 1000];
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
if 0 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end
if 0 % E perp par new
  hca = irf_panel('E perp par new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp_new.x.tlim(tint),mvaE?perp_new.y.tlim(tint),mvaE?perp_new.z.tlim(tint),mvaE?par_new},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end
if 1 % E perp par new
  hca = irf_panel('E new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?_new.x.tlim(tint),mvaE?_new.y.tlim(tint),mvaE?_new.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end
if 1 % V ExB LMN
  hca = irf_panel('V ExB mva new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaEVexB?_new.x.tlim(tint),mvaEVexB?_new.y.tlim(tint),mvaEVexB?_new.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E+v_exB','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);    
end
if 1 % gradPe/ne
  hca = irf_panel('Grad Pe/ne');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaOhmGradPe.x,mvaOhmGradPe.y,mvaOhmGradPe.z},'comp');
  hca.YLabel.String = {'\nabla\cdotP_e/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4 sc'},[0.06 0.9],'fontsize',12,'color',[0 0 0]);
end
if 1 % Ohm's law M
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE_new.y,mvaOhmGradPe.y,-mvaOhmVexB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4 sc'},[0.06 0.9],'fontsize',12,'color',[0 0 0]);
end
if 1 % Ohm's law N
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE_new.z,mvaOhmGradPe.z,-mvaOhmVexB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4 sc'},[0.06 0.9],'fontsize',12,'color',[0 0 0]);
end
if 0 % agyrotropy
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  if 1
    c_eval('irf_plot(hca,{Q?.^0.5},''comp'');',ic)
    hca.YLabel.String = '\sqrt{Q}';
  else
    irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
    hca.YLabel.String = 'Q';
  end
  hca.YLabel.String = 'Q';
end
if 1 % Epar new
  hca = irf_panel('E par new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?par_new.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};    
end
if 1 % E dot J
  hca = irf_panel('E dot J');
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval('irf_plot(hca,{EdotJ?_new,EdotJ?par_new.tlim(tint),EdotJ?perp_new,RedotJ?_new},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{'E \cdot J','E_{||}\cdot J_{||}','E_{\perp}\cdot J_{\perp}','(E+v_exB)\cdot J'},[0.98 0.9],'fontsize',12);      
  hca.YLabel.String = {'E \cdot J','(nW/m^3)'};  
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
%hca = irf_panel('e PA e32 deflux all low E mms4'); hca.YLim = [0 180];
%for  ii = 5:7; h(ii).CLim = [5 8.5]; end
%for  ii = 5:7; h(ii).YLim = [10 1000]; end
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

hca = irf_panel('B'); hca.YLim = [-15 15];

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 13;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Figure: reconnection signatures 2, mva fields
npanels = 7;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize')
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)

if 1 % B no abs
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-25 35];
end
if 0 % B with abs
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
  c_eval('irf_plot(hca,{ne?,ni?*0.96},''comp'');',ic)
  hca.YLabel.String = {'n','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'n_e','n_i'},[0.98 0.9],'fontsize',12);
  hca.YLim = [0 12];
end
if 0 % J
  hca = irf_panel('J');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);     
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
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta_B<' num2str(pas(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % eDist omni 64
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
if 1 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all low E');
  elim = [0 1000];
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
if 0 % ePDist pa 32 low E
  hca = irf_panel('e PA e32 deflux all low E mms4');
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.90],'fontsize',12,'color',[0 0 0]);
  hca.CLim = [7.3 8.2];
end
if 0 % ePDist pa 32
  hca = irf_panel('e PA e32 deflux all medium E mms4');
  elim = [400 1000];
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
if 0 % E
  hca = irf_panel('E perp par');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp.x.tlim(tint),mvaE?perp.y.tlim(tint),mvaE?perp.z.tlim(tint),mvaE?par},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end
if 0 % E perp par new
  hca = irf_panel('E perp par new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?perp_new.x.tlim(tint),mvaE?perp_new.y.tlim(tint),mvaE?perp_new.z.tlim(tint),mvaE?par_new},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','||'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end
if 1 % E perp par new
  hca = irf_panel('E new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?_new.x.tlim(tint),mvaE?_new.y.tlim(tint),mvaE?_new.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  irf_zoom(hca,'y')  
end
if 0 % V ExB LMN
  hca = irf_panel('V ExB mva new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaEVexB?_new.x.tlim(tint),mvaEVexB?_new.y.tlim(tint),mvaEVexB?_new.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E+v_exB','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);    
end
if 0 % gradPe/ne
  hca = irf_panel('Grad Pe/ne');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaOhmGradPe.x,mvaOhmGradPe.y,mvaOhmGradPe.z},'comp');
  hca.YLabel.String = {'\nabla\cdotP_e/ne','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4 sc'},[0.06 0.9],'fontsize',12,'color',[0 0 0]);
end
if 0 % Ohm's law M
  hca = irf_panel('Ohm M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE_new.y,mvaOhmGradPe.y,-mvaOhmVexB.y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4 sc'},[0.06 0.9],'fontsize',12,'color',[0 0 0]);
end
if 0 % Ohm's law N
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaAvE_new.z,mvaOhmGradPe.z,-mvaOhmVexB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'E','\nabla P_e','-v_e\times B'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4 sc'},[0.06 0.9],'fontsize',12,'color',[0 0 0]);
end
if 0 % agyrotropy
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  if 1
    c_eval('irf_plot(hca,{Q?.^0.5},''comp'');',ic)
    hca.YLabel.String = '\sqrt{Q}';
  else
    irf_plot(hca,{Q1.tlim(tint).abs,Q2.tlim(tint).abs,Q3.tlim(tint).abs,Q4.tlim(tint).abs},'comp');    
    hca.YLabel.String = 'Q';
  end
  hca.YLabel.String = 'Q';
end
if 1 % Epar new
  hca = irf_panel('E par new');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?par_new.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E_{||}','(mV/m)'};    
end
if 0 % E dot J
  hca = irf_panel('E dot J');
  set(hca,'ColorOrder',mms_colors('matlab'))  
  c_eval('irf_plot(hca,{EdotJ?_new,EdotJ?par_new.tlim(tint),EdotJ?perp_new,RedotJ?_new},''comp'');',ic)
  set(hca,'ColorOrder',mms_colors('matlab'))  
  irf_legend(hca,{'E \cdot J','E_{||}\cdot J_{||}','E_{\perp}\cdot J_{\perp}','(E+v_exB)\cdot J'},[0.98 0.9],'fontsize',12);      
  hca.YLabel.String = {'E \cdot J','(nW/m^3)'};  
end
if 1 % length scales
  hca = irf_panel('length scales');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{re?,Le?},''comp'');',ic)
  hca.YLabel.String = {'length scales','(km)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'r_e','L_e'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [0 12];
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
hca = irf_panel('e PA e32 deflux all low E'); hca.YLim = [0 180];


%hca = irf_panel('e PA e32 deflux all low E mms4'); hca.YLim = [0 180];
%for  ii = 5:7; h(ii).CLim = [5 8.5]; end
%for  ii = 5:7; h(ii).YLim = [10 1000]; end
irf_plot_axis_align

h(1).Title.String = {['MMS' num2str(ic)],''};

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legshift = 0; % the two sc configuration plots

hca = irf_panel('B'); hca.YLim = [-15 15];

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 13;
end

%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];
add_length_on_top(h(1),70,0.2)