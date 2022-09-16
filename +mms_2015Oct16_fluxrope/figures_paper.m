%% Overview
%% Figure 1A: Stripped version: LMN
iTint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z');
tint_zoom = irf.tint('2015-10-16T10:33:42.00Z/2015-10-16T10:33:52.00Z'); 
boundaryTint = EpochTT(['2015-10-16T10:33:26.20Z';...
                        '2015-10-16T10:33:27.00Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);
boundaryTint = EpochTT(['2015-10-16T10:33:27.30Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);                
npanels = 9;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end
if 1 % n
  iisub = iisub+1;
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?},''comp'');',ic)
  hca.YLabel.String = {'n_i','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % eDist omni 64
  iisub = iisub+1;
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
  hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
  hhleg.FontSize = 9;
end

hca = irf_panel('delete for space');

if 1 % B
  %iisub = iisub+1;
  hca = irf_panel('B zoom');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'}',[1.02 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 400];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 1 % eDist omni 64
  %iisub = iisub+1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap)   

end
if 0 % Ve all n
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-900 400];  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
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
if 0 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  labelsOutside = 1;
  labelFontSize = 10;
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_plot(hca,{mvaAvE.z.resample(mvaOhmVexB.time),1*mvaOhmVexB.z,mvaOhmGradPe.z,mvaOhmVexB.z+mvaOhmGradPe.z.resample(mvaOhmVexB.time)+mvaAvE.z.resample(mvaOhmVexB.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};  
  
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','v_{e}xB'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('z')); hl = irf_legend(hca,{'\nabla \cdot P_e/ne'},[1.01 0.6]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('1')); hl = irf_legend(hca,{'E+v_{e}xB +\nabla\cdot P_e/ne'},[1.01 0.3]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else      
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB'},[1.01 0.79],'color',mms_colors('y')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.59],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('1')); hl(1).VerticalAlignment = 'top';
    end
  else
    irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','v_{e}xB+\nabla \cdot P_e/ne+E'},[0.98 0.1],'fontsize',12);
  end
  
  %irf_legend(hca,{'4 sc average'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc average',ic),[0.06 0.9],'color','k','fontsize',11);    
end
if 0 % sqrt(Q)
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{sqrt(Q1),sqrt(Q2),sqrt(Q3),sqrt(Q4)},'comp');      
  hca.YLabel.Interpreter = 'tex';
  ylabel(hca,{'$$\sqrt{Q}$$',''},'interpreter','latex');
  %ylabel(hca,{'$$Q^{1/2}$$',''},'interpreter','latex');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.06 0.9],'fontsize',11);
end

if 1 % pressure
  hca = irf_panel('P');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{PB?.tlim(tint),gsePi?.trace.tlim(tint)/3,gsePe?.trace.tlim(tint)/3,PB?.tlim(tint).resample(gsePi?)+gsePi?.trace.tlim(tint)/3+gsePe?.trace.tlim(tint).resample(gsePi?)/3},''comp'');',ic)
  hca.YLabel.String = {'P','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'P_B','P_i','P_e','P_B+P_i+P_e'}',[1.02 0.9],'fontsize',12);
  %hca.YScale = 'lin';
  %hca.YLim = [-1200 1500];
end

%load('caa/cmap.mat');
for ii = [4 5]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  h(ii).YLabel.FontSize = 12;
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',iTint)
irf_zoom(h(iisub+2:npanels),'x',tint_zoom)
%irf_zoom(h(iisub+2:npanels),'y')
irf_zoom(h(5:7),'y')
%irf_zoom(h([1:3 6 10]),'y')


irf_pl_mark(h(1:iisub),tint_zoom.epochUnix')
%irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix')

hca = irf_panel('e DEF omni 64');  hca.YGrid = 'off'; hca.XGrid = 'off'; hca.YLabel.String = {'E_e','(eV)'}
hca.YLim = [10 1000]; hca.CLim = [5.8 8.9];
colormap(hca,irf_colormap('waterfall'))

hca = irf_panel('Vi'); hca.YLim = [-300 100];
hca = irf_panel('B'); hca.YLim = [-35 50];
hca = irf_panel('n'); hca.YLim = [0 29.9];
hca = irf_panel('J mom'); hca.YLim = [-999 699];
hca = irf_panel('P'); hca.YLabel.FontSize = 12;

if 0
hca = irf_panel('Te'); irf_zoom(hca,'y');
hca = irf_panel('Ve'); hca.YLim = [-800 700];
hca = irf_panel('Q'); hca.YLim = [0 0.09]; hca.YTick = [0 0.05];[0.04 0.08];
hca = irf_panel('Ohms law electrons: N'); hca.YLim = [-6 6];
hca = irf_panel('e DEF omni 64');  hca.YGrid = 'off'; hca.XGrid = 'off'; hca.YLabel.String = {'E_e','(eV)'};
hca.XLabel.String = '';
end

ihmark = irf_pl_mark(h(1:iisub),tint_zoom.epochUnix',[0.8 0.8 0.8]); 

hmark = irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark), hmark(ii).LineStyle = '-'; end

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

hca = irf_panel('delete for space'); 
hca.Visible = 'off';

if 0
% Add labels for the different regions
if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
if exist('hleg_outflow','var'); delete(hleg_outflow); end
if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
hca = h(iisub+2);
set(hca,'ColorOrder',mms_colors('11'))
hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.83 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';
end
%% Figure 1B: Stripped version: LMN
iTint = irf.tint('2015-10-16T10:32:55.00Z/2015-10-16T10:34:10.00Z');
tint_zoom = irf.tint('2015-10-16T10:33:42.00Z/2015-10-16T10:33:52.00Z'); 
boundaryTint = EpochTT(['2015-10-16T10:33:26.20Z';...
                        '2015-10-16T10:33:27.00Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);
boundaryTint = EpochTT(['2015-10-16T10:33:27.30Z';...
                        '2015-10-16T10:33:30.35Z';...
                        '2015-10-16T10:33:31.00Z']);                
npanels = 9;
cmap = 'jet';
h = irf_plot(npanels);
ic = 1;
iisub = 0;
if 1 % B
  iisub = iisub+1;
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end
if 1 % Vi
  iisub = iisub+1;
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint),gseVi?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVi?.x,mvaVi?.y,mvaVi?.z},''comp'');',ic)
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'V_i_x','V_i_y','V_i_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
  hca.YLim = [-300 220];
end
if 0 % iPDist omni 64
  iisub = iisub+1;
  hca = irf_panel('i DEF omni');  
  c_eval('irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
  colormap(hca,cmap)
end
if 0 % eDist omni 64
  iisub = iisub+1;
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
  hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
  hhleg.FontSize = 9;
end

hca = irf_panel('delete for space');

if 1 % B
  %iisub = iisub+1;
  hca = irf_panel('B zoom');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z,mvaB?.abs},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N','|B|'}',[1.02 0.9],'fontsize',12);
end
if 1 % J moments 
  hca = irf_panel('J mom');
  set(hca,'ColorOrder',mms_colors('xyz'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?.x.tlim(tint),mvaJ?.y.tlim(tint),mvaJ?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % J moments perp par
  hca = irf_panel('J mom perp par');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJ?perp.x.tlim(tint),mvaJ?perp.y.tlim(tint),mvaJ?perp.z.tlim(tint),mvaJ?par.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'J','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'\perp L','\perp M','\perp N','||'}',[1.02 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 1 % JxB moments 
  hca = irf_panel('JxB mom');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJxB?.x.tlim(tint),mvaJxB?.y.tlim(tint),mvaJxB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'JxB','(nA/m^2 nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
  hca.YScale = 'lin';
  hca.YLim = [-1200 1500];
end
if 0 % JxB moments perp par
  hca = irf_panel('JxB perp par mom');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaJxB?perp.x.tlim(tint),mvaJxB?perp.y.tlim(tint),mvaJxB?perp.z.tlim(tint),mvaJxB?par.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'JxB','(nA/m^2 nT)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'\perp L','\perp M','\perp N','||'}',[1.02 0.9],'fontsize',12);
  hca.YScale = 'lin';
  %hca.YLim = [-1200 1500];
end
if 1 % curvature
  hca = irf_panel('magnetic curvature');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  irf_plot(hca,{mvaCurvB.x.tlim(tint),mvaCurvB.y.tlim(tint),mvaCurvB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'\rho_B','(...)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
  hca.YScale = 'lin';
  %hca.YLim = [-1200 1500];
end

if 1 % E
  %iisub = iisub+1;
  hca = irf_panel('E zoom');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'L','M','N'}',[1.02 0.9],'fontsize',12);
end
if 0 % pressure
  hca = irf_panel('P');
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %c_eval('irf_plot(hca,{gseJcurl.x.tlim(tint),gseJcurl.y.tlim(tint),gseJcurl.z.tlim(tint),gseJcurl.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{PB?.tlim(tint),gsePi?.trace.tlim(tint)/3,gsePe?.trace.tlim(tint)/3,PB?.tlim(tint).resample(gsePi?)+gsePi?.trace.tlim(tint)/3+gsePe?.trace.tlim(tint).resample(gsePi?)/3},''comp'');',ic)
  hca.YLabel.String = {'P','(nPa/m^2)'};
  set(hca,'ColorOrder',mms_colors('xyz1'))
  %irf_legend(hca,{'J_x','J_y','J_z'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'P_B','P_i','P_e','P_B+P_i+P_e'}',[1.02 0.9],'fontsize',12);
  %hca.YScale = 'lin';
  %hca.YLim = [-1200 1500];
end
if 0 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('irf_plot(hca,{facTe?.xx.tlim(tint),(facTe?.yy+facTe?.zz)/2},''comp'');',ic)
  hca.YLabel.String = {'T_e','(eV)'};
  set(hca,'ColorOrder',mms_colors('12'))
  irf_legend(hca,{'T_{||}','T_{\perp}'}',[1.02 0.9],'fontsize',12);
  %hca.YScale = 'log'; %hca.YTick = [10:10:100 200:100:1000];
  hca.YLim = [10 50];
  %hca.YTick
  irf_zoom(hca,'y')
end
if 0 % eDist omni 64
  %iisub = iisub+1;
  hca = irf_panel('e DEF omni 64');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.omni.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0
  hold(hca,'on')
  c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
  lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
  hold(hca,'off')
  hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    hhleg.FontSize = 9;
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  %colormap(hout,irf_colormap('space')) 
  colormap(hca,cmap)   

end
if 0 % Ve all n
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint),gseVe?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'V_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-900 400];  
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'4 sc average'},[0.06 0.9],'fontsize',11,'color','k');
  set(hca,'ColorOrder',mms_colors('xyza')) 
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
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
if 0 % Ohms law electrons: N
  hca = irf_panel('Ohms law electrons: N');
  labelsOutside = 1;
  labelFontSize = 10;
  set(hca,'ColorOrder',mms_colors('xyz1'))
  irf_plot(hca,{mvaAvE.z.resample(mvaOhmVexB.time),1*mvaOhmVexB.z,mvaOhmGradPe.z,mvaOhmVexB.z+mvaOhmGradPe.z.resample(mvaOhmVexB.time)+mvaAvE.z.resample(mvaOhmVexB.time)},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};  
  
  if labelsOutside
    if 1
      set(hca,'ColorOrder',mms_colors('xy')); hl = irf_legend(hca,{'E','v_{e}xB'},[1.01 0.95]); hl(1).VerticalAlignment = 'top'; hl(2).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize; hl(2).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('z')); hl = irf_legend(hca,{'\nabla \cdot P_e/ne'},[1.01 0.6]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
      set(hca,'ColorOrder',mms_colors('1')); hl = irf_legend(hca,{'E+v_{e}xB +\nabla\cdot P_e/ne'},[1.01 0.3]); hl(1).VerticalAlignment = 'top';  hl(1).FontSize = labelFontSize;
    else      
      hl = irf_legend(hca,{'E'},[1.01 0.99],'color',mms_colors('x')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB'},[1.01 0.79],'color',mms_colors('y')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-\nabla \cdot P_e/ne'},[1.01 0.59],'color',mms_colors('z')); hl(1).VerticalAlignment = 'top';
      hl = irf_legend(hca,{'-v_{e}xB-\nabla \cdot P_e/ne'},[1.01 0.39],'color',mms_colors('1')); hl(1).VerticalAlignment = 'top';
    end
  else
    irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','v_{e}xB+\nabla \cdot P_e/ne+E'},[0.98 0.1],'fontsize',12);
  end
  
  %irf_legend(hca,{'4 sc average'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,irf_ssub('4 sc average',ic),[0.06 0.9],'color','k','fontsize',11);    
end
if 0 % sqrt(Q)
  hca = irf_panel('Q');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{sqrt(Q1),sqrt(Q2),sqrt(Q3),sqrt(Q4)},'comp');      
  hca.YLabel.Interpreter = 'tex';
  ylabel(hca,{'$$\sqrt{Q}$$',''},'interpreter','latex');
  %ylabel(hca,{'$$Q^{1/2}$$',''},'interpreter','latex');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.06 0.9],'fontsize',11);
end

if 0 % n  
  hca = irf_panel('n');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ni?},''comp'');',ic)
  hca.YLabel.String = {'n_i','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('1'))  
  hca.YLim = [0 30];
end

%load('caa/cmap.mat');
for ii = [4 5]
  %h(ii).CLim = [4 8];  
  %colormap(h(ii),irf_colormap('space'))
  colormap(h(ii),'jet')
end
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = [1:iisub iisub+2:npanels]  
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0])
  h(ii).YLabel.FontSize = 12;
  nInd = nInd + 1;
end

irf_zoom(h(1:iisub),'x',iTint)
irf_zoom(h(iisub+2:npanels),'x',tint_zoom)
%irf_zoom(h(iisub+2:npanels),'y')
irf_zoom(h(5:9),'y')
%irf_zoom(h([1:3 6 10]),'y')


irf_pl_mark(h(1:iisub),tint_zoom.epochUnix')
%irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix')

%hca = irf_panel('e DEF omni 64');  hca.YGrid = 'off'; hca.XGrid = 'off'; hca.YLabel.String = {'E_e','(eV)'}
%hca.YLim = [10 1000]; hca.CLim = [5.8 8.9];
%colormap(hca,irf_colormap('waterfall'))

hca = irf_panel('Vi'); hca.YLim = [-300 100];
hca = irf_panel('B'); hca.YLim = [-35 50];
hca = irf_panel('B zoom'); hca.YLim = [-29 35];
%hca = irf_panel('n'); hca.YLim = [0 29.9];
%hca = irf_panel('n'); hca.YLim = [0 29.9];
hca = irf_panel('J mom'); hca.YLim = [-999 699];
%hca = irf_panel('P'); hca.YLabel.FontSize = 12;

if 0
hca = irf_panel('Te'); irf_zoom(hca,'y');
hca = irf_panel('Ve'); hca.YLim = [-800 700];
hca = irf_panel('Q'); hca.YLim = [0 0.09]; hca.YTick = [0 0.05];[0.04 0.08];
hca = irf_panel('Ohms law electrons: N'); hca.YLim = [-6 6];
hca = irf_panel('e DEF omni 64');  hca.YGrid = 'off'; hca.XGrid = 'off'; hca.YLabel.String = {'E_e','(eV)'};
hca.XLabel.String = '';
end

ihmark = irf_pl_mark(h(1:iisub),tint_zoom.epochUnix',[0.8 0.8 0.8]); 

hmark = irf_pl_mark(h(iisub+2:npanels),boundaryTint.epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark), hmark(ii).LineStyle = '-'; end

irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);

irf_plot_zoomin_lines_between_panels(h(iisub),h(iisub+2))

hca = irf_panel('delete for space'); 
hca.Visible = 'off';

if 0
% Add labels for the different regions
if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
if exist('hleg_outflow','var'); delete(hleg_outflow); end
if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
hca = h(iisub+2);
set(hca,'ColorOrder',mms_colors('11'))
hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.83 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';
end

%% Electron pitchangle

% Eperp = mu*B
tint_ad = irf.tint('2015-10-16T10:33:44.40Z/2015-10-16T10:33:46.00Z');
fEperp = @(mu,B) mu*B;
fEpar = @(mu,B,E) E - mu*B;

tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosheath-magnetosphere
npanels = 4;
cmap = pic_colors('pasteljet');

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 1 % B
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
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
if 0 % B curv
  hca = irf_panel('Curv B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_plot(hca,{mvaCurvB.x.tlim(tint),mvaCurvB.y.tlim(tint),mvaCurvB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'Curv B','(1/km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-0.05 0.02];
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
if 0 % eDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.deflux.omni.tlim(tint).specrec,''log'');',ic)  
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
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF par');
  pas = [0 30];
  %c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)  
  c_eval('ePitch?lim = ePitch?par;',ic)   
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 1 % mu*B
    hold(hca,'on')
    c_eval('B1 = gseB?.abs.resample(tint_ad(1)).data;',ic)
    c_eval('Epar = fEpar(100/B1,gseB?.abs,300);',ic)
    aa = irf_plot(hca,Epar.tlim(tint_ad),'k--');
    aa.Color = [0 0 0]; aa.LineWidth = 1.0;
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
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF perp');
  pas = [75 115];
  %c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)  
  c_eval('ePitch?lim = ePitch?perp;',ic)  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 1 % mu*B
    hold(hca,'on')
    c_eval('Eperp = fEperp(70/18,gseB?.abs);',ic)
    aa = irf_plot(hca,Eperp.tlim(tint_ad),'k--');
    aa.Color = [0 0 0]; aa.LineWidth = 1.0;
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
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF apar');
  pas = [150 180];
  %c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)  
  c_eval('ePitch?lim = ePitch?apar;',ic)  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 1 % mu*B
    hold(hca,'on')
    c_eval('B1 = gseB?.abs.resample(tint_ad(1)).data;',ic)
    c_eval('Epar = fEpar(400/B1,gseB?.abs,600);',ic)
    aa = irf_plot(hca,Epar.tlim(tint_ad),'k--');
    aa.Color = [0 0 0]; aa.LineWidth = 1.0;
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

%% Ion pitchangle

% Eperp = mu*B


tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosheath-magnetosphere
tint = irf.tint('2015-10-16T10:32:54.00Z/2015-10-16T10:34:14.00Z'); % magnetosheath-magnetosphere
npanels = 5;
cmap = pic_colors('pasteljet');

ediv = 5000;

% Initialize figure
ic = 1;
h = irf_plot(npanels);

% Resize figure
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.5; figurePostition(4)=figurePostition(4)*0.9;
hcf = gcf; hcf.Position = figurePostition;

% Plot
if 0 % B
  hca = irf_panel('B gse');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  hca.YLim = [-32 40];
end
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  %hca.YLim = [-32 40];
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
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyz'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_{i}','(nT)'};  
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-500 500];  
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
if 1 % iDist omni 64
  hca = irf_panel('e DEF omni 64');    
  c_eval('[hout,hcb] = irf_spectrogram(hca,iPDist?.deflux.omni.tlim(tint).specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  hca.YLabel.String = {'E_i','(eV)'}; 
  colormap(hca,cmap)
  if 1
    hold(hca,'on')
    irf_plot(hca,irf.ts_scalar(iPDist1.time,ones(iPDist1.length,1)*ediv),'-k')
    hold(hca,'off')
  end
  %hca.YLim = [10 1000];
  %hca.CLim = [5 8.3];
end
if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all low E');  
  elim = [10 ediv];
  c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)  
  hold(hca,'on')
  hold(hca,'off')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  hca.YLabel.String = {'\theta_{PA}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all high E');
  elim = [ediv 34000];
  c_eval('irf_spectrogram(hca,iPitch?.tlim(tint).elim(elim).deflux.specrec(''pa''),''log'');',ic)  
  hold(hca,'on')
  hold(hca,'off')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  hca.YLabel.String = {'\theta_{PA}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
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
  %c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)  
  c_eval('ePitch?lim = ePitch?par;',ic)   
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 1 % mu*B
    hold(hca,'on')
    c_eval('B1 = gseB?.abs.resample(tint_ad(1)).data;',ic)
    c_eval('Epar = fEpar(100/B1,gseB?.abs,300);',ic)
    aa = irf_plot(hca,Epar.tlim(tint_ad),'k--');
    aa.Color = [0 0 0]; aa.LineWidth = 1.0;
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
  %c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)  
  c_eval('ePitch?lim = ePitch?perp;',ic)  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 1 % mu*B
    hold(hca,'on')
    c_eval('Eperp = fEperp(70/18,gseB?.abs);',ic)
    aa = irf_plot(hca,Eperp.tlim(tint_ad),'k--');
    aa.Color = [0 0 0]; aa.LineWidth = 1.0;
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
  %c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)  
  c_eval('ePitch?lim = ePitch?apar;',ic)  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePitch?lim.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  if 0 % sc pot
    hold(hca,'on')
    c_eval('lineScpot = irf_plot(hca,scPot?,''k'');',ic)  
    lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
    %hhleg=irf_legend(hca,'V_{sc}',[0.95 0.21],'color',[0 0 0]);
    %hhleg.FontSize = 9;
  end
  if 1 % mu*B
    hold(hca,'on')
    c_eval('B1 = gseB?.abs.resample(tint_ad(1)).data;',ic)
    c_eval('Epar = fEpar(400/B1,gseB?.abs,600);',ic)
    aa = irf_plot(hca,Epar.tlim(tint_ad),'k--');
    aa.Color = [0 0 0]; aa.LineWidth = 1.0;
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


