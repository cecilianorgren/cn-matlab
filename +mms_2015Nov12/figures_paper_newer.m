%% Timing 
dt0 = [0 0 0 0];
dt1 = [0.00     0.05     -0.10     -0.08]; v1 = 70*[-0.86,0.37,0.35]*lmn'; % good for second peak and BL
dt2 = [0.00     0.05     -0.05     -0.03]; v2 = 120*[-0.82 0.49 0.30]*lmn';% first part
n1 = irf_norm(v1);
n2 = irf_norm(v2);
itiming = 1;
c_eval('dt = dt?; vv = v?; vvn = vv/norm(vv); n = n?;',itiming)

pshift = 0;

%% Figure 2: Overview, started 12-12-2017
fastTint = irf.tint('2015-11-12T06:10:00.00Z/2015-11-12T07:30:00.00Z');
brstTint = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
tintBCS = irf.tint('2015-11-12T07:19:19.80Z/2015-11-12T07:19:22.70Z');
timeMP = irf_time('2015-11-12T07:19:17.80Z','utc>epochtt');
npanels = 11;
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
  irf_legend(hca,{'x','y','z','|B|'},[0.98 0.9],'fontsize',14);
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
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',14);  
  hca.YLim = [-490 900];  
end
if 1 % Vi
  hca = irf_panel('Vi');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVi?.x.tlim(tint),gseVi?.y.tlim(tint),gseVi?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',14);  
  %hca.YLim = [-1100 1100];  
end
if 1 % Ve  
  hca = irf_panel('Ve');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',14);  
  hca.YLim = [-910 780];  
  hca.YTick = [-1200:400:1200];
end
if 1 % Te par perp
  hca = irf_panel('Te');
  set(hca,'ColorOrder',mms_colors('1342'))
  refTi = 10;
  c_eval('irf_plot(hca,{facTe?.trace/3,facTe?.xx.tlim(brstTint),(facTe?.yy+facTe?.zz)/2,facTi?.trace/3/refTi},''comp'');',ic)
  hca.YLabel.String = {'T','(eV)'};
  set(hca,'ColorOrder',mms_colors('1342'))
  irf_legend(hca,{'T_e','T_{e,||}','T_{e,\perp}',['T_i/' num2str(refTi,'%.0f')]},[0.8 0.9],'fontsize',14);
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
if 1 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('123'))
  c_eval('irf_plot(hca,{betae?,betai?,beta?},''comp'');',ic)
  hca.YLabel.String = {'\beta'};
  hca.YScale = 'log';
  hca.YLim = [0.09 200];
  hca.YMinorTick = 'on';
  set(hca,'ColorOrder',mms_colors('123'))
  irf_legend(hca,{'\beta_e','\beta_i','\beta'},[0.98 0.9],'fontsize',14);  
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
if 0 % J x B
  hca = irf_panel('E + Ve x B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{gseJxB?.x.tlim(tint),gseJxB?.y.tlim(tint),gseJxB?.z.tlim(tint)},''comp'');',ic)
  %c_eval('irf_plot(hca,{gseVe?.x.tlim(tint),gseVe?.y.tlim(tint),gseVe?.z.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'J x B','()'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'x','y','z'},[0.98 0.9],'fontsize',12);  
  %hca.YLim = [-1100 1100];  
end
if 0 % diff P
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
  irf_legend(h(ii),legends{nInd},[0.01 0.9],'color',[0 0 0],'fontsize',14)
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

%% Figure 3: 4 sc, MVA, BL, BM, BN, curvBm, JL, JM, JN
tintZoom = irf.tint('2015-11-12T07:19:20.20Z/2015-11-12T07:19:22.20Z');
%tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.50Z');
npanels = 8;
pshift = 4;
h = irf_plot(npanels+pshift);


scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.2; figurePostition(4)=figurePostition(4)*0.5;
hcf = gcf; hcf.Position = figurePostition;

c_eval('hca = irf_panel(''delete1?''); hca.Visible = ''off'';',1:pshift)

dt = dt1; vv = v1; vvn = vv/norm(vv);

if 1 % BL
  hca = irf_panel('BL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.x.tlim(tint),mvaB2.x.tlim(tint),mvaB3.x.tlim(tint),mvaB4.x.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{L}','(nT)'};
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
end
if 1 % BM
  hca = irf_panel('BM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{M}','(nT)'};
end
if 1 % BM time shifted
  hca = irf_panel('BM dt');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.y.tlim(tint),mvaB2.y.tlim(tint),mvaB3.y.tlim(tint),mvaB4.y.tlim(tint)},'comp','dt',dt);
  hca.YLabel.String = {'B_{M}','(nT)'};
  set(hca,'ColorOrder',mms_colors('12341'))
  irf_legend(hca,{['\Delta t = [ ' num2str(dt(1))],num2str(dt(2)),num2str(dt(3)),num2str(dt(4)),'] s'},[0.95 0.9],'fontsize',12);
  %irf_legend(hca,{['v = ' num2str(norm(vv),'%.0f') '\times[' num2str(vvn(1),'%.2f') ', ' num2str(vvn(2),'%.2f') ', ' num2str(vvn(3),'%.2f') '] km/s']},[0.98 0.7],'fontsize',12);
  
end
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 1 % curv B
  hca = irf_panel('curv B');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaCurvB.x.tlim(tint),mvaCurvB.y.tlim(tint),mvaCurvB.z.tlim(tint)},'comp');
  hca.YLabel.String = {'\rho_c','(1/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 1 % JL
  %%
  hca = irf_panel('JL');
  if 1
    set(hca,'ColorOrder',mms_colors('b1234'))
    lines = irf_plot(hca,{mvaJcurl.x,mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x},'comp');    
    lines.Children(5).LineWidth = 2;
    hca.YLabel.String = {'J_L','(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234')) 
    irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
    set(hca,'ColorOrder',mms_colors('b')) 
    irf_legend(hca,{'curlometer'},[0.05 0.9],'fontsize',12);
  end
  %lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x},'comp');  
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.98 0.9],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % JM
  hca = irf_panel('JM');
  if 1
    set(hca,'ColorOrder',mms_colors('b1234'))
    lines = irf_plot(hca,{mvaJcurl.y,mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');    
    lines.Children(5).LineWidth = 2;
    hca.YLabel.String = {'J_M','(nA/m^2)'};
%     set(hca,'ColorOrder',mms_colors('1234')) 
%     irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
%     set(hca,'ColorOrder',mms_colors('b')) 
%     irf_legend(hca,{'Curlometer'},[0.05 0.9],'fontsize',12);
  end
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.98 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % JN
  hca = irf_panel('JN');
  if 1
    set(hca,'ColorOrder',mms_colors('b1234'))
    lines = irf_plot(hca,{mvaJcurl.z,mvaJ1.z,mvaJ2.z,mvaJ3.z,mvaJ4.z},'comp');    
    lines.Children(5).LineWidth = 2;
    hca.YLabel.String = {'J_N','(nA/m^2)'};
%     set(hca,'ColorOrder',mms_colors('1234')) 
%     irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);
%     set(hca,'ColorOrder',mms_colors('b')) 
%     irf_legend(hca,{'Curlometer'},[0.05 0.9],'fontsize',12);
  else
    set(hca,'ColorOrder',mms_colors('1234b'))
    lines = irf_plot(hca,{mvaJ1.z,mvaJ2.z,mvaJ3.z,mvaJ4.z,mvaJcurl.z},'comp');   
    %lines = irf_plot(hca,{mvaJ1.z,mvaJ2.z,mvaJ3.z,mvaJ4.z},'comp');   
    hca.YLabel.String = {'J_N','(nA/m^2)'};
    set(hca,'ColorOrder',mms_colors('1234b'))
    %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
    %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.98 0.9],'fontsize',12);
    %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
    %lines.Children(6).LineWidth = 1.5;
    %lines.Children(6).LineStyle = '--';
    %hca.YLim = [-1200 2000];
  end
end
if 0 % Ve  L
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,L}','(km/s)'},'interpreter','tex');
end
if 0 % Ve  M
  hca = irf_panel('Ve M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).y,mvaVe2.tlim(tint).y,mvaVe3.tlim(tint).y,mvaVe4.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,M}','(km/s)'},'interpreter','tex');
end
if 0 % Ve  N
  hca = irf_panel('Ve N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).z,mvaVe2.tlim(tint).z,mvaVe3.tlim(tint).z,mvaVe4.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{e,N}','(km/s)'},'interpreter','tex');
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');      
  ylabel(hca,{'v_{e,||}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp L
  hca = irf_panel('Ve perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).x,mvaVe2perp.tlim(tint).x,mvaVe3perp.tlim(tint).x,mvaVe4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,L}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp M
  hca = irf_panel('Ve perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1perp.tlim(tint).y,mvaVe2perp.tlim(tint).y,mvaVe3perp.tlim(tint).y,mvaVe4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{e,\perp,M}','(km/s)'},'interpreter','tex');
end
if 0 % Ve perp N
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

% 2 timings: Plot sc positions
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
%hca.XDir = 'reverse';
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
%hca.XDir = 'reverse';
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
%hca.XDir = 'reverse';
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
%irf_legend(h2(1),'<- MSH',[0.01 0.99],'fontsize',11)
%irf_legend(h2(1),'MSP ->',[0.99 0.99],'fontsize',11)

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

%% Figure 4: 1 sc more detailed
npanels = 8;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');
%tintZoom = irf.tint('2015-11-12T07:19:10.00Z/2015-11-12T07:19:40.00Z');
h = irf_plot(npanels);
ic = 1;
pshift = 0;
scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.3; figurePostition(4)=figurePostition(4)*0.7;
hcf = gcf; hcf.Position = figurePostition;

cmap = 'jet';

tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
c_eval('intNe? = irf_integrate(ne?-ne?.resample(tref).data,tref);',ic)
tref = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt');
c_eval('mvaIntE?N = irf_integrate(mvaE?.z,tref);',ic)

% prepare energy levels for adiabaticity
% reference time for mu in adiabatic inflow region
tref4 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength4 = [0 -0.30];
tref3 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength3 = [0 -0.35];
tref2 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength2 = [-0.50 0];
tref1 = irf_time('2015-11-12T07:19:22.00Z','utc>epochtt'); tlength1 = [0 -0.40];
c_eval('tref = tref?;',ic)
c_eval('tintB = tref+tlength?;',ic)
tintB = tintZoom;    
c_eval('Wrefperp = 0.5*(facTe?.yy.resample(tref).data+facTe?.zz.resample(tref).data);',ic)
c_eval('Wrefpar = facTe?.resample(tref).zz.data;',ic)
c_eval('Wreftot = 0.5*(Wrefpar+Wrefperp);',ic)
Wrefperp = 100;
Wrefpar = 100;
Wreftot = Wrefperp+Wrefpar;
c_eval('Bref = mvaB?.abs.resample(tref).data;',ic)
c_eval('muref = Wrefperp/mvaB?.abs.resample(tref).data;',ic)
c_eval('Wperp = irf.ts_scalar(mvaB?.tlim(tintB).time,mvaB?.tlim(tintB).abs.data*Wrefperp/Bref);',ic)
c_eval('Wpar = irf.ts_scalar(Wperp.time,Wreftot-Wperp.data);',ic)
Efactorpar = 1;
Efactorperp = 1;

isub = 0;
zoomE = [];
zoomPA = [];
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
if 0 % n
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
if 1 % ePDist pa low energies
  hca = irf_panel('e PA e32 deflux all low E'); isub = isub + 1;  zoomPA = [zoomPA isub];
  elim = [10 500];
  c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom+[-2 2]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
    
  if 1 % magnetic mirror angle
    hold(hca,'on')
    tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
    tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.30];
    tref2 = irf_time('2015-11-12T07:19:21.45Z','utc>epochtt'); tlength2 = [-0.40 0];
    tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
    tref1 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength1 = [-0.35 0];

    c_eval('tref = tref?;',ic)
    c_eval('tintB = tref+tlength?;',ic)
    
    c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
    
    c_eval('tref = tref?;',ic)    
    thetaref = -00;
    c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));',ic)    
    
    set(hca,'ColorOrder',mms_colors('11'))  
    alpha_shift = 40;
    irf_plot(hca,{alphaB-alpha_shift,180-alphaB-alpha_shift},'comp','-.'); 
    href = irf_plot(hca,{irf.ts_scalar(tref,90-alpha_shift)},'comp');
    href.Children(1).Marker = '*';
    href.Children(1).MarkerSize = 10;
    
    alpha_shift = 00;
    irf_plot(hca,{alphaB-alpha_shift,180-alphaB-alpha_shift},'comp','--');      
    %irf_pl_mark(hca,tref,'k')
    href = irf_plot(hca,{irf.ts_scalar(tref,90-alpha_shift)},'comp');
    href.Children(1).Marker = '*';
    href.Children(1).MarkerSize = 10;
    hca.XGrid = 'off';
    hca.YGrid = 'off';
    hold(hca,'off')
  end
  %hca.CLim = [7.4 8.0];
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  colormap(hca,cmap)
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 0 % ePDist pa low energies
  hca = irf_panel('e PA e32 dpflux all low E');
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength4 = [0 0.30];
  tref3 = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt'); tlength3 = [0 0.35];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.50 0];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 1000];
  if 0
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint).dpflux.elim(elim).specrec(''pa''),''log'');',ic)
  else
    c_eval('irf_spectrogram(hca,ePDist?.e64.tlim(tintZoom).pitchangles(dmpaB?,12).specrec(''pa''),''log'');',ic)  
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
if 1 % eDist omni 64 apar
  hca = irf_panel('e DEF 3'); isub = isub + 1; zoomE = [zoomE isub];
  pas = [150 180];
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom+[-2 2]).pitchangles(dmpaB?,pas);',ic)      
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
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    set(hca,'ColorOrder',mms_colors('1'))
    c_eval('lineWperp = irf_plot(hca,{Wpar*Efactorpar},''comp'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    %irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    %irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 perp
  hca = irf_panel('e DEF 2'); isub = isub + 1; zoomE = [zoomE isub];
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom+[-2 2]).pitchangles(dmpaB?,pas);',ic)    
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
  if 1 % plot perpendicular energy baed on magnetic moment conservation    
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('1'))
    c_eval('lineWperp = irf_plot(hca,{Wperp},''comp'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,cmap)   
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
end
if 1 % eDist omni 64 par
  hca = irf_panel('e DEF 1'); isub = isub + 1; zoomE = [zoomE isub];
  pas = [0 30]; 
  c_eval('ePitch?lim = ePDist?.tlim(tintZoom+[-2 2]).pitchangles(dmpaB?,pas);',ic)     
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
  if 1 % plot parallel energy based on magnetic moment conservation    
    hold(hca,'on')        
    set(hca,'ColorOrder',mms_colors('1'))
    c_eval('lineWperp = irf_plot(hca,{Wpar*Efactorpar},''comp'');',ic)  
    %lineScpot.Color = [0 0 0]; lineScpot.LineWidth = 1.5;
    %irf_legend(hca,{sprintf('%.0f E_{e,ref}',Efactorpar)},[0.77 0.2])
    hold(hca,'off')
    %irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';
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
if 0 % ePDist pa low energies
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
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x*1e3,mvaGradPe.y*1e3,mvaGradPe.z*1e3},'comp');
  hca.YLabel.String = {'\nabla \cdot P_e','(pPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);  
  irf_legend(hca,{'4 spacecraft'},[0.05 0.9],'fontsize',12,'color','k');
end
if 1 % E
  hca = irf_panel('E'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?.x.tlim(tint),mvaE?.y.tlim(tint),mvaE?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'E','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
  irf_zoom(hca,'y')
end
if 0 % Epar
  hca = irf_panel('E par'); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaE?par.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'E_{||}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'L_{\perp}','M_{\perp}','N_{\perp}','E_{||}'},[0.98 0.9],'fontsize',12);  
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
if 0 % Ohm's law N
  hca = irf_panel('Ohm N');
  set(hca,'ColorOrder',mms_colors('2134b'))
  irf_plot(hca,{mvaAvE.z,mvaOhmGradPe.z,-mvaOhmVexB.z,-mvaOhmVixB.z,mvaOhmJxB.z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('2134b'))
  irf_legend(hca,{'E ','\nabla P_e ','-v_e\times B ','-v_i\times B ','J\times B/ne '},[0.98 0.9],'fontsize',12);
  irf_legend(hca,{'4sc average'},[0.05 0.9],'color',[0 0 0],'fontsize',12);
end

irf_zoom(h,'x',tintZoom)
irf_zoom(h,'y')
irf_plot_axis_align

for  ii = zoomE; h(ii).CLim = [7 8.7]; end
for  ii = zoomE; h(ii).YLim = [11 1000]; end
for  ii = zoomPA; h(ii).CLim = [7.4 8.1]; end
for  ii = zoomPA; h(ii).YLim = [0 180]; end

%hca = irf_panel('B'); hca.YLim = [-14 14];
%hca = irf_panel('e PA e32 deflux all low E'); hca.CLim = [7.2 8.2];
%hca = irf_panel('e PA e32 deflux all mid E'); hca.CLim = [7.2 8.2];

h(1).Title.String = ['MMS' num2str(ic)];

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
legends_color = {'k','k','k','k','w','w','w','k','k','k'};
legshift = 0; % the two sc configuration plots

for ii = 1:npanels
  irf_legend(h(ii+pshift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0],'fontsize',14,'color',legends_color{ii+legshift});
  h(ii+pshift).FontSize = 14;  
  h(ii+pshift).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

%hca = irf_panel('gradPe'); hca.YLim = [-2.5 2.5];
%hca = irf_panel('E'); hca.YLim = [-3 3];
%irf_plot_axis_align
%for ii = 1:3; h(ii).Visible = 'off'; end

%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Figure 5: Crescents MMS3, crescents are 2x3

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

% plot crescents
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

% plot crescents
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

%% Figure 5: Crescents MMS3, crescents are 3x2

% Initialize figure
h1 = irf_plot(10);
cmap = 'jet';
ic = 3;
tintZoom = irf.tint('2015-11-12T07:19:19.40Z/2015-11-12T07:19:23.10Z');

% Plot
if 1 % B
  hca = irf_panel('B');
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{gseB?.x.tlim(tint),gseB?.y.tlim(tint),gseB?.z.tlim(tint),gseB?.abs.tlim(tint)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaB?.x.tlim(tint),mvaB?.y.tlim(tint),mvaB?.z.tlim(tint),mvaB?.abs.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N','|B|'},[0.98 0.9],'fontsize',12);
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
  elim = [10 1000];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tintZoom+[-1 1]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    c_eval('irf_spectrogram(hca,ePDist?.tlim(tintZoom+[-1 1]).pitchangles(dmpaB?,20).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  end
  %c_eval('irf_spectrogram(hca,ePDist?.e64.pitchangles(dmpaB?,20).elim([180 203]).deflux.specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};   
  hca.YTick = [45 90 135];   
  %hca.CLim = [7.4 8.1];%[7.4 8.1]
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
%for  ii = 6:8; h1(ii).CLim = [5 8.5]; end
%for  ii = 6:8; h1(ii).YLim = [10 1000]; end
irf_plot_axis_align
h1(1).Title.String = irf_ssub('MMS ?',ic);

h1 = h1(1:3);

ipanel = 0;
pshift = 0; legshift = 0;
for ii = 1:3
  ipanel = ipanel + 1;
  irf_legend(h1(ipanel),legends{ipanel},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h1(ipanel).FontSize = 14;  
  h1(ipanel).YLabel.FontSize = 14;
  %h(ii+pshift).Position(3) = h(ii+pshift).Position(3)*1.1;
end

% plot crescents
if exist('hmark'); delete(hmark); end
c_eval('times = ePDist?.time;',ic)
tind = 902:1:(960+9*2);
tind_mms1 = [903 905 907 909 911 913];
tind_mms3 = [900 901 902 907 908 909]; % sorted into rows
%tind_mms3 = [900 907 901 908 902 909]; % sorted into columns
legends = {'a)','b)','c)','d)','g)','e)','h)','f)','i)'};
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'};
c_eval('tind = tind_mms?;',ic)

% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  
nrows = 2;
ncols = 3;
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


if 0 % move around distribution axes, no adapated for 2x3 placement
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
end

%% Figure 6, TSeries Bmod, Emod 4 electrons, 3 planes E = 0
close
% Observed fields
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z')+-0.01;
tintObs = tintObs;
CS_normal_velocity = 70; % km/s
ic = 1;
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsCurvB = mvaCurvB.resample(obsB).tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsE = mvaE?_new.tlim(tintObs); obsE = obsE.resample(obsB);'... %'obsE = mvaE?_new.tlim(tintObs) - mvaEht.resample(mvaE?_new).tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp.resample(obsB);'...
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObsE = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;

% Colors
colors = mms_colors('xyz');

% Model fields
mms_2015Nov12.Bmodel;

% Initialize particles
colors = mms_colors('1234'); % particle colors
electron_energy = [170 40 170 40]; % eV
vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

% Initial positions and velocitites
x0 = [0 0 0 0]; 
y0 = [120 120 -40 -50]; % km
z0 = [30 30 -30 -30];

velocity_angle = [-20 -10 -10 -20];
velocity_angle_L = [-10 -20 20 10];
vx0 = vt.*cosd(velocity_angle_L); % km/s
vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);

x_init_all = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s

nParticles = numel(electron_energy);

% Integration
T = 1; % max integration time
limN = 32e3; % where integration stops

saveParticle = cell(1,nParticles);
for iParticle= 1:nParticles 
  % Initial positions and velocities                                   
  x_init = x_init_all(:,iParticle); 
  
  % Integrate trajectory
  stopfunction = @(t,y) eom.lim(t,y,limN);
  options = odeset('Events',stopfunction,'RelTol',1e-6);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
  EoM = @(ttt,xxx) eom.general(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
  [t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options    
  
  if 1 % interpolate to even time steps
    tlim = 1e-8;
    dtcut = 100*(t(end)-t(end-1));
    [x_sol, t] = resample(x_sol, t, 2/max(diff(t))); 
    x_sol = x_sol(t<(t(end)-dtcut),:); t = t(t<(t(end)-dtcut));
    x_sol = x_sol(t>(t(1)+dtcut),:); t = t(t>(t(1)+dtcut));
  end
  x = x_sol(:,1);
  y = x_sol(:,2);
  z = x_sol(:,3);
  vx = x_sol(:,4);
  vy = x_sol(:,5);
  vz = x_sol(:,6); 
  
  Bxyz = [Bx(x,y,z),By(x,y,z),Bz(x,y,z)]; normBxyz = irf_norm(Bxyz);
  Vxyz = [vx vy vz]; normVxyz = irf_norm(Vxyz);  
  pitchangle = acosd(normBxyz(:,1).*normVxyz(:,1) + ...
                     normBxyz(:,2).*normVxyz(:,2) + ...
                     normBxyz(:,3).*normVxyz(:,3));
      
  saveParticle{iParticle}.t = t;
  saveParticle{iParticle}.r = x_sol(:,1:3);
  saveParticle{iParticle}.r0 = [x0(iParticle),y0(iParticle),z0(iParticle)];
  saveParticle{iParticle}.v = x_sol(:,4:6);
  saveParticle{iParticle}.v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)];
  saveParticle{iParticle}.B = Bxyz;
  saveParticle{iParticle}.pa = pitchangle;
end

saveParticle_4 = saveParticle;
%% Figure 6, option 1, 4 panels, incl B, 3 particle panels
scrsz = get(groot,'ScreenSize');
figurePosition = scrsz;
%figurePosition(3) = 0.4*scrsz(3);
figurePosition(3) = 0.7*scrsz(3);
%figurePosition(4) = 0.5*scrsz(4);
figure('position',figurePosition)

colors = mms_colors('1234');
Bcolor = mms_colors('matlab'); Bcolor = Bcolor(3,:);

% Set up plot
nRows = 6;
nCols = 2;
units = irf_units;


clear h;

% left column
h(1) = subplot(nRows,nCols,1); 
h(2) = subplot(nRows,nCols,3); 
h(3) = subplot(nRows,nCols,5); 
h(4) = subplot(nRows,nCols,[7 9 11]); 
% right column
h(5) = subplot(nRows,nCols,[2 4 6]); 
h(6) = subplot(nRows,nCols,[8 10 12]); 
  

if 1 % rearrange positions
  h(1).Position(4) = h(1).Position(4)*1.2;  
  h(2).Position(4) = h(2).Position(4)*1.2;
  h(3).Position(4) = h(3).Position(4)*1.2;
  h(4).Position(4) = h(4).Position(4)*1.2;
  h(6).Position(4) = h(6).Position(4)*1.2;

  h(1).Position(2) = h(1).Position(2)+0.1;
  h(2).Position(2) = h(2).Position(2)+0.1;
  h(3).Position(2) = h(3).Position(2)+0.1;     

  h(1).Position(2) = h(1).Position(2)-h(1).Position(4);  
  h(2).Position(2) = h(1).Position(2)-h(2).Position(4);
  h(3).Position(2) = h(2).Position(2)-h(3).Position(4);  
  h(4).Position(2) = h(3).Position(2)-h(4).Position(4);  

  h(6).Position(2) = h(4).Position(2);
  h(5).Position(2) = h(3).Position(2);
  h(5).Position(4) = h(3).Position(4)*3;

  h(5).Position(3) = 0.4;h(5).Position(3)*1.2;
  h(6).Position(3) = 0.4;h(6).Position(3)*1.2;

end

isub = 1;
if 1 % BL, BM, BN, Babs 
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,50)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 0 % Curvature of magnetic field, obs+mod
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,obsCurvB.data,'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  modfac = 0.1;
  lineMod = plot(hca,zObs,modCurvB(:,1)*modfac,'--','color',B_colors(1,:));
  plot(hca,zObs,modCurvB(:,2)*modfac,'--','color',B_colors(2,:))
  plot(hca,zObs,modCurvB(:,3)*modfac,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'curb B (1/km)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30]; 
end
if 1 % EL, EM, EN, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data],'-');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.95],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'E','(mV/m)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3.9 3.9];
end  
if 0 % Time: ePDist
  hca = h(isub); isub = isub + 1;
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
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'x',tintObs)
end
if 1 % Distance: ePitch.deflux, test particle pitchangles to be added later
  hca = h(isub); isub = isub + 1;
  elim = [10 1000];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  %pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))  
  shading(hca,'flat')
  
  for iP = 1:numel(saveParticle) % Electron pitchangles
    % Pick out the data
    z = saveParticle{iP}.r(:,3); % N
    pa = saveParticle{iP}.pa; % pitch angles    
    colors = mms_colors('1234');
    step = 2;
    hold(hca,'on');
    plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
    hold(hca,'off') 
  end
  
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePitch.deflux, Lioville mapped
  hca = h(isub); isub = isub + 1;
  elim = [10 1000];
  
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end

maxx = 0;
for iP = 1:numel(saveParticle) % Electron trajectories
  isub = 4;
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  maxx = max([maxx, max(x)]);
  hca = h(isub);
  
  if 1
    for is = isub + [0 1 2];
      hca = h(is);
      hl(iP) =  plot3(hca,x*1e-3,y*1e-3,z*1e-3,'color',colors(iP,:));
      if is == 5; hleg(iP) = hl(iP); end
      if iP == 1; hold(hca,'on'); end
      plot3(hca,x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go',...
                x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'rx') % plot start and end in km's
      xlabel(hca,'L (km)')
      ylabel(hca,'M (km)')
      zlabel(hca,'N (km)');
      %hca.ZLim = [-30 30];
      %hca.YLim = [30 ceil(maxx*1e-3/10)*10];

      if iP == numel(saveParticle); hold(hca,'off'); end  
    end
  else 
    % NL
    hca = h(isub); isub = isub + 1;
    hl(iP) =  plot(hca,z*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,z(1)*1e-3,x(1)*1e-3,'go',...
             z(end)*1e-3,x(end)*1e-3,'rx') % plot start and end in km's
    ylabel(hca,'L (km)')
    xlabel(hca,'N (km)'); hca.XLim = [-30 30];
    hca.XLim = h(2).XLim;
    hca.YLim = [30 ceil(maxx*1e-3/10)*10];

    if iP == numel(saveParticle); hold(hca,'off'); end  

    % MN
    hca = h(isub); isub = isub + 1;
    hl(iP) =  plot(hca,-y*1e-3,-z*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,-y(1)*1e-3,-z(1)*1e-3,'go',...
             -y(end)*1e-3,-z(end)*1e-3,'rx') % plot start and end in km's
    ylabel(hca,'L (km)')
    xlabel(hca,'N (km)'); %hca.XLim = [-30 30];
    %hca.XLim = h(2).XLim;
    %hca.YLim = [30 ceil(maxx*1e-3/10)*10];

    if iP == numel(saveParticle); hold(hca,'off'); end  

    % ML
    hca = h(isub); isub = isub + 1;
               plot(hca,-y*1e-3,x*1e-3,'color',colors(iP,:));
    %hl(iP) =  plot(hca,-y*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,-y(1)*1e-3,x(1)*1e-3,'go',...
             -y(end)*1e-3,x(end)*1e-3,'rx') % plot start and end in km's
    ylabel(hca,'L (km)')
    xlabel(hca,'-M (km)');
    %hca.XLim = h(2).XLim;
    %hca.YLim = [30 ceil(maxx*1e-3/10)*10];    
  end
  
  
  if iP == numel(saveParticle); hold(hca,'off'); end    
  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
%legend(hl,labels{:},'location','west')
legend(hleg,labels{:},'location','northwest')

% add magnetic field
for iP = 4:6
  zBstart = -30*1e3;
  zBstop = 30*1e3;
  dr = 100;
  xyz = [0 0 zBstart];
  while xyz(end,3)<zBstop
    dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
  end
  if 1 % 3    
    nq = 20;
    qscale  = 0.2;
    
    toplot = fix(linspace(1,size(xyz,1),nq));
    toplot = [1 fix(size(xyz,1)*0.85)];
    hca = h(iP);
    hold(hca,'on')
    dL = 0;
    dM = 75; % km
    dN = 0;
    hB = plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
    quiver3(hca,...
       xyz(toplot,1)*1e-3 - dN,...
       xyz(toplot,2)*1e-3 - dM,...
       xyz(toplot,3)*1e-3 - dL,...
       Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       qscale,'color',Bcolor)
        
    dL = 0;
    dM = 10; % km
    dN = 0;
    plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
    quiver3(hca,...
       xyz(toplot,1)*1e-3 - dN,...
       xyz(toplot,2)*1e-3 - dM,...
       xyz(toplot,3)*1e-3 - dL,...
       Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
       qscale,'color',Bcolor)
     hold(hca,'off')
  else % 2D
    hca = h(iP);
    hold(hca,'on')
    yellow = mms_colors('matlab'); yellow = yellow(3,:);
    qscale = 0.5;
    plot(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3,'color',yellow);
    hold(hca,'on')
    nq = 30;
    toplot = fix(linspace(1,size(xyz,1),nq));
    quiver(hca,xyz(toplot,3)*1e-3,xyz(toplot,1)*1e-3,...
       Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale,'color',yellow)    
    dL = 20; % km
    plot(hca,xyz(:,3)*1e-3,xyz(:,1)*1e-3+dL,'color',yellow); 
    quiver(hca,xyz(toplot,3)*1e-3,xyz(toplot,1)*1e-3 + dL,...
       Bz(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,...
       Bx(xyz(toplot,1),-xyz(toplot,2),xyz(toplot,3))*1e9,qscale,'color',yellow)     
  end
end

h(4).ZLim = [-30 30];
c_eval('h(?).XTick = [-20 0 20];',1:2)
c_eval('h(?).XTick = [];',3)
c_eval('h(?).XLim = [-30 30];',1:3)
%c_eval('h(?).XLim = [0 115];',4:6)
%c_eval('h(?).YLim = [-50 110];',4:6)
%c_eval('h(?).ZLim = [-32 32];',4:6)
c_eval('h(?).XTick = -200:20:200;',4:6)
c_eval('h(?).YTick = -200:20:200;',4:6)
c_eval('h(?).ZTick = -200:20:200;',4:6)

irf_plot_axis_align(h(1:3))

c_eval('h(?).XGrid = ''on'';',1:6)
c_eval('h(?).YGrid = ''on'';',4:6)
c_eval('h(?).ZGrid = ''on'';',4:6)
c_eval('h(?).FontSize = 14;',1:6)
c_eval('h(?).Box = ''on'';',4:6)
c_eval('h(?).Position(3) = h(2).Position(3);',1:6);
h(5).Position(3) = 0.4;
h(6).Position(3) = 0.4;

h(5).YTickLabel = [];
%h(3).Position(2) = h(2).Position(2)-h(2).Position(4)-dy;
%h(4).Position(2) = h(3).Position(2)-h(3).Position(4)-dy;
%h(5).Position(2) = h(4).Position(2)-h(4).Position(4)-dy;

view(h(4),[0 -1 0]); h(4).ZDir = 'reverse'; camroll(h(4),90); h(4).XAxisLocation = 'top';
view(h(5),[1 0 0]); 
view(h(6),[0 0 -1]); camroll(h(6),90); h(6).XAxisLocation = 'top';

%% Figure 6A: tseries 
% load '/Users/cno062/Research/Events/2015-11-12_071854/saveParticle_finiteE_nP=32768'
%tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
t_center = irf_time('2015-11-12T07:19:21.175000012Z','utc>EpochTT'); 
tintObs = t_center + 1*[-1 1];
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'modB = modB.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPDistMap = tsFmap.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObsB = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zModB = (modB.time.epochUnix-mean(modB.time.epochUnix))*CS_normal_velocity;
zObsE = (obsE.time.epochUnix-mean(obsE.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;
zObsPDistMap = (obsPDistMap.time.epochUnix-mean(obsPDistMap.time.epochUnix))*CS_normal_velocity;

fig = figure(106);
fig.Position = [100   100   475   542];
% Set up plot
npanels = 7;
units = irf_units;
h = irf_plot(npanels);  

nlim = [-28 28];
elim = [30 350];
  
isub = 1;
if 1 % BL, BM, BN, Babs 
  hca = h(isub); isub = isub + 1; 
  B_colors = mms_colors('xyz1'); % Colors
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObsB,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N','|B|'},[0.95 0.2],'fontsize',14)
  
  hold(hca,'on')
  if 0
    zMod = linspace(nlim(1)-1,nlim(2)+1,100)*1e3;
    %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
    lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
    plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
    plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
    plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  else 
    linesMod = plot(hca,zModB,[modB.data modB.abs.data],'--');
    linesMod(1).Color = B_colors(1,:); linesMod(1).LineStyle = '--';
    linesMod(2).Color = B_colors(2,:); linesMod(2).LineStyle = '--';
    linesMod(3).Color = B_colors(3,:); linesMod(3).LineStyle = '--';
    linesMod(4).Color = B_colors(4,:); linesMod(4).LineStyle = '--';
  end
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 0 % Curvature of magnetic field, obs+mod
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,obsCurvB.data,'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  modfac = 0.1;
  lineMod = plot(hca,zObs,modCurvB(:,1)*modfac,'--','color',B_colors(1,:));
  plot(hca,zObs,modCurvB(:,2)*modfac,'--','color',B_colors(2,:))
  plot(hca,zObs,modCurvB(:,3)*modfac,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'curb B (1/km)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30]; 
end
if 1 % EL, EM, EN, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  
  linesObs = plot(hca,zObsE,[obsE.data],'-');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.95 0.95],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(nlim(1)-1,nlim(2)+1,100)*1e3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'E','(mV/m)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3.9 3.9];
end  
if 0 % Time: ePDist
  hca = h(isub); isub = isub + 1;
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
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'x',tintObs)
end
if 0 % Distance: ePitch.deflux, test particle pitchangles to be added later
  hca = h(isub); isub = isub + 1;
  elim = [10 1000];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  %pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))  
  shading(hca,'flat')
  
  for iP = 1:numel(saveParticle) % Electron pitchangles
    % Pick out the data
    z = saveParticle{iP}.r(:,3); % N
    pa = saveParticle{iP}.pa; % pitch angles    
    colors = mms_colors('1234');
    step = 2;
    hold(hca,'on');
    plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
    hold(hca,'off') 
  end
  
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 0 % Distance: ePitch.deflux, Liouville mapped, from left
  hca = h(isub); isub = isub + 1;
  plotPitch = tsFmap_left.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa')
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  shading(hca,'flat')     
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 0 % Distance: ePitch.deflux, Liouville mapped, from right
  hca = h(isub); isub = isub + 1;
  plotPitch = tsFmap_right.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa')
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  shading(hca,'flat')     
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 0 % Distance: ePitch.deflux, Liouville mapped, all
  hca = h(isub); isub = isub + 1;
  plotPitch = tsFmap_all.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa')
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  shading(hca,'flat')     
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePitch/psd, NO test particles
  hca = h(isub); isub = isub + 1;
  plotPitch = obsPDist.convertto('s^3/km^6').einterp.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  %pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))  
  shading(hca,'flat')
  
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePitch/psd, test particle pitchangles to be added later
  hca = h(isub); isub = isub + 1;
  plotPitch = obsPDist.convertto('s^3/km^6').einterp.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDist-0.5*diff(zObsPDist(1:2)); zObsPDist(end)+0.5*diff(zObsPDist(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  %pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))  
  shading(hca,'flat')
  
  for iP = 1:numel(saveParticle_4) % Electron pitchangles
    % Pick out the data
    z = saveParticle_4{iP}.r(:,3); % N
    pa = saveParticle_4{iP}.pa; % pitch angles    
    colors = mms_colors('1234');
    step = 2;
    hold(hca,'on');
    plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
    hold(hca,'off') 
  end
  
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePitch/psd, Liouville mapped, from left
  hca = h(isub); isub = isub + 1;
  plotPitch = tsFmap_left.convertto('s^3/km^6').pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDistMap-0.5*diff(zObsPDistMap(1:2)); zObsPDistMap(end)+0.5*diff(zObsPDistMap(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  shading(hca,'flat')     
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca); hcb_save = hcb;
  hcb.YLabel.String = plotPitch.p_label;
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePitch/psd, Liouville mapped, from right
  hca = h(isub); isub = isub + 1;
  plotPitch = tsFmap_right.convertto('s^3/km^6').pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDistMap-0.5*diff(zObsPDistMap(1:2)); zObsPDistMap(end)+0.5*diff(zObsPDistMap(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  shading(hca,'flat')     
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePitch/psd, Liouville mapped, all
  hca = h(isub); isub = isub + 1;
  plotPitch = tsFmap_all.convertto('s^3/km^6').pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa');
  edges_pa = [plotPitch.f-0.5*diff(plotPitch.f(1:2)) plotPitch.f(end)+0.5*diff(plotPitch.f(1:2))];
  edges_z = [zObsPDistMap-0.5*diff(zObsPDistMap(1:2)); zObsPDistMap(end)+0.5*diff(zObsPDistMap(1:2))]';
  surf_def = zeros(numel(edges_pa),numel(edges_z));
  surf_col = log10(plotPitch.p');
  surf(hca,edges_z,edges_pa,surf_def,surf_col)
  view(hca,[0 0 1])
  hca.YLim = [0 180];
  shading(hca,'flat')     
  hca.XGrid = 'off';
  hca.YGrid = 'off';   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = plotPitch.p_label;
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
if 0 % Distance: ePitch.deflux, Lioville mapped
  hca = h(isub); isub = isub + 1;
  elim = [10 1000];   
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log'); 
  hca.YLabel.String = {'Pitchangle','(\circ)'};
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',14,'color',[0 0 0]);  
  xlabel(hca,'N (km)')  
end
irf_plot_axis_align
h(3).Box = 'on';
c_eval('h(?).XLim = nlim;',1:npanels)
c_eval('h(?).XTickLabel = [];',1:npanels-1)
c_eval('h(?).XLabel.String = [];',1:npanels-1)
c_eval('h(?).FontSize = 14;',1:npanels)
c_eval('h(?).CLim = [3.0 3.9];',3:7)
colormap_str = 'jet';
c_eval('colormap(h(?),colormap_str);',3:7)

c_eval('h(?).YLabel.String = [];',[3 4 6 7])
c_eval('delete(colorbar(h(?)));',[3 4 6 7])
hcb_save.Position(2) = 0.10;% = [0.8042 0.3432    0.0421    0.1208];
hcb_save.Position(4) = 0.60;
hcb_save.Position(3) = 0.03;
%hcb_save.Position = [0.8042 0.3432    0.0421    0.1208];
labels = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'};
legends_color = {'k','k','k','k','k','w','w','k','k','k'};
c_eval('irf_legend(h(?),labels{?},[0.02 0.95],''fontsize'',16,''color'',legends_color{?})',1:npanels)

%% Figure 6B: trajectories, turned 3D plots
colors = mms_colors('1234');
Bcolor = mms_colors('matlab'); Bcolor = Bcolor(3,:);

nRows = 2;
nCols = 2;
h(1) = subplot(nRows,nCols,2); 
h(2) = subplot(nRows,nCols,3); 
h(3) = subplot(nRows,nCols,4); 


for iP = 1:numel(saveParticle) % Electron trajectories
  isub = 1; 
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  maxx = max([maxx, max(x)]);
  hca = h(isub);
  
  
  for is = isub + [0 1 2];
    hca = h(is);
    hl(iP) =  plot3(hca,x*1e-3,y*1e-3,z*1e-3,'color',colors(iP,:));
    if is == 5; hleg(iP) = hl(iP); end
    if iP == 1; hold(hca,'on'); end
    plot3(hca,x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go',...
              x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'rx') % plot start and end in km's
    xlabel(hca,'L (km)')
    ylabel(hca,'M (km)')
    zlabel(hca,'N (km)');
    %hca.ZLim = [-30 30];
    %hca.YLim = [30 ceil(maxx*1e-3/10)*10];

    if iP == numel(saveParticle); hold(hca,'off'); end  
  end  
  if iP == numel(saveParticle); hold(hca,'off'); end    
  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
%legend(hl,labels{:},'location','west')
legend(h(1),labels{:},'location','west')

% add magnetic field
for iP = 1:3
  zBstart = -30*1e3;
  zBstop = 30*1e3;
  dr = 100;
  xyz = [0 0 zBstart];
  while xyz(end,3)<zBstop
    dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
  end
  
  nq = 20;
  qscale  = 0.2;

  toplot = fix(linspace(1,size(xyz,1),nq));
  toplot = [1 fix(size(xyz,1)*0.85)];
  hca = h(iP);
  hold(hca,'on')
  dL = 0;
  dM = 75; % km
  dN = 0;
  hB = plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
  quiver3(hca,...
     xyz(toplot,1)*1e-3 - dN,...
     xyz(toplot,2)*1e-3 - dM,...
     xyz(toplot,3)*1e-3 - dL,...
     Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     qscale,'color',Bcolor)

  dL = 0;
  dM = 10; % km
  dN = 0;
  plot3(hca,xyz(:,1)*1e-3 - dL,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
  quiver3(hca,...
     xyz(toplot,1)*1e-3 - dN,...
     xyz(toplot,2)*1e-3 - dM,...
     xyz(toplot,3)*1e-3 - dL,...
     Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
     qscale,'color',Bcolor)
   hold(hca,'off')
end

% adjust axes
c_eval('axis(h(?),''equal'')',1:3)
view(h(2),[0 -1 0]); h(2).ZDir = 'reverse'; camroll(h(2),90); h(2).XAxisLocation = 'top';
view(h(1),[1 0 0]); 
view(h(3),[0 0 -1]); camroll(h(3),90); h(3).XAxisLocation = 'top';

%h(1).ZLim = [-30 30];
%c_eval('h(?).XLim = [0 115];',4:6)
%c_eval('h(?).YLim = [-50 110];',4:6)
%c_eval('h(?).ZLim = [-32 32];',4:6)
c_eval('h(?).XTick = -200:20:200;',1:3)
c_eval('h(?).YTick = -200:20:200;',1:3)
c_eval('h(?).ZTick = -200:20:200;',1:3)

c_eval('h(?).XGrid = ''on'';',1:3)
c_eval('h(?).YGrid = ''on'';',1:3)
c_eval('h(?).ZGrid = ''on'';',1:3)
c_eval('h(?).FontSize = 14;',1:3)
c_eval('h(?).Box = ''on'';',1:3)

l_width = h(3).Position(4);
m_width = h(3).Position(3);
n_width = h(1).Position(4);


% size
h(1).Position(3) = m_width;
h(1).Position(4) = n_width;
h(2).Position(3) = n_width;
h(2).Position(4) = l_width;
h(3).Position(3) = m_width;
h(3).Position(4) = l_width;
% position
h(1).Position(3) = m_width;
h(1).Position(4) = n_width;
h(2).Position(3) = n_width;
h(2).Position(4) = l_width;
h(3).Position(3) = m_width;
h(3).Position(4) = l_width;


%h(2).Position(3) = 0.4;
%h(6).Position(3) = 0.4;

h(1).YTickLabel = [];
h(1).YLabel.String = [];
h(3).XTickLabel = [];
h(3).XLabel.String = [];
%% Figure 6B: trajectories, 2D plots
close
fig = figure(10);
fig.Position = [100 100 1100 1100];
colors = mms_colors('1234');
Bcolor = mms_colors('matlab'); Bcolor = Bcolor(3,:).^.5;

if 0 % run this once to get the dimensions of the panels (m_width etc.), then run the 'else' option to properly position them
  nRows = 2;
  nCols = 2;
  h(1) = subplot(nRows,nCols,2); 
  h(2) = subplot(nRows,nCols,3); 
  h(3) = subplot(nRows,nCols,4); 
else
  edge_left = 0.8;
  edge_bot = 0.2;
  h(3) = axes('Position',[edge_left-m_width edge_bot m_width l_width]);
  h(2) = axes('Position',[edge_left-m_width-n_width edge_bot n_width l_width]);
  h(1) = axes('Position',[edge_left-m_width edge_bot+l_width m_width n_width]);  
end

lim_n = [-32 32];
lim_l = [-60 180];
lim_m = [-60 220];
maxx = 0;

addB1 = 1;
addB2 = 1;

if 1 % add magnetic field
  % get magnetic field line
  zBstart = -30*1e3;
  zBstop = 30*1e3;
  dr = 100;
  xyz = [0 0 zBstart];
  while xyz(end,3)<zBstop
    dx = dr*Bx(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dy = dr*By(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    dz = dr*Bz(xyz(end,1),xyz(end,2),xyz(end,3))/B(xyz(end,1),xyz(end,2),xyz(end,3));
    xyz(end+1,:) = xyz(end,:) + [dx dy dz];      
  end  

  nq = 20;
  qscale  = 0.2;

  toplot = fix(linspace(1,size(xyz,1),nq));
  toplot = [1 fix(size(xyz,1)*0.85)];
  quiver_width = 2;

  dL1 = -7;
  dL2 = -40;
  dM1 = 40;
  dM2 = 70;
  if 1 % MN
    hca = h(1);
    hold(hca,'on')
    if addB1
      dL = dL1;
      dM = dM1; % km
      dN = 0;  
      hB = plot(hca,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
      hq = quiver(hca,...
         xyz(toplot,2)*1e-3 - dM,...
         xyz(toplot,3)*1e-3 - dN,...          
         By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         qscale,'color',Bcolor);
      hq.LineWidth = quiver_width;
    end
    if addB2
      dL = dL1;
      dM = dM2; % km
      dN = 0;
      hB = plot(hca,xyz(:,2)*1e-3 - dM,xyz(:,3)*1e-3 - dN,'color',Bcolor,'LineWidth',2);      
      hq = quiver(hca,...
         xyz(toplot,2)*1e-3 - dM,...
         xyz(toplot,3)*1e-3 - dN,...          
         By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         qscale,'color',Bcolor);
      hq.LineWidth = quiver_width; 
      end
    hold(hca,'off')  
  end
  if 1 % NL
    hca = h(2);
    hold(hca,'on')
    if addB1
      dL = dL1;
      dM = dM1; % km
      dN = 0;  
      hB = plot(hca,xyz(:,3)*1e-3 - dN,xyz(:,1)*1e-3 - dL,'color',Bcolor,'LineWidth',2);      
      hq = quiver(hca,...
         xyz(toplot,3)*1e-3 - dN,...
         xyz(toplot,1)*1e-3 - dL,...          
         Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         qscale,'color',Bcolor);
       hq.LineWidth = quiver_width;
    end
    if addB2
      dL = dL2;
      dM = dM2; % km
      dN = 0;
      hB = plot(hca,xyz(:,3)*1e-3 - dN,xyz(:,1)*1e-3 - dL,'color',Bcolor,'LineWidth',2);      
      hq = quiver(hca,...
         xyz(toplot,3)*1e-3 - dN,...
         xyz(toplot,1)*1e-3 - dL,...          
         Bz(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         qscale,'color',Bcolor);
      hq.LineWidth = quiver_width;
    end
    hold(hca,'off')  
  end
  if 1 % ML
    hca = h(3);
    hold(hca,'on')
    if addB1
      dL = dL1;
      dM = dM1; % km
      dN = 0;  
      hB = plot(hca,xyz(:,2)*1e-3 - dM,xyz(:,1)*1e-3 - dL,'color',Bcolor,'LineWidth',2);      
      hq = quiver(hca,...
         xyz(toplot,2)*1e-3 - dM,...
         xyz(toplot,1)*1e-3 - dL,...          
         By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         qscale,'color',Bcolor);
       hq.LineWidth = quiver_width;
    end
    if addB2
      dL = dL2;
      dM = dM2; % km
      dN = 0;
      hB = plot(hca,xyz(:,2)*1e-3 - dM,xyz(:,1)*1e-3 - dL,'color',Bcolor,'LineWidth',2);      
      hq = quiver(hca,...
         xyz(toplot,2)*1e-3 - dM,...
         xyz(toplot,1)*1e-3 - dL,...          
         By(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         Bx(xyz(toplot,1),xyz(toplot,2),xyz(toplot,3))*1e9,...
         qscale,'color',Bcolor);
       hq.LineWidth = quiver_width;    
    end
     hold(hca,'off') 
  end
end

for iP = 1:numel(saveParticle) % Electron trajectories  
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  maxx = max([maxx, max(x)]);  
    
  % MN
  hca = h(1);  
  hold(hca,'on');
  hl(iP) =  plot(hca,y*1e-3,z*1e-3,'color',colors(iP,:));
  %if iP == 1; hold(hca,'on'); end  
  plot(hca,y(1)*1e-3,z(1)*1e-3,'go',...
           y(end)*1e-3,z(end)*1e-3,'rx') % plot start and end in km's  
  xlabel(hca,'M (km)'); 
  ylabel(hca,'N (km)')  
  
  if iP == numel(saveParticle); hold(hca,'off'); end  

  % NL
  hca = h(2);
  hold(hca,'on');
  plot(hca,z*1e-3,x*1e-3,'color',colors(iP,:));
  %if iP == 1; hold(hca,'on'); end
  plot(hca,z(1)*1e-3,x(1)*1e-3,'go',...
           z(end)*1e-3,x(end)*1e-3,'rx') % plot start and end in km's
  xlabel(hca,'N (km)');
  ylabel(hca,'L (km)')  
  

  if iP == numel(saveParticle); hold(hca,'off'); end  

  % ML
  hca = h(3);
  hold(hca,'on');
  plot(hca,y*1e-3,x*1e-3,'color',colors(iP,:));
  %hl(iP) =  plot(hca,-y*1e-3,x*1e-3,'color',colors(iP,:));
  %if iP == 1; hold(hca,'on'); end
  plot(hca,y(1)*1e-3,x(1)*1e-3,'go',...
           y(end)*1e-3,x(end)*1e-3,'rx') % plot start and end in km's
  ylabel(hca,'L (km)')
  xlabel(hca,'M (km)');
  %hca.XLim = h(2).XLim;
  %hca.YLim = [30 ceil(maxx*1e-3/10)*10];   
  
  %if iP == numel(saveParticle); hold(hca,'off'); end    
  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
%legend(hl,labels{:},'location','west')
hleg = legend(hl,labels{:},'location','west');
hleg.Position = [0.2228    0.5761    0.0927    0.0595];



% adjust axes
c_eval('axis(h(?),''equal'')',1:3)

h(1).XLim = lim_m;
h(1).YLim = lim_n;  
h(2).XLim = lim_n;
h(2).YLim = lim_l;
h(3).XLim = lim_m;
h(3).YLim = lim_l;
  
%h(1).ZLim = [-30 30maxx];
%c_eval('h(?).XLim = [0 115];',4:6)
%c_eval('h(?).YLim = [-50 110];',4:6)
%c_eval('h(?).ZLim = [-32 32];',4:6)
c_eval('h(?).XTick = -300:20:300;',1:3)
c_eval('h(?).YTick = -300:20:300;',1:3)
c_eval('h(?).ZTick = -300:20:300;',1:3)

c_eval('h(?).XGrid = ''on'';',1:3)
c_eval('h(?).YGrid = ''on'';',1:3)
c_eval('h(?).ZGrid = ''on'';',1:3)
c_eval('h(?).FontSize = 14;',1:3)
c_eval('h(?).Box = ''on'';',1:3)

l_width = h(3).Position(4);
m_width = l_width/diff(lim_l)*diff(lim_m);
n_width = l_width/diff(lim_l)*diff(lim_n);

%m_width = h(3).Position(3);
%n_width = h(2).Position(3);
set(h(1:3),'fontsize',18)

h(1).XTickLabel = [];
h(1).XLabel.String = [];
h(3).YTickLabel = [];
h(3).YLabel.String = [];

labels = {'h)','i)','j)'};
%labels = {'a)','b)','c)'};
c_eval('irf_legend(h(?),labels{?},[0.05 0.95],''fontsize'',18)',1:3)




