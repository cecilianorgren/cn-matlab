%% Figure 2: Plot figure with fields etc of 4 sc
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
npanels = 13;
h = irf_plot(npanels);

eint = [0 40000];

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
if 1 % BN
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
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaJcurl.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  irf_legend(hca,{'J_{curl}'},[0.98 0.9],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end
if 0 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');  
 
  % Add bars to indicate width of flow reversal
  hold(hca,'on')
  tintMMS1 = irf.tint('2015-10-16T10:33:30.25Z/2015-10-16T10:33:30.43Z');  
  tintMMS2 = irf.tint('2015-10-16T10:33:30.16Z/2015-10-16T10:33:30.37Z');
  tintMMS3 = irf.tint('2015-10-16T10:33:29.52Z/2015-10-16T10:33:30.25Z');
  tintMMS4 = irf.tint('2015-10-16T10:33:29.62Z/2015-10-16T10:33:30.40Z');
  irf_plot(hca,TSeries(tintMMS1,700*[1;1]),'linewidth',2,'color',mms_colors('1'))
  irf_plot(hca,TSeries(tintMMS2,650*[1;1]),'linewidth',2,'color',mms_colors('2'))
  irf_plot(hca,TSeries(tintMMS3,600*[1;1]),'linewidth',2,'color',mms_colors('3'))
  irf_plot(hca,TSeries(tintMMS4,550*[1;1]),'linewidth',2,'color',mms_colors('4'))
  hold(hca,'off')
  
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_e_L','(km/s)'},'interpreter','tex');
end
if 1 % Vi par
  hca = irf_panel('Vi par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVi1par.tlim(tint),mvaVi2par.tlim(tint),mvaVi3par.tlim(tint),mvaVi4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_{i,||}','(km/s)'},'interpreter','tex');
end
if 1 % Vi perp L
  hca = irf_panel('Vi perp L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVi1perp.tlim(tint).x,mvaVi2perp.tlim(tint).x,mvaVi3perp.tlim(tint).x,mvaVi4perp.tlim(tint).x},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'v_{i,\perp,L}','(km/s)'},'interpreter','tex');
end
if 1 % Vi perp M
  hca = irf_panel('Vi perp M');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVi1perp.tlim(tint).y,mvaVi2perp.tlim(tint).y,mvaVi3perp.tlim(tint).y,mvaVi4perp.tlim(tint).y},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,M}','(km/s)'};
  ylabel(hca,{'v_{i,\perp,M}','(km/s)'},'interpreter','tex');
end
if 1 % Vi perp N
  hca = irf_panel('Vi perp N');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVi1perp.tlim(tint).z,mvaVi2perp.tlim(tint).z,mvaVi3perp.tlim(tint).z,mvaVi4perp.tlim(tint).z},'comp');    
  %hca.YLabel.String = {'v_{e,\perp,N}','(km/s)'};
  ylabel(hca,{'v_{i,\perp,N}','(km/s)'},'interpreter','tex');
end
if 1 % ePDist omni 32 mms 1
  ic  = 1;
  hca = irf_panel('i DEF omni');
  c_eval('irf_spectrogram(hca,iPDist?.omni(''i'').specrec,''log'');',ic)
  hca.YLabel.String = {'E_i','(eV)'};  
  hca.YScale = 'log';
  hca.YTick = [1e1 1e2 1e3 1e4];
end
if 1 % ePDist pa 32 mms1
  ic = 1;
  hca = irf_panel(irf_ssub('i PA deflux, mms?',ic));
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,[]).deflux.elim(eint).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
  hca.YLabel.String = {'\theta_{PA,i}',irf_ssub('MMS ?',ic), '(\circ)'};     
end
if 1 % ePDist pa 32 mms2
  ic = 2;
  hca = irf_panel(irf_ssub('i PA deflux, mms?',ic));
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,[]).deflux.elim(eint).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'}; 
  hca.YTick = [45 90 135];   
  hca.YLabel.String = {'\theta_{PA,i}',irf_ssub('MMS ?',ic), '(\circ)'};   
end
if 1 % ePDist pa 32 mms3
  ic = 3; 
  hca = irf_panel(irf_ssub('i PA deflux, mms?',ic));
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,[]).deflux.elim(eint).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YLabel.String = {'\theta_{PA,i}','(\circ)'};   
  hca.YTick = [45 90 135];   
  hca.YLabel.String = {'\theta_{PA,i}',irf_ssub('MMS ?',ic), '(\circ)'};     
end
if 1 % ePDist pa 32 mms4
  ic = 4;
  hca = irf_panel(irf_ssub('i PA deflux, mms?',ic));
  c_eval('irf_spectrogram(hca,iPDist?.pitchangles(dmpaB?,[]).deflux.elim(eint).specrec(''pa''),''log'');',ic)
  %hca.YLabel.String = {'Pitchangle','(\circ)'};   
  hca.YTick = [45 90 135];   
  hca.YLabel.String = {'\theta_{PA,i}',irf_ssub('MMS ?',ic), '(\circ)'};   
end
if 0 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
%   irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
%                 mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
%                 mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
%                 mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
irf_plot(hca,{mvaEVexB1.z,mvaEVexB2.z,mvaEVexB3.z,mvaEVexB4.z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);  
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = RedJe?;',ic)
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP  
  
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h(:),'x',tintZoom)
irf_zoom(h(1:8),'y')
irf_plot_axis_align
%%
%h(9).YLim = [8 30000];
for ii = 9:npanels
  h(ii).CLim = [7.1 8.2];
  h(ii).YLim = [0 180];
  colormap(h(ii),'jet')
end
h(9).YLim = [8 3000];
h(9).CLim = [4 8.2];
%%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:npanels
  irf_legend(h(ii+3),legends{ii},[0.01 0.9],'color',[0 0 0])
  h(ii+3).FontSize = 12,
end

delete(h(1:3));
h = h(4:end);
h(7).YLim = 7*[-1 1];
h(8).YLim = 7*[-1 1];

% Plot sc positions
nrows = 5;
ncols = 3;
h2(1) = subplot(nrows,ncols,1);
h2(2) = subplot(nrows,ncols,2);

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
hca.XLim = 1.5*[-10 10];
hca.YLim = 1.5*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N';
hca.YLabel.String = 'L';

%plot(mvaRR1(3),mvaRR1(1))
hca = h2(2);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.2*[-10 10];
hca.YLim = 1.2*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N';
hca.YLabel.String = 'M';
hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside');
hleg.Box = 'off';
hleg.Position(1) = hleg.Position(1)+0.2;

% plot plane
v = [-0.90 0.26 0.36];
v = [-0.90 -0.28 -0.33];
lmnV = [-0.57 -0.03 -0.82];
lmnV = [-0.55 -0.08 -0.83];
lmnV = [-0.45 -0.05 -0.89];
lmnV = [-0.49 -0.05 -0.87];
v = [-0.9139   -0.4006    0.0653];
v = [-0.90 -0.28 -0.33];
v = [-0.88 -0.26 -0.40];
lmnV = v*[L' M' N'];
v = lmnV;
planeNormal = irf_norm(v);
planeVelocity = v;
planeSpeed = norm(v); % km/s
planeRadius = 10; % km

zFun = @(x,y) -(planeNormal(1)*x+planeNormal(2)*y)/planeNormal(3);
% (planeNormal(1)*x+planeNormal(2)*y)/planeNormal(3) = 0;
%  planeNormal(1)*x+planeNormal(2)*y = 0;
% x = - planeNormal(2)*y/planeNormal(1);
funM = @(y)  -planeNormal(3)*y/planeNormal(2);
funL = @(y)  -planeNormal(3)*y/planeNormal(1);

x = [-10 10];
hca = h2(1);
hold(hca,'on');
hplane = plot(hca,funL(x),x);
hold(hca,'off');

y = [-10 10];
hca = h2(2);
hold(hca,'on');
hplane = plot(hca,x,funM(x));
hold(hca,'off');

%% Figure 2: Different comparison plot of more field components
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
h = irf_plot(13);
%irf_panel('delete1');
%irf_panel('delete2');
%irf_panel('delete3');
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
if 1 % BN
  hca = irf_panel('BN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaB1.z.tlim(tint),mvaB2.z.tlim(tint),mvaB3.z.tlim(tint),mvaB4.z.tlim(tint)},'comp');
  hca.YLabel.String = {'B_{N}','(nT)'};
end
if 0 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{ne1.tlim(tint),ne2.tlim(tint),ne3.tlim(tint),ne4.tlim(tint)},'comp');
  hca.YLabel.String = {irf_ssub('n_{e}',ic),'(cm^{-3})'};
end
if 0 % VeL
  hca = irf_panel('Ve L');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1.tlim(tint).x,mvaVe2.tlim(tint).x,mvaVe3.tlim(tint).x,mvaVe4.tlim(tint).x},'comp');  
 
  % Add bars to indicate width of flow reversal
  hold(hca,'on')
  tintMMS1 = irf.tint('2015-10-16T10:33:30.25Z/2015-10-16T10:33:30.43Z');  
  tintMMS2 = irf.tint('2015-10-16T10:33:30.16Z/2015-10-16T10:33:30.37Z');
  tintMMS3 = irf.tint('2015-10-16T10:33:29.52Z/2015-10-16T10:33:30.25Z');
  tintMMS4 = irf.tint('2015-10-16T10:33:29.62Z/2015-10-16T10:33:30.40Z');
  irf_plot(hca,TSeries(tintMMS1,700*[1;1]),'linewidth',2,'color',mms_colors('1'))
  irf_plot(hca,TSeries(tintMMS2,650*[1;1]),'linewidth',2,'color',mms_colors('2'))
  irf_plot(hca,TSeries(tintMMS3,600*[1;1]),'linewidth',2,'color',mms_colors('3'))
  irf_plot(hca,TSeries(tintMMS4,550*[1;1]),'linewidth',2,'color',mms_colors('4'))
  hold(hca,'off')
  
  hca.YLabel.String = {'v_e_L','(km/s)'};
  ylabel(hca,{'v_e_L','(km/s)'},'interpreter','tex');
end
if 0 % Ve par
  hca = irf_panel('Ve par');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaVe1par.tlim(tint),mvaVe2par.tlim(tint),mvaVe3par.tlim(tint),mvaVe4par.tlim(tint)},'comp');    
  hca.YLabel.String = {'v_e_L','(km/s)'};
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
if 1 % EL
  hca = irf_panel('EL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).x,mvaE2.tlim(tint).x,mvaE3.tlim(tint).x,mvaE4.tlim(tint).x},'comp');
  hca.YLabel.String = {'E_L','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EM
  hca = irf_panel('EM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).y,mvaE2.tlim(tint).y,mvaE3.tlim(tint).y,mvaE4.tlim(tint).y},'comp');
  hca.YLabel.String = {'E_M','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, L 4sc
  hca = irf_panel('E + vexB L');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.x+1*mvaVexB1.resample(mvaE1.time).x,...
                mvaE2.x+1*mvaVexB2.resample(mvaE2.time).x,...
                mvaE3.x+1*mvaVexB3.resample(mvaE3.time).x,...
                mvaE4.x+1*mvaVexB4.resample(mvaE4.time).x},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_L','(mV/m)'};      
  hca.YLabel.String = {'R_e_L','(mV/m)'};      
end
if 1 % E + vexB, M 4sc
  hca = irf_panel('E + vexB M');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.y+1*mvaVexB1.resample(mvaE1.time).y,...
                mvaE2.y+1*mvaVexB2.resample(mvaE2.time).y,...
                mvaE3.y+1*mvaVexB3.resample(mvaE3.time).y,...
                mvaE4.y+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_M','(mV/m)'};
  hca.YLabel.String = {'R_e_M','(mV/m)'};            
end
if 1 % E + vexB, N 4sc
  hca = irf_panel('E + vexB N');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};   
  hca.YLabel.String = {'R_e_N','(mV/m)'};         
end
if 1 % JeL
  hca = irf_panel('JeL');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).x,mvaJe2.tlim(tint).x,mvaJe3.tlim(tint).x,mvaJe4.tlim(tint).x},'comp');
  hca.YLabel.String = {'J_e_L','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JeM
  hca = irf_panel('JeM');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).y,mvaJe2.tlim(tint).y,mvaJe3.tlim(tint).y,mvaJe4.tlim(tint).y},'comp');
  hca.YLabel.String = {'J_e_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % JeN
  hca = irf_panel('JeN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaJe1.tlim(tint).z,mvaJe2.tlim(tint).z,mvaJe3.tlim(tint).z,mvaJe4.tlim(tint).z},'comp');
  hca.YLabel.String = {'J_e_N','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 0 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaE1.z+1*mvaVexB1.resample(mvaE1.time).z,...
                mvaE2.z+1*mvaVexB2.resample(mvaE2.time).z,...
                mvaE3.z+1*mvaVexB3.resample(mvaE3.time).z,...
                mvaE4.z+1*mvaVexB4.resample(mvaE4.time).z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234'))
  %irf_legend(hca,{'mms 1','mms 2','mms 3','mms 4'},[0.98 0.9],'fontsize',12);  
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-10 10];
end
if 1 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = EdJe?;',ic)
  hca = irf_panel('(E + vexB) dot Je');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJe1,filtEdJe2,filtEdJe3,filtEdJe4},'comp'); 
  hca.YLabel.String = {'R_e\cdot J_e','(nW/m^3)'};  % (mV/m*mvaJe)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % curl(E + vexB)  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('curl(E + vexB)');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe.x,filtRotRe.y,filtRotRe.z},'comp'); 
  hca.YLabel.String = {'\nabla\times R_e','(mV/m/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % curl(E + vexB) dot Bhat  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  avB = (mvaB1+mvaB2.resample(mvaB1)+mvaB3.resample(mvaB1)+mvaB4.resample(mvaB1))/4;
  filtRotRe = filtRotRe.dot(avB.resample(filtRotRe))/avB.abs;
  hca = irf_panel('curl(E + vexB) dot Bhat');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{filtRotRe},'comp'); 
  hca.YLabel.String = {'(\nabla\times R_e)\cdot b','(mV/m/km)'};
  %set(hca,'ColorOrder',mms_colors('xyz'))   
  %irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Y: curl(E + vexB) proxy  
  filtRotRe = mvaRotRe.filt(0,50,[],5);
  hca = irf_panel('Y: curl(E + vexB) proxy');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,{Y},'comp'); 
  hca.YLabel.String = {'f_{ce}(\rho_e/L_P)^2','(Hz)'};  
end
if 0 % Y: rhoe/LgradP  
  
  hca = irf_panel('Y: rhoe/LgradP');
  set(hca,'ColorOrder',mms_colors('x'))
  irf_plot(hca,avRhoe/LgradP,'comp'); 
  hca.YLabel.String = {'\rho_e/L_P'};  
  ylabel(hca,'\rho_e/L_P','interpreter','tex')
end
if 0 % gradPe
  hca = irf_panel('gradPe');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaGradPe.x,mvaGradPe.y,mvaGradPe.z},'comp');
  hca.YLabel.String = {'\nabla \cdot P','(nPa/km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.9],'fontsize',12);
end
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end
if 0 % JN
  hca = irf_panel('JL');
  set(hca,'ColorOrder',mms_colors('1234ab'))
  lines = irf_plot(hca,{mvaJ1.x,mvaJ2.x,mvaJ3.x,mvaJ4.x,mvaAvJ.x,mvaJcurl.x},'comp');   
  hca.YLabel.String = {'J_L','(mvaJe)'};
  set(hca,'ColorOrder',mms_colors('1234ab'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  hca.YLim = [-1200 2000];
end

irf_zoom(h(1:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)'};
for ii = 1:9
  irf_legend(h(ii+3),legends{ii},[0.01 0.9],'color',[0 0 0])
end

%delete(h(1:3))
