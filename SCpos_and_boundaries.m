disp('Loading spacecraft position...')
R  = mms.get_data('R_gse',tint);
c_eval('gseR? = TSeries(R.time,R.gseR?,''to'',1);',ic);
c_eval('mvaR? = irf.ts_vec_xyz(gseR?.time,[gseR?.dot(L).data gseR?.dot(M).data gseR?.dot(N).data]);')
mvaR0 = (mvaR1.resample(mvaR1.time)+mvaR2.resample(mvaR1.time)+mvaR3.resample(mvaR1.time)+mvaR4.resample(mvaR1.time))/4;
c_eval('mvaRR? = mvaR?-mvaR0; mvaRR? = mvaRR?.resample(irf_time(''2015-10-16T10:33:30.000Z'',''utc>epochTT'')).data;',ic)


%% Figure 2: Plot figure with fields etc of 4 sc, for paper
tintZoom = irf.tint('2015-10-16T10:33:22.00Z/2015-10-16T10:33:33.00Z');
tintZoom = irf.tint('2015-10-16T10:33:24.00Z/2015-10-16T10:33:32.00Z');
npanels = 10;

boundaryTint = EpochTT(['2015-10-16T10:33:27.20Z';...
                        '2015-10-16T10:33:30.35Z']);                
pshift = 3;                      
h = irf_plot(npanels + pshift);
irf_panel('delete1');
irf_panel('delete2');
irf_panel('delete3');
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
if 1 % JM
  hca = irf_panel('JM');
  set(hca,'ColorOrder',mms_colors('1234b'))
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaAvJ.y,mvaJcurl.y},'comp');   
  %lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y,mvaJcurl.y},'comp');   
  lines = irf_plot(hca,{mvaJ1.y,mvaJ2.y,mvaJ3.y,mvaJ4.y},'comp');   
  hca.YLabel.String = {'J_M','(nA/m^2)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  %irf_legend(hca,{'E','v_{e}xB','\nabla \cdot P_e/ne','-v_{e}xB-\nabla \cdot P_e/ne'},[0.98 0.1],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','<J_{1234}>','J_{curl}'},[0.02 0.9],'fontsize',12);
  %irf_legend(hca,{'J_1','J_2','J_3','J_4','J_{curl}'},[0.02 0.2],'fontsize',12);
  %irf_legend(hca,{'J_{curl}'},[0.98 0.2],'fontsize',12,'color',mms_colors('b'));
  %lines.Children(6).LineWidth = 1.5;
  %lines.Children(6).LineStyle = '--';
  %hca.YLim = [-1200 2000];
end
if 1 % Ve par
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
  %hca.YLabel.String = {'v_{e,\perp,L}','(km/s)'};
  ylabel(hca,{'|v_{e,\perp}|','(km/s)'},'interpreter','tex');
end
if 1 % Ve perp L
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
if 1 % EN
  hca = irf_panel('EN');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{mvaE1.tlim(tint).z,mvaE2.tlim(tint).z,mvaE3.tlim(tint).z,mvaE4.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234'))
end
if 1 % E + vexB, 4sc
  hca = irf_panel('E + vexB');
  set(hca,'ColorOrder',mms_colors('1234b'))
 irf_plot(hca,{mvaEVexB1.z,mvaEVexB2.z,mvaEVexB3.z,mvaEVexB4.z},'comp'); 
  hca.YLabel.String = {'(E+v_{e}xB)_N','(mV/m)'};
  hca.YLabel.String = {'E''_N','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzba'))  
  set(hca,'ColorOrder',mms_colors('1234')) 
  hca.YLim = [-10 10];
end
if 1 % (E + vexB)*J, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  %c_eval('filtEdJ? = RedJ?;'); 
  hca = irf_panel('E'' dot J');
  %c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('(E + vexB) dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{RedJ1,RedJ2,RedJ3,RedJ4},'comp'); 
  hca.YLabel.String = {'E''\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  %hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  %hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('E dot J');
  %c_eval('filtEdJ? = EdJ?;'); hca = irf_panel('(E + vexB) dot J');
  set(hca,'ColorOrder',mms_colors('1234'))
  irf_plot(hca,{filtEdJ1,filtEdJ2,filtEdJ3,filtEdJ4},'comp'); 
  hca.YLabel.String = {'E''\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3'
  hca.YLabel.String = {'(E + v_e\timesB)\cdot J_e','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  hca.YLabel.String = {'E\cdot J','(nW/m^3)'};  % (mV/m*nA/m^2)*1e-3
  set(hca,'ColorOrder',mms_colors('1234'))   
end
if 0 % (E + vexB)*Je, 4sc
  %c_eval('filtEdJe? = EdJe?.filt(0,100,[],5);',ic)
  c_eval('filtEdJe? = RedJe?;')
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
if 0 % Length scales
  hca = irf_panel('Length scales');
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot({LgradP,rp1*1e-3,Lp1*1e-3},'comp')
  hca.YLabel.String = {'Length','(km)'};
  set(hca,'ColorOrder',mms_colors('xyz'))   
  irf_legend(hca,{'L_p','\rho_i','L_i'},[0.98 0.9],'fontsize',12);
end


irf_zoom(h(4:end),'x',tintZoom)
irf_zoom(h(:),'y')
irf_plot_axis_align

drawnow

% reset ylims 
hca = irf_panel('BN'); hca.YLim = [-7 7];
hca = irf_panel('BM'); hca.YTick = [-10 0 10];
hca = irf_panel('JM'); hca.YTick = [-1000 0 1000];
hca = irf_panel('E + vexB'); hca.YLim = [-4.5 6.5];
hca = irf_panel('EN'); hca.YLim = [-7.5 6.5];


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)'};
shift = 3; % three deleted panels
legshift = 2; % the two sc configuration plots
for ii = 1:npanels
  irf_legend(h(ii+shift),legends{ii+legshift},[0.01 0.9],'color',[0 0 0])
  h(ii+shift).FontSize = 12;  
  %if labelsOutside
  %  h(ii+shift).Position(3) = h(ii+shift).Position(3)*0.88;
  %end  
  h(ii+shift).YLabel.FontSize = 11;
end

for ii = 1:3; h(ii).Visible = 'off'; end

% Add labels for the different regions
if exist('hleg_mspsep','var'); delete(hleg_mspsep); end
if exist('hleg_outflow','var'); delete(hleg_outflow); end
if exist('hleg_mshsep','var'); delete(hleg_mshsep); end
hca = h(pshift+1);
set(hca,'ColorOrder',mms_colors('11'))
hleg_mspsep = irf_legend(hca,{{'Magnetosphere','inflow'}},[0.2 1],[0 0 0]); hleg_mspsep.VerticalAlignment = 'bottom'; hleg_mspsep.HorizontalAlignment = 'center';
hleg_outflow = irf_legend(hca,{{'Electron','outflow'}},[0.6 1],[0 0 0]); hleg_outflow.VerticalAlignment = 'bottom'; hleg_outflow.HorizontalAlignment = 'center';
hleg_mshsep = irf_legend(hca,{{'Magnetosheath','inflow'}},[0.9 1],[0 0 0]); hleg_mshsep.VerticalAlignment = 'bottom'; hleg_mshsep.HorizontalAlignment = 'center';

% Add borders/separtrices for the different regions
hmark1 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(1).epochUnix,[0.4 0.4 0.4]);
hmark2 = irf_pl_mark(h(pshift+1:npanels+pshift),boundaryTint(2).epochUnix,[0.4 0.4 0.4]);
for ii = 1:numel(hmark1), hmark1(ii).LineStyle = '-'; hmark1(ii).Color = [0.4 0.4 0.4]; end
for ii = 1:numel(hmark2), hmark2(ii).LineStyle = '-.'; hmark2(ii).Color = [0.4 0.4 0.4]; end




%h(10).YLim = 7*[-1 1];
%h(11).YLim = 7*[-1 1];

%% Plot sc positions and boundary orientations
if exist('h2'); delete(h2); end
nrows = 5;
ncols = 3;
h2(1) = subplot(nrows,ncols,1); %h2(1).Position(2) = h2(1).Position(2)+0.02;
h2(2) = subplot(nrows,ncols,2); %h2(2).Position(2) = h2(2).Position(2)+0.02;

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
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'L (km)';

%plot(mvaRR1(3),mvaRR1(1))
hca = h2(2);
hold(hca,'on');
for ic=1:4    
  % put Cluster markers
  %plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
  plot(hca,x{ic}(3),x{ic}(2),mms_marker{ic}{:});  
end
hold(hca,'off');
axis(hca,'square')
hca.XLim = 1.7*[-10 10];
hca.YLim = 1.7*[-10 10];
hca.XDir = 'reverse';
hca.XMinorGrid = 'on';
hca.YMinorGrid = 'on';
hca.Box = 'on';
grid(hca,'on');
hca.XLabel.String = 'N (km)';
hca.YLabel.String = 'M (km)';
hleg = legend({'MMS 1','MMS 2','MMS 3','MMS 4'},'location','EastOutside','fontsize',11);
hleg.Box = 'off';
hleg.Position(1) = hleg.Position(1)+0.2;
hleg.Position(2) = hleg.Position(2)+0.05;

% plot plane
v = [-0.90 0.26 0.36];
v = [-0.90 -0.28 -0.33];
lmnV = [-0.57 -0.03 -0.82];
lmnV = [-0.55 -0.08 -0.83];
lmnV = [-0.45 -0.05 -0.89];
lmnV = [-0.49 -0.05 -0.87];
v = [-0.9139   -0.4006    0.0653];
v = [-0.90 -0.28 -0.33];

mshV = 55*[-0.88 -0.26 -0.40]; % GSE

mspV= 31.8*[-0.88 -0.42 0.24]; % GSE
mshV = 55*[-0.90 -0.28 -0.33]; %GSE
mspVlmn = mspV*[L' M' N'];
mshVlmn = mshV*[L' M' N'];
mspN = irf_norm(mspVlmn);
mshN = irf_norm(mshVlmn);

x = 50*[-1 1];
y = 50*[-1 1];
z = 50*[-1 1];

%zFun = @(x,y,n) -(n(1)*x+n(2)*y)/n(3);
% (planeNormal(1)*x+planeNormal(2)*y)/planeNormal(3) = 0;
%  planeNormal(1)*x+planeNormal(2)*y = 0;
% x = - planeNormal(2)*y/planeNormal(1);
funX = @(y,z,n) (-n(2)*y-n(3)*z)/n(1);
funY = @(x,z,n) (-n(1)*x-n(3)*z)/n(2);
funZ = @(z,y,n) (-n(1)*x-n(2)*y)/n(3);



if exist('hmshN1'); delete(hmshN1); end
if exist('hmshN2'); delete(hmshN2); end
if exist('hmspN1'); delete(hmspN1); end
if exist('hmspN2'); delete(hmspN2); end

hca = h2(1); subplot(nrows,ncols,1);
hold(hca,'on');
%hmshN1 = plot3(hca,funX(y,z,mshN),funY(x,z,mshN),funZ(x,y,mshN),'r');
%hmspN1 = plot3(hca,funX(y,z,mspN),funY(x,z,mspN),funZ(x,y,mspN),'b');
hmshN1 = plot(hca,funZ(x,y,mshN),x+18,'k-.');
hmspN1 = plot(hca,funZ(x,y,mspN)-10,x,'k-');

ht = text(9.5,9,'MSH'); ht.Rotation = 55; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;
ht = text(-11.2,8,'MSP'); ht.Rotation = -80; ht.HorizontalAlignment = 'center'; ht.FontSize = 11;

hold(hca,'off');

y = [-10 10];
hca = h2(2);
hold(hca,'on');
%hmspN1 = plot3(hca,funX(y,z,mshN),funY(x,z,mshN),funZ(x,y,mshN),'r');
%hmspN2 = plot3(hca,funX(y,z,mspN),funY(x,z,mspN),funZ(x,y,mspN),'b');
%hmspN1 = plot(hca,funZ(x,y,mshN),funY(x,z,mshN),'r');
%hmspN2 = plot(hca,funZ(x,y,mshN),funY(x,z,mspN),'b');
hmspN1 = plot(hca,z+10,funY(x,z,mshN),'k-.');
hmspN2 = plot(hca,z-10,funY(x,z,mspN),'k-');
hold(hca,'off');
%h2(1).Position(2) = h2(1).Position(2)+0.02;
%h2(2).Position(2) = h2(2).Position(2)+0.02;
irf_legend(h2(1),'a)',[-0.4 1],'color',[0 0 0])
irf_legend(h2(2),'b)',[-0.4 1],'color',[0 0 0])
%
for ii = 1:2;
  h2(ii).YLabel.FontSize = 12;
  h2(ii).XLabel.FontSize = 12;
  h2(ii).FontSize = 12;
  h2(ii).Position(2) = h2(ii).Position(2)+0.05;
  h2(ii).Position(1) = h2(ii).Position(1)+0.05;
end