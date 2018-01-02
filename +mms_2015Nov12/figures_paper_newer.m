%% Timing 
dt0 = [0 0 0 0];
dt1 = [0.00     0.05     -0.10     -0.08]; v1 = 70*[-0.86,0.37,0.35]*lmn'; % good for second peak and BL
dt2 = [0.00     0.05     -0.05     -0.03]; v2 = 120*[-0.82 0.49 0.30]*lmn';% first part
n1 = irf_norm(v1);
n2 = irf_norm(v2);
itiming = 1;
c_eval('dt = dt?; vv = v?; vvn = vv/norm(vv); n = n?;',itiming)

pshift = 0;

%% Figure 2: 4 sc, MVA, BL, BM, BN, curvBm
tintZoom = irf.tint('2015-11-12T07:19:20.20Z/2015-11-12T07:19:22.20Z');
%tintZoom = irf.tint('2015-11-12T07:19:20.00Z/2015-11-12T07:19:23.50Z');
npanels = 8;
pshift = 4;
h = irf_plot(npanels+pshift);


scrsz = get(groot,'ScreenSize');
figurePostition = scrsz; figurePostition(3)=figurePostition(3)*0.4; figurePostition(4)=figurePostition(4)*0.9;
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
  hca.YLabel.String = {'curv B','(1/km)'};
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