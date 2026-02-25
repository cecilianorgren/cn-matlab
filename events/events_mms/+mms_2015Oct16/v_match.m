tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:33:38.00Z');
c_eval('B?=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint); dmpaB?brst = irf.ts_vec_xyz(dmpaB?brst.time,dmpaB?brst.data(:,1:3)); dmpaB1brst.units = ''nT'';',1:4);
c_eval('E?=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',1:4);
%c_eval('E? = dslE?brst.tlim(tint);',1:4)
%c_eval('B? = dmpaB?brst.tlim(tint);',1:4)
%mms_2015Oct16.mms4_v_gui(E1,E2,E3,E4,1);
dts = [0 -0.06 -0.18 -0.02];
velocity = mms_2015Oct16.mms4_v(dts+tint(1).epochUnix);
% velocity=61.7*[ -0.91 -0.27 -0.31]; % km/s gse
%%
h = irf_plot(8);

dtsE = [0.06, 0, -0.12, 0.06];
dtsBz = [0.06, 0, -0.12, 0.06];

dtsBz = [0.00     -0.02     -0.03      0.03]

hca = irf_panel('Bx');
irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bx timeshifted');
irf_plot(hca,{dmpaB1.x.tlim(tint),dmpaB2.x.tlim(tint),dmpaB3.x.tlim(tint),dmpaB4.x.tlim(tint)},'comp','dt',dtsBz);
hca.YLabel.String = {'B_{x,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('By');
irf_plot(hca,{dmpaB1.y.tlim(tint),dmpaB2.y.tlim(tint),dmpaB3.y.tlim(tint),dmpaB4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz');
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz timeshifted');
irf_plot(hca,{dmpaB1.z.tlim(tint),dmpaB2.z.tlim(tint),dmpaB3.z.tlim(tint),dmpaB4.z.tlim(tint)},'comp','dt',dtsBz);
hca.YLabel.String = {'B_{z,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
irf_legend(hca,{['dt = ' num2str(dtsBz(1))],num2str(dtsBz(2)),num2str(dtsBz(3)),num2str(dtsBz(4))},[0.95 0.05]);

hca = irf_panel('Ex');
irf_plot(hca,{dslE1.x.tlim(tint),dslE2.x.tlim(tint),dslE3.x.tlim(tint),dslE4.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

%
hca = irf_panel('Ex timeshifted');
irf_plot(hca,{dslE1.x.tlim(tint),dslE2.x.tlim(tint),dslE3.x.tlim(tint),dslE4.x.tlim(tint)},'comp','dt',dtsE);
hca.YLabel.String = {'E_{x,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
irf_legend(hca,{['dt = ' num2str(dtsE(1))],num2str(dtsE(2)),num2str(dtsE(3)),num2str(dtsE(4))},[0.95 0.05]);
%
hca = irf_panel('Ey');
irf_plot(hca,{dslE1.y.tlim(tint),dslE2.y.tlim(tint),dslE3.y.tlim(tint),dslE4.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);


% tint_zoom = irf.tint('2015-10-16T10:33:20.000Z',20);
tint_zoom = irf.tint('2015-10-16T10:33:25.500Z/2015-10-16T10:33:31.000Z');
irf_zoom(h,'x',tint_zoom)
irf_zoom(h,'y')

%%
tint = irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:34:00.000Z');
tintMark = irf.tint('2015-10-16T10:33:30.100Z/2015-10-16T10:33:30.500Z');

h = irf_plot(10);

dtsE = [0.06, 0, -0.12, 0.06];
dtsBz = [0.06, 0, -0.12, 0.06];

%dtsBz = [0.00     -0.02     -0.03      0.03];
%dtsE = dtsBz;

hca = irf_panel('Bx');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.x.tlim(tint),dmpaB2brst.x.tlim(tint),dmpaB3brst.x.tlim(tint),dmpaB4brst.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('Bx timeshifted');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.x.tlim(tint),dmpaB2brst.x.tlim(tint),dmpaB3brst.x.tlim(tint),dmpaB4brst.x.tlim(tint)},'comp','dt',dtsBz);
hca.YLabel.String = {'B_{x}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
irf_legend(hca,{['dt = ' num2str(dtsBz(1))],num2str(dtsBz(2)),num2str(dtsBz(3)),num2str(dtsBz(4))},[0.01 0.05]);
irf_legend(hca,{'shifted'},[0.90 0.95]);

hca = irf_panel('By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.y.tlim(tint),dmpaB2brst.y.tlim(tint),dmpaB3brst.y.tlim(tint),dmpaB4brst.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('shifted By');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.y.tlim(tint),dmpaB2brst.y.tlim(tint),dmpaB3brst.y.tlim(tint),dmpaB4brst.y.tlim(tint)},'comp','dt',dtsBz);
hca.YLabel.String = {'B_{y}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
%irf_legend(hca,{['dt = ' num2str(dtsBz(1))],num2str(dtsBz(2)),num2str(dtsBz(3)),num2str(dtsBz(4))},[0.95 0.05]);
irf_legend(hca,{'shifted'},[0.90 0.95]);

hca = irf_panel('Bz');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.z.tlim(tint),dmpaB2brst.z.tlim(tint),dmpaB3brst.z.tlim(tint),dmpaB4brst.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z}','[nT]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz timeshifted');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dmpaB1brst.z.tlim(tint),dmpaB2brst.z.tlim(tint),dmpaB3brst.z.tlim(tint),dmpaB4brst.z.tlim(tint)},'comp','dt',dtsBz);
hca.YLabel.String = {'B_{z}','[nT]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
%irf_legend(hca,{['dt = ' num2str(dtsBz(1))],num2str(dtsBz(2)),num2str(dtsBz(3)),num2str(dtsBz(4))},[0.95 0.05]);
irf_legend(hca,{'shifted'},[0.90 0.95]);

hca = irf_panel('Ex');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dslE1brst.x.tlim(tint),dslE2brst.x.tlim(tint),dslE3brst.x.tlim(tint),dslE4brst.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x}','[mV/m]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

%
hca = irf_panel('Ex timeshifted');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dslE1brst.x.tlim(tint),dslE2brst.x.tlim(tint),dslE3brst.x.tlim(tint),dslE4brst.x.tlim(tint)},'comp','dt',dtsE);
hca.YLabel.String = {'E_{x}','[mV/m]'};
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
%irf_legend(hca,{['dt = ' num2str(dtsE(1))],num2str(dtsE(2)),num2str(dtsE(3)),num2str(dtsE(4))},[0.95 0.05]);
irf_legend(hca,{'shifted'},[0.90 0.95]);

hca = irf_panel('Ey');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dslE1brst.y.tlim(tint),dslE2brst.y.tlim(tint),dslE3brst.y.tlim(tint),dslE4brst.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y}','[mV/m]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ey shifted');
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{dslE1brst.y.tlim(tint),dslE2brst.y.tlim(tint),dslE3brst.y.tlim(tint),dslE4brst.y.tlim(tint)},'comp','dt',dtsE);
hca.YLabel.String = {'E_{y}','[mV/m]'};
%irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);
set(hca,'ColorOrder',mms_colors('1234'))
%irf_legend(hca,{['dt = ' num2str(dtsE(1))],num2str(dtsE(2)),num2str(dtsE(3)),num2str(dtsE(4))},[0.95 0.05]);
irf_legend(hca,{'shifted'},[0.90 0.95]);


% tint_zoom = irf.tint('2015-10-16T10:33:20.000Z',20);
tint_zoom = irf.tint('2015-10-16T10:33:22.500Z/2015-10-16T10:33:34.000Z');
irf_zoom(h,'x',tint_zoom)
irf_zoom(h,'y')
for ii = 7:10; h(ii).YLim = 7*[-1 1]; end
irf_pl_mark(h([2:2:end]),tintMark.epochUnix','yellow')
irf_plot_axis_align
ax = add_length_on_top(h(1),norm(61),1);
ax.XLabel.String = 'km';
%% Plot sc location and the plane normal to velocity
ic = 1:4;
tint_loc = irf.tint('2015-10-16T10:33:30.02Z',0.3); tint_loc = tint+[-10 10];
c_eval('R?loc = mean(gseR?.resample(dmpaB1brst.tlim(tint_loc)).data,1);',ic);
c_eval('R?rel = R?loc-(R1loc+R2loc+R3loc+R4loc)/4;',ic);
Rrels = [R1rel;R2rel;R3rel;R4rel];

% make plane
planeNormal = irf_norm(velocity);
planeVelocity = velocity;
planeSpeed = norm(velocity); % km/s
planeRadius = 10; % km

zFun = @(x,y) -(planeNormal(1)*x+planeNormal(2)*y)/planeNormal(3);

angles = 0:10:360;
angx = acosd([1 0 0]*planeNormal');
angy = acosd([0 1 0]*planeNormal');
x = planeRadius*cosd(angles);%*sind(angy)*sind(angx);
y = planeRadius*sind(angles);%*sind(angy)*cosd(angx);

z = zFun(x,y);

% make plot
markersize = 10;

plh = plot3(R1rel(1),R1rel(:,2),R1rel(:,3)); hold on;
plh.Color = [0 0 0];
plh.Marker = 's';
plh.MarkerSize = markersize;
plh.LineWidth = 2;

plh = plot3(R2rel(1),R2rel(:,2),R2rel(:,3));
plh.Color = [1 0 0];
plh.Marker = 'd';
plh.MarkerSize = markersize;
plh.LineWidth = 2;

plh = plot3(R3rel(1),R3rel(:,2),R3rel(:,3));
plh.Color = [0 1 0];
plh.Marker = 'o';
plh.MarkerSize = markersize;
plh.LineWidth = 2;

plh = plot3(R4rel(1),R4rel(:,2),R4rel(:,3)); 
plh.Color = [0 0 1];
plh.Marker = 'v';
plh.MarkerSize = markersize;
plh.LineWidth = 2;

hplane = patch(x,y,z,[0 0.5 0.7]);
hplane.FaceAlpha = 0.1;
%axis square
axis equal
%patch(Rrels(:,1),Rrels(:,2),Rrels(:,3),'b')

hca=gca;
hca.XLabel.String = 'X (GSE)';
hca.YLabel.String = 'Y (GSE)';
hca.ZLabel.String = 'Z (GSE)';

qh = quiver3(hca,0,0,0,velocity(1),velocity(2),velocity(3),0.2);
qh.LineWidth = 2;
qh.Color = [0.4 0.4 0.7];
legend(hca,{'mms1','mms2','mms3','mms4'})

%hca.XLim = 10*[-1 1];
%hca.YLim = 10*[-1 1];
%hca.ZLim = 10*[-1 1];
%view(planeNormal)
hold off;


