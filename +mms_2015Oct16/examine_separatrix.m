tint = irf.tint('2015-10-16T10:33:29.50Z/2015-10-16T10:33:31.00Z');
c_eval('E? = dslE?brst.tlim(tint);',1:4)
c_eval('B? = dmpaB?brst.tlim(tint);',1:4)
%mms.db_init('local_file_db','/Volumes/Samsung/data');
%mms_2015Oct16.mms4_v_gui(E1,E2,E3,E4,1);
dts = [0 -0.06 -0.18 -0.02];
tmean = tint(1)+(tint(2)-tint(1))/2;
%v = mms_2015Oct16.mms4_v(dts+tmean.epochUnix);
velocity = 61.7 * [-0.91 -0.27 -0.31];

%% Integrate electric field * v
ic = 1:4;
c_eval('EdB? = irf_edb(irf.ts_vec_xy(dslE?brst.time,dslE?brst.data(:,1:2)),dmpaB?.resample(dslE?brst.time)); EdB?.units = ''mV/m'';',ic)
%c_eval('mvaEdb?=irf.ts_vec_xyz(EdB?.time,[EdB?.data*mva(1,:)'' EdB?.data*mva(2,:)'' EdB?.data*mva(3,:)'']);',ic)
c_eval('E? = EdB?;',ic);

c_eval(['E?filt = irf_filt(E?,0.05,50,[],3);',...
        'E?ts = irf.ts_vec_xyz(irf_time(E?mat(:,1),''epoch>utc''),E?mat(:,2:4));',...
        'intE?mat = irf_integrate([E?filt.time.epochUnix E?filt.data],tint(1).epochUnix);',...
        'intE? = irf.ts_vec_xyz(irf_time(intE?mat(:,1),''epoch>utc''),intE?mat(:,2:4));',...
        'Phi? = irf.ts_scalar(irf_time(intE?mat(:,1),''epoch>utc''),intE?mat(:,2:4)*(1)*velocity'');',...
        'Phi?x = irf.ts_scalar(irf_time(intE?mat(:,1),''epoch>utc''),intE?mat(:,2)*(1)*velocity(1)'');'],ic)

%% Plot sc location and the plane normal to velocity
ic = 1:4;
tint_loc = irf.tint('2015-10-16T10:33:30.02Z',0.3); tint_loc = tint+[-10 10];
c_eval('R?loc = mean(gseR?.resample(dmpaB1brst.tlim(tint_loc)).data,1);',ic);
c_eval('R?rel = R?loc-(R1loc+R2loc+R3loc+R4loc)/4;',ic);
Rrels = [R1rel;R2rel;R3rel;R4rel];

% make plane
planeNormal = irf_norm(v);
planeVelocity = v;
planeSpeed = norm(v); % km/s
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

%hca.XLim = 10*[-1 1];
%hca.YLim = 10*[-1 1];
%hca.ZLim = 10*[-1 1];
%view(planeNormal)
hold off;

%% Plot integrated electric fields
ic = 1;

tint = irf.tint('2015-10-16T10:33:29.50Z/2015-10-16T10:33:31.00Z') + 0.5*[-1 1];

h = irf_plot(6) ;
hca = irf_panel('B');
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B_{?,DMPA}',ic),'(nT)'};
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel('Ex');
irf_plot(hca,{dslE1brst.x,dslE2brst.tlim(tint).x,dslE3brst.tlim(tint).x,dslE4brst.tlim(tint).x},'comp');
hca.YLabel.String = {'E_{x}','(mV/m)'};
irf_legend(hca,{'mms1','mms2','mms3','mms3'},[0.95 0.95]);


hca = irf_panel('Ex filtered');
irf_plot(hca,{E1filt.tlim(tint).x,E2filt.tlim(tint).x,E3filt.tlim(tint).x,E4filt.tlim(tint).x},'comp');
hca.YLabel.String = {'E_{x,filtered}','(mV/m)'};
irf_legend(hca,{'mms1','mms2','mms3','mms3'},[0.95 0.95]);

hca = irf_panel('Ey');
irf_plot(hca,{dslE1brst.tlim(tint).y,dslE2brst.tlim(tint).y,dslE3brst.tlim(tint).y,dslE4brst.tlim(tint).y},'comp');
hca.YLabel.String = {'E_{y}','(mV/m)'};
irf_legend(hca,{'mms1','mms2','mms3','mms3'},[0.95 0.95]);


hca = irf_panel('Phi');
irf_plot(hca,{Phi1.tlim(tint),Phi2.tlim(tint),Phi3.tlim(tint),Phi4.tlim(tint)},'comp');
hca.YLabel.String = {'\phi','(V)'};
irf_legend(hca,{'mms1','mms2','mms3','mms3'},[0.95 0.95])

hca = irf_panel('Phi x');
irf_plot(hca,{Phi1x.tlim(tint),Phi2x.tlim(tint),Phi3x.tlim(tint),Phi4x.tlim(tint)},'comp');
hca.YLabel.String = {'\phi_x','(V)'};
irf_legend(hca,{'mms1','mms2','mms3','mms3'},[0.95 0.95]);

axtop = add_length_on_top(h(1),norm(velocity),0.25);
axtop.XLabel.String = ['km (v = ' num2str(norm(velocity),'%.1f') ' km/s)'];

if 0 % Ez
  hca = irf_panel('Ez');
  irf_plot(hca,{dslE1brst.tlim(tint).z,dslE2brst.tlim(tint).z,dslE3brst.tlim(tint).z,dslE4brst.tlim(tint).z},'comp');
  hca.YLabel.String = {'E_{z}','[mV/m]'};
  irf_legend(hca,{'mms1','mms2','mms3','mms3'},[0.95 0.95]);
end




irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h,'y')