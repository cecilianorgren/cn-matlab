%% Integrate trajectories for constant EB from simulation using gridden interpolant

% Integrate trajectories
vy = pictmp.viy;
vx = pictmp.vix;
vz = pictmp.viz;
F_vx = griddedInterpolant(X,Z,vx);
F_vy = griddedInterpolant(X,Z,vy);
F_vz = griddedInterpolant(X,Z,vz);


x0 = [20 21 22 23 24];
z0 = [3 3 3 3 3];

ip = 1;

v0 = [0 0.00 -0.05];
vx0 = F_vx(x0',z0');
vy0 = F_vy(x0',z0');
vz0 = F_vx(x0',z0');

tstart = 0;
tstop = 150;
m = 100;
q = 1;
%clear particles
for ipart = 4:numel(x0)
  r0 = [x0(ipart) 0 z0(ipart)];
  v0 = [vx0(ipart) vy0(ipart) vz0(ipart)];
  out = pic(it).integrate_trajectory_constant_EB(r0,v0,tstart,tstop,m,q);
  particles(ipart) = out;
end
%plot3(gca,out.x,out.y,out.z)
%%


%%


% See integrate_trajectories
no02m = PIC('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/fields.h5');
%trp =
%PICTraj('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/trajectories_paul.h5'); % not yet there

%% Make overview plot to pick starting locations
pic = no02m.twpelim(15000).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);
%pic.plot_map({'vex','vez','vtperp(4)'}','clim',{[-5 5],[-1 1],[0 5]},'A',0.5)
%pic.plot_map({'vex','vez','vtperp(4)','Ez','vey','vy([4 6])','ne','Ez.*Bx./Bx./Bx'}','clim',{[-5 5],[-1 1],[0 5],[-1 1],[-5 5],[-5 5],[0 0.3],[-5 5]},'A',0.5)
pic.plot_map({'Ez','ne','vey','Ez.*Bx./Bx./Bx','pxy([4 6])','divpy([4 6])','Ey'}','clim',{[-1 1],[0 0.3],[-5 5],[-5 5],0.003*[-1 1],0.3*[-1 1],0.3*[-1 1]},'A',0.5)

%%
traj = PICTraj('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/trajectories.h5');
ntr_pre = traj.ntr;

t0 = no02m.twpelim(16000).twci;
tspan = [t0-5 t0 t0+10];

r0 = [102,0,1];

vx = mean(pic.xlim(r0(1)+[-0.02 0.02]).zlim(r0(3)+[-0.02 0.02]).vx(4),'all');
vy = mean(pic.xlim(r0(1)+[-0.02 0.02]).zlim(r0(3)+[-0.02 0.02]).vy(4),'all') + 10*mean(pic.xlim(r0(1)+[-0.02 0.02]).zlim(r0(3)+[-0.02 0.02]).vtperp(4),'all');
vz = mean(pic.xlim(r0(1)+[-0.02 0.02]).zlim(r0(3)+[-0.02 0.02]).vz(4),'all') + 10*mean(pic.xlim(r0(1)+[-0.02 0.02]).zlim(r0(3)+[-0.02 0.02]).vtperp(4),'all');
v0 = [0, vy, vz];

m = 1/100;
q = -1;
tr_tmp = no02m.integrate_trajectory(r0,v0,tspan,m,q);


h5write_trajs('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/trajectories.h5',tr_tmp,'id',ntr_pre+1)
%% Plot
%traj = PICTraj('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/trajectories.h5');
itr = traj.ntr;

tr = traj(itr-1);


tr = tr.tlim([80 84.3]);
%tr_tmp = tr.tlim([81 85]);
tr_tmp = tr;

h = setup_subplots(2,2);
isub = 1;

if 0
  hca = h(isub); isub = isub + 1;
  scatter3(hca,tr.x,tr.y,tr.z,[],tr.t)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'time';
  colormap(hca,irf_colormap('waterfall'))
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter3(hca,tr.x,tr.y,tr.z,[],tr.Wy)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'W_y';
  colormap(hca,pic_colors('blue_gray_red'))
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
end

hca = h(isub); isub = isub + 1;
scatter(hca,tr.y,tr.z,[],tr.t)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'time';
colormap(hca,irf_colormap('waterfall'))
hca.XLabel.String = 'y';
hca.YLabel.String = 'z';

hca = h(isub); isub = isub + 1;
scatter(hca,tr.y,tr.z,[],tr.vy)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'vy';
colormap(hca,irf_colormap('waterfall'))
hca.XLabel.String = 'y';
hca.YLabel.String = 'z';

hca = h(isub); isub = isub + 1;
scatter(hca,tr.y,tr.z,[],tr.Wy)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'W_y';
colormap(hca,pic_colors('candy2'))
hca.XLabel.String = 'y';
hca.YLabel.String = 'z';
hca.CLim = 0.6*abs(max(hca.CLim))*[-1 1];

hca = h(isub); isub = isub + 1;
scatter(hca,tr.y,tr.z,[],tr.Bx)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'B_x';
colormap(hca,pic_colors('candy2'))
hca.XLabel.String = 'y';
hca.YLabel.String = 'z';
hca.CLim = 1.0*abs(max(hca.CLim))*[-1 1];


if 0
hca = h(isub); isub = isub + 1;
scatter(hca,tr.x,tr.z,[],tr.t)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'time';
colormap(hca,irf_colormap('waterfall'))
hca.XLabel.String = 'x';
hca.YLabel.String = 'z';
end
if 0
hca = h(isub); isub = isub + 1;
scatter(hca,tr.x,tr.z,[],tr.vy)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'vy';
colormap(hca,irf_colormap('waterfall'))
hca.XLabel.String = 'x';
hca.YLabel.String = 'z';
end
if 0
hca = h(isub); isub = isub + 1;
scatter(hca,tr.x,tr.z,[],tr.Wy)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'W_y';
colormap(hca,pic_colors('blue_gray_red'))
hca.XLabel.String = 'x';
hca.YLabel.String = 'z';
end
if 0 % (x,y,z)
  hca = h(isub); isub = isub + 1;
  plot3(hca,tr.x,tr.y,tr.z)
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hca.ZLabel.String = 'z';
end
if 0 %(y,z)
  hca = h(isub); isub = isub + 1;
  plot(hca,tr.y,tr.z)
  hca.XLabel.String = 'y';
  hca.YLabel.String = 'z';
end

if 0
hca = h(isub); isub = isub + 1;
tr_tmp = tr.tlim([73 85]);
plot(hca,tr_tmp.vy,tr_tmp.vz)
hcb = colorbar('peer',hca);
%hcb.YLabel.String = 'W_y';
%colormap(hca,pic_colors('blue_gray_red'))
hca.XLabel.String = 'v_y';
hca.YLabel.String = 'v_z';
end
if 0
hca = h(isub); isub = isub + 1;
scatter(hca,tr_tmp.vy,tr_tmp.vz,[],tr_tmp.t)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'time';
colormap(hca,irf_colormap('waterfall'))
hca.CLim = [tr.t(1) tr.t(end)];
hca.XLabel.String = 'v_y';
hca.YLabel.String = 'v_z';
end
if 0
hca = h(isub); isub = isub + 1;
scatter(hca,tr_tmp.vy,tr_tmp.vz,[],tr_tmp.Wy)
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'z';
colormap(hca,irf_colormap('waterfall'))
%hca.CLim = [tr.t(1) tr.t(end)];
hca.XLabel.String = 'v_y';
hca.YLabel.String = 'v_z';
end
%% Figure for paper
itr = traj.ntr;

tr = traj(itr-2);

%tr = tr.tlim([80 84.3]);
tr = tr.tlim([79.2 89]);
%tr_tmp = tr.tlim([81 85]);
tr_tmp = tr;

h = setup_subplots(1,1);
isub = 1;


hca = h(isub); isub = isub + 1;
hs = scatter(hca,tr.y*sqrt(100),tr.z*sqrt(100),30,tr.Wy);
hs.MarkerFaceColor = 'flat';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'dW = E_ydy';
%colormap(hca,pic_colors('blue_red'))
%colormap(hca,irf_colormap('bluered'))
hca.XLabel.String = 'y (d_e)';
hca.YLabel.String = 'z (d_e)';
%hca.CLim = 0.6*abs(max(hca.CLim))*[-1 1];


if 0
hca = h(isub); isub = isub + 1;
hs = scatter(hca,tr.x*sqrt(100),tr.z*sqrt(100),[],tr.Wy);
hs.MarkerFaceColor = 'flat';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'W_y';
%colormap(hca,pic_colors('blue_red'))
hca.XLabel.String = 'x';
hca.YLabel.String = 'z';
%hca.CLim = 0.6*abs(max(hca.CLim))*[-1 1];

hca = h(isub); isub = isub + 1;
hs = scatter(hca,tr.vy*sqrt(100),tr.vz*sqrt(100),[],tr.Wy);
hs.MarkerFaceColor = 'flat';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'W_y';
%colormap(hca,irf_colormap('waterfall'))
%colormap(hca,pic_colors('blue_red'))
%hca.CLim = [tr.t(1) tr.t(end)];
hca.XLabel.String = 'v_y';
hca.YLabel.String = 'v_z';
%hca.CLim = 0.9*abs(max(hca.CLim))*[-1 1];
end
hlinks = linkprop(h,{'CLim'});
colormap(irf_colormap('waterfall'))
colormap(pic_colors('blue_red'))
%colormap('parula')
h(1).CLim = max(tr.Wy)*[-1 1];

grid on
box on
fontsize = 18;
%hca.FontSize = 16;
hca.FontSize = fontsize;
hca.Position(2) = 0.16;

hcb.Location = 'north';
hcb.Label.String = 'negative \leftarrow dW = -E_ydy \rightarrow positive';
hcb.Position = [0.5269    0.1816    0.3239    0.0515];
hcb.Label.Position = [0 1.4 0];
hcb.YTick = [];
hcb.Label.FontSize = fontsize;

annotation('textarrow',[0.7 0.75],[0.9 0.9],'String',{'ExB drift in the inflow '},'fontsize',fontsize,'fontweight','light')
annotation('textarrow',[0.7 0.7],[0.55 0.65],'String',{'gradual acceleration by E_y'},'fontsize',fontsize,'fontweight','light')
annotation('textarrow',[0.47 0.46],[0.48 0.55],'String',{'increased acceleration by E_y','during meandering orbit'},'fontsize',fontsize,'fontweight','light','horizontalalignment','center')
annotation('textarrow',[0.35 0.3],[0.8 0.77],'String',{'trapped bouncing','in the outflow'},'fontsize',fontsize,'fontweight','light','horizontalalignment','center')

if 0 % just xy axis on the right
  ha = annotation('arrow',[0.82 0.92],[0.625 0.625]);
  ha = annotation('arrow',[0.82 0.82],[0.625 0.9]);
  text(17,-2,'y','fontsize',fontsize,'fontweight','light','horizontalalignment','left')
  text(7,12,'z','fontsize',fontsize,'fontweight','light','horizontalalignment','left')
else % entire z=0 axis
  ha = annotation('arrow',[0.17 0.82],[0.625 0.625]); % y
  ha = annotation('arrow',[0.17 0.17],[0.625 0.9]); % z
  text(6,-2,'y','fontsize',fontsize,'fontweight','light','horizontalalignment','left')
  text(-78,12,'z','fontsize',fontsize,'fontweight','light','horizontalalignment','left')
  
end
%hca.Visible = 'off';