tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');
tint_df = irf.tint('2017-06-22T03:01:20.00Z/2017-06-22T03:01:43.00Z');

ic = 1:3; % sc 4 not there?

%% Load data
% Magnetic field
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);

% Electric field
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);

% Particle distributions
c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)

% Densities
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% vExB
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s

%% Remove low counts
iDist = iPDist1.tlim(tint_df).elim(eint);
% first try: just remove everything below certain energy
iDist_rmE10 = iDist;
iDist_rmE10.data(:,1:10,:,:) = 0;
iDist_rmE12 = iDist;
iDist_rmE12.data(:,1:12,:,:) = 0;

%% Make reduced distributions
% 1D: vipar
% 2D: vix, vit, viz

eint = [000 40000];
vint = [-Inf Inf];



par = dmpaB1.resample(iDist).norm;

x = [1 0 0];
y = [0 1 0];
z = [0 0 1];

vg = [-2200:50:2200];


tic; if1D = iDist.reduce('1D',par,'vg',vg,'base','cart'); toc % reduced distribution along B
tic; if2D_xy = iDist.reduce('2D',x,y,'vg',vg,'base','cart'); toc
tic; if2D_xz = iDist.reduce('2D',x,z,'vg',vg,'base','cart'); toc
tic; if2D_yz = iDist.reduce('2D',y,z,'vg',vg,'base','cart'); toc

tic; if1D_rmE10 = iDist_rmE10.reduce('1D',par,'vg',vg,'base','cart'); toc % reduced distribution along B
tic; if2D_xy_rmE10 = iDist_rmE10.reduce('2D',x,y,'vg',vg,'base','cart'); toc
tic; if2D_xz_rmE10 = iDist_rmE10.reduce('2D',x,z,'vg',vg,'base','cart'); toc
tic; if2D_yz_rmE10 = iDist_rmE10.reduce('2D',y,z,'vg',vg,'base','cart'); toc

tic; if1D_rmE12 = iDist_rmE12.reduce('1D',par,'vg',vg,'base','cart'); toc % reduced distribution along B
tic; if2D_xy_rmE12 = iDist_rmE12.reduce('2D',x,y,'vg',vg,'base','cart'); toc
tic; if2D_xz_rmE12 = iDist_rmE12.reduce('2D',x,z,'vg',vg,'base','cart'); toc
tic; if2D_yz_rmE12 = iDist_rmE12.reduce('2D',y,z,'vg',vg,'base','cart'); toc

%% Figure: Timeseries of reduced distribution
npanels = 4;
h = irf_plot(npanels);
isub = 1; 

clim = [-3.5 0.3];
if 1 % B
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,gseB1);
end
if 1 % fi_red original
  hca = h(isub); isub = isub + 1;
  irf_spectrogram(hca,if1D.specrec);
  hca.CLim = clim;
end
if 1 % fi_red entire energy channels set to 0
  hca = h(isub); isub = isub + 1;
  irf_spectrogram(hca,if1D_rmE10.specrec);
  hca.CLim = clim;
end
if 1 % fi_red entire energy channels set to 0
  hca = h(isub); isub = isub + 1;
  irf_spectrogram(hca,if1D_rmE12.specrec);
  hca.CLim = clim;
end

colormap(pic_colors('candy'))
irf_plot_axis_align
irf_zoom(h,'x',tint_df)

%% Figure: Compare slices to reduced 2D distributions
nrows = 3;
ncols = 4;
npanels = nrows*ncols;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

vlim = [-2000 2000];
clim_red = [-10 -5];
clim_slice = [-28 -22];
time = irf_time('2017-06-22T03:01:32.30Z','utc>EpochTT');
tint_2d_plot = irf.tint('2017-06-22T03:01:32.00Z/2017-06-22T03:01:33.00Z');
tint_2d_plot = tint_2d_plot + [0.2 -0.2];
tint_2d_plot_str = tint_2d_plot.utc;
tint_2d_plot_str = tint_2d_plot_str(:,12:23);

if 1 % 2D reduced xy
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_xy.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = {'Reduced, E-channels 1-32',[tint_2d_plot_str(1,:), '-' ,tint_2d_plot_str(2,:)]};
end
if 1 % 2D reduced xz
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_xz.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced yz
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_yz.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced xy low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_xy_rmE10.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = {'Reduced, E-channels: 11-32',[tint_2d_plot_str(1,:), '-' ,tint_2d_plot_str(2,:)]};
end
if 1 % 2D reduced xz low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_xz_rmE10.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced yz low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_yz_rmE10.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced xy low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_xy_rmE12.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = {'Reduced, E-channels: 13-32',[tint_2d_plot_str(1,:), '-' ,tint_2d_plot_str(2,:)]};
end
if 1 % 2D reduced xz low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_xz_rmE12.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced yz low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_surf,h_axis,h_all] = if2D_yz_rmE12.tlim(tint_2d_plot).plot_plane(hca);
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced xy low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_axis,h_cb] = mms.plot_projection(hca,iDist.tlim(tint_2d_plot),'xyz',[x;y;z]);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
  hca.Title.String = {'Slice', hca.Title.String{:}};
end
if 1 % 2D reduced xz low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_axis,h_cb] = mms.plot_projection(hca,iDist.tlim(tint_2d_plot),'xyz',[x;z;-y]);
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_z';
end
if 1 % 2D reduced yz low energy bins zet to zero
  hca = h(isub); isub = isub + 1;
  [h_axis,h_cb] = mms.plot_projection(hca,iDist.tlim(tint_2d_plot),'xyz',[y;z;x]);
  hca.XLabel.String = 'v_y';
  hca.YLabel.String = 'v_z';
end


for ipanel = 1:9%npanels
  h(ipanel).XLim = vlim;
  h(ipanel).YLim = vlim;
  h(ipanel).CLim = clim_red;
  h(ipanel).Box = 'on';
  axis(h(ipanel),'square')
end
for ipanel = 10:npanels
  h(ipanel).XLim = vlim;
  h(ipanel).YLim = vlim;
  h(ipanel).CLim = clim_slice;
  h(ipanel).Box = 'on';
  axis(h(ipanel),'square')
end
for ipanel = 1:npanels
  h(ipanel).FontSize = 10;
end
hlink = linkprop(h,{'XLim','YLim'});
colormap(pic_colors('candy'))

%colormap(flipdim(pic_colors('blue_white'),1))




