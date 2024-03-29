%% Load data
tint = irf.tint('2017-07-25T22:04:44.80Z/2017-07-25T22:11:31.83Z');
units = irf_units;
ic = 1;

irf.log('critical') % suppress som putput
localuser = datastore('local','user');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   


c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);

c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);

c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)

c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic); % missing some ancillary data
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic); % missing some ancillary data

%% Prepare data
% Select distribution by time/index
dist = ePDist1(2000);
m = dist.mass;

% Get bulk speed in respective direction
v1d = dbcsVe1.resample(dist.time).dot(v1).data;
v2d = dbcsVe1.resample(dist.time).dot(v2).data;
v3d = dbcsVe1.resample(dist.time).dot(v3).data;

% Get pressure components
P = squeeze(gsePe1.resample(dist.time).data);

% Reduce distribution
v_grid = 80e3*linspace(-1,1,100);
v1 = [1 0 0];
v2 = [0 1 0];
v3 = [0 0 1];
f2D = dist.elim([100 40000]).reduce('2D',v1,v2,'vg',v_grid);


v_grid_center = f2D.depend{1}(1,:); % km/s
dv = v_grid_center(2) - v_grid_center(1); % km/s
[V1_base,V2_base] = ndgrid(v_grid_center,v_grid_center); % km/s
V1 = V1_base - v1d; % km/s
V2 = V2_base - v2d; % km/s

fvv = squeeze(f2D.data).*V1.*V2*1e3*1e3; % *1e-3*1e-3 to put v in m/s
f2D_vv = f2D.clone(f2D.time,reshape(fvv,[1 size(fvv)])); % s^2/m^5*(m/s)^2 = 1/m^3
p12 = m*sum(fvv(:))*dv*1e3*dv*1e3; % Pa
p12 = p12*1e9; % nPa

fvv11 = squeeze(f2D.data).*V1.*V1*1e3*1e3; % *1e-3*1e-3 to put v in m/s
p11 = m*sum(fvv11(:))*dv*1e3*dv*1e3; % Pa
p11 = p11*1e9; % nPa

% How to get the units right (we want nPa)
% Pa = kg/(ms^2)
% [fvv] = m^-3
% [fvv*dv^2] = m^-3*(m/s)^2 = 1/(m^1s^2)
% [fvv*dv^2*m] = kg/(m^1s^2) - ok, this is what we want

E = 0.5*m*sqrt((V1_base*1e3).^2 + (V2_base*1e3).^2);
E_frame = 0.5*m*sqrt((V1*1e3).^2 + (V2*1e3).^2);
fE = squeeze(f2D.data).*E;
fE_frame = squeeze(f2D.data).*E_frame;



% Plot data
nrows = 3;
ncols = 2;
h = setup_subplots(nrows,ncols,'vertical');
isub = 1;

if 1 % f
  hca = h(isub); isub = isub + 1;
  f2D.plot_plane(hca)
  colormap(hca,pic_colors('candy4'))
end
if 0
%hca = h(isub); isub = isub + 1;
%f2D_vv.plot_plane(hca)
end
if 1 % f*v*v
  hca = h(isub); isub = isub + 1;
  vscale = 1e-3;
  pcolor(hca,v_grid_center*vscale,v_grid_center*vscale,squeeze(fvv)')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = {'f(v1,v2)*(v_1-v_{1,bulk})*(v_2-v_{2,bulk})','(m^{-3})'};
  %hcb.YLabel.String = 'f(v1,v2)*(v_1-v_{1,bulk})*(v_2-v_{2,bulk})';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(fvv(:)))*[-1 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 'v_1 (10^3 km/s)';
  hca.YLabel.String = 'v_2 (10^3 km/s)';
  hold(hca,'on')
  plot(hca,v1d*vscale,v2d*vscale,'ko')
  hold(hca,'off')
  hca.Title.String = {sprintf('m*int f*v1*v2*dv*dv = %.g',p12),sprintf('m*int f*v1*v1*dv*dv = %.g',p11)};
end
if 1 % f*v*v, summed over each dimension
  hca = h(isub); isub = isub + 1;
  vscale = 1e-3;
  plot(hca,v_grid_center*vscale,sum(squeeze(fvv),2),v_grid_center*vscale,sum(squeeze(fvv),1))
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 'v (10^3 km/s)';
  hca.YLabel.String = 'int f v_{1}v_{2} dv_{1,2}(...)';
  legend(hca,{'fvv(v_1)','fvv(v_2)'},'box','off')
  %hca.Title.String = {sprintf('m*int f*v1*v2*dv*dv = %.g',p12),sprintf('m*int f*v1*v1*dv*dv = %.g',p11)};
  axis(hca,'square')
end
if 1 % |f*v*v|
  hca = h(isub); isub = isub + 1;
  vscale = 1e-3;
  pcolor(hca,v_grid_center*vscale,v_grid_center*vscale,abs(squeeze(fvv))')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = {'|f(v1,v2)*(v_1-v_{1,bulk})*(v_2-v_{2,bulk})|','(m^{-3})'};
  colormap(hca,pic_colors('candy4'))
  hca.CLim = max(abs(fvv(:)))*[0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 'v_1 (10^3 km/s)';
  hca.YLabel.String = 'v_2 (10^3 km/s)';
  hold(hca,'on')
  plot(hca,v1d*vscale,v2d*vscale,'ko')
  hold(hca,'off')
end
if 1 % f*E
  hca = h(isub); isub = isub + 1;
  vscale = 1e-3;
  pcolor(hca,v_grid_center*vscale,v_grid_center*vscale,abs(squeeze(fE))')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = {'|f(v1,v2)*E|','(...)'};
  colormap(hca,pic_colors('candy4'))
  %hca.CLim = max(abs(fvv(:)))*[0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 'v_1 (10^3 km/s)';
  hca.YLabel.String = 'v_2 (10^3 km/s)';
  hold(hca,'on')
  plot(hca,v1d*vscale,v2d*vscale,'ko')
  hold(hca,'off')
  hca.Title.String = 'f*E where E is in the still frame';
end
if 1 % f*E_frame
  hca = h(isub); isub = isub + 1;
  vscale = 1e-3;
  pcolor(hca,v_grid_center*vscale,v_grid_center*vscale,abs(squeeze(fE_frame))')
  shading(hca,'flat')
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = {'|f(v1,v2)*E|','(...)'};
  colormap(hca,pic_colors('candy4'))
  %hca.CLim = max(abs(fvv(:)))*[0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.XLabel.String = 'v_1 (10^3 km/s)';
  hca.YLabel.String = 'v_2 (10^3 km/s)';
  hold(hca,'on')
  plot(hca,v1d*vscale,v2d*vscale,'ko')
  hold(hca,'off')
  hca.Title.String = 'f*E where E is in the frame of the bulk flow';
end

i2D = [1 2 4 5 6];
hlinks = linkprop(h(i2D),{'YLim','XLim'});
h(1).XLim = v_grid([1 end])*vscale;
h(1).YLim = v_grid([1 end])*vscale;

for ip = i2D
  axis(h(ip),'square')
  h(ip).XTick = h(ip).YTick;
end

h(3).Position([1 3]) = h(2).Position([1 3]);
h(3).XLim = h(2).XLim;



