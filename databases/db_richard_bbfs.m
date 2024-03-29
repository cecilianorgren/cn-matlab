localuser = 'cecilia';
%time_table = readtable('/Users/cno062/Data/MMS/Richard_DB/mms_bbfsdb_2017-2021.csv');
time_table = readtable(['/Users/' localuser '/Data/MMS/Richard_DB/mms_bbfsdb_2017-2021.csv']);
t1 = irf_time(convertTo(time_table.Var1,'epochtime'),'epoch>EpochTT');
t2 = irf_time(convertTo(time_table.Var2,'epochtime'),'epoch>EpochTT');

time = ncread(['/Users/' localuser '/Data/MMS/Richard_DB/mms_b_gsm_2017-2021.nc'],'time');
vec_tot = ncread(['/Users/' localuser '/Data/MMS/Richard_DB/mms_b_gsm_2017-2021.nc'],'represent_vec_tot');
xarray = ncread(['/Users/' localuser '/Data/MMS/Richard_DB/mms_b_gsm_2017-2021.nc'],'__xarray_dataarray_variable__');



%finfo = ncinfo('/Users/cno062/Data/MMS/Richard_DB/mms_b_gsm_2017-2021.nc')

%%
%gseR1 = load('/Users/cno062/Data/MMS/DB_Lalti/Proba_full_mms1_pos.mat'); gseR1 = gseR1.dTmpR;
%gsmR1 = irf.geocentric_coordinate_transformation(gseR1,'gse>gsm');

units = irf_units;
path_nc_bffs = ['/Users/' localuser '/Data/MMS/Richard_DB/mms_bbfsdb_2017-2021.nc'];
finfo_bbfs = ncinfo(path_nc_bffs);

varstrs = {finfo_bbfs.Variables.Name};
for ivar = 1:numel(varstrs)
  varstr = varstrs{ivar};
  var = ncread(path_nc_bffs,varstr);
  eval(sprintf('%s = var;',varstr))
end
%%
path_nc_r_gsm = ['/Users/' localuser '/Data/MMS/Richard_DB/mms_r_gsm_betai_2017-2021.nc'];
finfo_betai = ncinfo(path_nc_r_gsm);
varstrs = {finfo_betai.Variables.Name};
for ivar = 1:numel(varstrs)
  varstr = varstrs{ivar};
  var = ncread(path_nc_r_gsm,varstr);
  eval(sprintf('%s = var;',varstr))
end

% finfo.Variables(1).Attributes(18).Value = 'microseconds since 2017-05-04 00:00:01.635116'
time_r =  EpochTT('2017-05-04T00:00:01.635116Z') + time*1e-3; % from ms to s
gsmR1 = irf.ts_scalar(time_r,r_gsm);
dt_r_db = 4.5;

%% Binning
db_tint = EpochTT(['2017-05-04T00:00:00.000Z';'2021-10-13T00:00:00.000Z']);
gsmR1 = gsmR1.tlim(db_tint);
dx = 1.; 
dy = 1.; 
dz = 1.;
xedges = -30:dx:0;  xcenter = xedges(1:end-1) - 0.5*dx; nx = numel(xcenter);
yedges = -30:dy:30; ycenter = yedges(1:end-1) - 0.5*dy; ny = numel(ycenter);
zedges = -10:dz:10; zcenter = zedges(1:end-1) - 0.5*dz; nz = numel(zcenter);


data_r_RE = gsmR1.data*1e3/units.RE;
data_r_RE = r_gsm'*1e3/units.RE;

dt_R = zeros(gsmR1.length,1);
t_R = gsmR1.time - gsmR1.time.start;
dt_R(2:end-1) = 0.5*(t_R(3:end) - t_R(1:end-2));
dt_R(1) = t_R(2) - t_R(1);
dt_R(end) = t_R(end) - t_R(end-1);
dt_R = repmat(dt_r_db,numel(data_r_RE),1);


[Txyz_mms_gse,edges,mid,loc] = histcn(data_r_RE,xedges,yedges,zedges,'AccumData',dt_R);
[Txy_mms_gse,edges,mid,loc] = histcn(data_r_RE(:,1:2),xedges,yedges,'AccumData',dt_R);
[Txz_mms_gse,edges,mid,loc] = histcn(data_r_RE(:,[1 3]),xedges,zedges,'AccumData',dt_R);

[Nxyz,edges,mid,loc] = histcn([x,y,z]*1e3/units.RE,xedges,yedges,zedges);
[Nxy,edges,mid,loc] = histcn([x,y]*1e3/units.RE,xedges,yedges);
[Nxz,edges,mid,loc] = histcn([x,z]*1e3/units.RE,xedges,zedges);
[Ny,edges,mid,loc] = histcn([y]*1e3/units.RE,yedges);

[Txyz,edges,mid,loc] = histcn([x,y,z]*1e3/units.RE,xedges,yedges,zedges,'AccumData',t);
[Txy,edges,mid,loc] = histcn([x,y]*1e3/units.RE,xedges,yedges,'AccumData',t);
[Txz,edges,mid,loc] = histcn([x,z]*1e3/units.RE,xedges,zedges,'AccumData',t);
[Ty,edges,mid,loc] = histcn([y]*1e3/units.RE,yedges,'AccumData',t);


[Nxy_df,edges,mid,loc] = histcn([x,y]*1e3/units.RE,xedges,yedges,'AccumData',is_df);
%[Txy_df,edges,mid,loc] = histcn([x,y]*1e3/units.RE,xedges,yedges,'AccumData',double(is_df).*double(t_df));

nrows = 2;
ncols = 2;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % T_mms(x,y)
  hca = h(isub); isub = isub + 1;
  zz = zeros(nx+1,ny+1);
  surf(hca,xedges,yedges,zz',Txy_mms_gse'/(60*60))
  view(hca,[0 0 1])
  shading(hca,'flat')
  hca.XLabel.String = 'x (R_E)';
  hca.YLabel.String = 'y (R_E)';
  colormap(hca,pic_colors('candy4'))
  hb = colorbar(hca);
  hb.YLabel.String = 'Time (hours)';
end
if 1 % T(x,y)
  hca = h(isub); isub = isub + 1;
  %pcolor(hca,xcenter,ycenter,squeeze(sum(Nxyz,3))')
  %[XED,ZED] = ndgrid(xedges,yedges);
  %surf(hca,XED,ZED,ZED*0,Nxy')
  zz = zeros(nx+1,ny+1);
  surf(hca,xedges,yedges,zz',Txy'/(60*60))
  view(hca,[0 0 1])
  shading(hca,'flat')
  hca.XLabel.String = 'x (R_E)';
  hca.YLabel.String = 'y (R_E)';
  colormap(hca,pic_colors('candy4'))
  hb = colorbar(hca);
  hb.YLabel.String = 'Time (hours)';
end
if 1 % T(x,y)/T_mms(x,y)
  hca = h(isub); isub = isub + 1;
  %pcolor(hca,xcenter,ycenter,squeeze(sum(Nxyz,3))')
  %[XED,ZED] = ndgrid(xedges,yedges);
  %surf(hca,XED,ZED,ZED*0,Nxy')
  zz = zeros(nx+1,ny+1);
  surf(hca,xedges,yedges,zz',Txy'./Txy_mms_gse')
  view(hca,[0 0 1])
  shading(hca,'flat')
  hca.CLim = [0 1];
  hca.XLabel.String = 'x (R_E)';
  hca.YLabel.String = 'y (R_E)';
  colormap(hca,pic_colors('candy4'))
  hb = colorbar(hca);
  hb.YLabel.String = 'Probability';
end
if 0 % N(x,z)
  hca = h(isub); isub = isub + 1;
  %pcolor(hca,xcenter,ycenter,squeeze(sum(Nxyz,3))')
  %pcolor(hca,xcenter,zcenter,Nxz')
  zz = zeros(nx+1,nz+1);
  surf(hca,xedges,zedges,zz',Nxz')
  view(hca,[0 0 1])
  shading(hca,'flat')
  hca.XLabel.String = 'x (R_E)';
  hca.YLabel.String = 'z (R_E)';
  colormap(hca,pic_colors('candy4'))
  hb = colorbar(hca);
  hb.YLabel.String = 'Counts';
end
if 1 % N(x,y)
  hca = h(isub); isub = isub + 1;
  %pcolor(hca,xcenter,ycenter,squeeze(sum(Nxyz,3))')
  %pcolor(hca,xcenter,ycenter,Nxy_df')
  zz = zeros(nx+1,ny+1);
  surf(hca,xedges,yedges,zz',Nxy_df')
  view(hca,[0 0 1])
  shading(hca,'flat')
  hca.XLabel.String = 'x (R_E)';
  hca.YLabel.String = 'y (R_E)';
  colormap(hca,pic_colors('candy4'))
  hb = colorbar(hca);
  hb.YLabel.String = 'Counts';
end