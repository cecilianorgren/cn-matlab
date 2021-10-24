%% Spacecraft, time interval, location
savePath = '/Users/cecilia/Research/energy_partitioning';
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
irf.log('critical')
units = irf_units;
ic = 1:4;
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

%% Run through files
data = [];
data_tmp = [];
nfiles = numel(files);

for ifile = 1:30%1:30%nfiles
  tic
  fileId = strsplit(files(ifile).name,'_'); fileId = fileId{5};
  tint = [files(ifile).start files(ifile).stop];
  strdisp = sprintf('#%g/%g, Time: %s - %s.',ifile,nfiles,tint(1).utc,tint(2).utc);
  disp(strdisp)    
  
  %% Load spacecraft location and check location criteria
  R = mms.get_data('R_gse',tint);
  
  if isfield(R,'gseR1') && size(R.gseR1,2) == 4
    c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
  elseif isfield(R,'gseR1')
    c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
  else
    disp('  No position data. Skipping time interval.')
    continue;
  end
  if all([isempty(gseR1),isempty(gseR2),isempty(gseR3),isempty(gseR4)])  
    disp('  No position data. Skipping time interval.')
    continue
  end
  %Rx = gsmR(1).x.data/(units.RE*1e-3);
  %Ry = gsmR(1).y.data/(units.RE*1e-3);
  %Rz = gsmR(1).z.data/(units.RE*1e-3);  
    
  %% Load plasma data
  tint_load = tint + [1 -1]; % avoid problems with overlapping files of different versions
  c_eval('ne?    = mms.get_data(''Ne_fpi_brst_l2'',tint_load,?);',ic);
  c_eval('ni?    = mms.get_data(''Ni_fpi_brst_l2'',tint_load,?);',ic);
  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint_load,?);',ic)
  c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint_load,?);',ic);    
  c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint_load,?);',ic) 
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint_load,?);',ic)
  c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint_load,?);',ic) 
  c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint_load,?);',ic);
  
  % If there's no data (assuming there's no data at all if there's no data
  % for ni*.)
  %c_eval('if isempty(ni?), gseR? = TSeries([]); end',ic)  

  % Downsample to 1s, and save all that data for post-analysis
  timeline = tint_load(1):10:tint_load(2);
  c_eval('ne?    =    ne?.resample(timeline);',ic)
  c_eval('ni?    =    ni?.resample(timeline);',ic)
  c_eval('gseVe? = gseVe?.resample(timeline);',ic)
  c_eval('gseVi? = gseVi?.resample(timeline);',ic)
  c_eval('gsePe? = gsePe?.resample(timeline);',ic)
  c_eval('gsePi? = gsePi?.resample(timeline);',ic)
  c_eval('gseTe? = gseTe?.resample(timeline);',ic)
  c_eval('gseTi? = gseTi?.resample(timeline);',ic)
  c_eval('gseR?  = gseR?.resample(timeline);',ic)
  
  c_eval('data_tmp.R?  = gseR?.data;',ic)
  c_eval('data_tmp.ne? = ne?.data;',ic)
  c_eval('data_tmp.ni? = ni?.data;',ic)
  c_eval('data_tmp.ve? = gseVe?.data;',ic)
  c_eval('data_tmp.vi? = gseVi?.data;',ic)
  c_eval('data_tmp.pe? = gsePe?.data;',ic)
  c_eval('data_tmp.pi? = gsePi?.data;',ic)
  c_eval('data_tmp.te? = gseTe?.data;',ic)
  c_eval('data_tmp.ti? = gseTi?.data;',ic)
  data = cat(1,data,data_tmp);
 
  %c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
  %c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
  %c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
  %c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic) 
  toc 
  %% Plot results
  if 0 % plot results
    %%    
    npanels = 4;
    h = irf_plot(npanels);
    
    if 1 % B
      hca = irf_panel('B');
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_plot(hca,{gsmB.x,gsmB.y,gsmB.z},'comp')
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_legend(hca,{'x','y','z'},[0.98 0.98])
      hca.YLabel.String = {'B_{GSM} (nT)'};
    end
    if not(isempty(gsmE)) % E
      hca = irf_panel('E');
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_plot(hca,{gsmE.x,gsmE.y,gsmE.z},'comp')
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_legend(hca,{'x','y','z'},[0.98 0.98])
      hca.YLabel.String = {'E_{GSM} (mV/m)'};
    end
    if not(isempty(ne)) % ne
      hca = irf_panel('n');
      set(hca,'ColorOrder',mms_colors('12'))
      irf_plot(hca,{ne},'comp')
      set(hca,'ColorOrder',mms_colors('12'))
      irf_legend(hca,{'n_e'},[0.98 0.98])
      hca.YLabel.String = {'n (cm^{-3})'};
    end
    if not(isempty(gsmVi)) % Vi
      hca = irf_panel('Vi');
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_plot(hca,{gsmVi.x,gsmVi.y,gsmVi.z},'comp')
      set(hca,'ColorOrder',mms_colors('xyz'))
      irf_legend(hca,{'x','y','z'},[0.98 0.98])
      hca.YLabel.String = {'V_{i,GSM} (km/s)'};
    end
    
    c_eval('irf_pl_mark(h(?),[df_save(!).t1 df_save(!).t2]);',1:npanels,1:numel(df_save))   
    %c_eval('irf_pl_mark(h(?),df_save(!).t1 + 0.5*(df_save(!).t2-df_save(!).t1),''color'',[0 0 0]);',1:npanels,1:numel(df_save))    
    irf_zoom(h,'x',tint)
    
    if not(isempty(df_save))
      str_print = sprintf('df_save_%s_nDF=%02.0f',df_save(1).data.fileId,numel(df_save));      
    else
      str_print = sprintf('df_save_%s_nDF=%02.0f',fileId,0);
    end
    cn.print(str_print,'path','/Users/cecilia/Research/dipolarization_fronts/initial_search')
  end
  
end
disp('Done for now.')

%% Collect results
% Go through data and 'delete' e.g. R2 is there is no fpi data for MMS2,
% etc.
data_orig = data;
fieldnames = fields(data);
nfields = numel(fieldnames);
doKeep = ones(size(data));
for ievent = 1:numel(data)
  for ifield = 1:nfields    
    tmp = data(ievent).(fieldnames{ifield});
    if isempty(tmp)
      data(ievent).(['R' fieldnames{ifield}(end)]) = [];
      continue
    end
  end
end
%%
R2 = cat(1,data.R2);
ti2 = cat(1,data.ti2);
te2 = cat(1,data.te2);

ti2_scalar = (ti2(:,1,1) + ti2(:,2,2) + ti2(:,3,3))/3;
te2_scalar = (te2(:,1,1) + te2(:,2,2) + te2(:,3,3))/3;
tite2 = ti2_scalar./te2_scalar;

x_edges = -25:1:10;
y_edges = -20:1:20;
[X_EDGES,Y_EDGES] = meshgrid(x_edges,y_edges);

[count,edges,mid,loc] = histcn(R2(:,1:2)*1e3/units.RE,x_edges,y_edges);
[accumTiTe,edges,mid,loc] = histcn(R2(:,1:2)*1e3/units.RE,x_edges,y_edges,'AccumData',tite2);

% Plot
nrows = 2;
ncols = 1;
h = setup_subplots(nrows,ncols);
isub = 1;

if 1 % Scatter plot of Ti/Te
  hca = h(isub); isub = isub + 1;
  scatter(hca,R2(:,1)*1e3/units.RE,R2(:,2)*1e3/units.RE,R2(:,1)*0+5,tite2)
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'T_i/T_e';
end
if 1 % Binned plot of Ti/Te
  hca = h(isub); isub = isub + 1;
  surf(hca,X_EDGES,Y_EDGES,X_EDGES*0,accumTiTe'./count')
  view(hca,[0 0 1]);
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'T_i/T_e';
  shading(hca,'flat')
end





