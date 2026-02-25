% Define time interval
tint = irf.tint('2015-12-02T01:14:14/2015-12-02T01:15:14Z');

%% Load data
ic = 1:4;

disp('Loading electric field...')
c_eval('tic; dslE?brst=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
%Eoffs1 = 2.05; Eoffs2 = 2.58; Eoffs3 = 2.53; Eoffs4 = 1.40;
%c_eval('dslE?brst.data(:,1)=dslE?brst.data(:,1)-Eoffs?;',ic)

disp('Loading magnetic field...')
c_eval('tic; dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint); toc',1:4);
c_eval('if iscell(dmpaB?brst); dmpaB?brst=dmpaB?brst{1}; end')
c_eval('tic; dmpaB?l2pre=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_dmpa'',tint); toc',1:4);
c_eval('tic; gseB?l2pre=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',tint); toc',1:4);

disp('Loading spacecraft potential...')
c_eval('P?brst=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',tint);',ic);

disp('Loading electron distribution data...')
%c_eval('tic; desDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint); toc',[1 2 3 4]);
load /Users/Cecilia/Data/MMS/2015Dec02/e-dist_Jan18.mat

disp('Loading ion distribution data...')
%c_eval('tic; disDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint); toc',[1 2 3 4]);
load /Users/Cecilia/Data/MMS/2015Dec02/i-dist_Jan18.mat

%disp('Constructing ion moments data...')
disp('Loading electron moments...')
c_eval('ne?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_numberDensity'',tint);',ic);

c_eval('vex?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkX'',tint);',ic);
c_eval('vey?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkY'',tint);',ic);
c_eval('vez?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkZ'',tint);',ic);
c_eval('ve?brst=irf.ts_vec_xyz(vex?brst.time,[vex?brst.data vey?brst.data vez?brst.data]);',ic)
c_eval('ve?brst.name = ''mms? ve brst'';')

%% Create coodinate system
% Define shorter time interval for MVA
tint = irf.tint('2015-12-02T01:14:14/2015-12-02T01:15:14Z');
irf_minvar_gui(gseB2l2pre.tlim(tint))
%%
tint = irf.tint('2015-12-02T01:14:53/2015-12-02T01:15:50Z');
[out,l,v] = irf_minvar(dmpaB2l2pre.tlim(tint));

%% Make psd plots
ic = 3;
tint = irf.tint('2015-12-02T01:14:56.6Z',0.04);
tint = tint+2*0.03;

% Get mean magnetic field direction
c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));
c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Projection coordinate system
x = hatB0;
y = hatExB0;
z = cross(x,y);


% Initialize figure
nrows = 3;
ncols = 2;
for ii = 1:(nrows*ncols); h(ii) = subplot(nrows,ncols,ii); end
isub = 1;

energies = [31 108];
c_eval('dist = desDist?',ic)
% Plot particle distribution for all directions but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,dist,'tint',tint,'energy',energies(1),'vectors',{hatB0,'B'},'flat');

% Plot particle distribution for all direction but a single energy
hca = h(isub); isub = isub + 1;
mms.plot_skymap(hca,dist,'tint',tint,'energy',energies(2),'vectors',{hatB0,'B'},'flat');

% Plot projection onto a plane perpendicular to B + the 2 other planes
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[x;y;z],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B';hatE0,'E'});
colormap(hca,'jet')
hca.CLim = [-1 4.5];

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[y;z;x],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B';hatE0,'E'});
colormap(hca,'jet')
hca.CLim = [-1 4.5];

hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,dist,'tint',tint,'xyz',[z;x;y],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B';hatE0,'E'});
colormap(hca,'jet')
hca.CLim = [-1 4.5];

% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint,''energies'',energies,''ylim'',[1e-2 1e6])',ic)

%% Plot several time steps, two planes for each timestep
ic = 3;
tsample =0.033;
tint = irf.tint('2015-12-02T01:14:56.81Z',tsample);
nTints = 5;
tints = EpochTT(tint(1)):(tsample):EpochTT(tint(1)+tsample*(nTints-1));


% Get mean magnetic field direction
c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(tint).time).data);',ic); 
hatB0 = double(irf_norm(B0));
c_eval('E0 = mean(dslE?brst.tlim(tint).data);',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Projection coordinate system
x = hatB0;
y = hatExB0;
z = cross(x,y);


% Initialize figure
nrows = 3;
ncols = nTints;
for ii = 1:(nrows*ncols); h(ii) = subplot(nrows,ncols,ii); end
isub = 1;

for ii = 1:nTints % Plot projection onto a plane perpendicular to B 
  % Get mean magnetic field direction
  c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(irf.tint(tints(ii),tsample)).time).data);',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = mean(dslE?brst.tlim(irf.tint(tints(ii),tsample)).data);',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
  % Projection coordinate system
  x = hatB0;
  y = hatExB0;
  z = cross(x,y);
  hca = h(isub); isub = isub + 1;
  mms.plot_projection(hca,dist,'tint',irf.tint(tints(ii),tsample),'xyz',[x;y;z],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B';hatE0,'E'});
  colormap(hca,'jet')
  hca.CLim = [-1 4.5];
end

for ii = 1:nTints % Plot projection onto a plane perpendicular to B 
  % Get mean magnetic field direction
  c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(irf.tint(tints(ii),tsample)).time).data);',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = mean(dslE?brst.tlim(irf.tint(tints(ii),tsample)).data);',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
  % Projection coordinate system
  x = hatB0;
  y = hatExB0;
  z = cross(x,y);
  hca = h(isub); isub = isub + 1;
  mms.plot_projection(hca,dist,'tint',irf.tint(tints(ii),tsample),'xyz',[y;z;x],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B';hatE0,'E'});
  colormap(hca,'jet')
  hca.CLim = [-1 4.5];
end
for ii = 1:nTints % Plot projection onto a plane perpendicular to ExB? 
  % Get mean magnetic field direction
  c_eval('B0 = mean(gseB?l2pre.resample(dslE?brst.tlim(irf.tint(tints(ii),tsample)).time).data);',ic); 
  hatB0 = double(irf_norm(B0));
  c_eval('E0 = mean(dslE?brst.tlim(irf.tint(tints(ii),tsample)).data);',ic); 
  hatE0 = double(irf_norm(E0));
  hatExB0 = cross(hatE0,hatB0);
  % Projection coordinate system
  x = hatB0;
  y = hatExB0;
  z = cross(x,y);
  hca = h(isub); isub = isub + 1;
  mms.plot_projection(hca,dist,'tint',irf.tint(tints(ii),tsample),'xyz',[z;x;y],'elevationlim',20,'vlim',15000,'vectors',{hatB0,'B';hatE0,'E'});
  colormap(hca,'jet')
  hca.CLim = [-1 4.5];
end
%%
% Plot particle distribution for pitch angles 0 90 and 180
hca = h(isub); isub = isub + 1;
c_eval('mms.plot_cross_section_psd(hca,dist,dmpaB?brst,''tint'',tint,''energies'',energies,''ylim'',[1e-2 1e6])',ic)
