%% Load data
ic = 1:4;
tint = irf.tint('2015-12-02T01:14:30.00Z/2015-12-02T01:15:15.00Z');
mms.db_init('local_file_db','/Volumes/Nexus/data');
eventPath = '/Users/Cecilia/Research/Events/2015-12-02_011414/';
units = irf_units;

c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);


try
R = mms.get_data('R_gse',tint);
if ~all([isfield(R,'gseR1') isfield(R,'gseR1') isfield(R,'gseR2') isfield(R,'gseR3') isfield(R,'gseR4')])  
  % not the right data fiels, try to load from irfu database instead
  db_info = datastore('mms_db');   
  mms.db_init('local_file_db','/data/mms');
  R = mms.get_data('R_gse',tint);
  mms.db_init('local_file_db',db_info.local_file_db_root);
end
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end
end


%% Make lmn vectors
ic_mva = 2;
c_eval('[out?,l?,v?] = irf_minvar(gseB?.tlim(irf.tint(''2015-12-02T01:14:50.030Z/2015-12-02T01:15:01.473Z'')));',ic_mva)
% GSE -> LMN
%L = [1 0 0];
%M = [0 1 0];
%N = [0 0 1];
%lmn = [L M N];
c_eval('lmn = v?;',ic_mva)

c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',ic)
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)
c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';',ic)
c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';',ic)

mvaR0 = (mvaR1.resample(mvaR1.time)+mvaR2.resample(mvaR1.time)+mvaR3.resample(mvaR1.time)+mvaR4.resample(mvaR1.time))/4;
c_eval('mvaRR? = mvaR?-mvaR0;',ic)

mvaB0 = (mvaB1.resample(mvaB1.time)+mvaB2.resample(mvaB1.time)+mvaB3.resample(mvaB1.time)+mvaB4.resample(mvaB1.time))/4;

%% Vector plot
ic = 2;      
tintQuivers = irf.tint('2015-12-02T01:14:54.00Z/2015-12-02T01:15:00.00Z');
tintQuivers = irf.tint('2015-12-02T01:14:55.00Z/2015-12-02T01:14:59.00Z');
c_eval('tsPlot = mvaVe?.tlim(tintQuivers);',ic)
c_eval('tsB = mvaB?.resample(tsPlot).tlim(tintQuivers);',ic)
c_eval('tsB? = mvaB?.resample(tsPlot).tlim(tintQuivers);',0)
tsB0 = tsB0.filt(0,1,[],5);
times = tsPlot.time;

% Get path xyz/lmn
vL = 300; % km/s, just took something
vM = 0;
vN = [];

posL = vL*(times-times(1)); % km
posM = (times-times(1))*0;

% normal distance is decided by BL
mpNtot = 100; % km, estimated thickness of magnetopause (I jsu took something), get this from timing or something
BL = tsB0.x.data;
dBtot = 60;
posN = BL/(dBtot/2)*mpNtot;
tsPos = irf.ts_vec_xyz(tsPlot.time,[posL posM posN]);

c_eval('posR? = mvaRR?.resample(times).data+tsPos.data;')
%c_eval('posV?perp = mvaVe?perp.resample(times).data;')
c_eval('posV? = mvaVe?.resample(times).data;')
c_eval('posB? = mvaB?.resample(times).data;')
c_eval('posB? = tsB?.resample(times).data;',0)

% Plot one quantity and time series
ic_plot = 1:4;
nrows = 2; ncols = 1; npanels = nrows*ncols; isub = 1;
for irow = 1:nrows
  for icol = 1:ncols
    h(irow,icol) = subplot(nrows,ncols,isub); isub = isub + 1;
  end
end

colors = mms_colors('matlab');

isub = 1;
if 0 % B
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  c_eval('plot(hca,posL,posB?(:,!),''color'',colors(!,:));',ic_plot,1:3)
  c_eval('hN = plot(hca,posL,posB?(:,!),''color'',[0 0 0]);',0,1)
  %c_eval('plot(hca,posL,posB?(:,!),''color'',mms_colors(''?''),''linestyle'',''.-'')',ic_plot,1:3)
  hold(hca,'off')
  %irf_zoom(hca,'x',tintQuivers)
  legend(hN,'used to get N position')
  
end
if 1 % V
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  c_eval('plot(hca,posL,posV?(:,!),''color'',colors(!,:))',ic_plot,1:3)
  %c_eval('plot(hca,posL,posB?(:,!),''color'',mms_colors(''?''),''linestyle'',''.-'')',ic_plot,1:3)
  hold(hca,'off')
  %irf_zoom(hca,'x',tintQuivers)
  hca.XLim = posL([1 end]);
end
if 1 % B arrows
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  S = 1; % arrow scale
  %c_eval('plot_quivers(hca,[posB?(:,1) posB?(:,2) posB?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  c_eval('quiver3(hca,posR?(:,1),posR?(:,2),posR?(:,3),posB?(:,1),posB?(:,2),posB?(:,3),S,''color'',mms_colors(''?''))',1:4)
  c_eval('hs = scatter3(hca,posR?(:,1),posR?(:,2),posR?(:,3),abs(posV?(:,2))*0.1,''MarkerEdgeColor'',mms_colors(''?''));',1:4)
  hold(hca,'off')
  hca.XLabel.String = 'L';
  hca.YLabel.String = 'M';
  hca.ZLabel.String = 'N';
  %axis(hca,'equal')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';  
  view(hca,[0 -1 0])
  %hca.ZLim = [-100 100];
  hca.ZDir = 'reverse';
  hca.XLim = posL([1 end]);
end

if 0 % V arrows 
  hca = h(isub); isub = isub + 1;
  hold(hca,'on')
  S = 1;
  %c_eval('plot_quivers(hca,[posV?(:,1) posV?(:,2) posV?(:,3)],[posR?(:,1) posR?(:,2) posR?(:,3)],mms_colors(''?''))')
  c_eval('quiver3(hca,posR?(:,1),posR?(:,2),posR?(:,3),posV?(:,1),posV?(:,2),posV?(:,3),S,''color'',mms_colors(''?''))',1:4)
  hold(hca,'off')
  hca.XLabel.String = 'L';
  hca.YLabel.String = 'M';
  hca.ZLabel.String = 'N';
  %axis(hca,'equal')
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.ZGrid = 'on';  
  view(hca,[0 -1 0])
  %hca.ZLim = [-100 100];
  hca.ZDir = 'reverse';
  hca.XLim = posL([1 end]);
end
