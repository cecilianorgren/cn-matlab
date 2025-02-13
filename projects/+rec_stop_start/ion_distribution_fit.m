%% Specify time and spacecraft
units = irf_units;
irf.log('critical')
ic = 1;

mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

% Time from time interval
%tint = irf.tint('2017-06-22T03:01:03.00Z/2017-06-22T03:01:43.00Z');

% Time from file name
fileId = '20170709171633'; %tint = [files(iFile-1).start files(iFile).stop] + [1 -1]; % also load the file before
%fileId = '20170709173053';
fileId = '20170611173913';

%fileId = '20170724125513';
fileId = '20170725220853';
%fileId = '20170726065803';
%fileId = '20170728200013'; %two front with 2 density decrease

iFile = find(cellfun(@(s) contains(s,fileId),{files.name}));
tint = [files(iFile-1).start files(iFile).stop] + [1 -1];



% Event path
eventPath = ['/Users/' localuser '/Research/Events/mms_' fileId '/']; % for saving figures
%matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/mms_' fileID '/'];
mkdir(eventPath)

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gsmVe? = c_coord_trans(''GSE'',''GSM'',gseVe?);',ic)
c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?);',ic)

c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
% Remove all one-count "noise"
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?.data(iPDist?.data < iPDistErr?.data*1.01) = 0;',ic)


disp('Preparing reduced distributions.')
vint = [-Inf Inf];
% Cutting off the edge of the distribution allows the maxwellian to go 
% beyond the edge. Otherwise, the obs red distribution is always forced to
% zero due to energy cutoff and spherical skymap.
vg = -1900:50:1900; 
c_eval('if1Dx? = iPDist?.reduce(''1D'',[1 0 0],''vint'',vint,''vg'',vg);',ic)
c_eval('if1Dy? = iPDist?.reduce(''1D'',[0 1 0],''vint'',vint,''vg'',vg);',ic)
c_eval('if1Dz? = iPDist?.reduce(''1D'',[0 0 1],''vint'',vint,''vg'',vg);',ic)

%% Calculate moments from vdf

%% Do fit
v = dbcsVi1.y;
ne = ne2;
ni = ni2;
vdf = if1Dy1(501:520);
vdf = if1Dy1(700:10:2000);
vdf = if1Dx1(250:20:500);

% particlemoments = mms.psd_moments(vdf,scPot2,'energyrange',[0 34000]);
X0 = [0.01e6 -100e3 100e3 0.01e6 1000e3 1000e3]; 
nPop = 2;
tic
% fitdata_nobg = fitdata; ts_nobg = ts;
[fitdata,ts] = funFitVDF(vdf,'nPop',nPop,'plot',1,'guessprevious',0,'X0',X0,'weight',[1 1 1]);
units = irf_units; ts.vt = irf.ts_scalar(ts.T.time,units.s*(1-1./(ts.T.data*units.e/(units.mp*units.c^2)+1).^2).^0.5*1e-3);
units = irf_units; ts.vt = irf.ts_scalar(ts.T.time,(2*ts.T.data*units.e/units.mp).^0.5*1e-3);
toc
vfitmaxtot = irf.ts_scalar(ts.n.time,sum(ts.n.data.*ts.vd.data,2)./sum(ts.n.data,2));
%% Plot results from funFitVDF()
h = irf_plot(8);
hca = irf_panel('vdf_obs');
irf_spectrogram(hca,vdf.specrec('velocity'),'log'); colormap(hca,irf_colormap('waterfall')) 
hca.YLabel.String = 'v^{obs} (km/s)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('vdf_fit');
irf_spectrogram(hca,ts.f.specrec('velocity'),'log'); colormap(hca,irf_colormap('waterfall'))
hold(hca,'on'); irf_plot(hca,ts.vd); hold(hca,'off'); 
hca.YLabel.String = 'v^{fit} (km/s)';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('n_fit');
irf_plot(hca,ts.n)
%hold(hca,'on'); hneobs = irf_plot(hca,ne/10); hold(hca,'off'); 
%hold(hca,'on'); hniobs = irf_plot(hca,ni/10); hold(hca,'off'); 
%irf_legend(hca,'ne obs',[0.98 0.98],'color',hneobs.Color)
%irf_legend(hca,'ni obs',[0.98 0.02],'color',hniobs.Color)
if 0
hold(hca,'on'); hniobs = irf_plot(hca,particlemoments.n_psd); hold(hca,'off'); 
irf_legend(hca,'ni obs mms.psdmoments',[0.98 0.98],'color',hniobs.Color)
hniobs.LineWidth = 0.5;
hniobs.LineStyle = ':';
end
hca.YLabel.String = 'n (cm^{-3})';
hca.YLabel.Interpreter = 'tex';

hca = irf_panel('vd_fit');
irf_plot(hca,ts.vd)
hold(hca,'on'); hvobs = irf_plot(hca,v); hold(hca,'off'); 
hold(hca,'on'); hvfittot = irf_plot(hca,vfitmaxtot); hold(hca,'off'); 
irf_legend(hca,{'obs'},[0.98 0.02],'color',hvobs.Color)
irf_legend(hca,{'sum(fit)'},[0.98 0.98],'color',hvfittot.Color)
hca.YLabel.String = 'v (km/s)';

hca = irf_panel('T_fit');
irf_plot(hca,ts.T*1e-3)
%hold(hca,'on'); hTobs = irf_plot(hca,T0.xx); hold(hca,'off'); 
%irf_legend(hca,'obs',[0.98 0.98],'color',hTobs.Color)
hca.YLabel.String = 'T (keV)';

hca = irf_panel('vt_fit');
irf_plot(hca,ts.vt)
%hold(hca,'on'); hTobs = irf_plot(hca,T0.xx); hold(hca,'off'); 
%irf_legend(hca,'obs',[0.98 0.98],'color',hTobs.Color)
hca.YLabel.String = 'v_t (km/s)';
hca.YLabel.Interpreter = 'tex';

if 1
  hca = irf_panel('vt_fit/vd_fit');
  irf_plot(hca,irf.ts_scalar(ts.vt.time,ts.vt.data./ts.vd.data))
  %hold(hca,'on'); hTobs = irf_plot(hca,T0.xx); hold(hca,'off'); 
  %irf_legend(hca,'obs',[0.98 0.98],'color',hTobs.Color)
  hca.YLabel.String = 'v_t/v_d';
  hca.YLabel.Interpreter = 'tex';
  
end

hca = irf_panel('cf_fit');
irf_plot(hca,ts.cf)

irf_plot_axis_align(h)
hlinks = linkprop([irf_panel('vdf_obs'),irf_panel('vdf_fit')],{'CLim','YLim'});
irf_zoom(h,'x',vdf.time([1 end]))



