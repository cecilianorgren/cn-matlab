% hwang2021_overview
%% Load data
% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   
localuser = datastore('local','user');
units = irf_units;

% Specify time intervals for waves
tint_waves{1} = irf.tint('2017-05-05T20:09:47.900Z/2017-05-05T20:09:47.960Z');
tint_waves{2} = irf.tint('2017-05-05T20:09:49.736Z/2017-05-05T20:09:49.746Z');
tint_waves{3} = irf.tint('2017-05-05T20:09:50.820Z/2017-05-05T20:09:50.839Z');

% Specify time interval for file
ic = 1:4;
tint_search = irf.tint('2017-05-05T20:09:20.00Z/2017-05-05T20:10:00.00Z');
files = mms.db_list_files('mms1_edp_brst_l2_scpot',tint_search);
tint = [files.start files.stop];

% Magnetic field
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);

% Electric field
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);

% Spacecraft potential
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);
c_eval('tic; psp?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_psp_brst_l2'',tint); toc;',ic);

c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);

c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);

%% Currents from moments, use ne also for Ji 
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
%c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
%c_eval('gseJ? = (gseJe?+gseJi?);',ic);

c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
%c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
%c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
%c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
%c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)


c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);

% Specify paths for figures and matlab files
fileName = gseB1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); 
eventPath = ['/Users/' localuser '/Research/Events/Hwang2021_probe_timing/'];  
matlabPath = ['/Users/' localuser '/MATLAB/cn-matlab/mms' fileNameSplit{5} '/'];

%% Plot probe data
ic = 1;
npanels = 5;
h = irf_plot(npanels);

hca = irf_panel('B');
c_eval('irf_plot(hca,gseB?);',ic)
irf_legend(hca,{'x','y','z'},[0.98 0.98]);
hca.YLabel.String = {'V GSE','(nT)'};

hca = irf_panel('E');
c_eval('irf_plot(hca,gseE?);',ic)
irf_legend(hca,{'x','y','z'},[0.98 0.98]);
hca.YLabel.String = {'E GSE','(mV/m)'};

hca = irf_panel('Vi');
c_eval('irf_plot(hca,gseVi?);',ic)
irf_legend(hca,{'x','y','z'},[0.98 0.98]);
hca.YLabel.String = {'Vi GSE','(km/s)'};

hca = irf_panel('Ve');
c_eval('irf_plot(hca,gseVe?);',ic)
irf_legend(hca,{'x','y','z'},[0.98 0.98]);
hca.YLabel.String = {'Ve GSE','(km/s)'};

hca = irf_panel('J');
c_eval('irf_plot(hca,gseJ?);',ic)
irf_legend(hca,{'x','y','z'},[0.98 0.98]);
hca.YLabel.String = {'J GSE','(nA/m2)'};

c_eval('irf_pl_mark(h(?),tint_waves{!},[0 0 0]);',1:npanels,1:numel(tint_waves))
tint_zoom = [tint_waves{1}(1) tint_waves{end}(2)] + [-2 2];
irf_zoom(h,'x',tint_zoom)

%% Check probe align times for spin plane probes
%mms.probe_align_times
%mms.sc_orient
ic = 1;
c_eval('Exyz = dslE?;',ic)
c_eval('Bxyz = dmpaB?;',ic)
c_eval('SCpot = dcv?;',ic)
c_eval('zphase = zphase?;',ic)
[starttime1,endtime1,starttime3,endtime3] = mms.probe_align_times(Exyz,Bxyz,SCpot,zphase,1);
h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('irf_pl_mark(h(?),tint_waves{!},[1 0 0]);',1:numel(h),1:numel(tint_waves))

%% Calculate frequency power spectrum and wave number
probecomb = 5;
trange = tint_waves{3};
c_eval('V6 = dcv?;',ic)
c_eval('Bxyz = dmpaB?;',ic)
c_eval('zphase = zphase?;',ic)
[fkpower,freq,wavenumber] = mms.fk_powerspectrum(probecomb,trange,V6,Bxyz);

hca = subplot(1,1,1);
pcolor(hca,wavenumber,freq,fkpower);
shading(hca,'flat');
hca.XLabel.String = 'k (m^{-1})';
hca.YLabel.String = 'f (Hz)';
hcb = colorbar('peer',hca);
hcb.YLabel.String = 'Power/max(Power)';
hca.Title.String = [trange(1).utc ' - ' trange(2).utc];
hca.XGrid = 'on';
hca.YGrid = 'on';

if 1 % tint 3
  hold(hca,'on')
  hca.YLim = [0 120];
  xfit = [0 0.2];
  yfit = [0 40];
  vph = yfit(2)*2*pi/xfit(2);
  plot(hca,xfit,yfit,'linewidth',1)
  hold(hca,'off')
end
if 0 % tint 1
  hold(hca,'on')
  hca.YLim = [0 200];
  xfit = [0 0.2];
  yfit = [0 70];
  vph = yfit(2)*2*pi/xfit(2);
  plot(hca,xfit,yfit,'linewidth',1)
  hold(hca,'off')
end