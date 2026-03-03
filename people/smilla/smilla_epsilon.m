%% Load data
mms.db_init('local_file_db','/Volumes/mms');
tint = irf.tint('2015-10-16T10:32:30.00Z/2015-10-16T10:34:10.00Z');
ic = 3;
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('dbcsTi? = mms.get_data(''Ti_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

%% Make model PDist
modelPDist = mms.make_model_dist(iPDist3,...
                                 dmpaB3.resample(iPDist3), ...
                                 scPot3.resample(iPDist3), ...
                                 ni3.resample(iPDist3), ...
                                 dbcsVi3.resample(iPDist3), ...
                                 dbcsTi3.resample(iPDist3));
modelPDist.data = modelPDist.data*modelPDist.siConversion/str2num(iPDist3.siConversion); % s^3/km^6 -> s^3/cm^6
modelPDist.units = iPDist3.units;

%% Calculate epsilon
epsilon = mms.calculate_epsilon(iPDist3,modelPDist,ni3,scPot3,'enchannels',1:32);
epsilon.name = 'epsilon';
%

%% Plot
h = irf_plot(3);

hca = irf_panel('obs dist');
irf_spectrogram(hca,iPDist3.omni.deflux.specrec);

hca = irf_panel('mod dist');
irf_spectrogram(hca,modelPDist.omni.deflux.specrec);

hca = irf_panel('epsilon');
irf_plot({epsilon});

linkprop(h(1:2),{'CLim'})

irf_plot_axis_align(h)
irf_zoom(h,'x',[iPDist3.time.start, iPDist3.time.stop])