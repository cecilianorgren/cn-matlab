% Load some example time interval
%savePath = '/Users/cecilia/Research/energy_partitioning';
localuser = datastore('local','user');
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   
irf.log('critical')
units = irf_units;

tint_all = irf.tint('2015-01-01T00:00:00.00Z/2018-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fgm_brst_l2',tint_all);

iFile = 2;
tint = [files(iFile).start files(iFile).stop];

%% Load data
ic = 1;
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);


tint_cs = irf.tint('2017-07-03T05:26:41.00Z/2017-07-03T05:26:44.00Z');

%results = ph_dev_energy_partition(ePDist1.tlim(tint_cs),dmpaB1.tlim(tint_cs),dbcsVe1.tlim(tint_cs),facTe1.tlim(tint_cs));