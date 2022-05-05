ic = 1;
tint_fast = irf.tint('2017-07-25T20:00:00.00Z/2017-07-25T24:00:00.00Z');

% Load datastore
localuser = datastore('local','user');
mms.db_init('local_file_db',['/Users/' localuser '/data']);
db_info = datastore('mms_db');   

%%

c_eval('gseB?_srvy = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',tint_fast);',ic);
c_eval('iPDist?_fast = mms.get_data(''PDi_fpi_fast_l2'',tint_fast,?);',ic) % missing some ancillary data
c_eval('gseVi?_fast = mms.get_data(''Vi_gse_fpi_fast_l2'',tint_fast,?);',ic);
c_eval('gseE?_fast = mms.get_data(''E_gse_edp_fast_l2'',tint_fast);',ic);

c_eval('gseVExB?_srvy = gseE?_srvy.resample(gseB?_srvy).cross(gseB?_srvy);',ic)
c_eval('gseVExB?_srvy = cross(gseE?_srvy.resample(gseB?_srvy.time),gseB?)/gseB?_srvy.abs2*1e3; gseVExB?_srvy.units = '''';',ic) % km/s

c_eval('nOp?_srvy = mms.get_data(''Noplus_hpca_srvy_l2'',tint_fast,?);',ic);
c_eval('nHp?_srvy = mms.get_data(''Nhplus_hpca_srvy_l2'',tint_fast,?);',ic);

c_eval('iPDist?_Op_fast = mms.get_data(''Omnifluxoplus_hpca_srvy_l2'',tint_fast,?);',ic) % missing some ancillary data 