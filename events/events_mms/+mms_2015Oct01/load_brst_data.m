%tint = irf.tint('2015-10-16T10:30:00.00Z/2015-10-16T10:30:40.00Z');
tint = irf.tint('2015-10-01T06:53:32.000Z/2015-10-01T06:53:48.000Z'); % the one with so different fields
ic = 1;
%% 
% Electric field
c_eval('dslE?brst=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('dmpaB?brst=mms.db_get_ts(''mms?_dfg_brst_ql'',''mms?_dfg_brst_dmpa'',tint);',1);

c_eval('P?brst=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',tint);',ic);

% Magnetic field (searchcoil)
%c_eval('load(''/Users/Cecilia/Data/MMS/2015Oct01/scmdata/mms?_scb_cleanup16_32_20151001_064900.mat''); B?sc = Bscm;  B?sc = irf.ts_vec_xyz(B?sc.time,B?sc.data);',ic);
c_eval('load(''/Users/Cecilia/Data/MMS/2015Oct01/scmdata/mms?_scb_cleanup16_32_20151001_065200.mat''); B?sc = Bscm;  B?sc = irf.ts_vec_xyz(B?sc.time,B?sc.data);',ic);

% Electron moments
c_eval('ne?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_numberDensity'',tint);',ic);

c_eval('vex?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkX'',tint);',ic);
c_eval('vey?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkY'',tint);',ic);
c_eval('vez?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_bulkZ'',tint);',ic);
c_eval('ve?brst=irf.ts_vec_xyz(vex?brst.time,[vex?brst.data vey?brst.data vez?brst.data]);',ic)

c_eval('Pexx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresXX'',tint);',ic);
c_eval('Peyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresYY'',tint);',ic);
c_eval('Pezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_PresZZ'',tint);',ic);
c_eval('Pe?brst=irf.ts_vec_xyz(Pexx?brst.time,[Pexx?brst.data Peyy?brst.data Pezz?brst.data]);',ic)

c_eval('Texx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempXX'',tint);',ic);
c_eval('Teyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempYY'',tint);',ic);
c_eval('Tezz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_des-moms'',''mms?_des_TempZZ'',tint);',ic);
c_eval('Te?brst=irf.ts_vec_xyz(Texx?brst.time,[Texx?brst.data Teyy?brst.data Tezz?brst.data]);',ic)

% Ion moments
c_eval('ni?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_numberDensity'',tint);',ic);

c_eval('vix?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkX'',tint);',ic);
c_eval('viy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkY'',tint);',ic);
c_eval('viz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_bulkZ'',tint);',ic);
c_eval('vi?brst=irf.ts_vec_xyz(vix?brst.time,[vix?brst.data viy?brst.data viz?brst.data]);',ic)

c_eval('Pixx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_PresXX'',tint);',ic);
c_eval('Piyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_PresYY'',tint);',ic);
c_eval('Pizz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_PresZZ'',tint);',ic);
c_eval('Pi?brst=irf.ts_vec_xyz(Pixx?brst.time,[Pixx?brst.data Piyy?brst.data Pizz?brst.data]);',ic)

c_eval('Tixx?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_TempXX'',tint);',ic);
c_eval('Tiyy?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_TempYY'',tint);',ic);
c_eval('Tizz?brst=mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_TempZZ'',tint);',ic);
c_eval('Ti?brst=irf.ts_vec_xyz(Tixx?brst.time,[Tixx?brst.data Tiyy?brst.data Tizz?brst.data]);',ic)

%% c_eval('ni?brst = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-moms'',''mms?_dis_numberDensity'',tint);',ic);
c_eval('desDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint);',ic);
c_eval('disDist? = mms.db_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint);',ic);
%c_eval('desDist? = mms.variable2ts(mms.db_get_variable(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint));',ic);
%%
%tmpDataObj = dataobj('/Volumes/SAMSUNG/data/mms1/fpi/brst/l1b/des-dist/2015/10/16/mms1_fpi_brst_l1b_des-dist_20151016103000_v0.2.0.cdf');
%dist = mms.variable2ts(get_variable(tmpDataObj,'mms1_des_brstSkyMap_dist'));
%disttemp = dist.tlim(tint);

%% Construct variables
% ExB velocity
c_eval('VExB? = cross(dslE?brst.resample(dmpaB?),dmpaB?)/dmpaB?.abs/dmpaB?.abs;',ic)
%%
e = 1.6022e-19;
c_eval('je? = -e*ne?brst*ve?brst*1e3*1e6*1e9; je?.units = ''nA/m^2''',ic);
c_eval('ji? = e*ni?brst*vi?brst*1e3*1e6*1e9; ji?.units = ''nA/m^2''',ic);
% Assuming GSE and DMPA are the same coordinate system.
%[j,divB,B,jxB,divTshear,divPb] = c_4_j('gseR?','dmpaB?');
%j.data = j.data*1e9; j.units = 'nAm^{-2}';