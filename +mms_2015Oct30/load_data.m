%% Load data
tint = irf.tint('2015-10-30T05:15:30.00Z/2015-10-30T05:16:20.00Z');
%tint = irf.tint('2015-10-30T05:15:44.00Z/2015-10-30T05:15:48.00Z');
ic = 1:4;

%% Load defatt, for coordinate tranformation
disp('Loading defatt...')
c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic);
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
  
%% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_dmpa'',tint); toc',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',tint); toc',ic);
%%
scmFilePath = '/Users/Cecilia/Data/MMS/2015Oct30/scm/';
scmFileList = dir(scmFilePath); scmFileList = scmFileList(3:end);
for ii = 1:numel(scmFileList)  
  load([scmFilePath scmFileList(ii).name]);
  c_eval('dmpaB?scm = Bscm;',ii)
  c_eval('gseB?scm = mms_dsl2gse(dmpaB?scm,defatt?);',ii)
end
%% Electric field
disp('Loading electric field...')
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
c_eval('dslE? = dslE?*1.2;',ic)
c_eval('gseE? = mms_dsl2gse(dslE?,defatt?);',ic)

%% Electron moments
disp('Loading electron moments...')
c_eval('ne? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',tint);',ic);

c_eval('vex?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkx_dbcs_brst'',tint);',ic);
c_eval('vey?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulky_dbcs_brst'',tint);',ic);
c_eval('vez?=mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_bulkz_dbcs_brst'',tint);',ic);
c_eval('dslVe?=irf.ts_vec_xyz(vex?.time,[vex?.data vey?.data vez?.data]);',ic)
c_eval('gseVe? = mms_dsl2gse(dslVe?,defatt?);',ic)
c_eval('gseVe?.name = ''mms? ve brst'';',ic)

c_eval('Pe? = mms.get_data(''Pe_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePe? = mms.rotate_tensor(Pe?,''gse'',?);',ic)
c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?);',ic)

c_eval('Te? = mms.get_data(''Te_fpi_brst_l2'',tint,?);',ic)
c_eval('gseTe? = mms.rotate_tensor(Te?,''gse'',?);',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)

%% Ion moments
disp('Loading ion moments...')
c_eval('ni? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_dbcs_brst'',tint);',ic);

%% Physical variables

% ExB velocity
c_eval('vExB?brst = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.abs.resample(gseE?.time)/gseB?.abs.resample(gseE?.time)*1e3;',ic) % km/s

% Frequencies
c_eval('Te? = gseTe?.resample(gseB?.time).trace.data/3;',ic)
c_eval('Ti? = Te?;',ic)
c_eval('flh? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Flh''); flh? = irf.ts_scalar(gseB?.time,flh?);',ic)
c_eval('fce? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Fce''); fce? = irf.ts_scalar(gseB?.time,fce?);',ic)
c_eval('fcp? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Fcp''); fcp? = irf.ts_scalar(gseB?.time,fcp?);',ic)
%c_eval('fpe? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Fpe''); fpe?brst = irf.ts_scalar(gseB?.time,fpe?brst);',ic)
%c_eval('fpp? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Fpp''); fpp?brst = irf.ts_scalar(gseB?.time,fpp?brst);',ic)

% Length scales
%c_eval('Lp? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Li''); Lp?brst = irf.ts_scalar(gseB?.time,Lp?brst);',ic)
%c_eval('Le? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Le''); Le?brst = irf.ts_scalar(gseB?.time,Le?brst);',ic)
%c_eval('Ld? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Ld''); Ld?brst = irf.ts_scalar(gseB?.time,Ld?brst);',ic)
c_eval('re? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Roe''); re? = irf.ts_scalar(gseB?.time,re?);',ic)
%c_eval('rp? = irf_plasma_calc(gseB?.abs.data,ne?.resample(gseB?.time).data,0,Te?,Ti?,''Rop''); rp?brst = irf.ts_scalar(gseB?.time,rp?brst);',ic)

