%% Load some survey data
sc=1:4;
c_eval('[scm?,dobj_scm?] = mms.cn_get_ts(''mms?_scm_fast_l2_scf'',''mms?_scm_scf_gse'',tint(1));',sc);
c_eval('[edp?,dobj_edp?] = mms.cn_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint(1));',sc);
c_eval('[P?,dobj_P?] = mms.cn_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_scpot'',tint(1));',sc);

c_eval('[dfg?,dobj_dfg?] = mms.cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint(1));',sc)
%% Particle moments
c_eval('[vix?fast,dobj_vi?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_X_DSC'',tint(1));',sc)
c_eval('[viy?fast,dobj_vi?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Y_DSC'',tint(1));',sc)
c_eval('[viz?fast,dobj_vi?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Z_DSC'',tint(1));',sc)
c_eval('[vex?fast,dobj_ve?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_X_DSC'',tint(1));',sc)
c_eval('[vey?fast,dobj_ve?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Y_DSC'',tint(1));',sc)
c_eval('[vez?fast,dobj_ve?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Z_DSC'',tint(1));',sc)
c_eval('ve?fast = irf.ts_vec_xyz(vey?fast.time,[vex?fast.data vey?fast.data vez?fast.data]);',sc)
c_eval('ve?fast.units = vez?fast.units;',sc)
c_eval('ve?fast.userData = vez?fast.userData;',sc)
c_eval('vi?fast = irf.ts_vec_xyz(viy?fast.time,[vix?fast.data viy?fast.data viz?fast.data]);',sc)
c_eval('vi?fast.units = viz?fast.units;',sc)
c_eval('vi?fast.userData = viz?fast.userData;',sc)
c_eval('[ni?fast,dobj_ni?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',tint(1));',sc)
c_eval('[ne?fast,dobj_ne?fast] = mms.cn_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint(1));',sc)


%% Burst data
% Electric field burst data
c_eval('[dslE?,dobjE] = mms.cn_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint(1));',sc)
c_eval('[P?,dobjP?] = cn_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',tint(1));',sc);
% Ion skymap
c_eval('[disDist?,dobjDisDist] = mms.cn_get_ts(''mms?_fpi_brst_l1b_dis-dist'',''mms?_dis_brstSkyMap_dist'',tint(1));',sc)

% Electron skymap
c_eval('[desDist?,dobjDesDist] = mms.cn_get_ts(''mms?_fpi_brst_l1b_des-dist'',''mms?_des_brstSkyMap_dist'',tint(1));',sc)

%% Electron moments
c_eval('[~,dobj] = cn_get_ts(''mms?_fpi_brst_l1b_des-moms'',[],tint(1));',sc)

c_eval('ne? = mms.variable2ts(get_variable(dobj,''mms?_des_numberDensity''));',sc)

c_eval('peXX = mms.variable2ts(get_variable(dobj,''mms?_des_PresXX'')); ',sc)
c_eval('peYY = mms.variable2ts(get_variable(dobj,''mms?_des_PresYY'')); ',sc)
c_eval('peZZ = mms.variable2ts(get_variable(dobj,''mms?_des_PresZZ'')); ',sc)
c_eval('pe? = irf.ts_scalar(peXX.time,(peXX.data + peYY.data + peZZ.data)/3);',sc)
c_eval('pe?.units = peXX.units;',sc)
c_eval('pe?.userData = peXX.userData;',sc)

c_eval('TeXX = mms.variable2ts(get_variable(dobj,''mms?_des_TempXX''));',sc)
c_eval('TeYY = mms.variable2ts(get_variable(dobj,''mms?_des_TempYY''));',sc)
c_eval('TeZZ = mms.variable2ts(get_variable(dobj,''mms?_des_TempZZ''));',sc)
c_eval('Te? = irf.ts_scalar(TeXX.time,(TeXX.data + TeYY.data + TeZZ.data)/3);',sc)
c_eval('Te?.units = TeXX.units;',sc)
c_eval('Te?.userData = TeXX.userData;',sc)

c_eval('veX = mms.variable2ts(get_variable(dobj,''mms?_des_bulkX''));',sc)
c_eval('veY = mms.variable2ts(get_variable(dobj,''mms?_des_bulkY''));',sc)
c_eval('veZ = mms.variable2ts(get_variable(dobj,''mms?_des_bulkZ''));',sc)
c_eval('ve? = irf.ts_vec_xyz(veX.time,[veX.data veY.data veZ.data]);',sc)
c_eval('ve?.units = veX.units;',sc)
c_eval('ve?.userData = veX.userData;',sc)

% Downsample electron moments
c_eval('fs = 1/(ne?.time(2)-ne?.time(1));',sc)
fny = fs/2;
c_eval('pe?_lowres = irf_filt(pe?,0,fny/2,fs,5);',sc)
c_eval('Te?_lowres = irf_filt(Te?,0,fny/2,fs,5);',sc)
c_eval('ne?_lowres = irf_filt(ne?,0,fny/2,fs,5);',sc)

% Ion moments
c_eval('[~,dobj] = cn_get_ts(''mms?_fpi_brst_l1b_dis-moms'',[],tint(1));',sc)

c_eval('ni? = mms.variable2ts(get_variable(dobj,''mms?_dis_numberDensity''));',sc)

c_eval('piXX = mms.variable2ts(get_variable(dobj,''mms?_dis_PresXX'')); ',sc)
c_eval('piYY = mms.variable2ts(get_variable(dobj,''mms?_dis_PresYY'')); ',sc)
c_eval('piZZ = mms.variable2ts(get_variable(dobj,''mms?_dis_PresZZ'')); ',sc)
c_eval('pi? = irf.ts_scalar(piXX.time,(piXX.data + piYY.data + piZZ.data)/3);',sc)
c_eval('pi?.units = piXX.units;',sc)
c_eval('pi?.userData = piXX.userData;',sc)

c_eval('TiXX = mms.variable2ts(get_variable(dobj,''mms?_dis_TempXX''));',sc)
c_eval('TiYY = mms.variable2ts(get_variable(dobj,''mms?_dis_TempYY''));',sc)
c_eval('TiZZ = mms.variable2ts(get_variable(dobj,''mms?_dis_TempZZ''));',sc)
c_eval('Ti? = irf.ts_scalar(TiXX.time,(TiXX.data + TiYY.data + TiZZ.data)/3);',sc)
c_eval('Ti?.units = TiXX.units;',sc)
c_eval('Ti?.userData = TiXX.userData;',sc)

c_eval('viX = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkX''));',sc)
c_eval('viY = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkY''));',sc)
c_eval('viZ = mms.variable2ts(get_variable(dobj,''mms?_dis_bulkZ''));',sc)
c_eval('vi? = irf.ts_vec_xyz(viX.time,[viX.data viY.data viZ.data]);',sc)
c_eval('vi?.units = viX.units;',sc)
c_eval('vi?.userData = viX.userData;',sc)

% Downsample electron moments
c_eval('fs = 1/(ni?.time(2)-ni?.time(1));',sc)
fny = fs/2;
c_eval('pi?_lowres = irf_filt(pi?,0,fny/2,fs,5);',sc)
c_eval('Ti?_lowres = irf_filt(Ti?,0,fny/2,fs,5);',sc)
c_eval('ni?_lowres = irf_filt(ni?,0,fny/2,fs,5);',sc)

