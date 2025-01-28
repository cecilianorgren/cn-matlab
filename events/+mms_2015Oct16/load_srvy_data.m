% Overview time interval
units = irf_units;
e = 1.6022e-19;
tintOverview = irf.tint('2015-10-16T05:30:00.00Z/2015-10-16T16:00:00.00Z');
tint = tintOverview;
ic = 1:4;

% OMNI database solar wind data
tintOMNI = [irf_time([2015 10 16 03 30 00.00]) irf_time([2015 10 16 16 00 00.00])]; 
omni_data2015Oct16 = irf_get_data_omni(tintOMNI,'n,V,Bx,By,Bz','omni_min');
tOMNI = omni_data2015Oct16(:,1); tOMNI = irf_time(tOMNI,'epoch>epochtt');
nOMNI = irf.ts_scalar(tOMNI,omni_data2015Oct16(:,2)); % cc
vOMNI = irf.ts_scalar(tOMNI,omni_data2015Oct16(:,3)); % km/s
BzOMNI = irf.ts_scalar(tOMNI,[omni_data2015Oct16(:,4) omni_data2015Oct16(:,5) omni_data2015Oct16(:,6)]); % nT

disp('Loading B, E, R and scPot...')
% Magnetic srvy field
c_eval('dmpaB?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);

% Spacecraft position
R  = mms.get_data('R_gse',tint);
c_eval('gseR? = TSeries(R.time,R.gseR?,''to'',1);',ic);

% Spacecraft potential
c_eval('P? = mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_scpot'',tint);',ic);

% Fast electric field data
c_eval('dslE?=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);

disp('Loading particle moments...')
% SITL electron moments
c_eval('ne?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);
c_eval('ni?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',tint);',ic);

c_eval('vex?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_X_DSC'',tint);',ic);
c_eval('vey?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Y_DSC'',tint);',ic);
c_eval('vez?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Z_DSC'',tint);',ic);
c_eval('vix?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_X_DSC'',tint);',ic);
c_eval('viy?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Y_DSC'',tint);',ic);
c_eval('viz?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Z_DSC'',tint);',ic);
c_eval('Tepar?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DEStempPara'',tint);',ic);
c_eval('Teper?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DEStempPerp'',tint);',ic);
c_eval('Tipar?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DIStempPara'',tint);',ic);
c_eval('Tiper?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DIStempPerp'',tint);',ic);
c_eval('Pexx?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_XX_DSC'',tint);',ic);
c_eval('Peyy?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_YY_DSC'',tint);',ic);
c_eval('Pezz?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_ZZ_DSC'',tint);',ic);
c_eval('j_fpi_x?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_X_DSC'',tint);',ic);
c_eval('j_fpi_y?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Y_DSC'',tint);',ic);
c_eval('j_fpi_z?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Z_DSC'',tint);',ic);

c_eval('ve?=irf.ts_vec_xyz(vex?.time,[vex?.data vey?.data vez?.data]);',ic)
c_eval('vi?=irf.ts_vec_xyz(vix?.time,[vix?.data viy?.data viz?.data]);',ic)
c_eval('Pe?=irf.ts_vec_xyz(Pexx?.time,[Pexx?.data Peyy?.data Pezz?.data]);',ic)
c_eval('j_fpi?=irf.ts_vec_xyz(j_fpi_x?.time,[j_fpi_x?.data,j_fpi_y?.data,j_fpi_z?.data]);',ic)
c_eval('Te?=irf.ts_vec_xy(Tepar?.time,[Tepar?.data Teper?.data]);',ic)
c_eval('clear vex? vey? vez? vix? viy? viz? Pexx? Peyy? Pezz? j_fpi_x? j_fpi_y? j_fpi_z?')

%% SITL omni differential energy flux
disp('Loading differential energy flux...')
c_eval('eEnSp_pX?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pX'',tint);',ic);%positiveX
c_eval('eEnSp_pY?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pY'',tint);',ic);
c_eval('eEnSp_pZ?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_pZ'',tint);',ic);
c_eval('eEnSp_mX?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mX'',tint);',ic);%negativeX
c_eval('eEnSp_mY?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mY'',tint);',ic);
c_eval('eEnSp_mZ?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_eEnergySpectr_mZ'',tint);',ic);

[~,energy] = hist([log10(10),log10(30e3)],32); energy = 10.^energy;
c_eval(['eEnSp?.t = irf_time(eEnSp_pX?.DEPEND_0.data,''ttns>epoch'');',...
        'eEnSp?.f = single(energy);',...
        'eEnSp?.p = log10(double((eEnSp_pX?.data(:,1:32)+eEnSp_mX?.data(:,1:32) + eEnSp_mY?.data(:,1:32)+eEnSp_mY?.data(:,1:32) + eEnSp_mZ?.data(:,1:32)+eEnSp_mZ?.data(:,1:32))/6));',...
        'eEnSp?.f_label=''E_e [eV]'';',...
        'eEnSp?.p_label={''log(dEF)''};'],ic)

c_eval('iEnSp_pX?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pX'',tint);',ic);%positiveX
c_eval('iEnSp_pY?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pY'',tint);',ic);
c_eval('iEnSp_pZ?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_pZ'',tint);',ic);
c_eval('iEnSp_mX?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mX'',tint);',ic);%negativeX
c_eval('iEnSp_mY?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mY'',tint);',ic);
c_eval('iEnSp_mZ?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_iEnergySpectr_mZ'',tint);',ic);

[~,energy] = hist([log10(10),log10(30e3)],32); energy = 10.^energy;
c_eval(['iEnSp?.t = irf_time(iEnSp_pX?.DEPEND_0.data,''ttns>epoch'');',...
        'iEnSp?.f = single(energy);',...
        'iEnSp?.p = log10(double((iEnSp_pX?.data(:,1:32)+iEnSp_mX?.data(:,1:32) + iEnSp_mY?.data(:,1:32)+iEnSp_mY?.data(:,1:32) + iEnSp_mZ?.data(:,1:32)+iEnSp_mZ?.data(:,1:32))/6));',...
        'iEnSp?.f_label=''E_i [eV]'';',...
        'iEnSp?.p_label={''log(dEF)''};'],ic)

%% Constructed variables
disp('Loading constructed variables...')
% ExB velocity
c_eval('vExB? = cross(dslE?.resample(dmpaB?),dmpaB?)/dmpaB?.abs/dmpaB?.abs*1e3;',ic) % km/s

% Current: assuming GSE and DMPA are the same coordinate system.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('gseR?','dmpaB?');
j = irf.ts_vec_xyz(j.time,j.data*1e9);
j.units = 'nAm^{-2}';

% Make hall electric field
c_eval('EJxB? = j.resample(dmpaB?).cross(dmpaB?)*1e-9*1e-9/e/(ni?*1e6)*1e3; EJxB?.units = ''mV/m'';',ic)

% Speeds
c_eval('vte? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Vte''); vte? = irf.ts_scalar(dmpaB?.time,vte?);',ic)
c_eval('vtp? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Vtp''); vtp? = irf.ts_scalar(dmpaB?.time,vtp?);',ic)
c_eval('vA? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Va''); vA? = irf.ts_scalar(dmpaB?.time,vA?);',ic)

% Frequencies
c_eval('flh? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Flh''); flh? = irf.ts_scalar(dmpaB?.time,flh?);',ic)
c_eval('fce? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Fce''); fce? = irf.ts_scalar(dmpaB?.time,fce?);',ic)
c_eval('fcp? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Fcp''); fcp? = irf.ts_scalar(dmpaB?.time,fcp?);',ic)
c_eval('fpe? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Fpe''); fpe? = irf.ts_scalar(dmpaB?.time,fpe?);',ic)
c_eval('fpp? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Fpp''); fpp? = irf.ts_scalar(dmpaB?.time,fpp?);',ic)

% Length scales
c_eval('Lp? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Li''); Lp? = irf.ts_scalar(dmpaB?.time,Lp?);',ic)
c_eval('Le? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Le''); Le? = irf.ts_scalar(dmpaB?.time,Le?);',ic)
c_eval('Ld? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Ld''); Ld? = irf.ts_scalar(dmpaB?.time,Ld?);',ic)
c_eval('re? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Roe''); re? = irf.ts_scalar(dmpaB?.time,re?);',ic)
c_eval('rp? = irf_plasma_calc(dmpaB?.abs.data,ne?.resample(dmpaB?.time).data,0,Teper?.resample(dmpaB?.time).data,Tiper?.resample(dmpaB?.time).data,''Rop''); rp? = irf.ts_scalar(dmpaB?.time,rp?);',ic)


