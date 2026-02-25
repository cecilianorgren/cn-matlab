% Overview time interval
units = irf_units;
e = 1.6022e-19;
tint = irf.tint('2015-10-16T05:30:00.00Z/2015-10-16T16:00:00.00Z');
ic = 1:4;

% Load srvy B field
c_eval('Bxyz?=mms.db_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint);',ic);
c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',[2:4]);
c_eval('Bxyz? = Bxyz?.tlim(tint);',ic)
c_eval('Bmag? = Bxyz?.abs.data;',ic)
c_eval('Bvec? = Bxyz?.data./([Bmag? Bmag? Bmag?]);',ic)
c_eval('Bvec? = TSeries(Bxyz?.time,Bvec?,''to'',1);',ic)
Bxyzav = (Bxyz1.data+Bxyz2.data+Bxyz3.data+Bxyz4.data)/4;
Bxyzav = TSeries(Bxyz1.time,Bxyzav,'to',1);

% Spacecraft position from magnetic field product
R  = mms.get_data('R_gse',tint);
c_eval('Rxyz? = TSeries(R.time,R.gseR?,''to'',1);',ic);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',ic);

% Load fast electric field data
c_eval('Exyz?=mms.db_get_ts(''mms?_edp_fast_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint);',ic);
c_eval('Exyz? = Exyz?.resample(Exyz1);',[2:4]);
c_eval('Exyz?.data(find(abs(Exyz?.data) > 100)) = NaN;',[1:4]); %Remove some questionable fields
Exyzav = (Exyz1.data+Exyz2.data+Exyz3.data+Exyz4.data)/4;
Exyzav = TSeries(Exyz1.time,Exyzav,'to',1);

% Calculate ExB velocity
c_eval('Exyz? = Exyz?.resample(Bxyz?);',ic)
c_eval('VExB? = cross(Exyz?,Bxyz?);',ic)
c_eval('VExB?.data = VExB?.data./[Bmag?.^2 Bmag?.^2 Bmag?.^2];',ic)
c_eval('VExBr? = VExB?; VExBr?.data = VExBr?.data*1e3;',ic)
c_eval('VExBav? = VExB?.abs.data*1e6;',ic)
c_eval('VExBav? = 0.5*1.673e-27*VExBav?.^2/1.6e-19;',ic)
c_eval('VExBav? = TSeries(VExB?.time,VExBav?,''to'',0);',ic)
%c_eval('VExBav? = mms.fftbandpass(VExBav?,0,0.5);',ic)

% Calculate current
R0 = (Rxyz1.data+Rxyz2.data+Rxyz3.data+Rxyz4.data)/4;
%c_eval('relR? = Rxyz?-R0;',sc)
% Assuming GSE and DMPA are the same coordinate system.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');
j.data = j.data*1e9; j.units = 'nAm^{-2}';




% load sitl electron moments
c_eval('ne?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);
c_eval('ne?=irf.ts_vec_xyz(ne?.time,[ne?.data ne?.data ne?.data]);',ic);
c_eval('ne? = ne?.resample(Bxyz1);',ic)
c_eval('ne?=irf.ts_scalar(ne?.time,ne?.data(:,1));',ic);
c_eval('ni?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',tint);',ic);
c_eval('ni?=irf.ts_vec_xyz(ni?.time,[ni?.data ni?.data ni?.data]);',ic);
c_eval('ni? = ni?.resample(Bxyz1);',ic)
c_eval('ni?=irf.ts_scalar(ni?.time,ni?.data(:,1));',ic);


% Make hall electric field
c_eval('Ehall? = j.cross(Bxyz?)*1e-9*1e-9/e/(ni?*1e6)*1e3; Ehall?.units = ''mV/m'';',ic)


% Spacecraft potential
c_eval('P? = mms.db_get_ts(''mms?_edp_fast_l2_scpot'',''mms?_edp_scpot'',tint);',ic);
c_eval('mP? = P?; mP?.data = -mP?.data;',ic)
%% SITL omni differential energy flux
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
%%
%fpi data
ic = 1:4;

c_eval('ni?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DISnumberDensity'',tint);',ic);
%c_eval('ne=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESnumberDensity'',tint);',ic);
c_eval('vex?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_X_DSC'',tint);',ic);
c_eval('vey?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Y_DSC'',tint);',ic);
c_eval('vez?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_eBulkV_Z_DSC'',tint);',ic);
c_eval('vix?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_X_DSC'',tint);',ic);
c_eval('viy?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Y_DSC'',tint);',ic);
c_eval('viz?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_iBulkV_Z_DSC'',tint);',ic);
c_eval('tepar?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DEStempPara'',tint);',ic);
c_eval('teper?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DEStempPerp'',tint);',ic);
c_eval('Pexx?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_XX_DSC'',tint);',ic);
c_eval('Peyy?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_YY_DSC'',tint);',ic);
c_eval('Pezz?=mms.db_get_ts(''mms?_fpi_fast_sitl'',''mms?_fpi_DESpress_ZZ_DSC'',tint);',ic);
c_eval('j_fpi_x?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_X_DSC'',tint);',ic);
c_eval('j_fpi_y?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Y_DSC'',tint);',ic);
c_eval('j_fpi_z?=mms.db_get_variable(''mms?_fpi_fast_sitl'',''mms?_fpi_Jplasma_Z_DSC'',tint);',ic);
%end

c_eval('ve?=irf.ts_vec_xyz(vex?.time,[vex?.data vey?.data vez?.data]);',ic)
c_eval('vi?=irf.ts_vec_xyz(vix?.time,[vix?.data viy?.data viz?.data]);',ic)
c_eval('Pev?=irf.ts_vec_xyz(Pexx?.time,[Pexx?.data Peyy?.data Pezz?.data]);',ic)
c_eval('j_fpi?=irf.ts_vec_xyz(j_fpi_x?.DEPEND_0.data,[j_fpi_x?.data,j_fpi_y?.data,j_fpi_z?.data]);',ic)
%c_eval('t_e?=irf.ts_vec_xy(tepar?.DEPEND_0.data,[tepar?.data teper?.data]);',ic)

%%
