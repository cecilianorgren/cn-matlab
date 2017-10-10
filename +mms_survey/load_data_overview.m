%% Load datastore
db_info = datastore('mms_db');   

    
%% Magnetic field
disp('Loading magnetic field...')
%c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc, dmpaB? = dmpaB?.clone(dmpaB?.time,dmpaB?.data(:,1:3));',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc, gseB? = gseB?.clone(gseB?.time,gseB?.data(:,1:3));',ic);

%% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);

%% Load spacecraft position
disp('Loading spacecraft position...')
tic
R = mms.get_data('R_gse',tint);
c_eval('gseR? = TSeries(R.time,R.gseR?,''to'',1);',1:4);
toc

%% Particle moments
% load dataobject first
% electrons
c_eval('edobj?=dataobj(mms.get_filepath(''mms?_fpi_brst_l2_des-moms'',tint(1)));',ic)

c_eval('ne? = mms.variable2ts(get_variable(edobj?,''mms?_des_numberdensity_dbcs_brst''));',ic);

c_eval('vex? = mms.variable2ts(get_variable(edobj?,''mms?_des_bulkx_dbcs_brst''));',ic);
c_eval('vey? = mms.variable2ts(get_variable(edobj?,''mms?_des_bulky_dbcs_brst''));',ic);
c_eval('vez? = mms.variable2ts(get_variable(edobj?,''mms?_des_bulkz_dbcs_brst''));',ic);
c_eval('dbcsVe?=irf.ts_vec_xyz(vex?.time,[vex?.data vey?.data vez?.data]);',ic)
c_eval('dbcsVe1.name = ''ve''; dbcsVe?.userData = vex?.userData;',ic)

c_eval('Pexx? = mms.variable2ts(get_variable(edobj?,''mms?_des_presxx_dbcs_brst''));',ic);
c_eval('Pexy? = mms.variable2ts(get_variable(edobj?,''mms?_des_presxy_dbcs_brst''));',ic);
c_eval('Pexz? = mms.variable2ts(get_variable(edobj?,''mms?_des_presxz_dbcs_brst''));',ic);  
c_eval('Peyy? = mms.variable2ts(get_variable(edobj?,''mms?_des_presyy_dbcs_brst''));',ic);
c_eval('Peyz? = mms.variable2ts(get_variable(edobj?,''mms?_des_presyz_dbcs_brst''));',ic);  
c_eval('Pezz? = mms.variable2ts(get_variable(edobj?,''mms?_des_preszz_dbcs_brst''));',ic);   
c_eval('DataPe = zeros(Pexx?.length,3,3); DataPe(:,1,1) = Pexx?.data; DataPe(:,1,2) = Pexy?.data; DataPe(:,1,3) = Pexz?.data; DataPe(:,2,2) = Peyy?.data; DataPe(:,2,3) = Peyz?.data; DataPe(:,3,3) = Pezz?.data;',ic);
c_eval('dbcsPe?=irf.ts_tensor_xyz(Pexx?.time,DataPe);',ic)
c_eval('dbcsPe?.name = ''Pe''; dbcsPe?.units = ''nPa''; dbcsPe?.userData = Pexx?.userData; dbcsPe?.coordinateSystem = ''DBCS'';',ic)

c_eval('Texx? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempxx_dbcs_brst''));',ic);
c_eval('Texy? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempxy_dbcs_brst''));',ic);
c_eval('Texz? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempxz_dbcs_brst''));',ic);  
c_eval('Teyy? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempyy_dbcs_brst''));',ic);
c_eval('Teyz? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempyz_dbcs_brst''));',ic);  
c_eval('Tezz? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempzz_dbcs_brst''));',ic);   
c_eval('DataTe = zeros(Texx?.length,3,3); DataTe(:,1,1) = Texx?.data; DataTe(:,1,2) = Texy?.data; DataTe(:,1,3) = Texz?.data; DataTe(:,2,2) = Teyy?.data; DataTe(:,2,3) = Teyz?.data; DataTe(:,3,3) = Tezz?.data;',ic);
c_eval('dbcsTe?=irf.ts_tensor_xyz(Texx?.time,DataTe);',ic)
c_eval('dbcsTe?.name = ''Te''; dbcsTe?.units = ''eV''; dbcsTe?.userData = Texx?.userData; dbcsTe?.coordinateSystem = ''DBCS'';',ic)

% ions
c_eval('idobj?=dataobj(mms.get_filepath(''mms?_fpi_brst_l2_dis-moms'',tint(1)));',ic)

c_eval('ni? = mms.variable2ts(get_variable(idobj?,''mms?_dis_numberdensity_dbcs_brst''));',ic);

c_eval('vix? = mms.variable2ts(get_variable(idobj?,''mms?_dis_bulkx_dbcs_brst''));',ic);
c_eval('viy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_bulky_dbcs_brst''));',ic);
c_eval('viz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_bulkz_dbcs_brst''));',ic);
c_eval('dbcsVi?=irf.ts_vec_xyz(vix?.time,[vix?.data viy?.data viz?.data]);',ic)
c_eval('dbcsVi?.name = ''vi''; dbcsVi?.userData = vix?.userData; dbcsVi?.coordinateSystem = ''DBCS'';',ic)

c_eval('Pixx? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presxx_dbcs_brst''));',ic);
c_eval('Pixy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presxy_dbcs_brst''));',ic);
c_eval('Pixz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presxz_dbcs_brst''));',ic);  
c_eval('Piyy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presyy_dbcs_brst''));',ic);
c_eval('Piyz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presyz_dbcs_brst''));',ic);  
c_eval('Pizz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_preszz_dbcs_brst''));',ic);   
c_eval('DataPi = zeros(Pixx?.length,3,3); DataPi(:,1,1) = Pixx?.data; DataPi(:,1,2) = Pixy?.data; DataPi(:,1,3) = Pixz?.data; DataPi(:,2,2) = Piyy?.data; DataPi(:,2,3) = Piyz?.data; DataPi(:,3,3) = Pizz?.data;',ic);
c_eval('dbcsPi?=irf.ts_tensor_xyz(Pixx?.time,DataPi);',ic)
c_eval('dbcsPi?.name = ''Pi''; dbcsPi?.units = ''nPa''; dbcsPi?.userData = Pixx?.userData; dbcsPi?.coordinateSystem = ''DBCS'';',ic)

c_eval('Tixx? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempxx_dbcs_brst''));',ic);
c_eval('Tixy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempxy_dbcs_brst''));',ic);
c_eval('Tixz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempxz_dbcs_brst''));',ic);  
c_eval('Tiyy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempyy_dbcs_brst''));',ic);
c_eval('Tiyz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempyz_dbcs_brst''));',ic);  
c_eval('Tizz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempzz_dbcs_brst''));',ic);   
c_eval('DataTi = zeros(Tixx?.length,3,3); DataTi(:,1,1) = Tixx?.data; DataTi(:,1,2) = Tixy?.data; DataTi(:,1,3) = Tixz?.data; DataTi(:,2,2) = Tiyy?.data; DataTi(:,2,3) = Tiyz?.data; DataTi(:,3,3) = Tizz?.data;',ic);
c_eval('dbcsTi?=irf.ts_tensor_xyz(Tixx?.time,DataTi);',ic)
c_eval('dbcsTi?.name = ''Ti''; dbcsTi?.units = ''eV''; dbcsTi?.userData = Tixx?.userData; dbcsTi?.coordinateSystem = ''DBCS'';',ic)

% Rotate into GSE coordinates
disp('Rotate into GSE...'); tic
c_eval('gsePe? = mms.rotate_tensor(dbcsPe?,''gse'',?); gsePe?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
c_eval('gsePi? = mms.rotate_tensor(dbcsPi?,''gse'',?); gsePi?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseTe? = mms.rotate_tensor(dbcsTe?,''gse'',?);',ic)
c_eval('gseTi? = mms.rotate_tensor(dbcsTi?,''gse'',?);',ic)
c_eval('gseVe? = mms_dsl2gse(dbcsVe?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseVi? = mms_dsl2gse(dbcsVi?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic); toc

% Rotate into FAC coordinates
disp('Rotate into FAC...'); tic
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?); facTe?.coordinateSystem = ''FAC'';',ic)
c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?); facTi?.coordinateSystem = ''FAC'';',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.coordinateSystem = ''FAC'';',ic); toc      

%% Calculate curlometer current
if 0
  %%
disp('Calculating ExB velocity, speeds, currents, frequencies, length scales...')
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
end

%% Calculate currents from moments
units = irf_units;
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);
try
  gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; 
end
%% Calculate gradient of pressure
if 0
try
  gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4);
  gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4);
end
end
%% Calculate some physical quantities
% ExB velocity
c_eval('gseVExB? = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.abs.resample(gseE?.time)/gseB?.abs.resample(gseE?.time)*1e3;',ic) % km/s
 
% Prepare input
if 0
  c_eval('matB? = gseB?.abs.data,ne?.resample(gseB?.time).data;',ic)
  c_eval('matParTe? = facTe?.xx.resample(gseB?.time).data;',ic)
  c_eval('matParTi? = facTi?.xx.resample(gseB?.time).data;',ic)
  c_eval('matPerTe? = (facTe?.yy.resample(gseB?.time).data + facTe?.zz.resample(gseB?.time).data)/2;',ic)
  c_eval('matPerTi? = (facTi?.yy.resample(gseB?.time).data + facTi?.zz.resample(gseB?.time).data)/2;',ic)
  c_eval('matTe? = facTe?.trace.resample(gseB?.time).data/3;',ic)
  c_eval('matTi? = facTi?.trace.resample(gseB?.time).data/3;',ic)
  c_eval('matNe? = ne?.resample(gseB?.time).data;',ic)

  % Speeds
  c_eval('vte?perp = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matTi?,''Vte''); vte?perp = irf.ts_scalar(gseB?.time,vte?perp)*1e-3; vte?.units = ''km/s'';',ic)
  c_eval('vte?par = irf_plasma_calc(matB?,matNe?,0,matParTe?,matTi?,''Vte''); vte?par = irf.ts_scalar(gseB?.time,vte?par)*1e-3; vte?.units = ''km/s'';',ic)
  c_eval('vte? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Vte''); vte? = irf.ts_scalar(gseB?.time,vte?)*1e-3; vte?.units = ''km/s'';',ic)
  c_eval('vtp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matPerTi?,''Vtp''); vtp? = irf.ts_scalar(gseB?.time,vtp?)*1e-3; vtp?.units = ''km/s'';',ic)
  c_eval('vA? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Va''); vA? = irf.ts_scalar(gseB?.time,vA?)*1e-3; vA?.units = ''km/s'';',ic)

  % Frequencies
  c_eval('flh? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Flh''); flh? = irf.ts_scalar(gseB?.time,flh?);',ic)
  c_eval('fce? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fce''); fce? = irf.ts_scalar(gseB?.time,fce?);',ic)
  c_eval('fcp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fcp''); fcp? = irf.ts_scalar(gseB?.time,fcp?);',ic)
  c_eval('fpe? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpe''); fpe? = irf.ts_scalar(gseB?.time,fpe?);',ic)
  c_eval('fpp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Fpp''); fpp? = irf.ts_scalar(gseB?.time,fpp?);',ic)

  % Length scales
  c_eval('Lp? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Li''); Lp? = irf.ts_scalar(gseB?.time,Lp?)*1e-3; Lp?.units = ''km'';',ic)
  c_eval('Le? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Le''); Le? = irf.ts_scalar(gseB?.time,Le?)*1e-3; Le?.units = ''km'';',ic)
  c_eval('Ld? = irf_plasma_calc(matB?,matNe?,0,matTe?,matTi?,''Ld''); Ld? = irf.ts_scalar(gseB?.time,Ld?)*1e-3; Ld?.units = ''km'';',ic)
  c_eval('re? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Roe''); re? = irf.ts_scalar(gseB?.time,re?)*1e-3; re?.units = ''km'';',ic)
  c_eval('rp? = irf_plasma_calc(matB?,matNe?,0,matPerTe?,matPerTi?,''Rop''); rp? = irf.ts_scalar(gseB?.time,rp?)*1e-3; rp?.units = ''km'';',ic)

  % Dimensionless
  c_eval('beta? = (re?/Le?).^2;',ic)
  c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
end