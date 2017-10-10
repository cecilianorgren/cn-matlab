tint = irf.tint('2015-11-30T00:21:44.00Z/2015-11-30T00:27:00.00Z'); % magnetosphere-magnetosheath-magnetosphere

ic = 1:4;

%% Load defatt, for coordinate tranformation
%disp('Loading defatt...')
%load /Users/Cecilia/Data/MMS/2015Oct16/defatt.mat

%% Magnetic field
disp('Loading magnetic field...')
%c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_dmpa'',tint); toc',ic);
%c_eval('tic; gseB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',tint); toc',ic);

%c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc, dmpaB? = dmpaB?.clone(dmpaB?.time,dmpaB?.data(:,1:3));',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc, gseB? = gseB?.clone(gseB?.time,gseB?.data(:,1:3));',ic);
%c_eval('tic; gsmB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc, gsmB? = gsmB?.clone(gsmB?.time,gsmB?.data(:,1:3));',ic);

c_eval('tic; gseB?scm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint); toc',ic);
%% Electric field
disp('Loading electric field...')
%Eoffs1 = 2.05; Eoffs2 = 2.58; Eoffs3 = 2.53; Eoffs4 = 1.40;
%c_eval('tic; dslE?dce2d=mms.db_get_ts(''mms?_edp_brst_ql_dce2d'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
%c_eval('tic; dslE?dce=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
%c_eval('dslE?dce2d.data(:,1)=dslE?dce2d.data(:,1)-Eoffs?;',ic)
%c_eval('dslE?dce.data(:,1)=dslE?dce.data(:,1)-Eoffs?;',ic)
%c_eval('dslE?=dslE?dce;',ic)
%c_eval('[dslEdB0?,d]=irf_edb(dslE1brst,dmpaB1,20);')
%c_eval('tic; dslE?l2pre=mms.db_get_ts(''mms?_edp_fast_l2a_dce2d'',''mms?_edp_dce_slow_l2a'',tint); toc',ic);

%c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l1b_hmfe'',''mms?_edp_hmfe_sensor'',tint); toc',ic);

% L2 data
%c_eval('tic; dslE?fast=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_dsl_fast_l2'',tint); toc',ic);
%c_eval('gseE?fast = mms_dsl2gse(dslE?fast,defatt?);',ic)

%% Load spacecraft position, not brst
disp('Loading spacecraft position...')
R  = mms.get_data('R_gse',tint);
c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',ic);

%% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);

%% Magnetic field (searchcoil)
%c_eval('load(''/Users/Cecilia/Data/MMS/2015Oct01/scmdata/mms?_scb_cleanup16_32_20151001_064900.mat''); B?sc = Bscm;  B?sc = irf.ts_vec_xyz(B?sc.time,B?sc.data);',ic);
%c_eval('load(''/Users/Cecilia/Data/MMS/2015Oct01/scmdata/mms?_scb_cleanup16_32_20151001_065200.mat''); B?sc = Bscm;  B?sc = irf.ts_vec_xyz(B?sc.time,B?sc.data);',ic);

%% Electron and ion data
% Particle distributions
disp('Loading particle distributions...')
if 1  
  db_info = datastore('mms_db');  
  c_eval('[ePDist? ePDistError?] = mms.make_pdist([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/10/16/mms?_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf'']);',ic)
  c_eval('[iPDist? iPDistError?] = mms.make_pdist([db_info.local_file_db_root ''/mms?/fpi/brst/l2/dis-dist/2015/10/16/mms?_fpi_brst_l2_dis-dist_20151016103254_v2.1.0.cdf'']);',ic)
else
  db_info = datastore('mms_db');   
  disp('Loading electron distribution...')
  c_eval('etmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l2/des-dist/2015/10/16/mms?_fpi_brst_l2_des-dist_20151016103254_v2.1.0.cdf'']);',ic);
  c_eval('desDist? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_dist_brst''));',ic);
  c_eval('eenergy0? = get_variable(etmpDataObj?,''mms?_des_energy0_brst'');',ic);
  c_eval('eenergy1? = get_variable(etmpDataObj?,''mms?_des_energy1_brst'');',ic);
  c_eval('ephi? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_phi_brst''));',ic);
  c_eval('etheta? = get_variable(etmpDataObj?,''mms?_des_theta_brst'');',ic);
  c_eval('estepTable? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_steptable_parity_brst''));',ic);
  c_eval('eEpoch_plus_var? = get_variable(etmpDataObj?,''Epoch_plus_var'');',ic);
  c_eval('eEpoch_minus_var? = get_variable(etmpDataObj?,''Epoch_minus_var'');',ic);

  % Make energytable from energy0, energy1 and energysteptable
  c_eval('energy? = repmat(torow(eenergy0?.data),numel(estepTable?.data),1);',ic)
  c_eval('energy?(estepTable?.data==1,:) = repmat(eenergy1?.data,sum(estepTable?.data),1);',ic)
  c_eval('dt_shift? = 0.5*(double(eEpoch_plus_var?.data)-double(eEpoch_minus_var?.data))*1e-3;',ic)
  c_eval('dt_minus? = double(eEpoch_minus_var?.data)*1e-3-dt_shift?;',ic)
  c_eval('dt_plus? = double(eEpoch_plus_var?.data)*1e-3-dt_shift?;',ic)
  c_eval('ePDist? = PDist(desDist?.time+dt_shift?,desDist?.data,''skymap'',energy?,ephi?.data,etheta?.data);',ic)
  c_eval('ePDist?.userData = desDist?.userData; ePDist?.name = desDist?.name; ePDist?.units = desDist?.units;',ic)
  c_eval('ePDist?.units = ''s^3/m^6''; ePDist?.species = ''electrons'';',ic)
  c_eval('ePDist?.ancillary.dt_minus = dt_minus?; ePDist?.ancillary.dt_plus = dt_plus?;',ic) 

  disp('Loading ion distribution...')
  c_eval('etmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l2/dis-dist/2015/10/16/mms?_fpi_brst_l2_dis-dist_20151016103254_v2.1.0.cdf'']);',ic);
  c_eval('disDist? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_dis_dist_brst''));',ic);
  c_eval('eenergy0? = get_variable(etmpDataObj?,''mms?_dis_energy0_brst'');',ic);
  c_eval('eenergy1? = get_variable(etmpDataObj?,''mms?_dis_energy1_brst'');',ic);
  c_eval('ephi? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_dis_phi_brst''));',ic);
  c_eval('etheta? = get_variable(etmpDataObj?,''mms?_dis_theta_brst'');',ic);
  c_eval('estepTable? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_dis_steptable_parity_brst''));',ic);
  c_eval('eEpoch_plus_var? = get_variable(etmpDataObj?,''Epoch_plus_var'');',ic);
  c_eval('eEpoch_minus_var? = get_variable(etmpDataObj?,''Epoch_minus_var'');',ic);

  % Make energytable from energy0, energy1 and energysteptable
  c_eval('energy? = repmat(torow(eenergy0?.data),numel(estepTable?.data),1);',ic)
  c_eval('energy?(estepTable?.data==1,:) = repmat(eenergy1?.data,sum(estepTable?.data),1);',ic)
  c_eval('dt_shift? = 0.5*(double(eEpoch_plus_var?.data)-double(eEpoch_minus_var?.data))*1e-3;',ic)
  c_eval('dt_minus? = double(eEpoch_minus_var?.data)*1e-3-dt_shift?;',ic)
  c_eval('dt_plus? = double(eEpoch_plus_var?.data)*1e-3-dt_shift?;',ic)

  c_eval('iPDist? = PDist(disDist?.time+dt_shift?,disDist?.data,''skymap'',energy?,ephi?.data,etheta?.data);',ic)
  c_eval('iPDist?.userData = disDist?.userData; iPDist?.name = disDist?.name; iPDist?.units = disDist?.units;',ic)
  c_eval('iPDist?.units = ''s^3/m^6''; iPDist?.species = ''ions'';',ic)
  c_eval('iPDist?.ancillary.dt_minus = dt_minus?; iPDist?.ancillary.dt_plus = dt_plus?;',ic) 
end
%
%c_eval('tic; desDist? = mms.db_get_ts(''mms?_fpi_brst_l2_des-dist'',''mms?_des_dist_brst'',tint); toc',ic);
%c_eval('tic; disDist? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-dist'',''mms?_dis_dist_brst'',tint); toc',ic);

% Pressure and temperature
disp('Loading pressure and temperature...')
%i_dt_shift = 0.075;
%e_dt_shift = 0.015;
c_eval('tic; dbcsPe? = mms.get_data(''Pe_fpi_brst_l2'',tint,?); toc;',ic) 
c_eval('tic; dbcsTe? = mms.get_data(''Te_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsPi? = mms.get_data(''Pi_fpi_brst_l2'',tint,?); toc;',ic) 
c_eval('tic; dbcsTi? = mms.get_data(''Ti_fpi_brst_l2'',tint,?); toc;',ic)

c_eval('gsePe? = mms.rotate_tensor(dbcsPe?,''gse'',?); gsePe?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseTe? = mms.rotate_tensor(dbcsTe?,''gse'',?);',ic)
c_eval('gsePi? = mms.rotate_tensor(dbcsPi?,''gse'',?); gsePi?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseTi? = mms.rotate_tensor(dbcsTi?,''gse'',?);',ic)

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...')
c_eval('tic; ne? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',tint); toc;',ic);
c_eval('tic; ni? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_dbcs_brst'',tint); toc;',ic);

% Velocity
disp('Loading bulk velocities...')
c_eval('tic; dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?); toc;',ic)

c_eval('gseVe? = mms_dsl2gse(dbcsVe?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseVi? = mms_dsl2gse(dbcsVi?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic)

if 0 % Load distribution help data
  db_info = datastore('mms_db');   
  % electrons
  disp('Loading electron distribution data and making moments...')
  c_eval('etmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l1b/des-dist/2015/10/16/mms?_fpi_brst_l1b_des-dist_20151016103254_v1.1.0.cdf'']);',ic);
  c_eval('desDist? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_brstSkyMap_dist''));',ic);
  c_eval('eenergy0? = get_variable(etmpDataObj?,''mms?_des_brstSkyMap_energy0'');',ic);
  c_eval('eenergy1? = get_variable(etmpDataObj?,''mms?_des_brstSkyMap_energy1'');',ic);
  c_eval('ephi? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_brstSkyMap_phi''));',ic);
  c_eval('etheta? = get_variable(etmpDataObj?,''mms?_des_brstSkyMap_theta'');',ic);
  c_eval('estepTable? = mms.variable2ts(get_variable(etmpDataObj?,''mms?_des_stepTable_parity''));',ic);
 

  %% ions
  disp('Loading ion distribution data and making moments...')
  c_eval('itmpDataObj? = dataobj([db_info.local_file_db_root ''/mms?/fpi/brst/l1b/dis-dist/2015/10/16/mms?_fpi_brst_l1b_dis-dist_20151016103254_v1.1.0.cdf'']);',ic);
  c_eval('disDist? = mms.variable2ts(get_variable(itmpDataObj?,''mms?_dis_brstSkyMap_dist''));',ic);
  c_eval('ienergy0? = get_variable(itmpDataObj?,''mms?_dis_brstSkyMap_energy0'');',ic);
  c_eval('ienergy1? = get_variable(itmpDataObj?,''mms?_dis_brstSkyMap_energy1'');',ic);
  c_eval('iphi? = mms.variable2ts(get_variable(itmpDataObj?,''mms?_dis_brstSkyMap_phi''));',ic);
  c_eval('itheta? = get_variable(itmpDataObj?,''mms?_dis_brstSkyMap_theta'');',ic);
  c_eval('istepTable? = mms.variable2ts(get_variable(itmpDataObj?,''mms?_dis_stepTable_parity''));',ic);
end  

%% Calculate curlometer current
disp('Calculating ExB velocity, speeds, currents, frequencies, length scales...')
c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';

%% Calculate currents from moments
units = irf_units;
c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
c_eval('gseJ? = (gseJe?+gseJi?);',ic);
gseAvJ = (gseJ1+gseJ2.resample(gseJ1.time)+gseJ3.resample(gseJ1.time)+gseJ4.resample(gseJ1.time))/4; 

%% Calculate gradient of pressure
gseGradPe = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePe1,gsePe2,gsePe3,gsePe4);
gseGradPi = mms_2015Oct16.gradP(gseR1,gseR2,gseR3,gseR4,gsePi1,gsePi2,gsePi3,gsePi4);

%% Calculate convective electric fields
c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)

%% Wave analysis
c_eval('wavE? = irf_wavelet(gseE?.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavE?.f_units = ''Hz''; wavE?.f_label = ''f [Hz]''; wavE?.p_label = {''log_{10} E^2'',''(mV/m)^2/Hz''};',ic)
c_eval('wavB? = irf_wavelet(gseB?scm.abs.tlim(tint),''wavelet_width'',5.36*2,''f'',[1 4000],''nf'',100);',ic)
c_eval('wavB?.f_units = ''nT''; wavB?.f_label = ''f [Hz]''; wavB?.p_label = {''log_{10} B^2'',''(nT)^2/Hz''};',ic)

tintPol = irf.tint('2015-10-16T10:33:40.00Z/2015-10-16T10:33:52.00Z');
c_eval('polarization? = irf_ebsp(gseE?.tlim(tintPol),gseB?scm.tlim(tintPol),gseB?.tlim(tintPol),gseB?.tlim(tintPol),gseR?brsttime.tlim(tintPol),[10 2200],''polarization'',''fac'');',1)

%% Calculate some additional parameters
% ExB velocity
c_eval('gseVExB? = cross(gseE?,gseB?.resample(gseE?.time))/gseB?.abs.resample(gseE?.time)/gseB?.abs.resample(gseE?.time)*1e3;',ic) % km/s

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?);',ic)
 
% Speeds
c_eval('matB? = gseB?.abs.data,ne?.resample(gseB?.time).data;',ic)
c_eval('matParTe? = facTe?.xx.resample(gseB?.time).data;',ic)
c_eval('matParTi? = facTi?.xx.resample(gseB?.time).data;',ic)
c_eval('matPerTe? = (facTe?.yy.resample(gseB?.time).data + facTe?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matPerTi? = (facTi?.yy.resample(gseB?.time).data + facTi?.zz.resample(gseB?.time).data)/2;',ic)
c_eval('matTe? = facTe?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matTi? = facTi?.trace.resample(gseB?.time).data/3;',ic)
c_eval('matNe? = ne?.resample(gseB?.time).data;',ic)

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

c_eval('beta? = (re?/Le?).^2;',ic)
c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)