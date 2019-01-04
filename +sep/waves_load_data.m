ic = 1:4;

%% Load data
units = irf_units;
if 0 % Particle distributions: electrons and ions
  disp('Loading particle distributions...')
  c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
  c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
  % Setting tint from ePDist1
  %c_eval('tint = ePDist?([1 ePDist?.length]).time;',ic)
end    
if 1 % Magnetic field
  disp('Loading magnetic field...')
  c_eval('dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',1:4);
  c_eval('gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
  c_eval('gseB?scm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint);',1:4);
end
if 1 % Electric field
  disp('Loading electric field...')
  c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',1:4);
  c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
end
if 1 % Load spacecraft position
  disp('Loading spacecraft position...')
  try
    R = mms.get_data('R_gse',tint);
    if ~all([isfield(R,'gseR1') isfield(R,'gseR1') isfield(R,'gseR2') isfield(R,'gseR3') isfield(R,'gseR4')])  
      % not the right data fielsd, try to load from irfu database instead
      db_info = datastore('mms_db');   
      mms.db_init('local_file_db','/data/mms');
      R = mms.get_data('R_gse',tint);
      mms.db_init('local_file_db',db_info.local_file_db_root);
    end
    if size(R.gseR1,2) == 4
      c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
    else
      c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
    end
  end
end
if 0 % Spacecraft potential
  disp('Loading spacecraft potential...')
  c_eval('scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
end
if 0 % Pressure and temperature
  disp('Loading pressure and temperature...')
  c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
  c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic)

  c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
  c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
  c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)
end
if 0 % Density
  disp('Loading density...')
  c_eval('tic; ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?); toc;',ic);
  c_eval('tic; ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?); toc;',ic);
end
if 0 % Velocity
  disp('Loading bulk velocities...')
  c_eval('tic; gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?); toc;',ic)
  c_eval('tic; gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?); toc;',ic)
  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
end

%% Prepeare data
if 0 % Current: use electron density to calculate current
  c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',ic);
  c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',ic);
  c_eval('gseJ? = (gseJe?+gseJi?);',ic);
end
if 0 % Current, curlometer
  if all([~isempty(gseB1) ~isempty(gseB2) ~isempty(gseB3) ~isempty(gseB4)])
    c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
    [Jcurl,divBbrst,Bbrst,JxBbrst,divTshearbrst,divPbbrst] = c_4_j('gseR?brsttime','gseB?');
    gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
    gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
    gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';
  else
    gseJcurl = irf.ts_vec_xyz(gseB1.time,gseB1.data*NaN);
  end
end
if 0 % Pitchangle distribution
  c_eval('ePitch? = ePDist?.pitchangles(dmpaB?,13); ePitch? = ePitch?.convertto(''s^3/km^6'');',ic)
end
if 0 % Motional electric fields, VxB
  c_eval('gseVexB? = gseVe?.cross(gseB?.resample(gseVe?))*1e-3; gseVexB?.units = ''mV/m'';',ic)
  c_eval('gseVixB? = gseVi?.cross(gseB?.resample(gseVi?))*1e-3; gseVixB?.units = ''mV/m'';',ic)
end
if 0 % E + V x B
  c_eval('gseEVexB? = gseE?.resample(gseVexB?.time)+gseVexB?; gseEVexB?.name = ''E+VexB'';',ic)
  c_eval('gseEVixB? = gseE?.resample(gseVixB?.time)+gseVixB?; gseEVixB?.name = ''E+VixB'';',ic)
end
if 1 % Epar, Eperp
  c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
end
if 0 % pressure and beta
  c_eval('PB? = gseB?.abs2/2/units.mu0*1e-9; PB?.name = ''Magnetic pressure''; PB?.units =''nPa'';',ic)
  c_eval('beta?e = gsePe?.trace/3/PB?.resample(gsePe?);',ic)
  c_eval('beta?i = gsePi?.trace/3/PB?.resample(gsePi?);',ic)
  c_eval('beta? = beta?i + beta?e.resample(beta?i);',ic)
end

disp('Ready.')