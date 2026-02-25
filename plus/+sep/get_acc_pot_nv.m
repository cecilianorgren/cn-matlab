%% Load datastore
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS/');
db_info = datastore('mms_db');   
localuser = datastore('local','user');

events = 1:20;
nevents = numel(events);
saveAccPotPath = ['/Users/' localuser '/MATLAB/cn-matlab/+sep/acc_potential/'];
printAccPotPath = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/acceleration_potential/'];

for ievent = 1:nevents
  event = events(ievent); % set event id
  sep.get_tints; % get tints for that event id
 
  %% Load data  
  ic = 1:4;
  units = irf_units;
  % Core data   
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
  if isempty(gseB4); ic = 1:3; end
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);  
  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)  
  c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)  
  c_eval('gseFlux? = ne?*gseVe?;',ic);
  c_eval('gseFlux?par = ne?*gseVe?par;',ic);
  c_eval('vte?par = (2*units.eV*facTe?.xx/units.me).^0.5;',ic)
  c_eval('vte?perp = (2*units.eV*(1/2)*(facTe?.yy+facTe?.zz)/units.me).^0.5;',ic)
  c_eval('vte? = (2*units.eV*(1/3)*(facTe?.xx+facTe?.yy+facTe?.zz)/units.me).^0.5;',ic)
  
  cad_resamp = 5;
  timeline = ne1.time(1:cad_resamp:end);
  c_eval('ne?_resamp = ne?.resample(timeline);',ic);
  c_eval('ve?_resamp = gseVe?par.resample(timeline);',ic);
  c_eval('flux?_resamp = gseFlux?par.resample(timeline);',ic);
  
  %% Pick out parameters in the lobe, sheet, and separatrix intervals  
  c_eval('n_lobe(?) = mean(ne?.tlim(tint_lobe).data,1);',ic)
  c_eval('n_sheet(?) = mean(ne?.tlim(tint_sheet).data,1);',ic)
  c_eval('n_sep(?) = mean(ne?.tlim(tint_sep).data,1);',ic)
  c_eval('n_sep_min(?) = min(ne?.tlim(tint_sep).data);',ic)
  c_eval('n_sep_min_resamp(?) = min(ne?_resamp.tlim(tint_sep).data);',ic)
  c_eval('ve_sep_max_resamp(?) = max(ve?_resamp.tlim(tint_sep).data);',ic)
  c_eval('flux_sep_max_resamp(?) = max(flux?_resamp.tlim(tint_sep).data);',ic)

  c_eval('Tepar_lobe(?) = mean(facTe?.tlim(tint_lobe).data(:,1,1),1);',ic)
  c_eval('Teperp_lobe(?)  = 0.5*mean(facTe?.tlim(tint_lobe).data(:,2,2)+facTe?.tlim(tint_lobe).data(:,3,3),1);',ic)

  c_eval('vtepar_lobe = sqrt(2*units.eV*Tepar_lobe./units.me)/1000;')
  c_eval('vteperp_lobe = sqrt(2*units.eV*Teperp_lobe./units.me)/1000;')
  
  colors = mms_colors('xzy');
  info_color = [0 0 0; colors; colors; colors; 0 0 0; 0 0 0; 0 0 0; colors];

  tint_day_utc = tint.utc; tint_day_utc = tint_day_utc(1,1:10);
  tint_lobe_utc = tint_lobe.utc; tint_lobe_utc = tint_lobe_utc(:,12:23);
  tint_sheet_utc = tint_sheet.utc; tint_sheet_utc = tint_sheet_utc(:,12:23);
  tint_sep_utc = tint_sep.utc; tint_sep_utc = tint_sep_utc(:,12:23);
  tint_phi_utc = tint_phi.utc; tint_phi_utc = tint_phi_utc(:,12:23); 
  
  n_sep_min_resamp./n_lobe
  n_sep_min./n_lobe
  n_sep./n_lobe
  ve_sep_max_resamp./vtepar_lobe
  flux_sep_max_resamp./vtepar_lobe./n_lobe
 
  %% Save data to file
  % .mat or .txt (together or separate?)  
  acc_nv_data.event_id = event;
  acc_nv_data.ic = ic;  
  acc_nv_data.n_lobe = n_lobe;
  acc_nv_data.n_sheet = n_sheet;
  acc_nv_data.n_sep = n_sep;
  acc_nv_data.n_sep_min = n_sep_min;
  acc_nv_data.n_sep_min_resamp = n_sep_min_resamp;
  acc_nv_data.ve_sep_max_resamp = ve_sep_max_resamp;
  acc_nv_data.flux_sep_max_resamp = flux_sep_max_resamp;
  acc_nv_data.Tepar_lobe = Tepar_lobe;
  acc_nv_data.Teperp_lobe = Teperp_lobe;
  acc_nv_data.vtepar_lobe = vtepar_lobe;
  acc_nv_data.vteperp_lobe = vteperp_lobe;
  
  save(sprintf('%s/acc_nv_data_event_%g',saveAccPotPath,event),'acc_nv_data')
  %fid = fopen([matlabPath 'acc_pot_statistics.txt'],'a+');
  % event_id
  %save_format = '%f';
  %fid = fclose(fid);
  %fprintf('printed to file %sesw_properties.txt: %s /n',matlabPath,print_str)
end