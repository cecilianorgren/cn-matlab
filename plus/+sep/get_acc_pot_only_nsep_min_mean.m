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
  c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
  %c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic);
  c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
  c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)
  
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)  
  %% Pick out parameters in the lobe, sheet, and separatrix intervals
  ic_orig = ic;
  ic = 1; 

  c_eval('n_lobe = mean(ne?.tlim(tint_lobe).data,1);',ic)
  c_eval('n_sheet = mean(ne?.tlim(tint_sheet).data,1);',ic)
  c_eval('n_sep = mean(ne?.tlim(tint_sep).data,1);',ic)
  c_eval('n_sep_min = min(ne?.tlim(tint_sep).data);',ic)


  colors = mms_colors('xzy');
  info_color = [0 0 0; colors; colors; colors; 0 0 0; 0 0 0; 0 0 0; colors];

  tint_day_utc = tint.utc; tint_day_utc = tint_day_utc(1,1:10);
  tint_lobe_utc = tint_lobe.utc; tint_lobe_utc = tint_lobe_utc(:,12:23);
  tint_sheet_utc = tint_sheet.utc; tint_sheet_utc = tint_sheet_utc(:,12:23);
  tint_sep_utc = tint_sep.utc; tint_sep_utc = tint_sep_utc(:,12:23);
  tint_phi_utc = tint_phi.utc; tint_phi_utc = tint_phi_utc(:,12:23);
  
  ic = ic_orig;     
  
 
  %% Save data to file
  % .mat or .txt (together or separate?)  
  acc_n_data.event_id = event;
  acc_n_data.ic = ic;  
  acc_n_data.n_lobe = n_lobe;
  acc_n_data.n_sheet = n_sheet;
  acc_n_data.n_sep = n_sep;
  acc_n_data.n_sep_min = n_sep_min;
  
  save(sprintf('%s/acc_n_data_event_%g',saveAccPotPath,event),'acc_n_data')
  %fid = fopen([matlabPath 'acc_pot_statistics.txt'],'a+');
  % event_id
  %save_format = '%f';
  %fid = fclose(fid);
  %fprintf('printed to file %sesw_properties.txt: %s /n',matlabPath,print_str)
end