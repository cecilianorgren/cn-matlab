function save_data(event,phi,B,n_lobe,n_sheet,n_sep,Tpar_lobe,Tpar_sheet,Tpar_sep,Tperp_lobe,Tperp_sheet,Tperp_sep,fce0,fpe0,beta_e_lobe,beta_e_sheet,beta_e_sep)
% saves data to separate text files, identified by their event id
localuser = datastore('local','user');
table_data_folder = datastore('acceleration','table_data_folder');
save_path = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/' table_data_folder '/'];

name_file = sprintf('event%g.txt',event);
fid = fopen([save_path name_file],'w');

%fprintf(fid,'%f %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f',...
fprintf(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f',...
  event,phi,B,n_lobe,n_sheet,n_sep,Tpar_lobe,Tpar_sheet,Tpar_sep,Tperp_lobe,Tperp_sheet,Tperp_sep,fce0,fpe0,beta_e_lobe,beta_e_sheet,beta_e_sep);
fid = fclose(fid);