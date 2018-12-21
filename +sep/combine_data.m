function table = combine_data()
localuser = datastore('local','user');
table_data_folder = datastore('acceleration','table_data_folder');
save_path = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/' table_data_folder '/'];
file_list = dir([save_path '*.txt']);
table = [];
n_events = numel(file_list);

for i_event = 1:n_events
  fid = fopen([save_path file_list(i_event).name],'r');
  tmp_data = textscan(fid,'%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f');
  for i_data = 1:numel(tmp_data)
    table(i_event,i_data) = tmp_data{i_data};
  end  
end