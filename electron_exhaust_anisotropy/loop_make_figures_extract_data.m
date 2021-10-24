ic = 1;
units = irf_units;
events = {'2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z'};
nEvents = numel(events);
userRoot = '/home/cecilia/';
userRootProject = '/home/cecilia/Research/ElectronExhaustAnisotropy/';


% Setup database
mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   
  
for iEvent = 1:nEvents
  tintEventStr = events{iEvent}; 
  tint = EpochTT(tintEventStr);  
  
  axh_anis.load_data;
  
  setup_dir(userRootProject,gseE1.userData.GlobalAttributes.Logical_file_id)
  
  %% Event path
  

  
end
