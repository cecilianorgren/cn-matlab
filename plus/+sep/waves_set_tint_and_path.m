%% Set datastore and savepath
%mms.db_init('local_file_db','/Volumes/Nexus/data');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   

%% Tint
sep.get_tints;
tint_brst = tint.utc; 
tint_day_utc = tint.utc; tint_day_utc = tint_day_utc(1,1:10);
tint_lobe_utc = tint_lobe.utc; tint_lobe_utc = tint_lobe_utc(:,12:23);
tint_sheet_utc = tint_sheet.utc; tint_sheet_utc = tint_sheet_utc(:,12:23);
tint_sep_utc = tint_sep.utc; tint_sep_utc = tint_sep_utc(:,12:23);
tint_phi_utc = tint_phi.utc; tint_phi_utc = tint_phi_utc(:,12:23);
%tint_waves_utc = tint_waves.utc; tint_waves_utc = tint_waves_utc(:,12:23);

%% Set paths
localuser = datastore('local','user');
pathRoot = ['/Users/' localuser '/GoogleDrive/Research/Separatrix_acceleration_events/']; 
pathRootWaves = [pathRoot 'saved_wave_data/'];
pathEvent = sprintf('waves_%s_%s_%s',tint_day_utc,tint_sep_utc(1,:),tint_sep_utc(2,:));
pathEvent(strfind(pathEvent,':'))=[];
pathFigures = [pathRootWaves pathEvent '/figures/'];
pathData = [pathRootWaves pathEvent '/wave_data/'];
mkdir([pathRootWaves pathEvent])
mkdir(pathData)
mkdir(pathFigures)

%% Set data format
data_format_write = '%s %s %s %s %s \n %s %s %s %s \n %f %f %f \n %f %f %f %f \n %f %f %f %f %f %f \n %f %f %f %f %f %f \n %f %f %f %f %f %f \n %f %f %f %f %f %f \n %f \n %f %f %f \n %f \n %f %f %f %f \n %f %f %f %f'; 
data_format_read =  '%s %s %s %s %s \n %s %s %s %s \n %f %f %f \n %f %f %f %f \n %f %f %f %f %f %f \n %f %f %f %f %f %f \n %f %f %f %f %f %f \n %f %f %f %f %f %f \n %f \n %f %f %f \n %f \n %f %f %f %f \n %f %f %f %f'; 
%cn.print(sprintf('sep_ov_%s_%s_%s',tint_day_utc,tint_sep_utc(1,:),tint_sep_utc(2,:)),'path',eventPath)