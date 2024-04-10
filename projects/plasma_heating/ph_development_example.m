%% Get list of magnetosheath time intervals
localuser = 'cecilia';
table_ubmshc = readtable('/Users/cecilia/MATLAB/cn-matlab/projects/plasma_heating/data/mms_unbiased_magnetosheath_campaign.txt');
region_prob = load(['/Users/' localuser '/Data/MMS/DB_Lalti/Proba_full.mat']);

%% Load data, unbiased magnetosheath campaign
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');

nIntervals = size(table_ubmshc,1);

files = struct([]);

for iInterval = 1:nIntervals
  t1 = table_ubmshc.time{iInterval}; t1 = [t1(1:10), 'T', t1(12:19), '.0Z'];
  t2 = table_ubmshc.xEnd{iInterval}; t2 = [t2(1:10), 'T', t2(12:19), '.0Z'];
  tint = EpochTT([t1; t2]);
  
  files_tmp = mms.db_list_files('mms1_fgm_brst_l2',tint);
  files = cat(1,files,files_tmp');
end

%% Go through burst intervals to identify current sheets

nFiles = size(files,1);
for iFile = 1%:nFiles
  tint = [files(iFile).start files(iFile).stop];
  
  % Load data  
  c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)
  c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
  c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
  [Jcurl,divB,gseB,JxB,gseCurvB,gseDivPb] = c_4_j('gseR?brsttime','gseB?');


  %% Find current sheets
  % Box average current to given window
  dT = 10;
  box_timeline = tint(1):dT:tint(2);
  Jcurl_box = Jcurl.resample(box_timeline);

  deltaJ = Jcurl + -1*Jcurl_box.resample(Jcurl);
  
  % Criteria
  Threshold = 10e-9;
  MinPeakProminence = 20e-9;
  [PKS,LOCS,W] = findpeaks(Jcurl.abs.data,'MinPeakProminence',MinPeakProminence);

  % Plot data
  
  colors = pic_colors('matlab');
  h = irf_plot(4);
  
  hca = irf_panel('B');
  set(hca,'colororder',mms_colors('xyza'))
  irf_plot(hca,{gseB.x,gseB.y,gseB.z,gseB.abs},'comp')
  
  hca = irf_panel('Jcurl');
  set(hca,'colororder',mms_colors('xyza'))
  irf_plot(hca,{Jcurl.x,Jcurl.y,Jcurl.z,Jcurl.abs},'comp')

  hca = irf_panel('Jabs');
  set(hca,'colororder',mms_colors('1234'))
  irf_plot(hca,{Jcurl.abs,Jcurl_box.abs},'comp')
  
  hca = irf_panel('delta J');
  set(hca,'colororder',mms_colors('1234'))
  irf_plot(hca,{deltaJ.abs},'comp')

  % Plot peaks
  time_pks = Jcurl.time(LOCS);  
  %tint_mark = [time_pks.epochUnix-0.5*W, time_pks.epochUnix+0.5*W];
  tint_mark = time_pks.epochUnix;
  for iTint = 1:numel(LOCS)
    irf_pl_mark(h,tint_mark(iTint,:),'k')
  end

end