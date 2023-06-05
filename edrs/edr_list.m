tint_all = irf.tint('2015-01-01T00:00:00.00Z/2024-01-01T00:00:00.00Z');


file_ids = ['20180827114043';...
            '20170703052613';...
            '20170703052703';...
            '20170810121733'];
% 20180827114043
% 20170703052613
% 20170703052703
% 20170703052703

for ifid = 1:size(file_ids,1)
  fid = file_ids(ifid,:);
  time_vector = [str2num(fid(1:4)),str2num(fid(5:6)),str2num(fid(7:8)),str2num(fid(9:10)),str2num(fid(11:12)),str2num(fid(13:14))];
  t1 = irf_time(time_vector,'vector6>utc');
  files = mms.db_list_files('mms1_fgm_brst_l2',EpochTT(t1) + [1 2]);
  tint_tmp = [files.start files.stop] + [1 -1];
  tints{ifid} = tint_tmp;
end


%% Load data
tint = tints{1};
fileId = fid;

ic = 1;
edr_load_data