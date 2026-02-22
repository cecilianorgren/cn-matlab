% generate 
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2023-01-01T00:00:00.00Z');
files = mms.db_list_files('mms1_fpi_brst_l2_des-dist',tint_all);


%
tstart = EpochTT(db_table_df.tstart_ttns);
tstop = EpochTT(db_table_df.tstop_ttns);

str_sync_cell = {};
str_sync_str = "";
count = 0;
for it = 1:tstart.length
  tint = [tstart(it) tstop(it)];
  file = mms.db_list_files('mms1_fpi_brst_l2_des-dist',tint);
  if not(isempty(file))
    count = count + 1;
    string_split = split(file.name,'_');
    str_sync_cell{count} = "--include='*" + "_" + string_split{6} + "*.cdf'\ ";
    str_sync_str = str_sync_str + "--include='*" + "_" + string_split{6} + "*.cdf'\ ";
  end
end
%time_epochtt = EpochTT(int64(db_table_df(33,end).t_df));

%file = mms.db_list_files('mms1_fpi_brst_l2_des-dist',time_epochtt);