%% Load DB table
%localuser = 'cecilianorgren';
localuser = 'cecilia';
%file = ['/Users/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_b_gsm_2017-2022.nc'];

file = ['/home/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.nc'];
%data = load(file);
file_csv = ['/home/' localuser '/Data/Databases/DB_Richard_2022_v2/mms_bbfs_db_2017-2021.csv'];
tlim = readtable(file_csv);
tstart = tlim.Var1; tstart.TimeZone = "UTCLeapSeconds";
tstop = tlim.Var2;  tstop.TimeZone = "UTCLeapSeconds";

tstart_ttns = convertTo(tstart, 'tt2000');
tstop_ttns = convertTo(tstop, 'tt2000');
t_duration_s = double((tstop_ttns - tstart_ttns))*1e-9;

%ncdisp(file)

info = ncinfo(file);
vars = {info.Variables.Name};
nvars = numel(vars);

clear db
for ivar = 1:nvars
  db.(vars{ivar}) = ncread(file,vars{ivar});
end

db_table_ff = struct2table(db);
% Add stop time of slow
db_table_ff = addvars(db_table_ff,tstart_ttns,tstop_ttns,t_duration_s,'After','time');


t0_ff = EpochTT('2017-05-05T19:41:44.790324');
% time: microseconds since 2017-05-05T19:41:44.790324
time_ff_ttns = t0_ff.ttns + db_table_ff.time*1e3;
db_table_ff.time = time_ff_ttns; % rewrite time in ttns


t0_df = EpochTT('2017-05-19T03:06:44.458185978');
% t_df: nanoseconds since  2017-05-19T03:06:44.458185978
time_df_ttns = double(t0_df.ttns + int64(db_table_ff.t_df));
time_df_ttns(~db_table_ff.is_df) = NaN;
db_table_ff.t_df = time_df_ttns; % rewrite time in ttns


db_table_df = db_table_ff(db_table_ff.is_df==1,:);
nDF = numel(db_table_df.time);

%% Generate string
tint_all = irf.tint('2017-01-01T00:00:00.00Z/2023-01-01T00:00:00.00Z');
%files = mms.db_list_files('mms1_fpi_brst_l2_des-dist',tint_all);


%
tstart = EpochTT(db_table_df.tstart_ttns);
tstop = EpochTT(db_table_df.tstop_ttns);

str_sync_cell = {};
str_sync_str = "";
count = 0;
for it = 1:tstart.length
  disp('it = %g',it)
  tint = [tstart(it) tstop(it)];
  %file = mms.db_list_files('mms1_fpi_brst_l2_des-dist',tint);
  try
    [filepath, filename] = mms.get_filepath(['mms1_fpi_brst_l2_dis-dist'], tint(1));
  catch
    disp(sprintft('%g: No file. Continuing.',it))
  end
  file.path = filepath;
  file.name = filename;
  if not(isempty(file))
    count = count + 1;
    string_split = split(file.name,'_');
    str_sync_cell{count} = "--include='*" + "_" + string_split{6} + "*.cdf'\ ";
    str_sync_str = str_sync_str + "--include='*" + "_" + string_split{6} + "*.cdf'\ ";
    if mod(count,10)
      disp(str_sync_str)
    end
  end
end
disp(str_sync_str)

%time_epochtt = EpochTT(int64(db_table_df(33,end).t_df));
%file = mms.db_list_files('mms1_fpi_brst_l2_des-dist',time_epochtt);