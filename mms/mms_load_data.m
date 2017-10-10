% mms_load_data

% Time representation
tint = irf.tint('2015-08-15T13:00:00Z/2015-08-15T13:30:00Z');
t1 = tint(1);
t1array = fromepoch(t1.epochUnix);
tmean = tint.start+(tint.stop-tint.start)/2
utc = tmean.toUtc;
t.year  = str2double(utc(1:4));
t.month = str2double(utc(6:7));
t.day   = str2double(utc(9:10));
t.hour  = str2double(utc(12:13));
t.min   = str2double(utc(15:16));
t.sec   = str2double(utc(18:end-1));
tdateNum = datenum(t.year,t.month,t.day,t.hour,t.min,t.sec);

% File information
scId = 2;
filePrefix = irf_ssub('mms?_edp_brst_l2_scpot',scId);
dataVar = irf_ssub('mms?_edp_psp',scId);

% Directory information
C = strsplit(filePrefix,'_');
varDir = '/data/mms';
for ix=1:length(C), varDir = [varDir filesep C{ix}]; end
fileDir = [utc(1:4) filesep utc(6:7) filesep utc(9:10) filesep];
totDir = [varDir filesep fileDir];

% List files in directory
listingD = dir([varDir filesep fileDir filesep filePrefix '*.cdf']);

% Find file that contains interval start time, t1
nFiles = numel(listingD);
for ii = 1:nFiles, 
    dateStr = listingD(ii).name(end-24:end-11);    
    dateNums(ii,1) = datenum(str2double(dateStr(1:4)),str2double(dateStr(5:6)),str2double(dateStr(7:8)),str2double(dateStr(9:10)),str2double(dateStr(11:12)),str2double(dateStr(13:14))); 
end
% plot(1:nFiles,dateNums,'-*',[1 nFiles],tdateNum*[1 1])
fileId = find(tdateNum>dateNums,1,'last');
listingD(fileId).name
utc
listingD(fileId+1).name

%% Load the data
tint = irf.tint('2015-08-15T13:00:00Z/2015-08-15T13:30:00Z');
t = tint.start+(tint.stop-tint.start)/2;
scId = 1:4;
c_eval('[dcvP?,dobjP?] = cn_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv'',t);',scId);
c_eval('[P?,dobjP?] = cn_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',t);',scId);
c_eval('[B?,dobjB?] = cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',t);',scId);
c_eval('[E?,dobjE?] = cn_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',t);',scId);
%c_eval('[P?,dobjP?] = cn_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot'',t);',scId);






