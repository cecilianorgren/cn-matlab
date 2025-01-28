%% 2015-12-28/29: overview
% SITL_20151228_214634-105644
tint = irf.tint('2015-12-28T22:00:00.00Z/2015-12-29T10:00:00.00Z');
mms_survey.load_data_overview

%% 2015-12-28/29: brst
time = irf_time('2015-12-28T22:19:05.00Z','utc>epochtt');
ic = 1;
mms_survey.load_data