% Script to load all data necessary in order to construct electric field
% and probe voltages from raw data.

storagePath = '/Users/Cecilia/Research/Faz/data/A';
cd(storagePath)
%event = '130522112527we.04';
%event = '130527000310we.04';
%event = '130414010159we.02';
%event = '130527000310we.04';
%event = '130527000629we.02';
scStr = event(end);
scNum = str2num(event(end));
dateiso = ['20' event(1:2) '-' event(3:4) '-' event(5:6) 'T' event(7:8) ':' event(9:10) ':' event(11:12) 'Z'];
%st = iso2epoch('2013-05-22T11:25:27Z');
st = iso2epoch(dateiso);

%% Check parameters in a file
% must mount the disc first by cmd + K in Finder and
% nfs://db.irfu.se/export/data/cluster
c_efw_burst_param(['/Volumes/cluster/burst/' event])

%% Load data and save to working directory
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'a')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'pburst')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'ibias')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'efwt')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'fdm')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'tmode')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'e')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'p')
getData(ClusterProc(storagePath),scNum,'whip')
getData(ClusterProc(storagePath),scNum,'sweep')
getData(ClusterProc(storagePath),scNum,'bdump')
getData(ClusterProc(storagePath),scNum,'dies')
getData(ClusterProc(storagePath),scNum,'badbias')
getData(ClusterProc(storagePath),scNum,'probesa')
getData(ClusterProc(storagePath),scNum,'p')
getData(ClusterProc(storagePath),scNum,'ps')
getData(ClusterProc(storagePath),scNum,'die')

%% Finally construct diE burst data.
getData(ClusterProc(storagePath),scNum,'dieburst')

%% Plot and compare
c_load('dibE4p1234')
c_load('diE4p1234')
h=irf_plot({dibE4p1234,diE4p1234},'comp');
%irf_plot(dibE4p1234)
