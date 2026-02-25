% Script to use caa_get_burst to construct electric field and probe data
% from raw data.

storagePath = '/Users/Cecilia/Research/Faz/data/Bb';
cd(storagePath)
event = '130522112527we.04';
%event = '130407060313we.02';
scStr = event(end);
scNum = str2num(event(end));
%st = iso2epoch('2013-05-22T11:25:27Z');

%% Check parameters in a file
% must mount the disc first by cmd + K in Finder and
% nfs://db.irfu.se/export/data/cluster
c_efw_burst_param(['/Volumes/cluster/burst/' event])

%% See if loading all these before caa_get_bursts help.
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

%% Call caa_get_burst
caa_get_bursts(event,2)

%% Plot and compare
%c_load('dibE4p1234')
%c_load('diE4p1234')
%irf_plot({dibE4p1234,diE4p1234},'comp')

