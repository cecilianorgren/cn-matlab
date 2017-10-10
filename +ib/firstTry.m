%% 
cd /Users/Cecilia/Research/Faz/data;
event = '130522112527we.04';
scStr = event(end);
scNum = str2num(event(end));
st = iso2epoch('2013-05-22T11:25:27Z');
storagePath = '/Users/Cecilia/Research/Faz/data';

%% Check parameters in a file
% must mount the disc first by cmd + K in Finder and
% nfs://db.irfu.se/export/data/cluster
c_efw_burst_param(['/Volumes/cluster/burst/' event])

%% Get spacecraft spin phase
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'a')
c_load(['Atwo' scStr])

%% Get probe data
% Fetches mEFWburstR containing P4kHz4p?
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'pburst')
%load mEFWburstR

%% Get e-field data
% Since I give a ClusterDB-object to the getData function, it will call the
% getData function/method which is in @ClusterDB folder
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'eburst')
% This should (right?) work, but doesn't, why?


%% 1
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'ibias')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'efwt')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'fdm')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'tmode')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'e')
getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'p')

%% 2
getData(ClusterProc(storagePath),scNum,'whip')
getData(ClusterProc(storagePath),scNum,'sweep')
getData(ClusterProc(storagePath),scNum,'bdump')
getData(ClusterProc(storagePath),scNum,'dies')
getData(ClusterProc(storagePath),scNum,'badbias')
getData(ClusterProc(storagePath),scNum,'probesa')
getData(ClusterProc(storagePath),scNum,'p')
getData(ClusterProc(storagePath),scNum,'ps')
%% 3
getData(ClusterProc(storagePath),scNum,'dies')

%% Get e-field data
% First run 1 2 3 above
% Since I give a ClusterProc-object to the getData function, it will call 
% the getData function/method whcCich is in @ClusterProc folder.
% The loaded data m____.mat files need to be in the storagePath directory.
% For example mA.mat which contains the phase.
getData(ClusterProc(storagePath),scNum,'dieburst')
irf_plot(wbE4p12)

%% Try to get bursts
dirOld = cd;
%cd /Volumes/cluster/burst/;
caa_get_bursts(event,2)
%cd(dirOld)

%% 
c_get_batch(st,60,scNum,'vars','e','sp','/Users/Cecilia/Research/Faz/data/data-c_get_batch/')