% study event 130527000629we.02
% 130527000310we.04 ? 9k ?4 ? ? 2 10 ?V1H, V2H, V3H, V4H
% 130527000330we.03 ? 9k ?4 ? ? 2 10 ?V3H, V2H, SCZ, SCY
% 130527000431we.01 ? 9k ?4 ? ? 2 10 ?V3H, V2H, SCZ, SCY
% 130527000629we.02 ? 9k ?4 ? ? 2 10 ?V3U, V2H, V3H, V4H

event = '130527000310we.04';
%event = '130522112527we.04';

%event = '130414010159we.02';
%event = '130407060313we.02';
%event = '130502031240we.02';
%event = '130511041127we.02';
%event = '130522112442we.02';
%event = '130527000629we.02';
%
disp(event)


% Make directory
storagePath = ['/Users/Cecilia/Research/Faz/data/' event(1:12)];
if ~exist(storagePath,'dir'); mkdir(storagePath); end
cd(storagePath)

% Start time etc.
scStr = event(end);
scNum = str2num(event(end));
isoDate = ['20' event(1:2) '-' event(3:4) '-' event(5:6) 'T' event(7:8) ':' event(9:10) ':' event(11:12) 'Z'];
st = iso2epoch(isoDate);
disp(isoDate)



% Check parameters in a file
% must mount the disc first by cmd + K in Finder and
% nfs://db.irfu.se/export/data/cluster
c_efw_burst_param(['/Volumes/cluster/burst/' event])

% Load data and save to working directory
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

% Finally construct diE burst data.
getData(ClusterProc(storagePath),scNum,'dieburst')

% Plot and compare
%c_load('dibE4p1234')
%c_load('diE4p1234')
%figure(51); irf_plot({dibE4p1234,diE4p1234},'comp')
figure(52); burstData = caa_get_bursts(event,2);
cn.print(event(1:12),'caa_get_bursts')
%% Zoom in to compare
tz1 = [2013 05 27 00 03 11.05]; 
tz2 = [2013 05 27 00 03 11.09];
tzoom = toepoch([tz1;tz2])';

figure(51); udata = get(gcf,'userdata'); h = udata.subplot_handles;
irf_zoom(h,'x',tzoom);
set(h,'ylim',[-200 200])
figure(52); h = get(gcf,'Children'); %h = udata.subplot_handles;
irf_zoom(h([4 7 8]),'x',tzoom);
set(h([4 7]),'ylim',[-200 200])

% Make additional plot with individual probe data

c_load('P4kHz4p?')
figure(53); h = irf_plot({P4kHz4p1,P4kHz4p2,P4kHz4p3,P4kHz4p4},'comp');
irf_zoom(h,'x',tzoom);

%% Try doing interferometry measurements
% use yk.my_interf
c_eval('p? = P4kHz4p?;',1:4)

yk.my_interf(p1,p2,p3,p4);



