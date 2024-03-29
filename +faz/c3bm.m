% Script to load all data necessary in order to construct electric field
% and probe voltages from raw data.

fid = fopen('/Users/Cecilia/Research/Faz/list_BM.txt');
list = textscan(fid,'%s%s%s');

listSC = list{1};
t1iso = list{2};
t2iso = list{3}; 
nEvent = numel(t1iso);

for ii = 1:nEvent
    listTint{ii} = [iso2epoch(t1iso{ii}) iso2epoch(t2iso{ii})];
    listDt{ii} = diff(listTint{ii});
end
%%

for ii = 1%:nEvent
    storagePath = ['/Users/Cecilia/Research/Faz/data/BM/' ...
                   irf_time(listTint{1}(1),'yyyymmddhhmmss')];
    if ~exist(storagePath,'dir'); mkdir(storagePath); end
    cd(storagePath)
    st = listTint{ii}(1);
    dt = listDt{ii};
    sc = str2num(listSC{ii}(2));
        
    % must mount the disc first by cmd + K in Finder and
    % nfs://db.irfu.se/export/data/cluster

    %% Load data and save to working directory
    if 0
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'a')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'pburst')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'ibias')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'efwt')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'fdm')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'tmode')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'e')
    getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'p')
    getData(ClusterProc(storagePath),sc,'whip')
    getData(ClusterProc(storagePath),sc,'sweep')
    %getData(ClusterProc(storagePath),sc,'bdump')
    %getData(ClusterProc(storagePath),scNum,'dies')
    getData(ClusterProc(storagePath),sc,'badbias')
    getData(ClusterProc(storagePath),sc,'probesa')
    getData(ClusterProc(storagePath),sc,'p')
    getData(ClusterProc(storagePath),sc,'ps')
    getData(ClusterProc(storagePath),sc,'die')
    
    %% Finally construct diE burst data.
    % These also look for normal burst mode, i.e. the 180 Hz filter.
    getData(ClusterProc(storagePath),sc,'dieburst')
    c_load('dibE4p1234')
    end
    %%
    if 1
        getData(ClusterDB('db:9','/Volumes/cluster'),st,dt,sc,'a')
        toLoad = {'tmode','fdm','efwt','ibias','p','e','a','sax',...
                  'r','v','b','vcis'};
        %c_get_batch(st-200+,dt,3:4,'vars',toLoad)        
        %c_get_batch(st,dt,3:4,'vars','tmode')        
        c_get_batch(st,dt,sc,'vars','e')        
        %c_eval('c_load(mTMode?)',sc)        
        %figure(3);c_eval('plot(mTMode?)',sc);
        c_load('diE3p1234')
        %figure(4);irf_plot(diE3p1234)
        %cn.f(diE4p1234)
        c_load('wE3p34')
        c_load('wE3p32')
        dwE3p24 = irf_add(1,wE3p34,-1,wE3p32);
        irf_plot({wE3p32,wE3p34,dwE3p24,diE3p1234})
    end
    
    
end
    

