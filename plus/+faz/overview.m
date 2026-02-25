% 2013-03-29T05:08:00.000Z 2013-03-29T05:05:53.000Z
%t1 = iso2epoch('2013-03-29T05:08:00.000Z');
%t2 = iso2epoch('2013-03-29T05:05:53.000Z');
%tint = [t1 t2];

fid = fopen('/Users/Cecilia/Research/Faz/list_BM.txt');
list = textscan(fid,'%s%s%s');



ntint=numel(list{:,1});
tint = cell(ntint,1);
tintiso = cell(ntint,2);
for ii = 1:ntint; 
    tint{ii} = [iso2epoch(list{:,2}{ii}) iso2epoch(list{:,3}{ii})]; 
    tintiso{ii,1} = list{:,2}{ii}; 
    tintiso{ii,2} = list{:,3}{ii};    
    scliststr{ii,1} = list{:,1}{ii};    
    sclist(ii,1) = str2num(scliststr{ii}(2));    
end

%% Load data from ISDAT
for ii = 6%1:ntint
    storagePath = ['/Users/Cecilia/Research/Faz/data/' tintiso{ii,1}];
    if ~exist(storagePath,'dir'); mkdir(storagePath); end
    cd(storagePath)
    scStr = scliststr{ii,1};
    scNum = sclist(ii,1);
    st = tint{ii}(1);
    dt = diff(tint{ii,1})
    
    %% Check parameters in a file
    % must mount the disc first by cmd + K in Finder and
    % nfs://db.irfu.se/export/data/cluster
    %c_efw_burst_param(['/Volumes/cluster/burst/' event])

    %% Load data and save to working directory
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'a')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-30,60,scNum,'pburst')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'ibias')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'efwt')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'fdm')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'tmode')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'e')
    getData(ClusterDB('db:9','/Volumes/cluster'),st-10,dt+20,scNum,'p')
    getData(ClusterProc(storagePath),scNum,'whip')
    getData(ClusterProc(storagePath),scNum,'sweep')
    getData(ClusterProc(storagePath),scNum,'bdump')
    getData(ClusterProc(storagePath),scNum,'dies')
    getData(ClusterProc(storagePath),scNum,'badbias')
    getData(ClusterProc(storagePath),scNum,'probesa')
    getData(ClusterProc(storagePath),scNum,'p')
    getData(ClusterProc(storagePath),scNum,'ps')
    getData(ClusterProc(storagePath),scNum,'die')
    
    load mEDSI;
    irf_plot(diELXs3p42)
    cn.print(tintiso{ii,1},'ov2')
end
