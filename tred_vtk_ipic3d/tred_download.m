function [to_download wanted existing] = tred_download(core)
unix('find /Users/Cecilia/TRED46/VirtualSatellite/* > /Users/Cecilia/TRED46/satellites.txt');
fid=fopen('/Users/Cecilia/TRED46/satellites.txt'); 
existing=textscan(fid,'/Users/Cecilia/TRED46/VirtualSatellite/%s');
existing=existing{1};
for k=1:max(size(core))
    eval(['wanted{k,1}=''VirtualSatelliteTraces',num2str(core(k)),'.txt'';'])
end
[~, i1]=intersect(wanted,existing);
[i2,~]=ismember(wanted,existing);
if isempty(i1)
    to_download=wanted;
else
    to_download=wanted(~i2,1);
end
if ~isempty(to_download)
    n_dl=max(size(to_download));
    download_string='scp -q cecilia@spis.irfu.se:"';
    for k=1:n_dl
        download_string=[download_string,...
            '/net/db/export/share/ReconnectionSimulations/TRED46_DATA/',...
            to_download{k,1},' '];
    end
    download_string=[download_string,'" /Users/Cecilia/TRED46/VirtualSatellite/']
    tic;
    eval(['unix(''',download_string,''');'])
    timeelapsed=toc;
    disp(['Time elapsed: ',num2str(timeelapsed)])
else
    disp('All wanted traces already exist!')
end
