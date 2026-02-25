function cores = tred_list_cores

unix('find /Users/Cecilia/TRED46/VirtualSatellite/* > /Users/Cecilia/TRED46/satellites.txt');
fid=fopen('/Users/Cecilia/TRED46/satellites.txt'); 
existing=textscan(fid,'/Users/Cecilia/TRED46/VirtualSatellite/%s');
existing=existing{1};
for k=1:max(size(existing))
    cores(k,1)=str2num(existing{k}(23:(end-4)));
end
cores=sort(cores);