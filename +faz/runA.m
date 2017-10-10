
fid = fopen('/Users/Cecilia/Research/Faz/list_IB.txt','r');
list = textscan(fid,'%s%*s%s%s%*s%*s%s%s%s%s%s%s');

events = list{1}(:);
nEvent = numel(events);

for ii=1:nEvent
    event = events{ii};
    dateiso = ['20' event(1:2) '-' event(3:4) '-' event(5:6) 'T' event(7:8) ':' event(9:10) ':' event(11:12) 'Z'];
    disp(dateiso)
    st = iso2epoch(dateiso);
    ib.B;
    %irf_zoom(h,'x',[st-20 st+20])
    %irf_pl_mark(h,[st-10 st+10])
    pause
   
    %
    %irf_pl_mark([st-10 st+10])
    cn.print(event(1:12),'faz/Bbb')
    clearvars -except events
    
end
    