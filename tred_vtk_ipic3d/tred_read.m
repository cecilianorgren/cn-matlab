format = repmat('%f',1,14*27); % 14*27 numbers
%core=[0 100 101 102 103 104 105 106 107 1000 1001 1002 1003 1004 1057 1058 1059 1060,...
%    1061 1062 1063 1064 1065 1066 1074 1075 1076 1077 1078 1079 1080];
core=[436 437 444 445 541 533 532 540];
for p=1:length(core)
   core(p)
eval(['fid=fopen(''/Users/Cecilia/TRED46/VirtualSatellite/VirtualSatelliteTraces',num2str(core(p)),'.txt'')'])
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', 0);
eval(['Position(',num2str(p),',:,:)=[Pos{1} Pos{2} Pos{3}];'])
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 0);
% every point has 14 own columns for each time step
ind=reshape(1:length(Data),length(Data)/27,27);
for k=1:27
    % Structure: Quantity(core,sc (1-27),time,xyz)
    eval(['B(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(1,k)),'} Data{1,',num2str(ind(2,k)),'} Data{1,',num2str(ind(3,k)),'}];']) 
    eval(['E(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(4,k)),'} Data{1,',num2str(ind(5,k)),'} Data{1,',num2str(ind(6,k)),'}];']) 
    eval(['Je(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(7,k)),'} Data{1,',num2str(ind(8,k)),'} Data{1,',num2str(ind(9,k)),'}];']) 
    eval(['Ji(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(10,k)),'} Data{1,',num2str(ind(11,k)),'} Data{1,',num2str(ind(12,k)),'}];']) 
    eval(['Ne(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(13,k)),'}];']) 
    eval(['Ni(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(14,k)),'}];']) 
end
end