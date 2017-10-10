format = repmat('%f',1,14*27); % 14*27 numbers

for k=1:length(Psi)
eval(['fid=fopen(''/Users/Cecilia/VTK/VirtualSatellite/VirtualSatelliteTraces',num2str(cores(Psi(k))),'.txt'');'])
Pos=textscan(fid, '%f%f%f', 27, 'delimiter', '\n', 'headerlines', sc(Psi(k))-1);
%eval(['Position(k,:,:)=[Pos{1} Pos{2} Pos{3}];'])
Data=textscan(fid, format, 'delimiter', '\n', 'headerlines', 27-sc(Psi)-1);
% every point has 14 own columns for each time step
ind=reshape(1:length(Data),length(Data)/27,27);
for k=1:27
    % Structure: Quantity(core,sc (1-27),time,xyz)
    eval(['B(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(1,k)),'} Data{1,',num2str(ind(2,k)),'} Data{1,',num2str(ind(3,k)),'}];']) 
    eval(['E(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(4,k)),'} Data{1,',num2str(ind(5,k)),'} Data{1,',num2str(ind(6,k)),'}];']) 
    eval(['Ji(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(7,k)),'} Data{1,',num2str(ind(8,k)),'} Data{1,',num2str(ind(9,k)),'}];']) 
    eval(['Je(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(10,k)),'} Data{1,',num2str(ind(11,k)),'} Data{1,',num2str(ind(12,k)),'}];']) 
    eval(['Ni(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(13,k)),'}];']) 
    eval(['Ne(',num2str(p),',',num2str(k),',:,:) = [Data{1,',num2str(ind(14,k)),'}];']) 
end
end