function data = tred_surr_sat(core,sat,thickness)
% TRED_SURR_SAT(core,satellite)
%
%   Finds the core-satellites surrounding the given satellite.
%

fid=fopen('/Users/Cecilia/TRED46/DATA/gridpoints.txt'); %add ,'a' to change file
for k=1:3; 
    gridpoints(:,k)=textscan(fid,'%f','headerlines',1); 
end; 
fclose(fid);

fid=fopen('/Users/Cecilia/TRED46/DATA/Positions.txt');
data=textscan(fid,['%*s %s %*s ' repmat('%f ',1,27*3)],...
    'headerlines',0,'Whitespace',' \b\t\n\r','CollectOutput',true);
x=unique(data{2}(:,1:3:end));
y=unique(data{2}(:,2:3:end));
z=unique(data{2}(:,3:3:end));
Positions=[reshape(data{2}(:,1:3:end)',1,27*size(data{2},1))',...
    reshape(data{2}(:,2:3:end)',1,27*size(data{2},1))',... 
    reshape(data{2}(:,3:3:end)',1,27*size(data{2},1))'];
for k=1:length(data{1}) % collect cores
    uniquecores(k)=sscanf(char(data{1}(k)),'VirtualSatelliteTraces%f.txt');
end
cores=reshape(repmat(uniquecores',1,27)',27*length(uniquecores),1);
satellites=repmat((1:27)',length(uniquecores/27),1);

% find position of given satellite
[~,satindex]=intersect([cores satellites],[core sat],'rows');
midpos=Positions(satindex,:);
xs=[x(find(x>midpos(1),thickness,'first'))' midpos(1) x(find(x<midpos(1),thickness,'last'))'];
ys=[y(find(y>midpos(2),thickness,'first'))' midpos(2) y(find(y<midpos(2),thickness,'last'))'];
zs=[z(find(z>midpos(3),thickness,'first'))' midpos(3) z(find(z<midpos(3),thickness,'last'))'];

% make combinations of vectors
combs=zeros(length(xs)*length(ys)*length(zs),3);
index=0;
for k=1:length(xs)
    for l=1:length(ys)
        for m=1:length(zs)
            index=index+1;
            combs(index,:)=[xs(k) ys(l) zs(m)];
        end
    end
end

[~,Psi,~]=intersect(Positions,combs,'rows');
data={cores(Psi),satellites(Psi),Positions(Psi,:)};