function data = tred_find_cores(xlims,ylims,zlims)
% TRED_FIND_CORES(xlim,ylim,zlim)
%   
%   Returns 1x3 cell array with 
%   cores - data{1}
%   satellites - data{2}
%   positions - data{3}
%
    
%% Find position
fid=fopen('/Users/Cecilia/TRED46/DATA/Positions.txt');
%hej=textscan(fid,['%*s %*s %*s ' repmat('%f %f %f',1,27) ' \n'],'headerlines',0);
hej=textscan(fid,['%*s %s %*s ' repmat('%f ',1,27*3)],'headerlines',0,'Whitespace',' \b\t\n\r','CollectOutput',true);
%core=textscan(sprintf,'VirtualSatelliteTraces%f.txt')
x=unique(hej{2}(:,1:3:end));
y=unique(hej{2}(:,2:3:end));
z=unique(hej{2}(:,3:3:end));
Ps=[reshape(hej{2}(:,1:3:end)',1,27*size(hej{2},1))',...
    reshape(hej{2}(:,2:3:end)',1,27*size(hej{2},1))',... 
    reshape(hej{2}(:,3:3:end)',1,27*size(hej{2},1))'];
for k=1:length(hej{1}) % collect cores
    core(k)=sscanf(char(hej{1}(k)),'VirtualSatelliteTraces%f.txt');
end
cores=reshape(repmat(core',1,27)',27*length(core),1);

fclose(fid);

%% Find core of positions between limits
xs=intersect(x(find(x>xlims(1))),x(find(x<xlims(2))));
ys=intersect(y(find(y>ylims(1))),y(find(y<ylims(2))));
zs=intersect(z(find(z>zlims(1))),z(find(z<zlims(2))));

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

[~,Psi,~]=intersect(Ps,combs,'rows');
sc=repmat(tocolumn(1:27),length(core),1);

%% Return data
data={cores(Psi), sc(Psi),Ps(Psi,:)};