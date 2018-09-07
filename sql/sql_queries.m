%% Perform queries
query = ['select X,Y,Z from MMS1a limit 20'];
 
[X,Y,Z] = mysql(query); %SymH assumed equal to DST
clear query;

%% Load all locations
query = ['select X,Y,Z from MMS1a'];
[X,Y,Z] = mysql(query); %SymH assumed equal to DST
clear query;
plot3(X,Y,Z,'o')
axis equal

%% Plot all Bmin at all locations
query = ['select EventId,DateStart,X,Y,Z,Bmin from MMS1a'];
[EventId,DateStart,X,Y,Z,Bmin] = mysql(query); %SymH assumed equal to DST
clear query;
ind = 1:numel(X);
ind = find(Bmin<1);
scatter3(X(ind),Y(ind),Z(ind),20,Bmin(ind))
axis equal
colorbar()

%% Plot Bshear Bmin at all locations
query = ['select EventId,DateStart,B1x,B1y,B1z,B2x,B2y,B2z,X,Y,Z,Bmin from MMS1']; % where Y<2 and Y>-2
[EventId,DateStart,B1x,B1y,B1z,B2x,B2y,B2z,X,Y,Z,Bmin] = mysql(query); %SymH assumed equal to DST
clear query;
B1abs = sqrt(B1x.^2+B1y.^2+B1z.^2);
B2abs = sqrt(B2x.^2+B2y.^2+B2z.^2);
Bshear = acosd((B1x.*B2x+B1y.*B2y+B1z.*B2z)./B1abs./B2abs);
dBx = B2x-B1x;
dBy = B2y-B1y;
dBz = B2z-B1z;

% plot
figure(41)
nrows = 2;
ncols = 5;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;
if 1
  hca = h(isub); isub = isub + 1;
  scatter(hca,B1y,B1z);
  hca.XLabel.String = 'B_{1y}';
  hca.YLabel.String = 'B_{1z}';
end
if 1
  hca = h(isub); isub = isub + 1;
  edgesB1y = linspace(min(B1y),max(B1y),21);
  edgesB1z = linspace(min(B1z),max(B1z),20);
  [count edges mid loc] = histcn([B1y,B1z], edgesB1y, edgesB1z);
  pcolor(hca,mid{1:2},count');
  hca.XLabel.String = 'B_{1y}';
  hca.YLabel.String = 'B_{1z}';
end
if 1
  hca = h(isub); isub = isub + 1;
  scatter(hca,B2y,B2z);
  hca.XLabel.String = 'B_{2y}';
  hca.YLabel.String = 'B_{2z}';    
end
if 1
  hca = h(isub); isub = isub + 1;
  edgesB2y = linspace(min(B2y),max(B2y),21);
  edgesB2z = linspace(min(B2z),max(B2z),20);
  [count edges mid loc] = histcn([B2y,B2z], edgesB2y, edgesB2z);
  pcolor(hca,mid{1:2},count');
  hca.XLabel.String = 'B_{2y}';
  hca.YLabel.String = 'B_{2z}';
end
if 1
  hca = h(isub); isub = isub + 1;
  scatter(hca,B2y-B1y,B2z-B1z);
  hca.XLabel.String = 'B_{2y}-B_{1y}';
  hca.YLabel.String = 'B_{2z}-B_{1z}';    
end
if 1
  hca = h(isub); isub = isub + 1;  
  edgesBy = linspace(min(B2y-B1y),max(B2y-B1y),10);
  edgesBz = linspace(min(B2z-B1z),max(B2z-B1z),11);
  edgesBy = linspace(-50,50,20);
  edgesBz = linspace(-50,50,21);
  [count edges mid loc] = histcn([B2y-B1y,B2z-B1z], edgesBy, edgesBz);
  pcolor(hca,mid{1:2},(count)');
  hca.XLabel.String = 'B_{2y}-B_{1y}';
  hca.YLabel.String = 'B_{2z}-B_{1z}';    
end

  
if 1
  hca = h(isub); isub = isub + 1;
  scatter(hca,Bmin,Bshear); 
  hca.XLabel.String = 'B_{min}';
  hca.YLabel.String = 'B_{shear}';    
end
if 1 % Bmin Bshear
  hca = h(isub); isub = isub + 1;
  edgesBy = linspace(min(Bmin),max(Bmin),21);
  edgesBz = linspace(min(Bshear),max(Bshear),20);
  [count edges mid loc] = histcn([Bmin,Bshear], edgesBy, edgesBz);
  pcolor(hca,mid{1:2},count');
  hca.XLabel.String = 'B_{min}';
  hca.YLabel.String = 'B_{shear}';    
end
if 1 % scatter Bmin, Bshear
  hca = h(isub); isub = isub + 1;
  scatter(hca,B1z,B2z,Bmin*10,Bshear); 
  hca.XLabel.String = 'B_{1z}';
  hca.YLabel.String = 'B_{2z}';
  hca.Title.String = 'B_{shear}';
  colorbar('peer',hca)
end
if 1 % scatter Bmin, Bshear
  hca = h(isub); isub = isub + 1;  
  ind = find(Bshear>170);
  quiver(hca,B1y(ind),B1z(ind),dBy(ind),dBz(ind),10)
  hca.XLabel.String = 'B_{y}';
  hca.YLabel.String = 'B_{z}'; 
  hca.Title.String = 'B_{2}-B_{1}';   
end



%% Get B profile
query = ['select B_profile from MMS1_profiles'];
[B_profile] = mysql(query); %SymH assumed equal to DST
clear query;

