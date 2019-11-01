%% Perform queries
query = ['select X,Y,Z from MMS1a limit 20'];

query = ['select X,Y,Z from MMS3 limit 20'];
 
[X,Y,Z] = mysql(query); %SymH assumed equal to DST
clear query;

%% Load all time intervals,
query = ['select X,Y,Z,Bmin,t_0,t_12,t_88 from MMS3'];
[X,Y,Z,Bmin,t_0,t_12,t_88] = mysql(query); %SymH assumed equal to DST
clear query;

% Make EpochTT list
t_0_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_0,'UniformOutput', false)));
t_12_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_12,'UniformOutput', false)));
t_88_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_88,'UniformOutput', false)));

dt = t_88_epochtt-t_12_epochtt;

hist(dt,0:1:max(dt))

%% Load ve LMN
query = ['select eventId,X,Y,Z,Bmin,t_0,t_12,t_88 from MMS3'];
[eventId,X,Y,Z,Bmin,t_0,t_12,t_88] = mysql(query); %SymH assumed equal to DST
clear query;

% Make EpochTT list
t_0_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_0,'UniformOutput', false)));
t_12_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_12,'UniformOutput', false)));
t_88_epochtt = EpochTT(cell2mat(cellfun(@(x) (cat(2,strrep(x,' ','T'),'Z')),t_88,'UniformOutput', false)));

dt = t_88_epochtt-t_12_epochtt;

hist(dt,0:1:max(dt))
%% Load B profiles, doesnt work
query = ['select EventId,B_profile,E_profile from MMS3_profiles'];
[eventId,B_profile,E_profile] = mysql(query); %SymH assumed equal to DST
clear query;

%% Separate data by flags

query = ['select eventId,X,Y,Z,Bmin,t_0,t_12,t_88,Lx,Ly,Lz,Mx,My,Mz,Nx,Ny,Nz,Flag,FlagStr from MMS4'];
[eventId,X,Y,Z,Bmin,t_0,t_12,t_88,Lx,Ly,Lz,Mx,My,Mz,Nx,Ny,Nz,Flag,FlagStr] = mysql(query); %SymH assumed equal to DST
clear query;

query = ['select eventId,Vex_tmaxJ,Vey_tmaxJ,Vez_tmaxJ from MMS3_elctrns'];
[eventId_elctrns,Vex_tmaxJ,Vey_tmaxJ,Vez_tmaxJ] = mysql(query); %SymH assumed equal to DST
% get lmn velocities

[C,ia,ib] = intersect(eventId,eventId_elctrns);
VeL_tmaxJ = Vex_tmaxJ(ib).*Lx(ia) + Vey_tmaxJ(ib).*Ly(ia) + Vez_tmaxJ(ib).*Lz(ia);
VeM_tmaxJ = Vex_tmaxJ(ib).*Mx(ia) + Vey_tmaxJ(ib).*My(ia) + Vez_tmaxJ(ib).*Mz(ia);
VeN_tmaxJ = Vex_tmaxJ(ib).*Nx(ia) + Vey_tmaxJ(ib).*Ny(ia) + Vez_tmaxJ(ib).*Nz(ia);


%% Load all locations
query = ['select X,Y,Z,Bmin from MMS3 where Bmin<10'];
[X,Y,Z,Bmin] = mysql(query); %SymH assumed equal to DST
clear query;
scatter3(X,Y,Z,Bmin,Bmin)
colorbar

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
query = ['select EventId,DateStart,B1x,B1y,B1z,B2x,B2y,B2z,X,Y,Z,Bmin from MMS3']; % where Y<2 and Y>-2
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

%%
query = ['select MMS1.X,MMS1.Y,MMS1.Z,MMS1.FlagStr from MMS1 where MMS1.FlagStr like(''%msh%'') limit 100'];
[X,Y,Z,FlagStr] = mysql(query); %SymH assumed equal to DST
clear query
%% Load Vimax at Bmin < 10
%query = ['select X,Y,Z,Vmax_tbrst,Bmin from MMS1 join MMS1_ions where MMS1_ions.Bmin<10'];
query = ['select MMS1.X,MMS1.Y,MMS1.Z,MMS1_ions.Vmax_tbrst,MMS1.Bmin from MMS1 join MMS1_ions on MMS1_ions.EventID = MMS1.EventId where MMS1_ions.Vmax_tbrst>100 and MMS1.FlagStr like(''mp%'')'];
[X,Y,Z,Vimax,Bmin] = mysql(query); %SymH assumed equal to DST
clear query;
%%
query = ['select MMS1.X,MMS1.Y,MMS1.Z,MMS1_ions.Vmax_tbrst,MMS1.Bmin from MMS1 join MMS1_ions on MMS1_ions.EventID = MMS1.EventId where MMS1.FlagStr like(''mp%'')'];
[mpX,mpY,mpZ,mpVimax,mpBmin] = mysql(query); %SymH assumed equal to DST
clear query
query = ['select MMS1.X,MMS1.Y,MMS1.Z,MMS1_ions.Vmax_tbrst,MMS1.Bmin from MMS1 join MMS1_ions on MMS1_ions.EventID = MMS1.EventId where MMS1.FlagStr like(''msh%'')'];
[mshX,mshY,mshZ,mshVimax,mshBmin] = mysql(query); %SymH assumed equal to DST
clear query

%% find science event in database
query = ['select MMS1.X,MMS1.Y,MMS1.Z,MMS1_ions.Vmax_tbrst,MMS1.Bmin,MMS1.DateStart from MMS1 join MMS1_ions on MMS1_ions.EventID = MMS1.EventId where MMS1.DateStart like(''%20151016_%'')'];
[X,Y,Z,Vimax,Bmin] = mysql(query); %SymH assumed equal to DST
clear query

%% find science event in database
query = ['select MMS1.X,MMS1.Y,MMS1.Z,MMS1_ions.Vmax_tbrst,MMS1.Bmin,MMS1.DateStart from MMS1 join MMS1_ions on MMS1_ions.EventID = MMS1.EventId where MMS1.DateStart like(''%20151016_%'')'];
[X,Y,Z,Vimax,Bmin] = mysql(query); %SymH assumed equal to DST
clear query


%% plot
colors = mms_colors('matlab');
figure(42)
nrows = 2;
ncols = 2;
npanels = nrows*ncols;
isub = 0;
for icol = 1:ncols
  for irow = 1:nrows  
    isub = isub + 1;         
    h(isub) = subplot(nrows,ncols,icol+(irow-1)*ncols);    
  end
end
isub = 1;
if 0
  hca = h(isub); isub = isub + 1;
  scatter(hca,Bmin,Vimax)
  hca.XLabel.String = 'B_{min}';
  hca.YLabel.String = 'V_{max}';
end
if 1
  hca = h(isub); isub = isub + 1;
  scatter3(hca,mpX,mpY,mpZ,mpBmin,colors(1,:))
  hold(hca,'on')
  scatter3(hca,mshX,mshY,mshZ,mshBmin,colors(2,:))
  hold(hca,'off')
  hca.XLabel.String = 'X';
  hca.YLabel.String = 'Y';
  hca.ZLabel.String = 'Z';
end
if 0
  hca = h(isub); isub = isub + 1;
  scatter3(hca,X,Y,Z,Bmin,Vimax)
  hca.XLabel.String = 'X';
  hca.YLabel.String = 'Y';
  hca.ZLabel.String = 'Z';
end