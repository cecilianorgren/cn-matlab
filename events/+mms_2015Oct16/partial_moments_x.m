% Choose spacecraft
ic = 1;

%% Time: 2015-10-16T10:33:30.266567000Z, dist = ePDist1(1209).time;
timeUTC = '2015-10-16T10:33:30.266567000Z';
time = irf_time(timeUTC,'utc>epochTT');

limPitchangleB = 90+40*[-1 1]; % degree
limPitchangleExB = [0 110]; % degree
limE = [60 40000];

% Pick out local values 
c_eval('dist = ePDist?.tlim(time+[-0.01  0.01]);',ic)
c_eval('scpot = scPot?.resample(time);',ic)
c_eval('n = ne?.resample(time);',ic)
c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
hatVe0 = double(irf_norm(Ve0));    
c_eval('Ve0perp = gseVe?perp.resample(time).data;',ic); 
c_eval('Ve0par = gseVe?par.resample(time).data;',ic); 
c_eval('dmpaB0 = dmpaB?.resample(time);',ic); 
hatB0 = double(irf_norm(dmpaB0.data));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Define coordinate system to plot in
x = hatB0;
y = hatExB0;
z = cross(x,y);
vlabels = {'B','ExB','Bx(ExB)'};

% Set up some projection parameters
vlim = 15*1e3;
elevlim = 30;
strCMap = 'jet';
projclim = [0 5];

% Plot projection in three planes
% Each row will have different components included: total, core, crescent 
% Initialize plot
nRows = 3; nCols = 3;
clear h;
isub = 1;
for ii = 1:nCols*nRows; h(isub) = subplot(nRows,nCols,ii); isub = isub + 1; end

vectors = {hatB0,'B';hatExB0,'ExB'};
vectors = {hatVe0,'v_0'};
isub = 1;
% Plot total distributions
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];

% Pick put partial distributions
energy =  dist.depend{1}; % energy
azimuthal = dist.depend{2}; % azimuthal angle
polar = dist.depend{3}; % polar angle
[EN,AZ,POL] = meshgrid(energy,azimuthal,polar);
EN = permute(EN,[2 1 3]);
AZ = permute(AZ,[2 1 3]);
POL = permute(POL,[2 1 3]);
dX = -sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
dY = -sind(POL).*sind(AZ);
dZ = -cosd(POL);
dataSize = size(dX);
xX = reshape(dX,prod(dataSize),1);
yY = reshape(dY,prod(dataSize),1);
zZ = reshape(dZ,prod(dataSize),1);

tmpPitchangleB = acosd([xX yY zZ]*hatB0');
tmpPitchangleExB = acosd([xX yY zZ]*hatExB0');
ePitchangleB = reshape(tmpPitchangleB,dataSize);
ePitchangleExB = reshape(tmpPitchangleExB,dataSize);
% surf(dX,dY,dZ,ePitchangleB)
% surf(dX,dY,dZ,ePitchangleExB)



%indB = [find(ePitchangleB>limPitchangleB(1)); find(ePitchangleB<limPitchangleB(2))];
matIncl = (ePitchangleB<limPitchangleB(2) & ...
           ePitchangleB>limPitchangleB(1) & ...
           ePitchangleExB>limPitchangleExB(1) & ...
           ePitchangleExB<limPitchangleExB(2) & ...
           EN>limE(1) & ...
           EN<limE(2));
tsMatIncl = irf.ts_scalar(time,single(reshape(matIncl,[1 size(matIncl)])));
tsMatExcl = irf.ts_scalar(time,single(not(reshape(matIncl,[1 size(matIncl)]))));

if 0 % Attempt to 3D visualize remaining distrbution
  X3D = -EN.*sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
  Y3D = -EN.*sind(POL).*sind(AZ);
  Z3D = -EN.*cosd(POL);
  X3D(matIncl) = NaN; X3D = reshape(X3D,1,prod(size(X3D)));
  Y3D(matIncl) = NaN; Y3D = reshape(Y3D,1,prod(size(Y3D)));
  Z3D(matIncl) = NaN; Z3D = reshape(Z3D,1,prod(size(Z3D)));
end

% Calculate moments
stepTable = TSeries(time,dist.ancillary.energyStepTable);
momsCrescent = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatIncl);
momsCore = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatExcl);

% Decompose velocities into parallel and perpendicular components
[dmpaVparCore,dmpaVperpCore] = irf_dec_parperp(dmpaB0,momsCrescent.V_psd);
[dmpaVparCrescent,dmpaVperpCrescent] = irf_dec_parperp(dmpaB0,momsCore.V_psd);

% DMPA LMN
%[out,l,v] = irf_minvar(dmpaB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
%dmpaL = v(1,:); dmpaM = v(2,:); dmpaN = v(3,:);


% Plot partial distributions
distCrescent = dist; distCrescent.data(find(~matIncl)) = NaN; 
distCore = dist; distCore.data(find((matIncl))) = NaN;

vectors = {irf_norm(momsCrescent.V_psd.data),'v_e';hatVe0,'v_0'};
%isub = 4;
% Crescent
nCrescent = momsCrescent.n_psd.data;
strTitle = ['n_e = ' num2str(nCrescent,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCrescent.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCrescent.abs.data,'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = strTitle;

% High pitch angles
hca = h(isub); isub = isub + 1; 
strTitle = ['n_e = ' num2str(momsCore.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCore.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCore.abs.data,'%.0f') ' km/s'];;
mms.plot_projection(hca,distCore,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]));       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]));       
hca.Title.String = strTitle;


% Fix figure
for ii = 1:nRows*nCols
  colormap(h(ii),strCMap);
  %h(ii).Title.String = ''; 
end
h(2).Title.String = {timeUTC(1:23),' ',h(2).Title.String};
%% Time: 2015-10-16T10:33:30.296567000Z, dist = ePDist1(1210).time;
timeUTC = '2015-10-16T10:33:30.296567000Z';
time = irf_time(timeUTC,'utc>epochTT');

limPitchangleB = 90+30*[-1 1]; % degree
limPitchangleExB = [0 100]; % degree
limE = [60 40000];

% Pick out local values 
c_eval('dist = ePDist?.tlim(time+[-0.01  0.01]);',ic)
c_eval('scpot = scPot?.resample(time);',ic)
c_eval('n = ne?.resample(time);',ic)
c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
hatVe0 = double(irf_norm(Ve0));    
c_eval('Ve0perp = gseVe?perp.resample(time).data;',ic); 
c_eval('Ve0par = gseVe?par.resample(time).data;',ic); 
c_eval('dmpaB0 = dmpaB?.resample(time);',ic); 
hatB0 = double(irf_norm(dmpaB0.data));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Define coordinate system to plot in
x = hatB0;
y = hatExB0;
z = cross(x,y);
vlabels = {'B','ExB','Bx(ExB)'};

% Set up some projection parameters
vlim = 15*1e3;
elevlim = 30;
strCMap = 'jet';
projclim = [0 5];

% Plot projection in three planes
% Each row will have different components included: total, core, crescent 
% Initialize plot
nRows = 3; nCols = 3;
clear h;
isub = 1;
for ii = 1:nCols*nRows; h(isub) = subplot(nRows,nCols,ii); isub = isub + 1; end

vectors = {hatB0,'B';hatExB0,'ExB'};
vectors = {hatVe0,'v_0'};
isub = 1;
% Plot total distributions
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];

% Pick put partial distributions
energy =  dist.depend{1}; % energy
azimuthal = dist.depend{2}; % azimuthal angle
polar = dist.depend{3}; % polar angle
[EN,AZ,POL] = meshgrid(energy,azimuthal,polar);
EN = permute(EN,[2 1 3]);
AZ = permute(AZ,[2 1 3]);
POL = permute(POL,[2 1 3]);
dX = -sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
dY = -sind(POL).*sind(AZ);
dZ = -cosd(POL);
dataSize = size(dX);
xX = reshape(dX,prod(dataSize),1);
yY = reshape(dY,prod(dataSize),1);
zZ = reshape(dZ,prod(dataSize),1);

tmpPitchangleB = acosd([xX yY zZ]*hatB0');
tmpPitchangleExB = acosd([xX yY zZ]*hatExB0');
ePitchangleB = reshape(tmpPitchangleB,dataSize);
ePitchangleExB = reshape(tmpPitchangleExB,dataSize);
% surf(dX,dY,dZ,ePitchangleB)
% surf(dX,dY,dZ,ePitchangleExB)

%indB = [find(ePitchangleB>limPitchangleB(1)); find(ePitchangleB<limPitchangleB(2))];
matIncl = (ePitchangleB<limPitchangleB(2) & ...
           ePitchangleB>limPitchangleB(1) & ...
           ePitchangleExB>limPitchangleExB(1) & ...
           ePitchangleExB<limPitchangleExB(2) & ...
           EN>limE(1) & ...
           EN<limE(2));
tsMatIncl = irf.ts_scalar(time,single(reshape(matIncl,[1 size(matIncl)])));
tsMatExcl = irf.ts_scalar(time,single(not(reshape(matIncl,[1 size(matIncl)]))));

if 0 % Attempt to 3D visualize remaining distrbution
  X3D = -EN.*sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
  Y3D = -EN.*sind(POL).*sind(AZ);
  Z3D = -EN.*cosd(POL);
  X3D(matIncl) = NaN; X3D = reshape(X3D,1,prod(size(X3D)));
  Y3D(matIncl) = NaN; Y3D = reshape(Y3D,1,prod(size(Y3D)));
  Z3D(matIncl) = NaN; Z3D = reshape(Z3D,1,prod(size(Z3D)));
end

% Calculate moments
stepTable = TSeries(time,dist.ancillary.energyStepTable);
momsCrescent = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatIncl);
momsCore = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatExcl);

% Decompose velocities into parallel and perpendicular components
[dmpaVparCore,dmpaVperpCore] = irf_dec_parperp(dmpaB0,momsCrescent.V_psd);
[dmpaVparCrescent,dmpaVperpCrescent] = irf_dec_parperp(dmpaB0,momsCore.V_psd);

% DMPA LMN
%[out,l,v] = irf_minvar(dmpaB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
%dmpaL = v(1,:); dmpaM = v(2,:); dmpaN = v(3,:);


% Plot partial distributions
distCrescent = dist; distCrescent.data(find(matIncl)) = NaN; 
distCore = dist; distCore.data(find((~matIncl))) = NaN;

vectors = {irf_norm(momsCrescent.V_psd.data),'v_e';hatVe0,'v_0'};
%isub = 4;
% Crescent
nCrescent = momsCrescent.n_psd.data;
strTitle = ['n_e = ' num2str(nCrescent,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCrescent.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCrescent.abs.data,'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = strTitle;

% High pitch angles
hca = h(isub); isub = isub + 1; 
strTitle = ['n_e = ' num2str(momsCore.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCore.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCore.abs.data,'%.0f') ' km/s'];;
mms.plot_projection(hca,distCore,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]));       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]));       
hca.Title.String = strTitle;


% Fix figure
for ii = 1:nRows*nCols
  colormap(h(ii),strCMap);
  %h(ii).Title.String = ''; 
end
h(2).Title.String = {timeUTC(1:23),' ',h(2).Title.String};
%% Time: 2015-10-16T10:33:30.386567000Z, dist = ePDist1(1211).time;
timeUTC = '2015-10-16T10:33:30.386567000Z';
time = irf_time(timeUTC,'utc>epochTT');

% Pick out local values 
c_eval('dist = ePDist?.tlim(time+[-0.01  0.01]);',ic)
c_eval('scpot = scPot?.resample(time);',ic)
c_eval('n = ne?.resample(time);',ic)
c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
hatVe0 = double(irf_norm(Ve0));    
c_eval('Ve0perp = gseVe?perp.resample(time).data;',ic); 
c_eval('Ve0par = gseVe?par.resample(time).data;',ic); 
c_eval('dmpaB0 = dmpaB?.resample(time);',ic); 
hatB0 = double(irf_norm(dmpaB0.data));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Define coordinate system to plot in
x = hatB0;
y = hatExB0;
z = cross(x,y);
vlabels = {'B','ExB','Bx(ExB)'};

% Set up some projection parameters
vlim = 15*1e3;
elevlim = 30;
strCMap = 'jet';
projclim = [0 5];

% Plot projection in three planes
% Each row will have different components included: total, core, crescent 
% Initialize plot
nRows = 3; nCols = 3;
clear h;
isub = 1;
for ii = 1:nCols*nRows; h(isub) = subplot(nRows,nCols,ii); isub = isub + 1; end

vectors = {hatB0,'B';hatExB0,'ExB'};
vectors = {hatVe0,'v_0'};
isub = 1;
% Plot total distributions
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];

% Pick put partial distributions
energy =  dist.depend{1}; % energy
azimuthal = dist.depend{2}; % azimuthal angle
polar = dist.depend{3}; % polar angle
[EN,AZ,POL] = meshgrid(energy,azimuthal,polar);
EN = permute(EN,[2 1 3]);
AZ = permute(AZ,[2 1 3]);
POL = permute(POL,[2 1 3]);
dX = -sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
dY = -sind(POL).*sind(AZ);
dZ = -cosd(POL);
dataSize = size(dX);
xX = reshape(dX,prod(dataSize),1);
yY = reshape(dY,prod(dataSize),1);
zZ = reshape(dZ,prod(dataSize),1);

tmpPitchangleB = acosd([xX yY zZ]*hatB0');
tmpPitchangleExB = acosd([xX yY zZ]*hatExB0');
ePitchangleB = reshape(tmpPitchangleB,dataSize);
ePitchangleExB = reshape(tmpPitchangleExB,dataSize);
% surf(dX,dY,dZ,ePitchangleB)
% surf(dX,dY,dZ,ePitchangleExB)

limPitchangleB = 90+30*[-1 1]; % degree
limPitchangleExB = [0 110]; % degree
limE = [60 40000];

%indB = [find(ePitchangleB>limPitchangleB(1)); find(ePitchangleB<limPitchangleB(2))];
matIncl = (ePitchangleB<limPitchangleB(2) & ...
           ePitchangleB>limPitchangleB(1) & ...
           ePitchangleExB>limPitchangleExB(1) & ...
           ePitchangleExB<limPitchangleExB(2) & ...
           EN>limE(1) & ...
           EN<limE(2));
tsMatIncl = irf.ts_scalar(time,single(reshape(matIncl,[1 size(matIncl)])));
tsMatExcl = irf.ts_scalar(time,single(not(reshape(matIncl,[1 size(matIncl)]))));

if 0 % Attempt to 3D visualize remaining distrbution
  X3D = -EN.*sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
  Y3D = -EN.*sind(POL).*sind(AZ);
  Z3D = -EN.*cosd(POL);
  X3D(matIncl) = NaN; X3D = reshape(X3D,1,prod(size(X3D)));
  Y3D(matIncl) = NaN; Y3D = reshape(Y3D,1,prod(size(Y3D)));
  Z3D(matIncl) = NaN; Z3D = reshape(Z3D,1,prod(size(Z3D)));
end

% Calculate moments
stepTable = TSeries(time,dist.ancillary.energyStepTable);
momsCrescent = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatIncl);
momsCore = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatExcl);

% Decompose velocities into parallel and perpendicular components
[dmpaVparCore,dmpaVperpCore] = irf_dec_parperp(dmpaB0,momsCrescent.V_psd);
[dmpaVparCrescent,dmpaVperpCrescent] = irf_dec_parperp(dmpaB0,momsCore.V_psd);

% DMPA LMN
%[out,l,v] = irf_minvar(dmpaB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
%dmpaL = v(1,:); dmpaM = v(2,:); dmpaN = v(3,:);


% Plot partial distributions
distCrescent = dist; distCrescent.data(find(~matIncl)) = NaN; 
distCore = dist; distCore.data(find((matIncl))) = NaN;

vectors = {irf_norm(momsCrescent.V_psd.data),'v_e';hatVe0,'v_0'};
%isub = 4;
% Crescent
nCrescent = momsCrescent.n_psd.data;
strTitle = ['n_e = ' num2str(nCrescent,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCrescent.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCrescent.abs.data,'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = strTitle;

% High pitch angles
hca = h(isub); isub = isub + 1; 
strTitle = ['n_e = ' num2str(momsCore.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCore.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCore.abs.data,'%.0f') ' km/s'];;
mms.plot_projection(hca,distCore,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]));       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]));       
hca.Title.String = strTitle;


% Fix figure
for ii = 1:nRows*nCols
  colormap(h(ii),strCMap);
  %h(ii).Title.String = ''; 
end
h(2).Title.String = {timeUTC(1:23),' ',h(2).Title.String};
%% Time: 2015-10-16T10:33:30.356567000Z, dist = ePDist1(1212).time;
timeUTC = '2015-10-16T10:33:30.356567000Z';
time = irf_time(timeUTC,'utc>epochTT');

% Pick out local values 
c_eval('dist = ePDist?.tlim(time+[-0.01  0.01]);',ic)
c_eval('scpot = scPot?.resample(time);',ic)
c_eval('n = ne?.resample(time);',ic)
c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
hatVe0 = double(irf_norm(Ve0));    
c_eval('Ve0perp = gseVe?perp.resample(time).data;',ic); 
c_eval('Ve0par = gseVe?par.resample(time).data;',ic); 
c_eval('dmpaB0 = dmpaB?.resample(time);',ic); 
hatB0 = double(irf_norm(dmpaB0.data));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Define coordinate system to plot in
x = hatB0;
y = hatExB0;
z = cross(x,y);
vlabels = {'B','ExB','Bx(ExB)'};

% Set up some projection parameters
vlim = 15*1e3;
elevlim = 30;
strCMap = 'jet';
projclim = [0 5];

% Plot projection in three planes
% Each row will have different components included: total, core, crescent 
% Initialize plot
nRows = 3; nCols = 3;
clear h;
isub = 1;
for ii = 1:nCols*nRows; h(isub) = subplot(nRows,nCols,ii); isub = isub + 1; end

vectors = {hatB0,'B';hatExB0,'ExB'};
vectors = {hatVe0,'v_0'};
isub = 1;
% Plot total distributions
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];

% Pick put partial distributions
energy =  dist.depend{1}; % energy
azimuthal = dist.depend{2}; % azimuthal angle
polar = dist.depend{3}; % polar angle
[EN,AZ,POL] = meshgrid(energy,azimuthal,polar);
EN = permute(EN,[2 1 3]);
AZ = permute(AZ,[2 1 3]);
POL = permute(POL,[2 1 3]);
dX = -sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
dY = -sind(POL).*sind(AZ);
dZ = -cosd(POL);
dataSize = size(dX);
xX = reshape(dX,prod(dataSize),1);
yY = reshape(dY,prod(dataSize),1);
zZ = reshape(dZ,prod(dataSize),1);

tmpPitchangleB = acosd([xX yY zZ]*hatB0');
tmpPitchangleExB = acosd([xX yY zZ]*hatExB0');
ePitchangleB = reshape(tmpPitchangleB,dataSize);
ePitchangleExB = reshape(tmpPitchangleExB,dataSize);
% surf(dX,dY,dZ,ePitchangleB)
% surf(dX,dY,dZ,ePitchangleExB)

limPitchangleB = 90+30*[-1 1]; % degree
limPitchangleExB = [0 110]; % degree
limE = [60 40000];

%indB = [find(ePitchangleB>limPitchangleB(1)); find(ePitchangleB<limPitchangleB(2))];
matIncl = (ePitchangleB<limPitchangleB(2) & ...
           ePitchangleB>limPitchangleB(1) & ...
           ePitchangleExB>limPitchangleExB(1) & ...
           ePitchangleExB<limPitchangleExB(2) & ...
           EN>limE(1) & ...
           EN<limE(2));
tsMatIncl = irf.ts_scalar(time,single(reshape(matIncl,[1 size(matIncl)])));
tsMatExcl = irf.ts_scalar(time,single(not(reshape(matIncl,[1 size(matIncl)]))));

if 0 % Attempt to 3D visualize remaining distrbution
  X3D = -EN.*sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
  Y3D = -EN.*sind(POL).*sind(AZ);
  Z3D = -EN.*cosd(POL);
  X3D(matIncl) = NaN; X3D = reshape(X3D,1,prod(size(X3D)));
  Y3D(matIncl) = NaN; Y3D = reshape(Y3D,1,prod(size(Y3D)));
  Z3D(matIncl) = NaN; Z3D = reshape(Z3D,1,prod(size(Z3D)));
end

% Calculate moments
stepTable = TSeries(time,dist.ancillary.energyStepTable);
momsCrescent = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatIncl);
momsCore = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatExcl);

% Decompose velocities into parallel and perpendicular components
[dmpaVparCore,dmpaVperpCore] = irf_dec_parperp(dmpaB0,momsCrescent.V_psd);
[dmpaVparCrescent,dmpaVperpCrescent] = irf_dec_parperp(dmpaB0,momsCore.V_psd);

% DMPA LMN
%[out,l,v] = irf_minvar(dmpaB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
%dmpaL = v(1,:); dmpaM = v(2,:); dmpaN = v(3,:);


% Plot partial distributions
distCrescent = dist; distCrescent.data(find(~matIncl)) = NaN; 
distCore = dist; distCore.data(find((matIncl))) = NaN;

vectors = {irf_norm(momsCrescent.V_psd.data),'v_e';hatVe0,'v_0'};
%isub = 4;
% Crescent
nCrescent = momsCrescent.n_psd.data;
strTitle = ['n_e = ' num2str(nCrescent,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCrescent.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCrescent.abs.data,'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = strTitle;

% High pitch angles
hca = h(isub); isub = isub + 1; 
strTitle = ['n_e = ' num2str(momsCore.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCore.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCore.abs.data,'%.0f') ' km/s'];;
mms.plot_projection(hca,distCore,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]));       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]));       
hca.Title.String = strTitle;


% Fix figure
for ii = 1:nRows*nCols
  colormap(h(ii),strCMap);
  %h(ii).Title.String = ''; 
end
h(2).Title.String = {timeUTC(1:23),' ',h(2).Title.String};
%% Time: 2015-10-16T10:33:30.566567000Z, dist = ePDist1(1219).time
timeUTC = '2015-10-16T10:33:30.566567000Z';
time = irf_time(timeUTC,'utc>epochTT');

% Pick out local values 
c_eval('dist = ePDist?.tlim(time+[-0.01  0.01]);',ic)
c_eval('scpot = scPot?.resample(time);',ic)
c_eval('n = ne?.resample(time);',ic)
c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
hatVe0 = double(irf_norm(Ve0));    
c_eval('Ve0perp = gseVe?perp.resample(time).data;',ic); 
c_eval('Ve0par = gseVe?par.resample(time).data;',ic); 
c_eval('dmpaB0 = dmpaB?.resample(time);',ic); 
hatB0 = double(irf_norm(dmpaB0.data));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

% Define coordinate system to plot in
x = hatB0;
y = hatExB0;
z = cross(x,y);
vlabels = {'B','ExB','Bx(ExB)'};

% Set up some projection parameters
vlim = 15*1e3;
elevlim = 30;
strCMap = 'jet';
projclim = [0 5];

% Plot projection in three planes
% Each row will have different components included: total, core, crescent 
% Initialize plot
nRows = 3; nCols = 3;
clear h;
isub = 1;
for ii = 1:nCols*nRows; h(isub) = subplot(nRows,nCols,ii); isub = isub + 1; end

vectors = {hatB0,'B';hatExB0,'ExB'};
vectors = {hatVe0,'v_0'};
isub = 1;
% Plot total distributions
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,dist,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = ['n_e = ' num2str(n.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(norm(Ve0perp),'%.0f') ' km/s, v_{e,||} = ' num2str(norm(Ve0par),'%.0f') ' km/s'];

% Pick put partial distributions
energy =  dist.depend{1}; % energy
azimuthal = dist.depend{2}; % azimuthal angle
polar = dist.depend{3}; % polar angle
[EN,AZ,POL] = meshgrid(energy,azimuthal,polar);
EN = permute(EN,[2 1 3]);
AZ = permute(AZ,[2 1 3]);
POL = permute(POL,[2 1 3]);
dX = -sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
dY = -sind(POL).*sind(AZ);
dZ = -cosd(POL);
dataSize = size(dX);
xX = reshape(dX,prod(dataSize),1);
yY = reshape(dY,prod(dataSize),1);
zZ = reshape(dZ,prod(dataSize),1);

tmpPitchangleB = acosd([xX yY zZ]*hatB0');
tmpPitchangleExB = acosd([xX yY zZ]*hatExB0');
ePitchangleB = reshape(tmpPitchangleB,dataSize);
ePitchangleExB = reshape(tmpPitchangleExB,dataSize);
% surf(dX,dY,dZ,ePitchangleB)
% surf(dX,dY,dZ,ePitchangleExB)

limPitchangleB = 90+30*[-1 1]; % degree
limPitchangleExB = [0 110]; % degree
limE = [60 40000];

%indB = [find(ePitchangleB>limPitchangleB(1)); find(ePitchangleB<limPitchangleB(2))];
matIncl = (ePitchangleB<limPitchangleB(2) & ...
           ePitchangleB>limPitchangleB(1) & ...
           ePitchangleExB>limPitchangleExB(1) & ...
           ePitchangleExB<limPitchangleExB(2) & ...
           EN>limE(1) & ...
           EN<limE(2));
tsMatIncl = irf.ts_scalar(time,single(reshape(matIncl,[1 size(matIncl)])));
tsMatExcl = irf.ts_scalar(time,single(not(reshape(matIncl,[1 size(matIncl)]))));

if 0 % Attempt to 3D visualize remaining distrbution
  X3D = -EN.*sind(POL).*cosd(AZ); % '-' because the data shows which direction the particles were coming from
  Y3D = -EN.*sind(POL).*sind(AZ);
  Z3D = -EN.*cosd(POL);
  X3D(matIncl) = NaN; X3D = reshape(X3D,1,prod(size(X3D)));
  Y3D(matIncl) = NaN; Y3D = reshape(Y3D,1,prod(size(Y3D)));
  Z3D(matIncl) = NaN; Z3D = reshape(Z3D,1,prod(size(Z3D)));
end

% Calculate moments
stepTable = TSeries(time,dist.ancillary.energyStepTable);
momsCrescent = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatIncl);
momsCore = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatExcl);

% Decompose velocities into parallel and perpendicular components
[dmpaVparCore,dmpaVperpCore] = irf_dec_parperp(dmpaB0,momsCrescent.V_psd);
[dmpaVparCrescent,dmpaVperpCrescent] = irf_dec_parperp(dmpaB0,momsCore.V_psd);

% DMPA LMN
%[out,l,v] = irf_minvar(dmpaB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
%dmpaL = v(1,:); dmpaM = v(2,:); dmpaN = v(3,:);


% Plot partial distributions
distCrescent = dist; distCrescent.data(find(~matIncl)) = NaN; 
distCore = dist; distCore.data(find((matIncl))) = NaN;

vectors = {irf_norm(momsCrescent.V_psd.data),'v_e';hatVe0,'v_0'};
%isub = 4;
% Crescent
nCrescent = momsCrescent.n_psd.data;
strTitle = ['n_e = ' num2str(nCrescent,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCrescent.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCrescent.abs.data,'%.0f') ' km/s'];
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
hca.Title.String = strTitle;

% High pitch angles
hca = h(isub); isub = isub + 1; 
strTitle = ['n_e = ' num2str(momsCore.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCore.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCore.abs.data,'%.0f') ' km/s'];;
mms.plot_projection(hca,distCore,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels);       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]));       
hca.Title.String = strTitle;
hca = h(isub); isub = isub + 1; 
mms.plot_projection(hca,distCore,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]));       
hca.Title.String = strTitle;


% Fix figure
for ii = 1:nRows*nCols
  colormap(h(ii),strCMap);
  %h(ii).Title.String = ''; 
end
h(2).Title.String = {timeUTC(1:23),' ',h(2).Title.String};
