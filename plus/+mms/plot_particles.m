% Example time interval
tint = EpochTT(irf_time(toepoch([2015 08 15 12 30 00]),'epoch>utc'));
tint2 = EpochTT(irf_time(toepoch([2015 08 15 12 30 00;2015 08 15 12 33 00])','epoch>utc'));

%% Load data
[desSkyMap,dobj] = cn_get_ts('mms1_fpi_brst_l1b_des-dist','mms1_des_brstSkyMap_dist',tint(1));
[desEnergy,dobj] = cn_get_ts('mms1_fpi_brst_l1b_des-dist','mms1_des_energy_index',tint(1));
% desSkyMap: [time energy angle angle] ?

%% Set up coordinates for skymap plot
r = 1;
phi_edges = linspace(0,2*pi,size(desSkyMap.data,3)+1);
theta_edges = linspace(0,pi,size(desSkyMap.data,4)+1);
[PHI,THETA] = meshgrid(phi_edges,theta_edges);
X = r*sin(THETA).*cos(PHI);
Y = r*sin(THETA).*sin(PHI);
Z = r*cos(THETA);

%% Plot one energy level skymap mean of all times
energyLevel = 10;
C = squeeze(nanmean(desSkyMap.data(:,energyLevel,:,:),1))';
%C = squeeze(nanmean(desSkyMap.data(:,:,energyLevel,:),1))';
%C = squeeze(nanmean(desSkyMap.data(:,:,:,energyLevel),1))';
hs = surf(X,Y,Z,C);
axis square
hc = colorbar;
hca=gca;
hca.XLabel.String = 'X';
hca.YLabel.String = 'Y';
hca.ZLabel.String = 'Z';

%% Plot a few energy levels on a skymap
nrows = 3;
ncols = 4;
nEnergyLevels = size(desSkyMap.data,2);
levelStep = ceil(nEnergyLevels/nrows/ncols);
nPlots = ceil(nEnergyLevels/levelStep);
CLims = [];
for k = 1:nPlots
    energyLevel = 1+(k-1)*levelStep;
    hca = subplot(nrows,ncols,k);
    h(k) = hca;
    C = squeeze(nanmean(desSkyMap.data(100:150,energyLevel,:,:),1))';
    hs = surf(hca,X,Y,Z,C);
    hc = colorbar('peer',hca);
    hca.XLabel.String = 'X';
    hca.YLabel.String = 'Y';
    hca.ZLabel.String = 'Z';
    hca.Title.String = ['Energy level = ' num2str(energyLevel)];
    %view([-1 1 1])
    CLims = [CLims; hca.CLim];
    shading flat;
end

%% Plot a few energy levels on a 'flat skymap'
nrows = 3;
ncols = 4;
nEnergyLevels = size(desSkyMap.data,2);
maxEnergyLevel = 23;
levelStep = ceil((maxEnergyLevel)/nrows/ncols);
nPlots = floor(maxEnergyLevel/levelStep);
CLims = [];
for k = 1:nPlots
    energyLevel = 1+(k-1)*levelStep;
    hca = subplot(nrows,ncols,k);
    h(k) = hca;
    C = squeeze(nanmean(desSkyMap.data(:,energyLevel,:,:),1))';
    phi = 1:size(desSkyMap.data,3);
    theta = 1:size(desSkyMap.data,4);
    %hs = pcolor(hca,theta*180/pi,phi*180/pi,C');
    hs = surf(hca,PHI*180/pi,THETA*180/pi,THETA*0,C);
    view([0 0 1])
    hca.YLim = [0 pi]*180/pi;
    hca.XLim = [0 2*pi]*180/pi;
    hca.XTick = [0:60:360];
    hca.YTick = [0:30:180];
    hca.Box = 'on';
    hca.XLabel.String = 'Azimuthal angle (phi)';
    hca.YLabel.String = 'Polar angle (theta)';    
    hc = colorbar('peer',hca);
    hca.Title.String = ['Energy level = ' num2str(energyLevel)];
    %view([-1 1 1])
    CLims = [CLims; hca.CLim];
    shading flat;
end

%% Plot a timeseries of one energylevel sum over phi
h = irf_plot(3);

energyLevel = 10;
skyMap = squeeze(desSkyMap.data(:,energyLevel,:,:)); size(skyMap)
skyMap = squeeze(nansum(skyMap,2)); size(skyMap)
%skyMap = nanmean(skyMap,3); size(skyMap)
CC = [desSkyMap.time.epochUnix skyMap];
specrec.t = desSkyMap.time.epochUnix;
specrec.p = skyMap;
specrec.p_label = 'Electron data';
specrec.f_label = {'Polar','angle'};
specrec.f = 1:16;
irf_spectrogram(h(1),specrec)

energyLevel = 10;
skyMap = squeeze(desSkyMap.data(:,energyLevel,:,:)); size(skyMap)
skyMap = squeeze(nansum(skyMap,3)); size(skyMap)
%skyMap = nanmean(skyMap,3); size(skyMap)
CC = [desSkyMap.time.epochUnix skyMap];
specrec.t = desSkyMap.time.epochUnix;
specrec.p = skyMap;
specrec.p_label = 'Electron data';
specrec.f_label = {'Azimuthal','angle'};
specrec.f = 1:size(skyMap,2);
irf_spectrogram(h(2),specrec)

skyMap = desSkyMap.data;
skyMap = squeeze(nansum(skyMap,3)); size(skyMap)
skyMap = squeeze(nansum(skyMap,2)); size(skyMap)
%skyMap = nanmean(skyMap,3); size(skyMap)
CC = [desSkyMap.time.epochUnix skyMap];
specrec.t = desSkyMap.time.epochUnix;
specrec.p = skyMap;
specrec.p_label = 'Electron data';
specrec.f_label = {'Energy','level'};
specrec.f = 1:size(skyMap,2);
irf_spectrogram(h(3),specrec)





