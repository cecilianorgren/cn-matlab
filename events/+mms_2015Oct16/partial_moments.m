% Choose spacecraft
ic = 1;
doPlot = 0;
doPrint = 0;
clear PCore PCrescent TCore TCrescent vCore vCrescent nCore nCrescent nBG momsMinE momsMaxE
tUTC = [];
momsCrescent = {};
momsCore = {};
momsBG = {};
momsLimE = [];

iInd = 0;
timeInds = 1205:1213;1206:1220;1213;1209;1206:1225;1213;
for indT = timeInds % 1213% 1210;
iInd = iInd + 1;
doVePitch = 1; % Choose 2nd pitch angle from Veperp instead of ExB
% Choose time
indTimes = {1205:1213,...
            [],...
            [],...
            []};
%indT = 1225;
time = ePDist1(indT).time;
timeUTC = time.utc;

switch ic
  case 1
    switch indT
      case 1204
      case 1205
        limPitchangleB = 90+[-40 35]; % degree
        limPitchangleExB = [0 110]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [60 40000];
      case 1206
        limPitchangleB = 90+[-40 35]; % degree
        limPitchangleExB = [0 110]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [60 40000];
      case 1207
        limPitchangleB = 90+[-45 35]; % degree
        limPitchangleExB = [0 110]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [60 40000];
      case 1208
        limPitchangleB = 90+[-40 30]; % degree
        limPitchangleExB = [0 110]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [60 40000];
      case 1209
        limPitchangleB = 90+[-35 25]; % degree
        limPitchangleExB = [130 180]; % degree
        limPitchangleVeperp = [0 130]; % degree
        limE = [60 40000];
        %doVePitch = 0;
      case 1210
        limPitchangleB = 90+30*[-1 1]; % degree
        limPitchangleExB = [100 180]; % degree
        limPitchangleVeperp = [110 180]; % degree
        limE = [60 40000];
      case 1211
        limPitchangleB = 90+25*[-1 1]; % degree
        limPitchangleExB = [0 70]; % degree
        limPitchangleVeperp = [0 70]; % degree
        limE = [60 40000];
      case 1212
        limPitchangleB = 90+20*[-1 1]; % degree
        limPitchangleExB = [0 70]; % degree
        limPitchangleVeperp = [0 60]; % degree
        limE = [80 40000];
      case 1213
        limPitchangleB = 90+[-25 20]; % degree
        limPitchangleExB = [0 90]; % degree
        limPitchangleVeperp = [0 90]; % degree
        limE = [100 40000];
        %doVePitch = 0;
      case 1214
        limPitchangleB = 90+20*[-1 1]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [110 40000];
      case 1215
        limPitchangleB = 90+25*[-1 1]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [110 40000];
      case 1216
        limPitchangleB = 90+[-10 30]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [80 40000];
      case 1217
        limPitchangleB = 90+[-00 30]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [80 40000];
      case 1218
        limPitchangleB = 90+[-00 30]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [80 40000];
      case 1219
        limPitchangleB = 90+[-00 30]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [100 40000];
      case 1220
        limPitchangleB = 90+[-00 30]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [100 40000];
      case 1225
        limPitchangleB = 90+[-00 30]; % degree
        limPitchangleExB = [0 60]; % degree
        limPitchangleVeperp = [0 180]; % degree
        limE = [100 40000];
      otherwise
        limPitchangleB = 90+40*[-1 1]; % degree
        limPitchangleExB = [0 110]; % degree
        limE = [60 40000];
      end
    case 2
    case 3 
    case 4
    
end
limPitchangleBGperp = [160 180];
limPitchangleVeperp = [0 180]; % degree


% Pick out local values 
c_eval('dist = ePDist?.tlim(time+[-0.01  0.01]);',ic)
c_eval('scpot = scPot?.resample(time);',ic)
c_eval('n = ne?.resample(time);',ic)
c_eval('Ve0 = gseVe?.resample(time).data;',ic); 
hatVe0 = double(irf_norm(Ve0));    
c_eval('Ve0perp = gseVe?perp.resample(time).data;',ic); 
hatVe0perp = irf_norm(Ve0perp);
c_eval('Ve0par = gseVe?par.resample(time).data;',ic); 
c_eval('dmpaB0 = dmpaB?.resample(time);',ic); 
hatB0 = double(irf_norm(dmpaB0.data));
c_eval('E0 = dslE?.resample(time).data;',ic); 
hatE0 = double(irf_norm(E0));
hatExB0 = cross(hatE0,hatB0);

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
tmpPitchangleVeperp = acosd([xX yY zZ]*hatVe0perp');
ePitchangleB = reshape(tmpPitchangleB,dataSize);
ePitchangleExB = reshape(tmpPitchangleExB,dataSize);
ePitchangleVeperp = reshape(tmpPitchangleVeperp,dataSize);
% surf(dX,dY,dZ,ePitchangleB)
% surf(dX,dY,dZ,ePitchangleExB)



%indB = [find(ePitchangleB>limPitchangleB(1)); find(ePitchangleB<limPitchangleB(2))];
if doVePitch
  matIncl = (ePitchangleB<limPitchangleB(2) & ...
             ePitchangleB>limPitchangleB(1) & ...
             ePitchangleVeperp>limPitchangleVeperp(1) & ...
             ePitchangleVeperp<limPitchangleVeperp(2) & ...
             EN>limE(1) & ...
             EN<limE(2));
else
  matIncl = (ePitchangleB<limPitchangleB(2) & ...
             ePitchangleB>limPitchangleB(1) & ...
             ePitchangleExB>limPitchangleExB(1) & ...
             ePitchangleExB<limPitchangleExB(2) & ...
             EN>limE(1) & ...
             EN<limE(2));
end

matBG = (ePitchangleB<limPitchangleB(2) & ...
         ePitchangleB>limPitchangleB(1) & ...
         ePitchangleVeperp>limPitchangleBGperp(1) & ...
         ePitchangleVeperp<limPitchangleBGperp(2) & ...
         EN>limE(1) & ...
         EN<limE(2));
           
       
tsMatIncl = irf.ts_scalar(time,single(reshape(matIncl,[1 size(matIncl)])));
tsMatExcl = irf.ts_scalar(time,single(not(reshape(matIncl,[1 size(matIncl)]))));
tsMatBG = irf.ts_scalar(time,single((reshape(matBG,[1 size(matBG)]))));

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
momsBG = mms.psd_moments(dist,TSeries(time,azimuthal),polar,stepTable,dist.ancillary.energy0,dist.ancillary.energy1,scpot,'electron','partialmoms',tsMatBG);

%nBG = momsBG.n_psd*180/abs(diff(limPitchangleBGperp));

% Decompose velocities into parallel and perpendicular components
[dmpaVparCore,dmpaVperpCore] = irf_dec_parperp(dmpaB0,momsCore.V_psd);
[dmpaVparCrescent,dmpaVperpCrescent] = irf_dec_parperp(dmpaB0,momsCrescent.V_psd);

% DMPA LMN
%[out,l,v] = irf_minvar(dmpaB4.tlim(irf.tint('2015-10-16T10:33:20.000Z/2015-10-16T10:33:38.000Z')));
%dmpaL = v(1,:); dmpaM = v(2,:); dmpaN = v(3,:);


% Plot partial distributions
distCrescent = dist; distCrescent.data(find(not(matIncl))) = NaN; 
distCore = dist; distCore.data(find((matIncl))) = NaN;
  
minE = energy(find(energy<limE(1),1,'last'));
maxE = energy(find(energy<limE(2),1,'last'));
  
  
if doPlot
  % Define coordinate system to plot in
  x = hatB0;
  y = hatExB0;
  z = cross(x,y);
  vlabels = {'B','ExB','Bx(ExB)'};

  % Set up some projection parameters
  vlim = 15*1e3;
  elevlim = 20;
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
  vectors = {hatVe0,'v_e';irf_norm(Ve0perp),'v_{e,\perp }'};
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

  vectors = {irf_norm(momsCrescent.V_psd.data),'v_e'};  
  % Crescent
  strTitle = ['n_e = ' num2str(momsCrescent.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCrescent.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCrescent.abs.data,'%.0f') ' km/s'];
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);
  hca.Title.String = strTitle;
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
  hca.Title.String = strTitle;
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,distCrescent,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
  hca.Title.String = strTitle;

  % Core
  vectors = {irf_norm(momsCore.V_psd.data),'v_e'};
  hca = h(isub); isub = isub + 1; 
  strTitle = ['n_e = ' num2str(momsCore.n_psd.data,'%.2f') ' cc, v_{e,\perp} = ' num2str(dmpaVperpCore.abs.data,'%.0f') ' km/s, v_{e,||} = ' num2str(dmpaVparCore.abs.data,'%.0f') ' km/s'];;
  mms.plot_projection(hca,distCore,'tint',time,'xyz',[x;y;z],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);       
  hca.Title.String = strTitle;
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,distCore,'tint',time,'xyz',[y;z;x],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([2,3,1]),'vectors',vectors);       
  hca.Title.String = strTitle;
  hca = h(isub); isub = isub + 1; 
  mms.plot_projection(hca,distCore,'tint',time,'xyz',[z;x;y],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels([3,1,2]),'vectors',vectors);       
  hca.Title.String = strTitle;


  % Fix figure
  for ii = 1:nRows*nCols
    colormap(h(ii),strCMap);
    %h(ii).Title.String = ''; 
  end
  strE = ['E = [' num2str(minE,'%g') ' ' num2str(maxE,'%g') ']'];
  strB = ['\theta_{B} = [' num2str(limPitchangleB(1),'%g') ' ' num2str(limPitchangleB(2),'%g') ']'];
  if doVePitch
    strV = ['\theta_{v_{e,\perp}} = [' num2str(limPitchangleVeperp(1),'%g') ' ' num2str(limPitchangleVeperp(2),'%g') ']'];
  else
    strV = ['\theta_{ExB} = [' num2str(limPitchangleExB(1),'%g') ' ' num2str(limPitchangleExB(2),'%g') ']'];
  end

  extraStrTitle = [timeUTC(1:23) ': ' strE ', ' strB ', ' strV];
  h(2).Title.String = {extraStrTitle,' ',h(2).Title.String};
  if doPrint 
    cn.print(['crescent_partial_moments_' irf_ssub('mms?',ic) '_' timeUTC([1:4 6:7 9:10 11 12:13 15:16 18:19]) 'ms' timeUTC(21:end-1)])
  end
end
%disp(['n_crescent = ' num2str(momsCrescent.n_psd.data)])
nBG(iInd) = momsBG.n_psd.data/diff(limPitchangleBGperp)*360;
nCore(iInd) = momsCore.n_psd.data;
nCrescent(iInd) = momsCrescent.n_psd.data; %nCrescent
vCore(iInd,:) = momsCore.V_psd.data;
vCrescent(iInd,:) = momsCrescent.V_psd.data;

%PCore(iInd,:) = momsCore.P_psd.data;
%PCrescent(iInd,:,:) = momsCrescent.P_psd.data;
PCore(iInd,:,:) = [momsCore.P_psd.data(1) momsCore.P_psd.data(2) momsCore.P_psd.data(3); momsCore.P_psd.data(2) momsCore.P_psd.data(4) momsCore.P_psd.data(5); momsCore.P_psd.data(3) momsCore.P_psd.data(5) momsCore.P_psd.data(6)];
PCrescent(iInd,:,:) = [momsCrescent.P_psd.data(1) momsCrescent.P_psd.data(2) momsCrescent.P_psd.data(3); momsCrescent.P_psd.data(2) momsCrescent.P_psd.data(4) momsCrescent.P_psd.data(5); momsCrescent.P_psd.data(3) momsCrescent.P_psd.data(5) momsCrescent.P_psd.data(6)];
%TCore(iInd,:) = momsCore.T_psd.data;
%TCrescent(iInd,:) = momsCrescent.T_psd.data;
TCore(iInd,:,:) = [momsCore.T_psd.data(1) momsCore.T_psd.data(2) momsCore.T_psd.data(3); momsCore.T_psd.data(2) momsCore.T_psd.data(4) momsCore.T_psd.data(5); momsCore.T_psd.data(3) momsCore.T_psd.data(5) momsCore.T_psd.data(6)];
TCrescent(iInd,:,:) = [momsCrescent.T_psd.data(1) momsCrescent.T_psd.data(2) momsCrescent.T_psd.data(3); momsCrescent.T_psd.data(2) momsCrescent.T_psd.data(4) momsCrescent.T_psd.data(5); momsCrescent.T_psd.data(3) momsCrescent.T_psd.data(5) momsCrescent.T_psd.data(6)];

momsCrescent(iInd) = momsCrescent;
momsCore(iInd) = momsCore;
momsBG(iInd) = momsBG;

momsMinE(iInd,1) = minE(1);
momsLimE(iInd,1) = limE(1);
end
%
c_eval('timeEpochTT = ePDist?(timeInds).time;',ic)
c_eval('limE? = irf.ts_scalar(timeEpochTT,tocolumn(momsLimE)); limE?.units =''eV'';',ic)
c_eval('minE? = irf.ts_scalar(timeEpochTT,tocolumn(momsMinE)); mimnE?.units =''eV'';',ic)
c_eval('limV? = sqrt(limE?*units.eV*2/units.me)/1000; limV?.units =''km/s'';',ic)
c_eval('minV? = sqrt(minE?*units.eV*2/units.me)/1000; minV?.units =''km/s'';',ic)

c_eval('ne?BG = irf.ts_scalar(timeEpochTT,tocolumn(nBG));',ic)
c_eval('ne?Core = irf.ts_scalar(timeEpochTT,tocolumn(nCore));',ic)
c_eval('ne?Crescent = irf.ts_scalar(timeEpochTT,tocolumn(nCrescent));',ic)
c_eval('ve?Core = irf.ts_vec_xyz(timeEpochTT,vCore);',ic)
c_eval('ve?Crescent = irf.ts_vec_xyz(timeEpochTT,vCrescent);',ic)

c_eval('Pe?Core = irf.ts_tensor_xyz(timeEpochTT,PCore);',ic)
c_eval('Pe?Crescent = irf.ts_tensor_xyz(timeEpochTT,PCrescent);',ic)

c_eval('Te?Core = irf.ts_tensor_xyz(timeEpochTT,TCore);',ic)
c_eval('Te?Crescent = irf.ts_tensor_xyz(timeEpochTT,TCrescent);',ic)

%% Transform to gse coordinates
c_eval('gseVe?Core = mms_dsl2gse(ve?Core,defatt?);',ic)
c_eval('gseVe?Crescent = mms_dsl2gse(ve?Crescent,defatt?);',ic)
c_eval('gsePe?Core = mms_dsl2gse(Pe?Core,defatt?);',ic)
c_eval('gsePe?Crescent = mms_dsl2gse(Pe?Crescent,defatt?);',ic)
c_eval('gseTe?Core = mms_dsl2gse(Te?Core,defatt?);',ic)
c_eval('gseTe?Crescent = mms_dsl2gse(Te?Crescent,defatt?);',ic)
%%
c_eval('mvaVe?Core = irf.ts_vec_xyz(gseVe?Core.time,[gseVe?Core.dot(L).data gseVe?Core.dot(M).data gseVe?Core.dot(N).data]); mvaVe?Core.name = ''Ve core LMN'';',ic)
c_eval('mvaVe?Crescent = irf.ts_vec_xyz(gseVe?Crescent.time,[gseVe?Crescent.dot(L).data gseVe?Crescent.dot(M).data gseVe?Crescent.dot(N).data]); mvaVe?Crescent.name = ''Ve cres LMN'';',ic)

c_eval('facPe?Core = mms.rotate_tensor(Pe?Core,''fac'',dmpaB?); facPe?Core = irf.ts_tensor_xyz(facPe?Core.time,facPe?Core.data);',ic)
c_eval('facPe?ppCore = mms.rotate_tensor(gsePe?Core,''fac'',gseB?,''pp''); facPe?pp = irf.ts_tensor_xyz(facPe?ppCore.time,facPe?Core.data);',ic)
c_eval('facPe?qqCore = mms.rotate_tensor(gsePe?Core,''fac'',gseB?,''qq''); facPe?pp = irf.ts_tensor_xyz(facPe?qqCore.time,facPe?Core.data);',ic)
c_eval('facTe?Core = mms.rotate_tensor(gseTe?Core,''fac'',gseB?); facTe?Core = irf.ts_tensor_xyz(facTe?Core.time,facTe?Core.data);',ic)

c_eval('facPe?Crescent = mms.rotate_tensor(Pe?Crescent,''fac'',dmpaB?); facPe?Crescent = irf.ts_tensor_xyz(facPe?Crescent.time,facPe?Crescent.data);',ic)
c_eval('facPe?ppCrescent = mms.rotate_tensor(gsePe?Crescent,''fac'',gseB?,''pp''); facPe?pp = irf.ts_tensor_xyz(facPe?ppCrescent.time,facPe?Crescent.data);',ic)
c_eval('facPe?qqCrescent = mms.rotate_tensor(gsePe?Crescent,''fac'',gseB?,''qq''); facPe?pp = irf.ts_tensor_xyz(facPe?qqCrescent.time,facPe?Crescent.data);',ic)
c_eval('facTe?Crescent = mms.rotate_tensor(gseTe?Crescent,''fac'',gseB?); facTe?Crescent = irf.ts_tensor_xyz(facTe?Crescent.time,facTe?Crescent.data);',ic)


% Agyrotropy and anisotropy
c_eval('Q?Core = irf.ts_scalar(facPe?ppCore.time,(facPe?ppCore.data(:,1,2).^2+facPe?ppCore.data(:,1,3).^2+facPe?ppCore.data(:,2,3).^2)./(facPe?ppCore.data(:,2,2).^2+2*facPe?ppCore.data(:,2,2).*facPe?ppCore.data(:,1,1))); Q?Core.name =''Q core'';',ic);
c_eval('Dng?Core = irf.ts_scalar(facPe?ppCore.time,sqrt(8*(facPe?ppCore.data(:,1,2).^2+facPe?ppCore.data(:,1,3).^2+facPe?ppCore.data(:,2,3).^2))./(facPe?ppCore.data(:,1,1)+2*facPe?ppCore.data(:,2,2))); Dng?Core.name =''Dng core'';',ic);
c_eval('Q?Crescent = irf.ts_scalar(facPe?ppCrescent.time,(facPe?ppCrescent.data(:,1,2).^2+facPe?ppCrescent.data(:,1,3).^2+facPe?ppCrescent.data(:,2,3).^2)./(facPe?ppCrescent.data(:,2,2).^2+2*facPe?ppCrescent.data(:,2,2).*facPe?ppCrescent.data(:,1,1))); Q?Crescent.name =''Q crescent'';',ic);
c_eval('Dng?Crescent = irf.ts_scalar(facPe?ppCrescent.time,sqrt(8*(facPe?ppCrescent.data(:,1,2).^2+facPe?ppCrescent.data(:,1,3).^2+facPe?ppCrescent.data(:,2,3).^2))./(facPe?ppCrescent.data(:,1,1)+2*facPe?ppCrescent.data(:,2,2))); Dng?Crescent.name =''Dng crescent'';',ic);

c_eval('T?ratioCore = irf.ts_scalar(facTe?Core.time,facPe?ppCore.data(:,1,1)./(facPe?ppCore.data(:,2,2))); T?ratioCore.name = ''Tpar/Tperp Core'';',ic);
c_eval('T?ratioCrescent = irf.ts_scalar(facTe?Crescent.time,facPe?ppCrescent.data(:,1,1)./(facPe?ppCrescent.data(:,2,2))); T?ratioCrescent.name = ''Tpar/Tperp Crescent'';',ic);

c_eval('Ao?Core = irf.ts_scalar(facPe?Core.time,2*abs(facPe?qqCore.data(:,2,2)-facPe?qqCore.data(:,3,3))./(facPe?qqCore.data(:,2,2)+facPe?qqCore.data(:,3,3))); Ao?Core.name = ''Ao Core'';',ic);
c_eval('Ao?Crescent = irf.ts_scalar(facPe?Crescent.time,2*abs(facPe?qqCrescent.data(:,2,2)-facPe?qqCrescent.data(:,3,3))./(facPe?qqCrescent.data(:,2,2)+facPe?qqCrescent.data(:,3,3))); Ao?Crescent.name = ''Ao Crescent'';',ic);

%irf_plot({mvaB1,mvaVe1perp,mvaVe1par,ne1,Dng1,Dng1Core,Dng1Crescent,Dng1Crescent/Dng1.resample(Dng1Crescent)})