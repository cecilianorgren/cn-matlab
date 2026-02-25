% check_B_contribution

% Load data
doLoad = 1;
if doLoad
loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/';
matToLoad = '20140508T085034-2007-08-31-200V-3509485';
varToLoad = {'rSurf','zSurf','sumMVXxyz','sumMVYxyz','sumMVZxyz',...
             'xGrid','yGrid','zGrid','sumMxyz',...
             'xSurf','ySurf','zSurf',...
             'lr','lz','B0','nr','nz','rlim','phi0','n',...
             'zlim','meanVrz','zGrid','rGrid','Tper'};
for kk = 1:numel(varToLoad)
    load([loadPath,matToLoad],varToLoad{kk});    
end
end
%%
% Subtract noise level from simulated ExB-drift in xyz-coordinates
% Take noise to be the mean between 10 and 20 and 80 and 90%.
limits = [0.02 0.25 0.75 0.98];
zlimits = -zlim+zlim*limits*2;
int1 = find(zSurf>(zlimits(1)) & zSurf<(zlimits(2)));
int2 = find(zSurf>(zlimits(3)) & zSurf<(zlimits(4)));
VXxyz = sumMVXxyz./sumMxyz;
VYxyz = sumMVYxyz./sumMxyz;
VZxyz = sumMVZxyz./sumMxyz;
VXnoise = repmat(mean(VXxyz(:,:,[int1, int2]),3),1,1,nz);
VYnoise = repmat(mean(VYxyz(:,:,[int1, int2]),3),1,1,nz);
VZnoise = repmat(mean(VZxyz(:,:,[int1, int2]),3),1,1,nz);
VXfix = VXxyz-VXnoise;
VYfix = VYxyz-VYnoise;
VZfix = VZxyz-VZnoise;


% Calculate magnetic field with noise subtracted
rbin = 34;
thbin= 23;
zbin = 34;
for sim = [ 1 ]; ExB.calculateB_rz; end
% Magnetic field
simB = sqrt(simBx.^2+simBy.^2+simBz.^2);
simBp = sqrt(simBx.^2+simBy.^2);
modB = sqrt(modBx.^2+modBy.^2+modBz.^2);
modBp = sqrt(modBx.^2+modBy.^2);
fixB = sqrt(fixBx.^2+fixBy.^2+fixBz.^2);
fixBp = sqrt(fixBx.^2+fixBy.^2);

% Theoretical dB at the origin
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

% Plot the three different fields
fig = figure(29); 
set(fig,'position',[1 1 400 800]);
nPanels = 4;
for k = 1:nPanels; h(k) = subplot(nPanels,1,k); end
isub = 1;
irf_colormap('poynting');

if 1 % Plot model Bz in rz-plane         
    hca = h(isub); isub = isub + 1;         
    modBz_rz = squeeze(mean(modBz,2))*1e9; % nT
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',modBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[0 1],'ylim',zlim*[-1 1]) 
    title(hca,'modBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    colorbar('peer',hca);
    caxis(hca,tdB*1*[-1 1])
    grid(hca,'off')
    box(hca,'on')
end
if 1 % Plot simulation Bz in rz-plane         
    hca = h(isub); isub = isub + 1;         
    simBz_rz = squeeze(mean(simBz,2))*1e9; % nT
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',simBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[0 1],'ylim',zlim*[-1 1]) 
    title(hca,'simBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    colorbar('peer',hca);
    caxis(hca,tdB*1*[-1 1])
    grid(hca,'off')
    box(hca,'on')
end
if 1 % Plot noise-fixed simulation Bz in rz-plane         
    hca = h(isub); isub = isub + 1;       
    fixBz_rz = squeeze(mean(fixBz,2))*1e9; % nT
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',fixBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[0 1],'ylim',zlim*[-1 1]) 
    title(hca,'fixBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    colorbar('peer',hca);
    caxis(hca,tdB*1*[-1 1])
    grid(hca,'off')
    box(hca,'on')
end
if 1 % Plot difference between model and noise-fixed field
    hca = h(isub); isub = isub + 1;         
    diffBz_rz = squeeze(mean(fixBz,2))*1e9-squeeze(mean(modBz,2))*1e9; % nT
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',diffBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[0 1],'ylim',zlim*[-1 1]) 
    title(hca,'fixBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    colorbar('peer',hca);
    caxis(hca,tdB*1*[-1 1])
    grid(hca,'off')
    box(hca,'on')
end