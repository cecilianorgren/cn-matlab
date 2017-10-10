% Make figure for article
% It will have two subplots, depicting the theoretical velocity field (or
% current?) with magnetic field lines, and the simulated magnetic
% field,also with magnetic field lines.
plotAll = 0;
plotOne = 1;
doLoad = 0;

loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/';
allMats = ls([loadPath 'Spis/*.mat']);
indexDot = strfind(allMats,'.');
indexStart = [1 indexDot(1:(end-1))+5];
indexEnd = indexDot+3;
indexFileSep = strfind(allMats,'/');
nFiles = numel(indexDot);
vecToLoad = cell(1,nFiles);
vecFileName = cell(1,nFiles);
for kk = 1:nFiles, vecToLoad{kk} = allMats(indexStart(kk):indexEnd(kk)); end
for pp = 1:nFiles
    son = find(indexFileSep<indexEnd(pp),1,'last'); 
    vecFileName{pp} = allMats(indexFileSep(son)+1:indexEnd(pp)); 
end
nLast = 2;
if plotAll; plotThis = 1:nFiles; 
elseif plotOne; plotThis = nFiles-nLast+1;
else plotThis = (nFiles-nLast+1):nFiles; 
end 
    

for ii =plotThis, 
    try
toLoad = vecToLoad{ii};
fig = figure(13);
for k = 1:2; h(k) = subplot(1,2,k); end
set(fig,'position',[839   200   700   350]); % not full screen
irf_colormap('poynting')

%
%toLoad = 'Spis/20140424T075537-2007-08-31-wider-less-dense-2193397.mat'; 
varToLoad = {'modBz','rg','zg','rbin','zbin','rSurf','zSurf','phi0',...
             'lr','lz','fR','fZ','n','B0','modBr','simBr','simBz',...
             'nParticles','name'};
for kk = 1:numel(varToLoad)
    %load([loadPath toLoad],varToLoad{kk});
    load(vecToLoad{ii},varToLoad{kk});
end

mu0 = 1.2566e-6; e = 1.6022e-19;
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

isub = 1;

if 1 % 2007-08-31 model
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(modBz,2))*1e9; % nT
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',modBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])    
    title(hca,'Model magnetic field')
    xlabel(hca,'r')
    ylabel(hca,'z') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contour(hca,RS,ZS,PHI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(modBz,2))*1e9; % nT    
    modBr_rz = squeeze(mean(modBr,2))*1e9; % nT
    RR=squeeze(fR(:,1,:));
    ZZ=squeeze(fZ(:,1,:));
    st_r = 2;
    st_z = 2;
    quiver(hca,RR(1:st_r:end,1:st_z:end),ZZ(1:st_r:end,1:st_z:end),modBr_rz(1:st_r:end,1:st_z:end),modBz_rz(1:st_r:end,1:st_z:end),'k')  
    hold(hca,'off')
    box(hca,'on')
    axis(hca,'equal')
end
if 1 % 2007-08-31 simulation
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simBz_rz = squeeze(mean(simBz,2))*1e9; % nT
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',simBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])    
    title(hca,'Simulation magnetic field')
    xlabel(hca,'r')
    set(hca,'ytick',[])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contour(hca,RS,ZS,PHI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simBz_rz = squeeze(mean(simBz,2))*1e9; % nT    
    simBr_rz = squeeze(mean(simBr,2))*1e9; % nT
    RR=squeeze(fR(:,1,:));
    ZZ=squeeze(fZ(:,1,:));
    st_r = 2;
    st_z = 2;
    quiver(hca,RR(1:st_r:end,1:st_z:end),ZZ(1:st_r:end,1:st_z:end),simBr_rz(1:st_r:end,1:st_z:end),simBz_rz(1:st_r:end,1:st_z:end),'k')      
    %quiver(hca,RR,ZZ,simBr_rz,simBz_rz,'k')  
    hold(hca,'off')
    box(hca,'on')
    axis(hca,'equal')
end

ch=colorbar('peer',h(2));
ylabel(ch,'B [nT]');
set(h(1),'position',[0.20 0.20 0.30 0.65]);
set(h(2),'position',[0.50 0.20 0.30 0.65]);

for kk=1:2;
    set(h(kk),'clim',tdB*[-1 1],'xlim',[0 4*lr],'ylim',lz*6*[-1 1],'xtick',[0 20 40 60 80 100]) 
end
%cn.print(['B_' vecFileName{ii}])
    catch err
    end
end