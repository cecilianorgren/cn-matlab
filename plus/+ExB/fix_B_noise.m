% Calculates B from simulation, model simulation, and noised-fixed
% simulation. Must define toSim = [], doSave, doPlot, doLoad (and also matToLoad)
%

% Load data
%doLoad = 0;
[~, hostname] = system('hostname');
if findstr(hostname,'spis')
    loadPath = '/home/cecilia/Research/ElectronHoleSimulation/SavedData/';
    savePath = '/home/cecilia/Research/ElectronHoleSimulation/SavedData/B/';
elseif any([findstr(hostname,'mother') findstr(hostname,'phasespace')]);
    loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/';
    savePath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/B/';
end
if doLoad     
    %matToLoad = '20140508T085034-2007-08-31-200V-3509485';
    % kolla ExB.see_list f?r lista p? alla runs...
    
    matToLoad = '20140611T144721-2007-08-31-500V_lr12_lz5_2000eV-1585950.mat';
    matToLoad = '20140516T103255-2007-08-31-200V_lr0_lz5_16eV-2175774.mat';
    matToLoad = '20140517T080958-2007-08-31-200V_lr0_lz5_1300eV-3260494.mat';
    matToLoad = '20140517T192323-2007-08-31-200V_lr0_lz5_2500eV-3261101.mat';
    matToLoad = '20140626T010016-2007-08-31-500V_lr12_lz5_2000eV-3961255.mat';
    varToLoad = {'rSurf','zSurf','sumMVXxyz','sumMVYxyz','sumMVZxyz',...
                 'xGrid','yGrid','zGrid','sumMxyz',...
                 'xSurf','ySurf','zSurf',...
                 'lr','lz','B0','nr','nz','rlim','phi0','n',...
                 'zlim','meanVrz','zGrid','rGrid','Tper'};
    for kk = 1:numel(varToLoad)
        load([loadPath,matToLoad],varToLoad{kk});    
    end
end

%% Subtract noise level from simulated ExB-drift in xyz-coordinates
% Take noise to be the mean between 10 and 20 and 80 and 90%.
limits = [0.02 0.20 0.80 0.98];
zlimits = -zlim+zlim*limits*2;
zlimits = [-zlim -11 11 zlim];
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


%% Calculate magnetic field with noise subtracted
%rbin = 34;
%thbin= 23;
%zbin = 34;
%rbin = 11;
%thbin= 11;
%zbin = 11;
for sim = toSim; ExB.calculateB_rz; end
% Magnetic field
simB = sqrt(simBx.^2+simBy.^2+simBz.^2);
simBp = sqrt(simBx.^2+simBy.^2);
modB = sqrt(modBx.^2+modBy.^2+modBz.^2);
modBp = sqrt(modBx.^2+modBy.^2);
fixB = sqrt(fixBx.^2+fixBy.^2+fixBz.^2);
fixBp = sqrt(fixBx.^2+fixBy.^2);

% Theoretical dB at the origin
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

%% Save the magnetic fields
if doSave
    saveFilePath = [savePath datestr(now,'yyyymmddTHHMMSS') '-' name '-B-' num2str(nParticles)];
    save(saveFilePath)
    disp(['Saved ' saveFilepath])
end

%% Plot the three different fields
if doPlot
    fig = figure(29); 
    set(fig,'position',[1 1 1000 800]);
    nPanels = 9;
    nRows = 3;
    for k = 1:nPanels; h(k) = subplot(ceil(nPanels/nRows),nRows,k); end
    isub = 1;
    irf_colormap('poynting');
    Vind = [];
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
    if 1
        hca = h(isub); Vind = [Vind isub]; isub = isub + 1;
        pcolor(hca,squeeze(VXmod(:,:,fix(size(VXmod,3)/2))))
        colorbar('peer',hca)
        caxmax = [-1 1]*max(get(hca,'clim'));
        title(hca,'before trim')
    end
    if 1
        hca = h(isub); Vind = [Vind isub]; isub = isub + 1;
        pcolor(hca,squeeze(modVX(:,:,fix(size(modVX,3)/2))))
        colorbar('peer',hca)
        %caxmax = [-1 1]*max(get(hca,'clim'));        
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
    if 1
        hca = h(isub); Vind = [Vind isub]; isub = isub + 1;
        pcolor(hca,squeeze(VXsim(:,:,fix(size(VXsim,3)/2))))
        colorbar('peer',hca)
        title(hca,'before trim')
    end
    if 1
        hca = h(isub); Vind = [Vind isub]; isub = isub + 1;
        pcolor(hca,squeeze(simVX(:,:,fix(size(simVX,3)/2))))
        colorbar('peer',hca)
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
    if 1
        hca = h(isub); Vind = [Vind isub]; isub = isub + 1;
        pcolor(hca,squeeze(fixVX(:,:,fix(size(fixVX,3)/2))))
        colorbar('peer',hca)
        title(hca,'before trim')
    end
    if 1
        hca = h(isub); Vind = [Vind isub]; isub = isub + 1;
        pcolor(hca,squeeze(fixVX(:,:,fix(size(fixVX,3)/2))))
        colorbar('peer',hca)
    end
    if 0 % Plot difference between model and noise-fixed field
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
    
    
    
    for  k = Vind
        shading(h(k),'flat')
        caxis(h(k),caxmax)
    end
end