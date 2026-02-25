% Make figure for article
% It will have four subplots, depicting the theoretical velocity field (or
% current?) with magnetic field lines, and below it the simulated magnetic
% field,a lso with magnetic field lines. This should be done for two
% different models, one that give magnetic field and one that doesnt. Maybe
% a third line should also be seen with the initial particle distributions
% that have a theoretical cure fitted to it?

fig = figure(13);
for k = 1:4; h(k) = subplot(2,2,k); end
set(fig,'position',[839   200   700   700]); % not full screen
set(h(1),'position',[0.15 0.50 0.35 0.35])
set(h(2),'position',[0.50 0.50 0.35 0.35])
set(h(3),'position',[0.15 0.15 0.35 0.35])
set(h(4),'position',[0.50 0.15 0.35 0.35])

isub = 1;

if 1 % 2007-08-31 model
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(vecModBz{1},2))*1e9; % nT
    rg = vecRG{1};
    zg = vecZG{1};
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',modBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rg([1 end]),'ylim',zg([1 end])) 
    title(hca,'modBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modelnr = 1; ExB.model;
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contour(hca,RS,ZS,PHI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(vecModBz{1},2))*1e9; % nT    
    modBr_rz = squeeze(mean(vecModBr{1},2))*1e9; % nT
    RR=squeeze(vecModfR{1}(:,1,:));
    ZZ=squeeze(vecModfZ{1}(:,1,:));
    quiver(hca,RR,ZZ,modBr_rz,modBz_rz,'k')  
    hold(hca,'off')
end
if 1 % 2007-08-31 simulation
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simBz_rz = squeeze(mean(vecSimBz{1},2))*1e9; % nT
    rg = vecRG{1};
    zg = vecZG{1};
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',simBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rg([1 end]),'ylim',zg([1 end])) 
    title(hca,'simBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modelnr = 1; ExB.model;
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contour(hca,RS,ZS,PHI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simBz_rz = squeeze(mean(vecSimBz{1},2))*1e9; % nT    
    simBr_rz = squeeze(mean(vecSimBr{1},2))*1e9; % nT
    RR=squeeze(vecModfR{1}(:,1,:));
    ZZ=squeeze(vecModfZ{1}(:,1,:));
    quiver(hca,RR,ZZ,simBr_rz,simBz_rz,'k')  
    hold(hca,'off')
end
if 1 % Tao2011 model
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(vecModBz{2},2))*1e9; % nT
    rg = vecRG{2};
    zg = vecZG{2};
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',modBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rg([1 end]),'ylim',zg([1 end])) 
    title(hca,'modBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modelnr = 2; ExB.model;
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contour(hca,RS,ZS,PHI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(vecModBz{2},2))*1e9; % nT    
    modBr_rz = squeeze(mean(vecModBr{2},2))*1e9; % nT
    RR=squeeze(vecModfR{2}(:,1,:));
    ZZ=squeeze(vecModfZ{2}(:,1,:));
    quiver(hca,RR,ZZ,modBr_rz,modBz_rz,'k')  
    hold(hca,'off')
end
if 1 % Tao2011 simulation
    hca = h(isub); isub = isub + 1;
    hold(hca,'on')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simBz_rz = squeeze(mean(vecSimBz{2},2))*1e9; % nT
    rg = vecRG{2};
    zg = vecZG{2};
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',simBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rg([1 end]),'ylim',zg([1 end])) 
    title(hca,'simBz')
    xlabel(hca,'r')
    zlabel(hca,'z') 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modelnr = 2; ExB.model;
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contour(hca,RS,ZS,PHI)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    simBz_rz = squeeze(mean(vecSimBz{2},2))*1e9; % nT    
    simBr_rz = squeeze(mean(vecSimBr{2},2))*1e9; % nT
    RR=squeeze(vecSimfR{2}(:,1,:));
    ZZ=squeeze(vecSimfZ{2}(:,1,:));
    quiver(hca,RR,ZZ,simBr_rz,simBz_rz,'k')  
    hold(hca,'off')
    caxis(hca,ceil(tdB/10^round(log10(tdB)))*10^round(log10(tdB))*[-1 1])%1*[-1 1]*tdB)
end

irf_colormap('poynting')
for ii=1:4
      pos=get(h(ii), 'Position');
      %set(h(ii), 'Position', [0.9*pos(1) pos(2) 0.82*pos(3) pos(4)]);
      caxis(h(ii),ceil(tdB/10^round(log10(tdB)))*10^round(log10(tdB))*[-1 1])%1*[-1 1]*tdB)
end
ch=colorbar('peer',h(1));
ylabel(ch,'B [nT]')
set(ch, 'Position', [.8714 .15 .0381 .70])