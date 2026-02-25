% Make figure for article
% It will have four subplots, depicting the theoretical velocity field (or
% current?) with magnetic field lines, and below it the simulated magnetic
% field,a lso with magnetic field lines. This should be done for two
% different models, one that give magnetic field and one that doesnt. Maybe
% a third line should also be seen with the initial particle distributions
% that have a theoretical cure fitted to it?

fig = figure(13);
for k = 1:6; h(k) = subplot(3,2,k); end
set(fig,'position',[839   515   838   441]); % not full screen
isub = 1;

if 1 % Plot of particle distributions
    hca = h(isub); isub = isub + 1;
    
end
if 1
    hca = h(isub); isub = isub + 1;
    
end
if 1 % Plot of model current
    hca = h(isub); isub = isub + 1;
    [RS,ZRS] = meshgrid(rSurf,zSurf);
    Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZRS/lz).^2)*1e-3; % V    
    ExB_az = - Er / (B0*1e-9); % km/s 
    e = 1.6022e-19;
    J_az = -n*e*ExB_az*1e9; % nA
    surfa = zeros(nz+1,nr+1);
    surf(hca,rGrid,zGrid,surfa,J_az)
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'nA')
    vlim = max(max(abs(J_az)));
    caxis(hca,vlim*[-1 1]*2)
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
    %title(hca,'ExB-drift')
end

if 1 % Plot of particle distributions
    hca = h(isub); isub = isub + 1;
    
end
if 1 % Plot of simulation current   
    hca = h(isub); isub = isub + 1;
    meanVrz = sumMVrz./sumMrz;    
    V_az = - meanVrz / (B0*1e-9); % km/s 
    e = 1.6022e-19;
    J_az = -n*e*V_az; % nA
    surfa = zeros(nz+1,nr+1);
    surf(hca,rGrid,zGrid,surfa,J_az')
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'nA')
    %vlim = max(max(abs(J_az)));
    caxis(hca,vlim*[-1 1]*2)
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
    %title(hca,'ExB-drift')
end
if 1
    hca = h(isub); isub = isub + 1;
    
end