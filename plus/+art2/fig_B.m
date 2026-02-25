% Assume all relevant data is loaded.
mu0 = 1.2566e-6; e = 1.6022e-19;
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

fig = figure(13);
for k = 1:2; h(k) = subplot(1,2,k); end
set(fig,'position',[839   200   700   350]); % not full screen
irf_colormap('poynting')
isub = 1;
doOutliers=1;
if 1 % 2007-08-31 model
    hca = h(isub); isub = isub + 1;
    set(gcf,'CurrentAxes',hca)
    hold(hca,'on')    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    modBz_rz = squeeze(mean(modBz,2))*1e9; % nT
    if doOutliers
    %[outMax outRow outCol] = cn.Nmax(fixBz_rz,0);
    outRow = [19 19 19];
    outCol = [13 18 23];
    for kk=outRow
        for pp=outCol
            meanof4 = nanmean([modBz_rz(kk-1,pp-1) modBz_rz(kk-1,pp+1) modBz_rz(kk+1,pp-1) modBz_rz(kk+1,pp+1)]);            
            modBz_rz(kk,pp) = meanof4;
        end
    end
    end
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',modBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])    
    title(hca,'Analytical magnetic field ')
    xlabel(hca,'r [km]')
    ylabel(hca,'z [km]') 
    text(2,25,'a)')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contours=phi0*0.01*[90 80 70 60 50 40 30 20 10];
    contour(hca,RS,ZS,PHI,contours(1:2:end))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
    modBz_rz = squeeze(mean(modBz,2))*1e9; % nT    
    modBr_rz = squeeze(mean(modBr,2))*1e9; % nT
    RR=squeeze(fR(:,1,:));
    ZZ=squeeze(fZ(:,1,:));
    st_r = 2;
    st_z = 2;
    quiver(hca,RR(1:st_r:end,1:st_z:end),ZZ(1:st_r:end,1:st_z:end),modBr_rz(1:st_r:end,1:st_z:end),modBz_rz(1:st_r:end,1:st_z:end),'k')  
    end
    hold(hca,'off')
    box(hca,'on')
    axis(hca,'equal')
end
if 1 % 2007-08-31 simulation
    hca = h(isub); isub = isub + 1;
    set(gcf,'CurrentAxes',hca)
    hold(hca,'on')    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot colorplot of Bz
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fixBz_rz = squeeze(mean(fixBz,2))*1e9; % nT
    % take away real outliers
    if doOutliers
    [outMax outRow outCol] = cn.Nmax(fixBz_rz,5);
    for kk=outRow
        for pp=outCol
            meanof4 = nanmean([fixBz_rz(kk-1,pp-1) fixBz_rz(kk-1,pp+1) fixBz_rz(kk+1,pp-1) fixBz_rz(kk+1,pp+1)]);
            %if abs(fixBz_rz(kk,pp)) > 2*meanof4
                fixBz_rz(kk,pp) = meanof4;
            %end
        end
    end
    end
            
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',fixBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])    
    title(hca,'Simulation magnetic field ')
    xlabel(hca,'r [km]')
    set(hca,'ytick',[])
    text(2,25,'b)')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot equipotential contours of phi
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    [RS,ZS] = meshgrid(rSurf,zSurf);
    PHI = phi0*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2); % V
    contours=phi0*0.01*[90 80 70 60 50 40 30 20 10];
    contour(hca,RS,ZS,PHI,contours(1:2:end))
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot magnetic field direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if 0
    simBz_rz = squeeze(mean(simBz,2))*1e9; % nT    
    simBr_rz = squeeze(mean(simBr,2))*1e9; % nT
    RR=squeeze(fR(:,1,:));
    ZZ=squeeze(fZ(:,1,:));
    st_r = 2;
    st_z = 2;
    quiver(hca,RR(1:st_r:end,1:st_z:end),ZZ(1:st_r:end,1:st_z:end),fixBr_rz(1:st_r:end,1:st_z:end),fixBz_rz(1:st_r:end,1:st_z:end),'k')      
    %quiver(hca,RR,ZZ,simBr_rz,simBz_rz,'k')  
    end
    hold(hca,'off')
    box(hca,'on')
    axis(hca,'equal')
end

ch=colorbar('peer',h(2));
ylabel(ch,'\delta B_{||} [nT]');
set(h(1),'position',[0.20 0.20 0.30 0.65]);
set(h(2),'position',[0.50 0.20 0.30 0.65]);

for kk=1:2;
    set(h(kk),'clim',tdB*[-1 1],'xlim',[0 4.5*lr],'ylim',lz*6*[-1 1],'xtick',[0 20 40 60 80 100]) 
end


%% Plot ratio
if 0
%%
h=subplot(1,1,1);
hca=h;
Bratio=fixBz_rz./modBz_rz;
surf(hca,rg,zg,zeros(rbin+1,zbin+1)',Bratio');
shading(hca,'flat')
view(hca,[0 0 1])    
title(hca,'fixB/modB ')
xlabel(hca,'r')
ylabel(hca,'z') 
colorbar('peer',hca)
%set(hca,'clim',[-1 1])

set(hca,'clim',[-1 1],'xlim',[0 4*lr],'ylim',lz*6*[-1 1],'xtick',[0 20 40 60 80 100]) 

    
end