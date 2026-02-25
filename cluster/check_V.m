

% Plot the three different fields
fig = figure(30); 
set(fig,'position',[1 1 800 800]);
nPanels = 6;
for k = 1:nPanels; h(k) = subplot(3,2,k); end
isub = 1;
irf_colormap('poynting');

if 1 % 
    [XS,YS,ZS] = cn.meshgrid(xSurf,ySurf,zSurf);
    Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
    Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    sVX = Ey/(B0*1e-9); % m/s 
    sVY = -Ex/(B0*1e-9); % m/s 
    pV2x = squeeze(sVX(:,:,ceil(nz/2)))*1e-3;
    pV2y = squeeze(sVY(:,:,ceil(nz/2)))*1e-3;
    
    hca = h(isub); isub = isub + 1;      
    surf(hca,xGrid,yGrid,zeros(nx+1,ny+1)',pV2x');
    title(hca,'modVX')
    hca = h(isub); isub = isub + 1;  
    surf(hca,xGrid,yGrid,zeros(nx+1,ny+1)',pV2y');
    title(hca,'modVY')
end
if 1 % 
    sVX = sumMVXxyz./sumMxyz; % m/s
    sVY = sumMVYxyz./sumMxyz; % m/s
    sVZ = (sumMVZxyz./sumMxyz)*0; % m/s
    pV1x = squeeze(sVX(:,:,ceil(nz/2)))*1e-3;
    pV1y = squeeze(sVY(:,:,ceil(nz/2)))*1e-3;    
    
    hca = h(isub); isub = isub + 1;      
    surf(hca,xGrid,yGrid,zeros(nx+1,ny+1)',pV1x');
    title(hca,'simVX')
    hca = h(isub); isub = isub + 1;  
    surf(hca,xGrid,yGrid,zeros(nx+1,ny+1)',pV1y');
    title(hca,'simVY')
end
if 1 % 
    limits = [0.05 0.25 0.75 0.95];
    zlimits = -zlim+zlim*limits*2;
    int1 = find(zSurf>(zlimits(1)) & zSurf<(zlimits(2)));
    int2 = find(zSurf>(zlimits(3)) & zSurf<(zlimits(4)));
    VXnoise = repmat(mean(meanVXxyz(:,:,[int1, int2]),3),1,1,nz);
    VYnoise = repmat(mean(meanVYxyz(:,:,[int1, int2]),3),1,1,nz);
    VZnoise = repmat(mean(meanVZxyz(:,:,[int1, int2]),3),1,1,nz);
    VXfix = meanVXxyz-VXnoise;
    VYfix = meanVYxyz-VYnoise;
    VZfix = meanVZxyz-VZnoise;       
    sVX = VXfix./sumMxyz;
    sVY = VYfix./sumMxyz;
    sVZ = VZfix./sumMxyz;
    pV1x = squeeze(sVX(:,:,ceil(nz/2)))*1e-3;
    pV1y = squeeze(sVY(:,:,ceil(nz/2)))*1e-3;    
    
    hca = h(isub); isub = isub + 1;      
    surf(hca,xGrid,yGrid,zeros(nx+1,ny+1)',pV1x');
    title(hca,'fixVX')
    hca = h(isub); isub = isub + 1;  
    surf(hca,xGrid,yGrid,zeros(nx+1,ny+1)',pV1y');
    title(hca,'fixVY')
end

for ii = 1:nPanels
    hca = h(ii);
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])     
    xlabel(hca,'x')
    ylabel(hca,'y') 
    colorbar('peer',hca);
    caxis(hca,3000*[-1 1])
    grid(hca,'off')
    box(hca,'on')
end