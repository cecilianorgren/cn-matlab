% check how model v checks out with simulation v

%% Make new grid with different cells in each direction to check plots

nx = 120;
ny = 120;
nz = 120; 
nr = 60; 
% centers of bins
xGrid = linspace(-rlim,rlim,nx+1); 
yGrid = linspace(-rlim,rlim,ny+1); 
zGrid = linspace(-zlim,zlim,nz+1); 
rGrid = linspace(0,rlim,nr+1); 
% edges between bins
xSurf = xGrid(1:end-1)+diff(xGrid(1:2))/2; 
ySurf = yGrid(1:end-1)+diff(yGrid(1:2))/2;  
zSurf = zGrid(1:end-1)+diff(zGrid(1:2))/2; 
rSurf = rGrid(1:end-1)+diff(rGrid(1:2))/2;

%% Make position and electric field matrices from model
[XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
E = sqrt(Ex.^2+Ey.^2+Ez.^2);
ExB_y = - Ex/(B0*1e-9); % m/s 
ExB_x = + Ey/(B0*1e-9); % m/s 
ExB_y = ExB_y*1e-3; % km/s
ExB_x = ExB_x*1e-3; % km/s
ExB = sqrt(ExB_x.^2+ExB_y.^2);

%% Make position and electric field matrices from simulation
VAZ = sumMVxyz./sumMxyz*1e-3; % km/s
VX = sumMVXxyz./sumMxyz*1e-3; % km/s
VY = sumMVYxyz./sumMxyz*1e-3; % km/s
%% Make plot
fig = figure(12);
for k = 1:6; h(k) = subplot(2,3,k); end
set(fig,'position',[839   515   838   441]); % not full screen
isub = 1;
irf_colormap('poynting'); 
vlim = max(max(max(abs(ExB))))*1.2; % common colorbar axis
irange = -8:8;
% From model
if 1 % Plot ExB_y electric field in xz-plane 
    meanind = fix(size(ExB_y,2)/2);
    m_ExB_y = squeeze(mean(ExB_y(meanind+irange,:,:),1));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_ExB_y')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')
    caxis(hca,vlim*[-1 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'ExB_y')
end 
if 1 % Plot ExB_y electric field in yz-plane 
    meanind = fix(size(ExB_y,2)/2);
    m_ExB_y = squeeze(mean(ExB_y(:,meanind+irange,:),2));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,ny+1);
    surf(hca,yGrid,zGrid,surfa,m_ExB_y')
    shading(hca,'flat')
    xlabel(hca,'y'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')
    caxis(hca,vlim*[-1 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'ExB_y')
end 
if 1 % Plot ExB_y in xy-plane 
    meanind = fix(size(ExB_y,3)/2);
    m_ExB_y = squeeze(mean(ExB_y(:,:,meanind+irange),3));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(ny+1,nx+1);
    surf(hca,xGrid,yGrid,surfa,m_ExB_y)
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'y')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')    
    caxis(hca,vlim*[-1 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'ExB_y')
end
% From simulation
if 1 % Plot ExB_y electric field in xz-plane
    meanind = fix(size(VY,2)/2);
    m_VY = squeeze(mean(VY(:,meanind+irange,:),2));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_VY')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')
    caxis(hca,vlim*[-1 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'ExB_y')
end 
if 1 % Plot ExB_y electric field in yz-plane
    meanind = fix(size(VY,2)/2);
    m_VY = squeeze(mean(VY(meanind+irange,:,:),1));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,ny+1);
    surf(hca,yGrid,zGrid,surfa,m_VY')
    shading(hca,'flat')
    xlabel(hca,'y'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')
    caxis(hca,vlim*[-1 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'ExB_y')
end 
if 1 % Plot ExB_y in xy-plane
    meanind = fix(size(VY,3)/2);
    m_VY = squeeze(mean(VY(:,:,meanind+irange),3));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(ny+1,nx+1);
    surf(hca,xGrid,yGrid,surfa,m_VY')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'y')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')    
    caxis(hca,vlim*[-1 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'ExB_y')
end
