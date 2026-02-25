%% Make position and electric field matrices
[XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
E = sqrt(Ex.^2+Ey.^2+Ez.^2);
ExB_y = - Ex/(B0*1e-9); % m/s 
ExB_x = + Ey/(B0*1e-9); % m/s 
ExB_y = ExB_y*1e-3; % km/s
ExB_x = ExB_x*1e-3; % km/s
%% Plot quivers
fig = figure(11);
for k = 1:6; h(k) = subplot(2,3,k); end
set(fig,'position',[1 1 1400 820]); % not full screen
isub = 1;
irf_colormap('poynting'); 
    tpx = 5:10:numel(xSurf);
    tpy = 5:10:numel(ySurf);
    tpz = 5:10:numel(zSurf);
if 0 % Plot total electric field in arrows
    [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
    hca = h(isub); isub = isub + 1;        
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz),Ey(tpx,tpy,tpz),Ez(tpx,tpy,tpz))
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E field')
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end
if 1 % Plot x-electric field in arrows
    hca = h(isub); isub = isub + 1;      
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz),Ey(tpx,tpy,tpz)*0,Ez(tpx,tpy,tpz)*0)
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E_x')
    view(hca,[0 1 0])
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end   
if 1 % Plot y-electric field in arrows
    hca = h(isub); isub = isub + 1;      
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz)*0,Ey(tpx,tpy,tpz),Ez(tpx,tpy,tpz)*0)
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E_y')
    view(hca,[0 1 0])
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end  


%% Plot surface
%fig = figure(11);
%for k = 1:3; h(k) = subplot(1,3,k); end
%set(fig,'position',[1 1 1400 820]); % not full screen
%isub = 1;
irf_colormap('poynting');

if 0 % Plot total electric field in arrows
    [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
    hca = h(isub); isub = isub + 1;        
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz),Ey(tpx,tpy,tpz),Ez(tpx,tpy,tpz))
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E field')
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end
if 1 % Plot the x electric field in xz-plane            
    meanind = fix(size(Ex,2)/2);
    m_Ex = squeeze(mean(Ex(:,meanind+(-1:1),:),2));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_Ex')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    %vlim = max(max(abs(ExBdrift)));
    %caxis(hca,vlim*[-1 1]*2)
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'E_x')
end         
if 1 % Plot the y electric field in xz-plane            
    meanind = fix(size(Ey,2)/2);
    m_Ey = squeeze(mean(Ey(:,meanind+(-1:1),:),2));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_Ey')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    %vlim = max(max(abs(ExBdrift)));
    %caxis(hca,vlim*[-1 1]*2)
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'E_y')
end
if 1 % Plot the y electric field in xy-plane            
    meanind = fix(size(Ey,2)/2);
    m_Ey = squeeze(mean(Ey(:,:,meanind+(-1:1)),3));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(ny+1,nx+1);
    surf(hca,xGrid,yGrid,surfa,m_Ey')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    %vlim = max(max(abs(ExBdrift)));
    %caxis(hca,vlim*[-1 1]*2)
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'E_y')
end
if 1 % Plot the x electric field in xy-plane            
    meanind = fix(size(Ex,2)/2);
    m_Ex = squeeze(mean(Ex(:,:,meanind+(-1:1)),3));
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(ny+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_Ey')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    %vlim = max(max(abs(ExBdrift)));
    %caxis(hca,vlim*[-1 1]*2)
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,'E_y')
end


%% Plot E's
fig = figure(11);
for k = 1:6; h(k) = subplot(2,3,k); end
set(fig,'position',[1 1 1400 820]); % not full screen
isub = 1;
irf_colormap('poynting'); 

if 1 % Plot x-electric field in arrows
    hca = h(isub); isub = isub + 1;      
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz),Ey(tpx,tpy,tpz)*0,Ez(tpx,tpy,tpz)*0)
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E_x')
    view(hca,[0 1 0])
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end   
if 1 % Plot y-electric field in arrows
    hca = h(isub); isub = isub + 1;      
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz)*0,Ey(tpx,tpy,tpz),Ez(tpx,tpy,tpz)*0)
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E_y')
    view(hca,[0 1 0])
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end
if 1 % Plot y-electric field in arrows
    hca = h(isub); isub = isub + 1;      
    quiver3(hca,XS(tpx,tpy,tpz),YS(tpx,tpy,tpz),ZS(tpx,tpy,tpz),...
                Ex(tpx,tpy,tpz)*0,Ey(tpx,tpy,tpz),Ez(tpx,tpy,tpz)*0)
    xlabel(hca,'x')
    ylabel(hca,'y')
    zlabel(hca,'z')
    title(hca,'E_y')
    view(hca,[0 0 1])
    axis(hca,'equal')
    grid(hca,'off')            
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1],'zlim',zlim*[-1 1])    
end  
if 1 % Plot the y electric field in xz-plane            
    hca = h(isub); isub = isub + 1;
    ind = fix(size(Ey,2)*1/2);
    m_Ey = squeeze(Ey(:,ind,:));    
    surfa = zeros(nz+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_Ey')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    elim = max(max(max(E)));
    caxis(hca,elim*[-1 1]);
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,['E_y at y = ' num2str(yGrid(ind))])
end
if 1 % Plot the y electric field in xy-plane            
    hca = h(isub); isub = isub + 1;
    ind = fix(size(Ey,2)*3/4);
    m_Ey = squeeze(Ey(:,ind,:));    
    surfa = zeros(nz+1,nx+1);
    surf(hca,xGrid,zGrid,surfa,m_Ey')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    elim = max(max(max(E)));
    caxis(hca,elim*[-1 1]);
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,['E_y at y = ' num2str(yGrid(ind))])
end
if 1 % Plot the y electric field in xy-plane            
    hca = h(isub); isub = isub + 1;
    ind = fix(size(Ey,3)*2/4); % halfway
    m_Ey = squeeze(Ey(:,:,ind));    
    surfa = zeros(ny+1,nx+1);
    surf(hca,xGrid,yGrid,surfa,m_Ey')
    shading(hca,'flat')
    xlabel(hca,'x'); ylabel(hca,'y')
    cha = colorbar('peer',hca); ylabel(cha,'V')
    elim = max(max(max(E)));
    caxis(hca,elim*[-1 1]);
    view(hca,[0 0 1])
    axis(hca,'equal')
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',[-zlim zlim]) 
    title(hca,['E_y at z = ' num2str(zGrid(ind))])
end

%% Look at XS YS ZS
[tocolumn(squeeze(XS(:,1,1))) tocolumn(squeeze(XS(1,:,1))) tocolumn(squeeze(XS(1,1,:))),...
 tocolumn(squeeze(YS(:,1,1))) tocolumn(squeeze(YS(1,:,1))) tocolumn(squeeze(YS(1,1,:))),...
 tocolumn(squeeze(ZS(:,1,1))) tocolumn(squeeze(ZS(1,:,1))) tocolumn(squeeze(ZS(1,1,:)))]