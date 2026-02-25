% Plotting options
overrideVlim = 0;
fig = figure(27); 
set(fig,'position',[1 1 1676 955]);
for k = 1:15; h(k) = subplot(4,4,k); end
isub = 1;
irf_colormap('poynting');

% See how many particles bounced
maxBounce = max(nBounce);
[nnn,~] = histc(nBounce,0:maxBounce);    
%bar(hca,0:maxBounce,nnn);
fractionBounce = sum(nnn(2:end))/sum(nnn);

% Info to be sued in plots
units = irf_units;
me=units.me;
oce = e*25e-9/me ; % rad/s
fce = oce/2/pi; % Hz 
re = vtper/fce/2/pi;
% Add info     
if exist('toLoad','var'); toLoadStr = toLoad; else toLoadStr = ' '; end
strTitle = {['Load data at: ' parseunderline(toLoadStr)],...
            ['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V,  '],...    
            ['l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km'],...
            ['N_{Box} = ' num2str(nParticlesBox) ',  N_{Flow} = ' num2str(nParticlesFlow) ',  n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')],...
            ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV'],...
            ['v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
            };
annotation(gcf,'textbox',[0.75 0.10 0.17 0.2],'string',strTitle)

% Electric field and ExB-drift
% rz
[RS,ZS] = cn.meshgrid(rSurf,zSurf);
Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
ExBdrift = - Er/(B0*1e-9)*1e-3; % km/s
% xyz
[XS,YS,ZS] = cn.meshgrid(zSurf,ySurf,zSurf);
Ex = sqrt(XS.^2+YS.^2)./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ey = sqrt(XS.^2+YS.^2)./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ez = ZS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
ExBdriftXYZ = - sqrt(Ex.^2+Ey.^2)/(B0*1e-9)*1e-3; % km/s

vlim = max(max(abs(ExBdrift)));
if vlim == 0; vlim = vtper/100; end
if overrideVlim; vlim = 484.3; end

% Magnetic field
simB = sqrt(simBx.^2+simBy.^2+simBz.^2);
simBp = sqrt(simBx.^2+simBy.^2);
modB = sqrt(modBx.^2+modBy.^2+modBz.^2);
modBp = sqrt(modBx.^2+modBy.^2);
% Theoretical dB at the origin
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

% Starting position
r0 = sqrt(x0.^2+y0.^2);
%rc = sqrt(xc.^2+yc.^2);

% Azimuthal velocity in xy plane
iz50 = ceil(nz*2/4)+[-2 3]; meanVxyz_050 = mean(meanVxyz(:,:,iz50),3);
iz75 = ceil(nz*3/4)+[-2 3]; meanVxyz_075 = mean(meanVxyz(:,:,iz75),3);

% Azimuthal velocity in xz plane
if mod(ny+1,2)==0 % uneven number of ybins
    nyc = ceil(ny/2); yind = (nyc-2):(nyc+3);
else % even number of bins
    nyc = ceil(ny/2); yind = (nyc-1):(nyc+3);
end

if 1 % Illustrate the expected ExB drift           
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);
    surf(hca,rGrid,zGrid,surfa',ExBdrift')
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')    
    caxis(hca,2*vlim*[-1 1]);
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) ;
    title(hca,'ExB-drift')
    box(hca,'on')
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);    
    surf(hca,rGrid,zGrid,surfa',meanVrz'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'km/s');    
    caxis(hca,2*vlim*[-1 1])           
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
    title(hca,'Average azimuthal velocity')
    box(hca,'on')
end
if 1 % Plot average azimuthal velocity, xz-plane
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nx+1);    
    surf(hca,xGrid,zGrid,surfa,squeeze(nanmean(vazmat2(:,yind,:),2))'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'x')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'10^3 km/s')
    vmax = abs(max(max(vazmat2xzCenter)));
    if ~isnan(vmax); caxis(hca,vlim*2*[-1 1]); end
    title(hca,['Azimuthal velocity, z = [' num2str(yGrid(yind(1)),'%.1f') ' ' num2str(yGrid(yind(end)),'%.1f') '] km'])
    %caxis(hca,20*[-1 1]*1e3)
    %caxis(hca,vlim*[-1 1]) 
    %axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
end
if 1 % Plot mass coverage
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nr+1,nz+1);   
    normnumbrz = repmat(rSurf',1,numel(zSurf));
    surf(hca,rGrid,zGrid,surfa',sumMrz')%sumMrz'./normnumbrz')    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    chca = colorbar('peer',hca);
    ylabel(chca,'mass');         
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
    title(hca,'Mass distribution')
    box(hca,'on')
end
if 1 % Starting velocity histograms
    hca = h(isub); isub = isub + 1;
    hold(hca,'on');
    
    % v_perp from all
    hist(hca,sqrt(vx0.^2+vy0.^2)*1e-3,30);
    set(hca,'xtick',[0 0.5 1 1.5]*1e2)
        
    % v_z from flow
    iFlow = (nParticlesBox+1):nParticles;
    hist(hca,vz0(iFlow)*1e-3,30); %hold(hca,'on');
    
    % v_z from box
    iBox = 1:nParticlesBox;
    hist(hca,vz0(iBox)*1e-3,30); 
        
    % change transparency
    hp = findobj(hca,'Type','patch');
    set(hp(1),'FaceColor','r','EdgeColor','w','facealpha',0.5)
    set(hp(2),'FaceColor','b','EdgeColor','w','facealpha',0.5)
    set(hp(3),'FaceColor','g','EdgeColor','w','facealpha',0.5)
    
    % add fit from distribution
    histxlim = [-1 1]*1e2;
    v = linspace(histxlim(1),histxlim(2),100);
    f_flow = @(v,vt,veh) abs((v/vt/sqrt(pi)).*exp(-((v-veh)/vt).^2)); 
    
    %f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
    f_box = @(v,vt,veh) (1/vt/sqrt(pi)).*exp(-((v-veh)/vt).^2);    
    
    % flow    
    [bincountsFlow,binpositionsFlow] = hist(hca,vz0(iFlow),30); %hold(hca,'on');
    binwidth = binpositionsFlow(2) - binpositionsFlow(1);
    histarea = binwidth*sum(bincountsFlow);    
    plot(hca,v,f_flow(v*1e6,vtpar*1e3,-veh*1e3)*histarea*1e-4,'b'); 
    % box
    [bincountsBox,binpositionsBox] = hist(hca,vz0(iBox),30); %hold(hca,'on');
    binwidth = binpositionsBox(2) - binpositionsBox(1);
    histarea = binwidth*sum(bincountsBox);
    %v = linspace(binpositions(1),binpositions(end),100)*vtpar*1e-3-veh*1e-3;
    plot(hca,v,f_box(v*1e6,vtpar*1e3,-veh*1e3)*histarea*1e3,'r'); %hold(hca,'off');        
    
    % General information
    title(hca,'Starting velocities')
    xlabel(hca,'v [10^3 km/s]'); ylabel(hca,'# of particles');
    legend(hca,'v_{\perp}','v_{z,flow}','v_{z,box}')
    set(hca,'xtick',[-1:0.2:1]*1e2,'xlim',[-1 1]*1e2)   
    
    box(hca,'on');
    hold(hca,'off');
    
    %vt = 50e3; % km/s
    %plot(hca,v,nParticles/10*f(v,vt*1e-3));
     %hold off
end
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
if 1 % Plot simulation Bz in xy-plane         
    hca = h(isub); isub = isub + 1;         
    simBz_rth = squeeze(mean(simBz,3))*1e9; % nT
    surf(hca,fX(:,:,1),fY(:,:,1),simBz_rth);
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[0 1],'ylim',rlim*[-1 1]) 
    title(hca,'simBz')
    xlabel(hca,'x')
    zlabel(hca,'y') 
    colorbar('peer',hca);
    caxis(hca,tdB*1*[-1 1])
    grid(hca,'off')
    box(hca,'on')
end
if 1 % histogram of starting radius
    hca = h(isub); isub = isub + 1;
    hist(hca,sqrt(x0.^2+y0.^2)*1e-3,30);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Perpendicular starting positions')
    xlabel(hca,'r [km]'); ylabel(hca,'# of particles');
end 
if 1 % histogram of initial mass distribution
    hca = h(isub); isub = isub + 1;
    hist(hca,m,30);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Mass')
    xlabel(hca,'r [km]'); ylabel(hca,'# of particles');
end 
if 1 % Plot average azimuthal velocity, xy-plane, center about z = 0
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nx+1,ny+1);    
    surf(hca,xGrid,yGrid,surfa,meanVxyz_050'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'x')
    ylabel(hca,'y')
    chca = colorbar('peer',hca);
    ylabel(chca,'v_{\phi} [km/s]')        
    title(hca,['Azimuthal velocity, z = [' num2str(zGrid(iz50(1)),'%.1f') ' ' num2str(zGrid(iz50(end)+1),'%.1f') '] km'])    
    caxis(hca,4*vlim*[-1 1])  
    axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    box(hca,'on')
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
end
if 1 % Plot average azimuthal velocity, xy-plane, center about z = higher up
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nx+1,ny+1);    
    surf(hca,xGrid,yGrid,surfa,meanVxyz_075'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'x')
    ylabel(hca,'y')
    chca = colorbar('peer',hca);
    ylabel(chca,'v_{\phi} [km/s]')        
    title(hca,['Azimuthal velocity, z = [' num2str(zGrid(iz75(1))) ' ' num2str(zGrid(iz75(end)+1)) '] km'])    
    caxis(hca,4*vlim*[-1 1])  
    axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    box(hca,'on')
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
end
        