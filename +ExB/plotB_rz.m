% plot magnetic field in rz binning system
%% Compare the magnetic fields
if 1
% Plot magnetic field in red
quiver3(fX(rind,thind,zind)*1e-3,fY(rind,thind,zind)*1e-3,fZ(rind,thind,zind)*1e-3,...
        modBx(rind,thind,zind)*1,modBy(rind,thind,zind)*1,modBz(rind,thind,zind),'r'); hold on;
%%
    % Plot magnetic field in blue
quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
        simBx(xind,yind,zind)*1,simBy(xind,yind,zind)*1,simBz(xind,yind,zind),'b'); hold off;
end
%% Plot difference
if 1
diffBx = modBx-simBx; diffBy = modBy-simBy; diffBz = modBz-simBz;

quiver3(fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,fZ(xind,yind,zind)*1e-3,...
        diffBx(xind,yind,zind)*1,diffBy(xind,yind,zind)*1,diffBz(xind,yind,zind),'b'); hold off;
end
%% Plot surface plots of B next to each other with same colorbar

fig = figure(13);
for k = 1:6; h(k) = subplot(2,3,k); end
set(fig,'position',[839   515   838   441]); % not full screen
isub = 1;
%irf_colormap('poynting'); 
simB = sqrt(simBx.^2+simBy.^2+simBz.^2);
simBp = sqrt(simBx.^2+simBy.^2);
modB = sqrt(modBx.^2+modBy.^2+modBz.^2);
modBp = sqrt(modBx.^2+modBy.^2);

% Theoretical dB at the origin
tdB = e*phi0*n*1e6*mu0/(B0*1e-9)*art2.g(0.999999*lr/lz)*1e9; % nT

irf_colormap('poynting')
% From simulation
if 1 % Plot Bz in rz-plane     
    hca = h(isub); isub = isub + 1;         
    modBz_rz = squeeze(mean(modBz,2))*1e9; % nT
    %surfarea = zeros(numel(r_rz),numel(z_rz));
    %dr = diff(rSurf(1:2)); r_rz = [rSurf(rind)-dr/2 rSurf(rind(end))+dr/2];
    %dz = diff(zSurf(1:2)); z_rz = [zSurf(zind)-dz/2 zSurf(zind(end))+dz/2];
    %z_rz = zGrid(zind);
    %pcolor(hca,r_rz,z_rz,modBz_rz');
    surf(hca,rg,zg,zeros(rbin+1,zbin+1)',modBz_rz');
    shading(hca,'flat')
    view(hca,[0 0 1])
    zzlim = z_rz([1 end]);
    rrlim = r_rz([1 end]);
    set(hca,'xlim',rrlim,'ylim',zzlim) 
    title(hca,'simBz')
    xlabel('r')
    zlabel('z')    
end
if 0 % Plot Bz in rth-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,squeeze(fX(xind,yind,zind))*1e-3,squeeze(fY(xind,yind,zind))*1e-3,squeeze(simBz(xind,yind,zind)))
    modBz_rth = squeeze(mean(Bz(rind,thind,zind),3))*1e9; % nT
    surfarea = zeros(numel(r_rz),numel(th_rz));
    dr = diff(rSurf(1:2)); r_rz = [rSurf(rind)-dr/2 rSurf(rind(end))+dr/2];
    dth = diff(thSurf(1:2)); th_rz = [thSurf(thind)-dth/2 thSurf(thind(end))+dth/2];
    %z_rz = zGrid(zind);
    %pcolor(hca,r_rz,z_rz,modBz_rz');
    xxx = squeeze(mean(fX(rind,thind,zind),3))*1e-3;
    yyy = squeeze(mean(fY(rind,thind,zind),3))*1e-3;
    surf(hca,xxx,yyy,surfarea,modBz_rth);
    shading(hca,'flat')
    view(hca,[0 0 1])
    zzlim = z_rz([1 end]);
    rrlim = r_rz([1 end]);
    %set(hca,'xlim',rrlim,'ylim',rrlim) 
    title(hca,'simBz')
    xlabel('r')
    zlabel('th')
    %ch(isub-1)=colorbar('peer',hca);
end
if 0 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simB(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simB(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simB')
end
% From model
if 0 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBz(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBz(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'modBz')
end
if 0 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBp(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modBp(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'modBperp')
end
if 0 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modB(xind,yind,zind)))
    pcolor(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(modB(xind,yind,zind))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'modB')
end

ch=colorbar;
ylabel(ch,'B [nT]')
set(ch, 'Position', [.8514 .11 .0381 .8150])
if 1
for ii=1:6
      pos=get(h(ii), 'Position');
      set(h(ii), 'Position', [0.9*pos(1) pos(2) 0.82*pos(3) pos(4)]);
      caxis(h(ii),1*[-1 1]*tdB)
end
end
%ylabel(h(1),'y')
%ylabel(h(4),'y')
%xlabel(h(4),'x')
%xlabel(h(5),'x')
%xlabel(h(6),'x')