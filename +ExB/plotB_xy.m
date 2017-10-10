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
if 1 % Plot Bz in xy-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,squeeze(fX(xind,yind,zind))*1e-3,squeeze(fY(xind,yind,zind))*1e-3,squeeze(simBz(xind,yind,zind)))
    pcolor(hca,squeeze(fX(:,:,1)),squeeze(fY(:,:,1)),squeeze(simBz(:,:,ceil(zbin/2)))*1e9);
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simBz')
    %ch(isub-1)=colorbar('peer',hca);
end
if 1 % Plot Bperp in xy-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simBp(xind,yind,zind)))
    pcolor(hca,fX(:,:,1),fY(:,:,1),squeeze(simBp(:,:,ceil(zbin/2)))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simBperp')
end
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simB(xind,yind,zind)))
    pcolor(hca,fX(:,:,1),fY(:,:,1),squeeze(simB(:,:,ceil(zbin/2)))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simB')
end
% From model
if 1 % Plot Bz in xy-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,squeeze(fX(xind,yind,zind))*1e-3,squeeze(fY(xind,yind,zind))*1e-3,squeeze(simBz(xind,yind,zind)))
    pcolor(hca,squeeze(fX(:,:,1)),squeeze(fY(:,:,1)),squeeze(modBz(:,:,ceil(zbin/2)))*1e9);
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simBz')
    %ch(isub-1)=colorbar('peer',hca);
end
if 1 % Plot Bperp in xy-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simBp(xind,yind,zind)))
    pcolor(hca,fX(:,:,1),fY(:,:,1),squeeze(modBp(:,:,ceil(zbin/2)))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simBperp')
end
if 1 % Plot ExB_y electric field in xz-plane     
    hca = h(isub); isub = isub + 1;     
    %surf(hca,fX(xind,yind,zind)*1e-3,fY(xind,yind,zind)*1e-3,squeeze(simB(xind,yind,zind)))
    pcolor(hca,fX(:,:,1),fY(:,:,1),squeeze(modB(:,:,ceil(zbin/2)))*1e9)
    shading(hca,'flat')
    view(hca,[0 0 1])
    set(hca,'xlim',rlimproc*rlim*[-1 1],'ylim',rlimproc*rlim*[-1 1]) 
    title(hca,'simB')
end

ch=colorbar;
ylabel(ch,'B [nT]')
set(ch, 'Position', [.8514 .11 .0381 .8150])
for ii=1:6
      pos=get(h(ii), 'Position');
      set(h(ii), 'Position', [0.9*pos(1) pos(2) 0.82*pos(3) pos(4)]);
      caxis(h(ii),1*[-1 1]*tdB)
end
ylabel(h(1),'y')
ylabel(h(4),'y')
xlabel(h(4),'x')
xlabel(h(5),'x')
xlabel(h(6),'x')