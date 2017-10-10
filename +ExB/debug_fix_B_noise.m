nPanels = 3;
for k = 1:nPanels;
    h(k) = subplot(nPanels,1,k);
end

isub = 1;
if 1
    hca = h(isub); isub = isub + 1;
    pcolor(hca,squeeze(VXfix(:,:,60)))
    colorbar('peer',hca)
end
if 1
    hca = h(isub); isub = isub + 1;
    pcolor(hca,squeeze(VXsim(:,:,60)))
    colorbar('peer',hca)
end
if 1
    hca = h(isub); isub = isub + 1;
    pcolor(hca,squeeze(VXmod(:,:,60)))
    colorbar('peer',hca)
    caxmax = [-1 1]*max(get(hca,'clim'));
end



for  k = 1:nPanels
    shading(h(k),'flat')
    caxis(h(k),'clim',caxmax)
end