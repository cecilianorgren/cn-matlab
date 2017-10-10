

tints = toepoch([2007 08 31 10 15 50;...
                 2007 08 31 10 17 10;...
                 2007 08 31 10 17 30;...
                 2007 08 31 10 17 50;...
                 2007 08 31 10 19 25]);

hiaPSD4 = c_caa_distribution_data('C3_CP_CIS-HIA_HS_MAG_IONS_PSD');
%%
nPanels = numel(tints);
for oo = 1:nPanels
    h(oo) = subplot(2,3,oo);
    hca = h(oo);
    c_caa_plot_distribution_function(hca,'tint',tints(oo)+[0 12],'polar',hiaPSD4);
    caxis(hca,[4.1 6.4]);
    irf_colormap('space')
end

%%
nPanels = 20;
h = irf_plot(nPanels,'newfigure');
for oo = 1:nPanels
    hca = h(oo);    
    pcolor(hca,log(hiaPSD4.p(:,:,32-oo))');    
    hc(oo) = colorbar('peer',hca);
    shading(hca,'flat')
    caxis(hca,[4.1 6.4]);
    irf_colormap('space')
end
        
