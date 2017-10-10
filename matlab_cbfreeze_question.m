%% Plot
n_plots = 3;
x = -10:10; y = -10:10; [X,Y] = meshgrid(x,y);
for k = 1:n_plots; 
    h(k) = subplot(3,1,k); 
    pcolor(h(k),X,Y,X+k); 
    ch(k) = colorbar('peer',h(k));
    ylabel(ch(k),['x+' num2str(k)]);
    c_lims(k,:) = get(h(k),'clim'); 
end
c_lim = [min(min(c_lims)) max(max(c_lims))];
for k = 1:n_plots; set(h(k),'clim',c_lim); end

%
doIRF = 1;
% Change colormap of first panel
if doIRF
    irf_colormap(h(1:2),'poynting')
else    
    apply_cmap_to = h(1);

    all_plot = findall(gcf,'type','axes','tag','');
    all_cbar = findobj(gcf,'tag','Colorbar');
    active_plot = apply_cmap_to;
    active_cbar = findobj(gcf, 'Type', 'axes', 'Tag', 'Colorbar','Axes', active_plot);
    nonactive_plot = all_plot; nonactive_plot(find(nonactive_plot==active_plot)) = [];
    nonactive_cbar = all_cbar; nonactive_cbar(find(nonactive_cbar==active_cbar)) = [];

    for ii = 1:numel(nonactive_plot); freezeColors(nonactive_plot(ii)); end   
    for pp = 1:numel(nonactive_cbar), % workaround cbfreeze bug that cbfreeze removes cblabel
        hcb = nonactive_cbar(pp);
        hy=get(hcb,'ylabel');
        ylabel_string=get(hy,'string');
        ylabel_fontsize=get(hy,'fontsize');
        cbar_position = get(hcb,'position');
        new_hcb = cbfreeze(hcb);
        new_hy=get(new_hcb,'ylabel');        
        set(new_hy,'string',ylabel_string,'fontsize',ylabel_fontsize);
        set(new_hcb,'position',cbar_position);
    end
    colormap(active_plot,'hot');
end
%% Change colormap of second panel
if doIRF
    irf_colormap(h(1),'default')
else    
    apply_cmap_to = h(2);

    all_plot = findall(gcf,'type','axes','tag','');
    all_cbar = findobj(gcf,'tag','Colorbar');
    active_plot = apply_cmap_to;
    active_cbar = findobj(gcf, 'Type', 'axes', 'Tag', 'Colorbar','Axes', active_plot)
    nonactive_plot = all_plot; nonactive_plot(find(nonactive_plot==active_plot)) = [];
    nonactive_cbar = all_cbar; nonactive_cbar(find(nonactive_cbar==active_cbar)) = [];

    for ii = 1:numel(nonactive_plot); freezeColors(nonactive_plot(ii)); end   
    for pp = 1:numel(nonactive_cbar), % workaround cbfreeze bug that cbfreeze removes cblabel
        hcb = nonactive_cbar(pp);
        hy=get(hcb,'ylabel');
        ylabel_string=get(hy,'string');
        ylabel_fontsize=get(hy,'fontsize');
        cbar_position = get(hcb,'position');
        new_hcb = cbfreeze(hcb);
        new_hy=get(new_hcb,'ylabel');        
        set(new_hy,'string',ylabel_string,'fontsize',ylabel_fontsize);
        set(new_hcb,'position',cbar_position);
    end
    colormap(active_plot,'hot');
end