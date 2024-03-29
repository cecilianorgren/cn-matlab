% Reduced distribution along separatrix
twpe = 24000;
pic = no02m.twpelim(twpe);
ds = ds100.twpelim(twpe).update_inds({430:522});
ds = ds100.twpelim(24000).update_inds({358:418});
xdist1 = (ds.xi1{1}+ds.xi2{1})/2;
zdist1 = (ds.zi1{1}+ds.zi2{1})/2;
Bx1_ = pic.Bx;
By1_ = pic.By;
Bz1_ = pic.Bz;
Bx1 = interpfield(pic.xi,pic.zi,Bx1_,xdist1,zdist1); 
By1 = interpfield(pic.xi,pic.zi,By1_,xdist1,zdist1); 
Bz1 = interpfield(pic.xi,pic.zi,Bz1_,xdist1,zdist1);
fred35_A75_1 = ds.reduce_1d_new('x',[3 5],[],'vpar',{Bx1,By1,Bz1},'pitch',{Bx1,By1,Bz1});
fred3_A75_1 = ds.reduce_1d_new('x',[3],[],'vpar',{Bx1,By1,Bz1},'pitch',{Bx1,By1,Bz1});
fred46_A75_1 = ds.reduce_1d_new('x',[4 6],[],'vpar',{Bx1,By1,Bz1},'pitch',{Bx1,By1,Bz1});
fred4_A75_1 = ds.reduce_1d_new('x',[4],[],'vpar',{Bx1,By1,Bz1},'pitch',{Bx1,By1,Bz1});
arclength1 = [0 cumsum(sqrt(diff(xdist1).^2 + diff(zdist1).^2))];
arc01 = arclength1(find(abs(zdist1)==min(abs(zdist1))));
arclength1 = arclength1 - arc01;
darcs1 = diff(arclength1);
arcedges1 = [arclength1(1)-0.5*darcs1(1) arclength1(1:end-1)+0.5*darcs1 arclength1(end)+0.5*darcs1(end)];




%% Figure
colors = pic_colors('matlab');

nrows = 2;
ncols = 1;
%h = setup_subplots(nrows,ncols,'vertical');
h = setup_subplots(nrows,ncols,'horizontal');
%[h,h2] = initialize_combined_plot(nrows,2,2,0.4,'vertical')
isub = 1;
doE = 0; colorE = [0.0 0.0 0.0];
doV = 0; colorV = 0*[1 1 1];
doN = 0; colorN = [0 0 0];
doExB = 0; colorExB = 0*[1 1 1]+0.5;
isMap = [];
xlimMap = mean(pic.xi) + [-30 30];
zlimMap = [-10 10];


if 1 % Epar
  hca = h(isub); isub = isub + 1;
  isMap(end + 1) = isub - 1;
  pic.xlim(xlimMap).zlim(zlimMap).plot_map(hca,{'Epar'});
  shading(hca,'flat')
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'z (d_i)';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_{||}';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-0.5 0.5];  
  hca.XGrid = 'on';
  %hca.YGrid = 'on';
  hca.Layer = 'top';
  
  hold(hca,'on')
  ds.plot_boxes(hca)
  hold(hca,'off')
end

if 1 % log 10 fi46(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred46_A75_1;  
  %surf(hca,arcedges1,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',log10(fred.fvpar)')
  dx = diff(fred.x);
  x_edges = [fred.x(1:end-1)-0.5*dx; fred.x(end)-0.5*dx(end); fred.x(end)+0.5*dx(end)];
  surf(hca,x_edges,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'x (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{ec}'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
 % hca.YLim = [-2.5 2.5];
  %hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 0 % fi4/fi46(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred46_A75_1;
  fred = fred4_A75_1;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges1,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength1,fred.vpar_center,log10(fredall.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  %hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end

if 0 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75_1;  
  surf(hca,arcedges1,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 0 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred35_A75_1;
  fred = fred3_A75_1;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges1,fred.vpar_edges,zeros(numel(arcedges1),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength1,fred.vpar_center,log10(fredall.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end

if 0 % Epar
  hca = h(isub); isub = isub + 1;
  pcolor(hca,ss,times,smooth2(stE,2))
  shading(hca,'flat')
  hca.XLabel.String = 's';
  hca.YLabel.String = 'twpe';
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'E_{||}';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = [-0.5 0.5];  
  hca.XGrid = 'on';
  %hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YDir = 'reverse';
end

if 0 % log 10 fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fred = fred35_A75_2;  
  surf(hca,arcedges2,fred.vpar_edges,zeros(numel(arcedges2),numel(fred.vpar_edges))',log10(fred.fvpar)')
  view(hca,[0 0 1]); 
  shading(hca,'flat')
  
  if 0 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength2,fred.vpar_center,log10(fred.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 'arclength (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('candy4'))  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = ['log_{10}f_{ic}'];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2.5 2.5];
  hca.CLim = [-6 -1];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)    
    %plot(hca,arclength_interp,smooth(Epar,100)*max(abs(hca.YLim))/max(abs(smooth(Epar,100))),'color',colorE+0.2)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
if 0 % fi3/fi35(vpar)
  hca = h(isub); isub = isub + 1;
  fredall = fred35_A75_2;
  fred = fred3_A75_2;
  fplot = (fred.fvpar./fredall.fvpar);
  %pcolor(hca,arclength,fred.vpar_center,fred.fvpar')
  surf(hca,arcedges2,fred.vpar_edges,zeros(numel(arcedges2),numel(fred.vpar_edges))',fplot')
  view(hca,[0 0 1]);   
  shading(hca,'flat')
  if 1 % contour level for reference
    hold(hca,'on')
    clim = hca.CLim;
    contour(hca,arclength2,fred.vpar_center,log10(fredall.fvpar)',[-2 -2],'k')
    hca.CLim = clim;
    hold(hca,'off')
  end
  hca.XLabel.String = 's_{||} (d_i)';
  hca.YLabel.String = 'v_{||}';
  colormap(hca,pic_colors('pasteljet'))
  hcb = colorbar('peer',hca);  
  hcb.YLabel.String = ['f_{ic}^{top}/f_{ic}^{tot}'];
  %hca.CLim(2) = prctile(fred.fvpar(:),99);
  %hca.CLim = fi_clim;
  hca.CLim = [0 1];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  hca.Layer = 'top';
  hca.YLim = [-2 2];
  if doE
    hold(hca,'on')
    %plot(hca,arclength,Epar*max(abs(hca.YLim))/max(abs(Epar)),'color',colorE)
    plot(hca,arclength_interp,Epar_*max(abs(hca.YLim))/max(abs(Epar_)),'color',colorE+0.2)
    hold(hca,'off')
  end
end
compact_panels(0.01)

hlinksx = linkprop(h,{'XLim'});
%hlinksy = linkprop(h([1 2 4 5]),{'YLim'});