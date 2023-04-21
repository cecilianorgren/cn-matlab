vtx = 1;
vty = 1;
vtz = 1;
vscale = max([vtx,vty,vtz]);

np = 600000;
vx = vtx*randn(np,1);
vy = vty*randn(np,1);
vz = vtz*randn(np,1);


vxbin = vscale*linspace(-3,3,50);
vybin = vscale*linspace(-3,3,51);
vzbin = vscale*linspace(-3,3,52);

[N,edges,mid,loc] = histcn([vx vy vz],vxbin,vybin,vzbin);


ipanel = 0;
nrows = 1;
ncols = 1;
h = gobjects([nrows,ncols]);
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;


if 1
  hca = h(isub); isub = isub + 1;
  pcolor(mid{1:2},squeeze(sum(N,3))')
  shading(hca,'flat')
  colormap(hca,irf_colormap('waterfall'))
  hca.XLabel.String = 'v_x';
  hca.YLabel.String = 'v_y';
end


%%
vtx = 1;
vty = 1;
vtz = 1;
vscale = max([vtx,vty,vtz]);

vxbin = vscale*linspace(-3,3,50);
vybin = vscale*linspace(-3,3,51);
vzbin = vscale*linspace(-3,3,52);

nps = [1e2 1e3 1e4 1e5 1e6];
nps = [1e2 1e3 1e4 1e5 ];


ipanel = 0;
nrows = numel(nps);
ncols = 3;
h = gobjects([nrows,ncols]);
for icol = 1:ncols, for irow = 1:nrows, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;


for np = nps
  vx = vtx*randn(np,1);
  vy = vty*randn(np,1);
  vz = vtz*randn(np,1);
  [N,edges,mid,loc] = histcn([vx vy vz],vxbin,vybin,vzbin);

  if 1 % f(vx,vy)
    hca = h(isub); isub = isub + 1;
    plot3(hca,vx,vy,vz,'.')
    %hs = scatter3(hca,vx,vy,vz,1,'marker','o','MarkerFaceColor','k','MarkerFaceAlpha',0.2,'MarkerEdgeAlpha',0.2);
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
    hca.ZLabel.String = 'v_z';
    hca.Box = 'on';
  end
  if 1 % f(vx,vy)
    hca = h(isub); isub = isub + 1;
    pcolor(hca,mid{1:2},squeeze(sum(N,3))')
    shading(hca,'flat')
    %colormap(hca,irf_colormap('waterfall'))
    cmap = pic_colors('blue_red');
    colormap(hca,flip(cmap(1:fix(end/2),:),1));
    hca.XLabel.String = 'v_x';
    hca.YLabel.String = 'v_y';
    hca.Title.String = sprintf('N_p = %g',np);
  end
  if 1 % f(vx)
    hca = h(isub); isub = isub + 1;
    plot(hca,mid{1},squeeze(sum(sum(N,3),2))')    
    hca.XLabel.String = 'v_x';    
  end
end

