nrows = 1;
ncols = 3;
for ip = 1:nrows*ncols
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;

  
clim = [0 5];

%%%
hca = h(isub); isub = isub + 1;

f_ = f01;
f_ = f_.convertto('s^3/km^6');
f = f_;
flim = 0.5e3;
f.data(f.data < flim) = NaN;
f.data(f.data > 1e5) = NaN;
[vx,vy,vz] = f.v('squeeze');
step = 1;

plx = vx(1:step:end)*1e-3;
ply = vy(1:step:end)*1e-3;
plz = vz(1:step:end)*1e-3;
plc = log10(f.data(1:step:end));
pls = plc*0+60;

scatter3(hca,plx,ply,plz,pls,plc,'filled'); 
axis(hca,'equal')
hca.XLabel.String = 'X';
hca.YLabel.String = 'Y';
hca.ZLabel.String = 'Z';
hca.XLim = [-10 10];
hca.YLim = [-10 10];
hca.ZLim = [-10 10];

%%%

hca = h(isub); isub = isub + 1;

f = f_;
[x,y,z] = f.xyz('squeeze');
step = 1;
energylevel = 10;
velocity = sqrt(f.depend{1}(energylevel)*units.eV*2/units.me)/1000;

plx = x(1:step:end)*1e-3*velocity;
ply = y(1:step:end)*1e-3*velocity;
plz = z(1:step:end)*1e-3*velocity;
plc = log10(squeeze(f.data(1,energylevel,:)));
pls = plc*0+40;

plX = x*1e-3*velocity;
plY = y*1e-3*velocity;
plZ = z*1e-3*velocity;
plC = log10(squeeze(f_.data(1,energylevel,:,:)));


scatter3(hca,plx,ply,plz,pls,plc,'filled'); 
%surf(hca,plX,plY,plZ,plC);  shading(hca,'flat');
axis(hca,'equal')
hca.XLabel.String = 'X';
hca.YLabel.String = 'Y';
hca.ZLabel.String = 'Z';
hca.XLim = [-10 10];
hca.YLim = [-10 10];
hca.ZLim = [-10 10];

hcb = colorbar('peer',hca);
hca.CLim = clim;
colormap(hca,'jet')


%%%
hca = h(isub); isub = isub + 1;
mms.plot_projection(hca,f,'clim',clim);
hca.XLim = [-10 10];
hca.YLim = [-10 10];
colormap(hca,'jet')

%%
hca = h(isub); isub = isub + 1;
fisovalue = 4e3;
axes(hca)
fv = isosurface(vx,vy,vz,squeeze(f.data),fisovalue,'noshare');
[F,V] = isosurface(vx,vy,vz,(squeeze(f.data)),fisovalue);
h = trisurf(F,V(:,1),V(:,2),V(:,3),'FaceColor',[0.7,0.7,0.7], 'EdgeColor', 'none');
axis(hca,'equal')
%patch(hca,isosurface(vx,vy,vz,squeeze(f.data),fisovalue)) 
% fv = isosurface(V,isovalue)
% fvc = isosurface(...,colors)
% fv = isosurface(...,'noshare')
% fv = isosurface(...,'verbose')
% [f,v] = isosurface(...)
% [f,v,c] = isosurface(...)
% isosurface(...)

%% Example for PDist.v

f = ePDist(100).convertto('s^3/km^6');
f.data(f.data < 2e3) = NaN;
[vx,vy,vz] = f.v('squeeze');
dotsize = 50;
scatter3(vx(:)*1e-3,vy(:)*1e-3,vz(:)*1e-3,f.data(:)*0+dotsize,log10(f.data(:)),'filled'); 
axis equal; colorbar;
vlim = [-5 5]; clim = [3 5];
set(gca,'clim',clim,'xlim',vlim,'ylim',vlim,'zlim',vlim)

%% Example for PDist.v

f = ePDist(100).convertto('s^3/km^6');
f.data(f.data < 2e3) = NaN;
[vx,vy,vz] = f.v('squeeze');
dotsize = 50;
scatter3(vx(:)*1e-3,vy(:)*1e-3,vz(:)*1e-3,f.data(:)*0+dotsize,log10(f.data(:)),'filled'); 
axis equal; colorbar;
vlim = [-5 5]; clim = [3 5];
set(gca,'clim',clim,'xlim',vlim,'ylim',vlim,'zlim',vlim)


