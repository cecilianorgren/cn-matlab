lon = 0:1:360;
colat = 60:1:90;
RE = 6371; % km

[LON,COLAT] = ndgrid(lon,colat);


X = RE*cosd(LON).*cosd(COLAT);
Y = RE*sind(LON).*cosd(COLAT);
RHO = sqrt(X.^2 + Y.^2);

x = linspace(-3000,3000,100);
y = linspace(-3000,3000,100);

[X,Y] = ndgrid(x,y);

scale = RE*sind(10);
phi0 = 30; % keV
pot_conv = @(x,y) -1.0*phi0*(y/scale).*exp(-( (x.^2 + y.^2)/scale^2));
pot_corot = @(x,y) .30*(x.^2 + y.^2)/scale^2;

POT_conv = pot_conv(X,Y);
POT_corot = pot_corot(X,Y);
POT = POT_conv + POT_corot;

[X,Y] = ndgrid(x,y);
[EX,EY] = gradient(-POT',X',Y');
EX = EX';
EY = EY';

EX(EX==inf) = nan;

% Find saddle point  
[nx,ny,nz] = surfnorm(X,Y,POT);
[az,el,rho] = cart2sph(nx,ny,nz);   % find azimuth and elevation
[~,ix] = max(el(:));  
[~,ind]=max(abs(nz(:)));
hold(hca,'on')


fontsize = 12;

nrows = 3;
ncols = 1;
h = gobjects([nrows,ncols]);
ipanel = 1;
for irow = 1:nrows
  for icol = 1:ncols
    h(irow,icol) = subplot(nrows,ncols,ipanel); ipanel = ipanel + 1;
  end 
end
isub = 1;

philev = [-100:1:100]-mod(POT(ind),1)+0.5;

hca = h(isub); isub = isub + 1;
contourf(hca,Y,X,POT_conv',philev)
colormap(hca,pic_colors('blue_red'))
hcb = colorbar(hca);
hcb.Label.String = 'Convection potential';
hcb.Label.Rotation = 270;

hca = h(isub); isub = isub + 1;
contourf(hca,Y,X,POT_corot',philev)
colormap(hca,pic_colors('blue_red'))
hcb = colorbar(hca);
hcb.Label.String = 'Corotation potential';
hcb.Label.Rotation = 270;

if 0
hca = h(isub); isub = isub + 1;
contourf(hca,Y,X,POT',philev)
colormap(hca,pic_colors('blue_red'))
end
if 1
hca = h(isub); isub = isub + 1;
contourf(hca,Y,X,POT',philev)
colormap(hca,pic_colors('blue_red'))
hca.CLim = [-13 13]+POT(ind);
hcb = colorbar(hca);
hcb.Label.String = 'Total potential';
hcb.Label.Rotation = 270;

if 0
hold(hca,'on')
contour(X,Y,POT,POT(ind)*[1 1],'k:','linewidth',2)
hold(hca,'off')
end
end

if 0
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Y,EX)
colormap(hca,pic_colors('blue_red'))
end

%hca.CLim = [-13 13];
if 0
hca = h(isub); isub = isub + 1;
pcolor(hca,X,Y,EY*1e3)
colormap(hca,pic_colors('blue_red'))
end


c_eval('shading(h(?),''flat'')',1:numel(h))
%c_eval('hb(?) = colorbar(h(?));',1:numel(h))
linkprop(h(1:2),{'CLim'});
h(1).CLim = [-13 13];

c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).XLabel.String = ''x'';',1:numel(h))
c_eval('h(?).YLabel.String = ''y'';',1:numel(h))
c_eval('h(?).YLabel.Rotation = 270;',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('axis(h(?),''square'');',1:numel(h))
c_eval('h(?).Position(1) = 0.15;',1:numel(h))

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
c_eval('hb(?).FontSize = fontsize;',1:numel(h))
%hca.XDir = 'reverse';