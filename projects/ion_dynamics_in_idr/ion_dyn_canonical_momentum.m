%% Canonical momentum for Hall magnetic field
units = irf_units;

% Should lz be a function of x?
lz = @(x,lx,lz0) -lz0*(x/lx).^1;
%lz = @(x,lx,lz0) X*0 + lz0;
Ax = @(x,z,lx,lz0,B0) (B0./lz(x,lx,lz0)).*exp(-z.^2./lz(x,lx,lz0).^2);
vx = @(x,z,lx,lz0,B0,q,m) - (q/m)*Ax(x,z,lx,lz0,B0);
vx = @(x,z,lx,lz0,B0) vx(x,z,lx,lz0,B0,units.e,units.mp);

xvec = linspace(-2000,2000,500);
zvec = linspace(-300,300,400);
[X,Z] = ndgrid(xvec,zvec);


B0 = 5;
lx = 1000;
lz0 = 100;
LZ = lz(X,lx,lz0);
AX = Ax(X,Z,lx,lz0,B0); 
[dxAx,dzAx] = gradient(AX,X,Z); % By = dzAx - dxAz
By = dzAx;

nRows = 3;
nCols = 1;
h = setup_subplots(nRows,nCols);
isub = 1;


if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,LZ)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'l_z (km)';
end

if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,AX)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'A_x (...)';
end

if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,X,Z,BY)
  shading(hca,'flat')
  hcb = colorbar(hca);
  hcb.YLabel.String = 'B_y (...)';
end

hlinks = linkprop(h,{'XLim','YLim'});
colormap(pic_colors('blue_red'))