x = linspace(0,10,10);
y = linspace(0,10,10);
dx = x(2) - x(1);
dy = y(2) - y(1);
a = @(x,y) x.^2;
b = @(x,y) y.^2;

[X,Y] = meshgrid(x,y);
A = a(X,Y);
B = b(X,Y);

[Agx,Agy] = gradient(A,dx,dy);
[Bgx,Bgy] = gradient(B,dx,dy);

Bfield = cross([Agx(:), Agy(:), 0*Agy(:)],[Bgx(:), Bgy(:), 0*Agy(:)]);
Bx = reshape(Bfield(:,1),size(A));
By = reshape(Bfield(:,2),size(A));
Bz = reshape(Bfield(:,3),size(A));


nrows = 4;
ncols = 2;
isub = 1;

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,A)
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\alpha(x,y)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,B)
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\beta(x,y)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,Agx)
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\partial_x\alpha(x,y)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,Bgx)
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\partial_x\beta(x,y)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,Agy)
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\partial_y\alpha(x,y)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
pcolor(hca,X,Y,Bgy)
shading(hca,'flat')
hcb = colorbar('peer',hca);
hcb.YLabel.String = '\partial_y\beta(x,y)';

hca = subplot(nrows,ncols,isub); isub = isub + 1;
quiver3(hca,X,Y,X*0,Bx,By,Bz)
