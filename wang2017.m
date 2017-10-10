%

b0 = 2;
bn = 0.04;
bg = 0.4;
n0 = -0.15;

AL = @(l,m,n) 0.5*bg*(n-n0).^2 - 0.5*bg*n0^2;
AM = @(l,m,n) 0.5*b0*n.^2 - bn*l;
AN = @(l,m,n) 0;
A = @(l,m,n) sqrt(AL(l,m,n).^2 + AM(l,m,n).^2 + AN(l,m,n).^2);



l = linspace(0,10,100);
n = linspace(-2,2,100);
[L,N] = meshgrid(l,n);
M = 0*L;

nrows = 3;
ncols = 1;
npanels = nrows*ncols;
for ip = 1:npanels;
  h(ip) = subplot(nrows,ncols,ip);
end
isub = 1;  

hca = h(isub); isub = isub + 1;
[c,h_] = contourf(hca,L,N,AL(L,M,N));
clabel(c,h_)
hca.XLabel.String = 'L';
hca.YLabel.String = 'N';
hca.Title.String = 'A_L';

hca = h(isub); isub = isub + 1;
[c,h_] = contourf(hca,L,N,AM(L,M,N));
clabel(c,h_)
hca.XLabel.String = 'L';
hca.YLabel.String = 'N';
hca.Title.String = 'A_M';

hca = h(isub); isub = isub + 1;
%pcolor(hca,L,N,AM(L,M,N)); shading(hca,'flat')
contourf(hca,L,N,AM(L,M,N));
hca.XLabel.String = 'L';
hca.YLabel.String = 'N';
hca.Title.String = '(A_L^2+A_M^2)^{1/2}';

%% 3D sphere quadrant
N=100;
[oX,oY,oZ] = sphere(N);
X = oX;
Y = oY;
Z = oZ;

h = surf(X,Y,Z,Z*0+1);
shading flat
h.FaceAlpha = 0.5;
axis square
hca = gca;
hca.XLim = [0 1];
hca.YLim = [0 1];
hca.ZLim = [0 1];
