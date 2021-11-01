a = 2.8;
b = 1.01;
m = 2;
n = 1;
B0 = 1;

w = 1;
t = 0;
k = 0.1;
l = 2*pi/k;
l = 2*pi;

x = linspace(0,a,11);
y = linspace(0,b,6);
z = linspace(0,2*l,11);
z = 1.5;

[X,Y,Z] = meshgrid(x,y,z);

% vecB = vecB0(x,y)*sin(kz-wt) =
%      =    xhat*Bx0(x,y)*sin(kz-wt) 
%         + yhat*By0(x,y)*sin(kz-wt) 
%         + zhat*Bz0(x,y)*sin(kz-wt)
%
%
% Bz = B0z(x,y)*sin(kx-wt)
Bz0 = @(x,y) B0*cos(m*pi*x/a).*cos(n*pi*y/b);
Bz = @(x,y,z,t) Bz0(x,y).*sin(k*z-w*t);

Bx = @(x,y,z,t) x.*0;
By = @(x,y,z,t) x.*0;

BX = Bx(X,Y,Z,t);
BY = By(X,Y,Z,t);
BZ = Bz(X,Y,Z,t);

hca = subplot(1,1,1);

quiver3(hca,X,Y,Z,BX,BY,BZ)
hca.XLabel.String = 'x';
hca.YLabel.String = 'y';
hca.ZLabel.String = 'z';

hold(hca,'on')
S = abs(BZ(:))*50;
S(S==0) = nan;
scatter3(X(:),Y(:),Z(:),S)
hold(hca,'off')

hca.Title.String = sprintf('m = %g, n = %g',m,n);



