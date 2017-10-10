% fluid
nx = 30; ny = 30;
xl = 20;
yl = 20
x = linspace(-xl,xl,nx);
y = linspace(-0,yl,ny);
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);
THETA = atan(Y./X);
tp = 1:2:nx;
%theta = linspace(0,2*pi,100);
%r = linspace(0,10,20);
%[R,THETA] = meshgrid(r,theta);
sf1 = @(A,x,y) A*x.*y;      vp1 = @(A,x,y) .5*A*(x.^2-y.^2);
sf2 = @(m,theta) m*theta;   vp2 = @(m,r) m*log(r);
sf12 = @(A,x,y,m) A*x.*y + m*atan(y./x);
A = 1;
m = 7;
SF = sf2(m,THETA) + sf1(A,X,Y);
VP = vp2(m,R)     + vp1(A,X,Y);
%[VX,VY] = gradient(VP,X,Y);
VX=A*X;
VY=-A*Y;
quiver(X(tp,tp),Y(tp,tp),VX(tp,tp),VY(tp,tp)); hold on;
%level = [-50:2:,-7 8:2:40]; % m=5
%level = [-95:4:-11 11:4:95]; % m=7
contour(X,Y,SF,40)
contour(X,Y,VP,20); 
axis equal
set(gca,'ylim',[0 yl],'xlim',[-xl xl])
title(['Fluid falling on a hump, m=' num2str(m) ' A=' num2str(A)])
hold off;

xzero = x;
yzero = x*0;
for k=1:nx
   yzero(k) = fzero(@(y) sf1(A,xzero(k),y)+sf2(m,atan(xzero(k)/y)),3);
end
hold on;
plot(xzero,yzero)

%%
nx = 30; ny = 30;
xl = 0.1;
yl = 0.03;
x = linspace(-xl,xl,nx);
y = linspace(-0,yl,ny);
[X,Y] = meshgrid(x,y);
R = sqrt(X.^2+Y.^2);
THETA = atan(Y./X);
tp = 1:2:nx;
m=-0.01;
v0=1.2;

theta= linspace(0,pi,20);
r = linspace(0,0.1,20);
[R,THETA]=meshgrid(r,theta);
X=R.*cos(THETA);
Y=R.*sin(THETA);

sp = -m/v0;

sf1 = @(v0,y,m,theta)(v0*y+m*theta);
vx = @(v0,m,r,theta) v0 + m*cos(theta)./r;  
vy = @(v0,m,r,theta) m*sin(theta)./r;  
VX = vx(v0,m,R,THETA);
VY = vy(v0,m,R,THETA);
VX(find(abs(VX)>5*v0))=nan;
VY(find(abs(VY)>5*v0))=nan;

%vx = @(v0,m,x,y) m*x/(x.^2+y.^2) + v0;
%vy = @(v0,m,x,y) m*y/(x.^2+y.^2);
%VX = vx(v0,m,X,Y);
%VY = vy(v0,m,X,Y);
plot(0,0)
subplot(2,1,1)
contour(X,Y,sf1(v0,Y,m,THETA),40); hold on;
quiver(X,Y,VX,VY,1)
plot(sp,0,'kx','linewidth',2,'markersize',15)
plot(0,0,'ro','linewidth',2,'markersize',10)
axis equal
set(gca,'ylim',[0 yl],'xlim',[-xl xl])
title(['Fluid sucked into sink, m=' num2str(m) ' v0=' num2str(v0)])
hold off
hold on
sl = @(y,v0,m) y./tan(-v0*y/m);
y=linspace(0,-pi*m/v0,200);
plot(sl(y,v0,m),y,'linewidth',2)
hold off
%%
hca=subplot(1,1,1);

contourf(X,Y,sf1(v0,Y,m,THETA),70); hold on;
quiver(X,Y,VX,VY,1,'k')
plot(sl(y,v0,m),y,'k','linewidth',2)
%axis equal
title(['Fluid sucked into sink, m=' num2str(m) ' v0=' num2str(v0)])
ch = colorbar;
ch.YLabel.String = 'Stream function';
hca.XLabel.String = 'x';
hca.YLabel.String = 'y';
axis equal
set(gca,'clim',0.05*[-1 1])
set(gca,'ylim',[0 yl],'xlim',[-xl*0.9 xl*0.9])
plot(sp,0,'kx','linewidth',2,'markersize',15)
plot(0,0,'ro','linewidth',2,'markersize',10)
colormap(hca,cn.cmap('bluered'))
hca.FontSize = 14;
hold off

%%
nx = 30; ny = 30;
xl = 10;
yl = 1  0;
x = linspace(-xl,xl,nx);
y = linspace(-yl,yl,ny);
[X,Y] = meshgrid(x,y);

a=3;
v0=1;
m=3;

sf1 = @(x,y,a) atan(y./(x+a)) -atan(y./(x-a));

sf1 = @(x,y,a) atan((x+a)./y) -atan((x-a)./y);
sf2 = @(x,y,a,m) m*atan(2*a*y./(a^2+x.^2+y.^2)) +v0*y;
SF1 = sf1(X,Y,a);
SF2 = sf2(X,Y,a,m);

vx  = @(x,y,m,v0) m*(x+a)./((x+a).^2+y.^2)-m*(x-a)./((x-a).^2+y.^2)+v0;
vy  = @(x,y,m,v0) m*y./((x+a).^2+y.^2)-m*y./((x-a).^2+y.^2);
VX  = vx(X,Y,m,v0);
VY  = vy(X,Y,m);

%h1=subplot(1,3,1);
contour(X,Y,SF1,30); hold on;
quiver(X,Y,VX,VY)
hc1=colorbar;

%h2=subplot(1,3,2);
%contour(X,Y,SF2,30)
%hc2=colorbar;

%h3=subplot(1,3,3);
%contour(X,Y,SF2-SF1)
%hc3=colorbar;

%set(h1,'clim',get(h2,'clim'))
hold off;
