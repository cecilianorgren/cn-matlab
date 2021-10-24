% EDI

n = 32;
ithB = 1;
iphiB = 26;
theta = (-n:2:n)/n*pi;
phi = (-n:2:n)'/n*pi/2;
cosphi = cos(phi); cosphi(1) = 0; cosphi(n+1) = 0;
sintheta = sin(theta); sintheta(1) = 0; sintheta(n+1) = 0;

X = cosphi*cos(theta);
Y = cosphi*sintheta;
Z = sin(phi)*ones(1,n+1);
isnan = X>0.01;
isnan = Z<-0;
X(isnan)=NaN;
Y(isnan)=NaN;
Z(isnan)=NaN;
C = Z*0;
C(iphiB,:) = 1;

c_eval('C(iphiB,ithB+?) = 0.2+0.1*?;',0:3);

%N = 12;
%[X,Y,Z] = sphere(N);
hs = surf(X,Y,Z,C);
hs.FaceAlpha = 0.8;

hca = hs.Parent;

hold(hca,'on')
phiB = phi(iphiB)*180/pi+180/32/2; % rad
thB = theta(ithB)*180/pi+180/32; % rad
%thB = 185; % deg
rB = 1.9;
xB = rB*cosd(phiB)*cosd(thB);
yB = rB*cosd(phiB)*sind(thB);
zB = rB*sind(phiB);
%plot3(hca,[0 xB],[0 yB],[0 zB],'color',[0 0 0])
hq = quiver3(hca,0,0,0,xB,yB,zB,'color',[0 0 0]);
hq.LineWidth = 6;
hold(hca,'off')

%view([0 0 1])

cmap = pic_colors('thermal');
colormap('gray');
colormap(cmap(20:end-20,:));
axis equal
hca.Visible = 'off';