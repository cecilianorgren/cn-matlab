a = 4; % a and b defines the aspect ratio of the magnetic field
b = 1;

% Used for plotting magnetic field lines with A 
xvec = a*linspace(-10,10,700);
zvec = b*linspace(-10,10,100);
yvec = linspace(-0,0,1);
[X,Y,Z] = ndgrid(xvec,yvec,zvec);
Ay = @(x,y,z) (x/a).^2 - (z/b).^2;
AY0 = Ay(X,Y,Z);
%xlim = [-20 20];
%zlim = [-5 5];
S = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[-105:10:110]); % B field lines
%SX = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[0 0]-.01);  % X line  
SX = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',(-40:40:40)-0.01);  % X line  
  

colors = pic_colors('matlab');
color_sep = colors(1,:);
color_sep = [0.8 0.8 0.8];
color_sep = [0.8 0.7 0.6];
color_sep = [0.9 0.7 0.3];
%color_sep = colors(3,:);
color_box = colors(3,:);
color_box = colors(2,:);
color_box = [0.6 0.3 0.8];

color_sep = [0.6 0.3 0.8];
color_box = [0.9 0.7 0.3];
facealpha_sep = 0.95;
facealpha_box = 0.7;

hca = subplot(1,1,1);

yy = [0 20];
xx = [-40 40];
zz = xx*b/a;

plot3(hca,[0 0],yy,[0 0],'k--','linewidth',2)
hold(hca,'on')


% Draw xz box
xdata = [xx(1) xx(2) xx(2) xx(1)];
ydata = [yy(1) yy(1) yy(2) yy(2)];
zdata = [zz(1) zz(1) zz(1) zz(1)];
hp = patch(hca,xdata,ydata,zdata,'r');
hp.FaceColor = color_box;
hp.FaceAlpha = facealpha_box;

xdata = [xx(2) xx(2) xx(2) xx(2)];
ydata = [yy(1) yy(2) yy(2) yy(1)];
zdata = [zz(1) zz(1) zz(2) zz(2)];
hp = patch(hca,xdata,ydata,zdata,'r');
hp.FaceColor = color_box;
hp.FaceAlpha = facealpha_box;


xdata = [xx(1) xx(2) xx(2) xx(1)];
ydata = [yy(1) yy(1) yy(2) yy(2)];
zdata = [zz(2) zz(2) zz(2) zz(2)];
hp = patch(hca,xdata,ydata,zdata,'r');
hp.FaceColor = color_box;
hp.FaceAlpha = facealpha_box;

xdata = [xx(1) xx(1) xx(1) xx(1)];
ydata = [yy(1) yy(2) yy(2) yy(1)];
zdata = [zz(1) zz(1) zz(2) zz(2)];
hp = patch(hca,xdata,ydata,zdata,'r');
hp.FaceColor = color_box;
hp.FaceAlpha = facealpha_box;

% Draw separatrices
for isep = 1:numel(SX)
  xdata = [SX(isep).X SX(isep).X(end:-1:1)];
  ydata = [repmat(yy(1),[1,numel(SX(isep).X)]),repmat(yy(2),[1,numel(SX(isep).X)])];
  zdata = [SX(isep).Y SX(isep).Y(end:-1:1)];
  hp = patch(hca,xdata,ydata,zdata,'r');
  hp.FaceColor = color_sep;
  hp.FaceAlpha = facealpha_sep;
end

axis equal

hca.XTick = [];
hca.YTick = [];
hca.ZTick = [];
hca.Visible = 'off';
hca.XDir = 'reverse';
hca.View = [-33.1000   14.8000];

hold(hca,'off')

hp = findobj(gcf,'type','patch'); hp = hp(end:-1:1);
c_eval('hp(?).LineWidth = 1.5;',1:numel(hp))
% export transparent background with export_fig from matlab central user library 
