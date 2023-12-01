%% One layer
a = 4;
b = 1;
x = a*linspace(-10,10,200);
y = b*linspace(-10,10,100);
z = 0;
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%dz = z(2) - z(1);

Ay = @(x,y) (x/a).^2 - (y/b).^2;
AY = Ay(X,Y);

[FX,FY] = gradient(AY,dx,dy);
Bx = -FX;
By = FY;
AYlev = linspace(min(AY(:)),0.99*max(AY(:)),10);

S = contourcs(x,y,AY,AYlev+2);

colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];
color = [0 0 0];
color_sep = [0 0 0];

hca = subplot(1,1,1);

x_xline = x;
y_xline = x*b/a;

linewidth = 2;
linewidth_sep = 2;

plot3(hca,0,0,0)
%plot3(hca,x_xline,x_xline*0,y_xline,'linewidth',linewidth_sep,'color',color_sep)
hold(hca,'on')
%plot3(hca,x_xline,x_xline*0,-y_xline,'linewidth',linewidth_sep,'color',color_sep)

%for iz = 1:numel(z)
  for is = 1:numel(S)
    %plot3(hca,S(is).X,S(is).Y*0+z(iz),S(is).Y,'color',color,'linewidth',linewidth)    
    xp = [S(is).X, S(is).X(end:-1:1)];
    yp = [S(is).Y, S(is).Y(end:-1:1)];
    zp = [S(is).Y*0, S(is).Y(end:-1:1)*2];
    patch(hca,xp,yp,zp)
  end
%end
hold(hca,'off')
axis(hca,'equal')
%axis(hca,'off')


