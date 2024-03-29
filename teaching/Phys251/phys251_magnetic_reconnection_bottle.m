a = 5;
b = 1;
x = a*linspace(-15,15,200);
y = b*linspace(-15,15,100);
z = linspace(-10,10,5);
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%dz = z(2) - z(1);
x_xline = x;
y_xline = x*b/a;

Ay = @(x,y) (x/a).^2 - (y/b).^2;
AY0 = Ay(X,Y);

%[FX,FY] = gradient(AY,dx,dy);
%Bx = -FX;
%By = FY;

colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];

t = 0;
it = 1;
Astep = 20;
dA = Astep/numel(t);
AYlev0 = -200:Astep:(200 + Astep);

hca = subplot(1,1,1);

% Draw separatrix
  plot(hca,x_xline,y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  hold(hca,'on')
  plot(hca,x_xline,-y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
 
  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.5:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.5:numel(S(is).Y));%S(is).Y;
    plot(hca,sx,sy,'color',0.8+[0 0 0],'linewidth',1)   
  end
  
  i0 = 15;
  sx1 = [S(i0).X S(i0+2).X(end:-1:1)];
  sy1 = [S(i0).Y S(i0+2).Y(end:-1:1)];
  patch(hca,sx1,sy1,[0.8500    0.3250    0.0980])
  
  i0 = 23;
  sx2 = [S(i0).X S(i0+2).X(end:-1:1)];
  sy2 = [S(i0).Y S(i0+2).Y(end:-1:1)];
  patch(hca,sx2,sy2,[     0    0.4470    0.7410])  
  
  pause(0.1)
  drawnow
  hold(hca,'off')
  hca.XLim = x([1 end]);
  hca.YLim = y([1 end]);
  hca.Visible = 'off';
  hca.Position = [0 0 1 1];
