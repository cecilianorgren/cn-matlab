% illustrate py

xvec = linspace(-2,2,20);
yvec = linspace(-2,2,20);
Ay = @(x,y) x + 0*y;
Bx = 1;

x = @(t,T) cos(t*2*pi/T);
y = @(t,T) -sin(t*2*pi/T);
vy = @(t,T) -cos(t*2*pi/T);
vx = @(t,T) -sin(t*2*pi/T);

x = @(t,T) sin(t*2*pi/T);
y = @(t,T) cos(t*2*pi/T);
vy = @(t,T) -sin(t*2*pi/T);
vx = @(t,T) sin(t*2*pi/T);

[X,Y] = ndgrid(xvec,yvec);

T = 1;
nt = 40;
t = linspace(0,1*T,nt);

fontsize = 12; 
vidObj = VideoWriter([printpath 'illustration_Ay_gyration'],'MPEG-4');
open(vidObj);

fig = figure(100);
fig.Color = [1 1 1];
h(1) = subplot(1,2,1);
h(2) = subplot(1,2,2);
for it = 1:nt
  hca = h(1);
  hp = plot(hca,x(t(it),T),y(t(it),T),'*k');
  hold(hca,'on')
  contourf(hca,X,Y,Ay(X,Y),10,'k')
  colormap(hca,irf_colormap('waterfall'))
  hold(hca,'off')
  hca.XLim = xvec([1,end]);
  hca.YLim = xvec([1,end]);
  axis(hca,'square')
  h(1).Children = h(1).Children(end:-1:1);
  hca.XLabel.String = 'x';
  hca.YLabel.String = 'y';
  hb = colorbar(hca);
  hb.YLabel.String = 'A_y';
  hb.Location = 'manual';
  hca.FontSize = fontsize;
  hca.Title.String = 'A_y = x, B_z = \partial_x A_y = 1';
  drawnow
  
  if 0
  hca = h(2);
  plot(hca,t(it),Ay(x(t(it),T),y(t(it),T)),'*k',t(1:it),Ay(x(t(1:it),T),y(t(1:it),T)),'-')
  hold(hca,'on')
  plot(hca,t(it),vy(t(it),T),'*k',t(1:it),vy(t(1:it),T),'-')
  hold(hca,'off')
  hca.XLim = t([1 end]);
  hca.YLim = [-1 1];
  irf_legend(hca,{'A_y','v_y'},[0.2 0.5],'fontsize',fontsize)
  hca.XLabel.String = 'time';
  hca.FontSize = fontsize;
  
  hca.Position(2) = h(1).Position(2);
  hca.Position(4) = h(1).Position(4);
  
  hb.Position(1) = 0.4;  
  h(1).Position(1) = 0.08;
  end
  if 1
    hca = h(2);
    plot(hca,x(t(it),T),vy(t(it),T),'*k',t(1:it),x(t(1:it),vy(t(1:it),T)),'-')
    hca.XLim = [-2 2];
    hca.YLim = [-2 2];
  end
  pause(0.1)
  currFrame = getframe(gcf);
  writeVideo(vidObj,currFrame);
end
  
% Close the file.
close(vidObj);