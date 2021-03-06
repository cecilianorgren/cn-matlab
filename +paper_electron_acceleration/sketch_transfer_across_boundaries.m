% sketch related to discussion of momentum/energy transfer across boundary
h = setup_subplots(1,2);
isub = 1;

if 0 % parallel phase space, z,vz 
  hca = h(isub); isub = isub + 1;
  
end
if 0 % real space x,y with particla trajectories phase mixed
  hca = h(isub); isub = isub + 1;
  ntheta = 200;
  theta = linspace(0,2*pi,ntheta);
  fun_x = @(r,x0,theta) x0+r*cos(theta);
  fun_y = @(r,y0,theta) y0+r*sin(theta);
  fun_z = @(r,z0,theta) z0 + theta*0;
  r1 = 1;
  r2 = 0.1;
  x1c = r1;
  x2c = -r2;
  y1c = 0;
  y2c = 0;
  z1c = 0;
  z2c = 1;
  x1 = r1*cos(theta);
  y1 = r1*sin(theta);
  x2 = r2*cos(theta);
  y2 = r2*sin(theta);
  
  xx1 = x1c+x1;
  xx2 = x2c+x2;
  yy1 = y1c+y1;
  yy2 = y2c+y2;
  zz1 = x1*0+z1c;
  zz2 = x2*0+z2c; 
  
  xx1 = x1c+x1;
  xx2 = x2c+x2;
  yy1 = y1c+y1;
  yy2 = y2c+y2;
  zz1 = x1*0+z1c;
  zz2 = x2*0+z2c;
  plot3(hca,xx1,yy1,zz1,xx2,yy2,zz2)
  view(hca,[0,1,1])
end

if 1 % parallel phase space, z,vz 
  hca = h(isub); isub = isub + 1;
  % spatial coordinates throughout parallel motion
  x1 = -0.5;
  x2 = 0.5;
  xturn = 0;  
  fun_x = @(r,x0,theta) x0+r*cos(theta);
  fun_z = @(r,x0,theta) x0+r*sin(theta);
  
  % initial (and final) parallel velocities
  z1c = 0;
  z2c = 0.5;
  z0 = 0.5*(z1c+z2c);
  
  % plot plasma sheet electron
  %theta = 
  plot(hca,[x1,xturn],[z1c z1c],'--','LineWidth',1.5,'Color',colors(1,:))
  hold(hca,'on')
  plot(hca,[x1,xturn],[z2c z2c],'-','LineWidth',1.5,'Color',colors(1,:))
  
  
  plot(hca,[x2,xturn],[z2c z2c],'--','LineWidth',1.5,'Color',colors(2,:))  
  plot(hca,[x2,xturn],[z1c z1c],'-','LineWidth',1.5,'Color',colors(2,:))
  
  theta1 = linspace(0,pi,ntheta)+pi/2;
  theta2 = linspace(pi,2*pi,ntheta)-pi/2;
  plot(hca,fun_x(z0,0,theta1),fun_z(z0,z0,theta1),':','LineWidth',1.5,'Color',colors(2,:))
  plot(hca,-fun_x(z0,0,theta2),fun_z(z0,z0,theta2),':','LineWidth',1.5,'Color',colors(1,:))
  
  hold(hca,'off')
end
if 1 % real space x,y with particla trajectories phase mixed
  hca = h(isub); isub = isub + 1;
  
  colors = pic_colors('matlab');
  
  ntheta = 100; 
  theta1 = linspace(0,pi,ntheta);
  theta2 = linspace(pi,2*pi,ntheta);
  fun_x = @(r,x0,theta) x0+r*cos(theta);
  fun_y = @(r,y0,theta) y0+r*sin(theta);
  fun_z = @(r,z0,theta) z0 + theta*0;
  r1 = 1;
  r2 = 0.1;
  % spatial coordinates throughout gyromotion
  x1c = r1;
  x2c = -r2*1.1;
  y1c = 0;
  y2c = 0;
  % initial (and final) parallel velocities
  z1c = 0;
  z2c = 0.5;
  
  hlines = plot3(hca,fun_x(r1,x1c,theta1),fun_y(r1,y1c,theta1),fun_z(r1,z1c,theta1),...
                     fun_x(r1,x1c,theta2),fun_y(r1,y1c,theta2),fun_z(r1,z2c,theta2));
  hlines(1).Color = colors(2,:);
  hlines(2).Color = colors(2,:);
  hlines(1).LineStyle = '-';
  hlines(2).LineStyle = '--';
  hlines(1).LineWidth = 1.5;
  hlines(2).LineWidth = 1.5;
  
  hold(hca,'on')
  angle = pi;
  plot3(hca,fun_x(r1,x1c,[angle angle]),fun_y(r1,y1c,[angle angle]),[fun_z(r1,z1c,angle),fun_z(r1,z2c,angle)],':','color',colors(2,:),'LineWidth',1.5)
  
  
  hlines = plot3(hca,fun_x(r2,x2c,theta2),fun_y(r2,y2c,theta2),fun_z(r2,z2c,theta2),...
                     fun_x(r2,x2c,theta1),fun_y(r2,y2c,theta1),fun_z(r2,z1c,theta1));
  hlines(1).Color = colors(1,:);
  hlines(2).Color = colors(1,:);
  hlines(1).LineStyle = '-';
  hlines(2).LineStyle = '--';
  hlines(1).LineWidth = 1.5;
  hlines(2).LineWidth = 1.5;
  
  angle = 0;
  plot3(hca,fun_x(r2,x2c,[angle angle]),fun_y(r2,y2c,[angle angle]),[fun_z(r2,z1c,angle),fun_z(r2,z2c,angle)],':','color',colors(1,:),'LineWidth',1.5)
  
  hold(hca,'off')
  view(hca,[.2,1,1])
end
h(1).Box = 'off';
h(2).Box = 'off';

axis(h(1),'off')
axis(h(2),'off')