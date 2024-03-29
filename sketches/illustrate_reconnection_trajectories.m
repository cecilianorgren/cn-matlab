% illustrate_reconnection_trajectories

%% Define magnetic field
a = 5; % a and b defines the aspect ratio of the magnetic field
b = 1;

% Used for plotting magnetic field lines with A 
xvec = a*linspace(-10.1,10,500);
zvec = b*linspace(-10,10,100);
yvec = linspace(-0,0,1);
[X,Y,Z] = ndgrid(xvec,yvec,zvec);
Ay = @(x,y,z) (x/a).^2 - (z/b).^2;
AY0 = Ay(X,Y,Z);

% Integrate orbits
m = 2;
q = 1;
Bx = @(x,y,z) 2*z/b^2;
By = @(x,y,z) x*0;
Bz = @(x,y,z) 2*x/a^2;
Ex = @(x,y,z) x*0;
Ey = @(x,y,z) x*0 + 0.1;
Ez = @(x,y,z) x*0;

options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,1,[-40 40]),...
                 'AbsTol',1e-16);

EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 200;
% Good for y, but do not cross at the same z.
x_init_all = [-1 5 4 0 0 0;
              -18.1 0 -3.5 1 0 0
              18.1 -3 3.9 -1.18 0 0];
x_init_all = [-17.81 0 -3.5 1 0 0;
              -1 7 4 0 0 0;              
              18.1 -1.8 3.9 -1.18 0 0];
if 0 % other particles
x_init_all = [-1 5 4 0 0 0;
              -17.5 0 3.5 1 0 0
              18.1 -3 3.9 -1.2 0 0];
            
x_init_all = [-2 0 0.1 -1 -2 0;
              -2 0 0.1 -1 2 0
              -2 0 0.1 1 -2 0];
end

x_init_all = [0.1 7 4 0 1 0;              
              2 7 4 0 1 0;              
              4 7 4 0 1 0;              
              ];

clear p
for ip = 1:size(x_init_all,1)
  x_init = x_init_all(ip,:);
  [t,x_sol] = ode45(EoM,[tstart tstop],x_init,options); % 
  p(ip).t = t;
  p(ip).x = x_sol(:,1);
  p(ip).y = x_sol(:,2);
  p(ip).z = x_sol(:,3);
  p(ip).vx = x_sol(:,4);
  p(ip).vy = x_sol(:,5);
  p(ip).vz = x_sol(:,6);
end

colors = pic_colors('matlab');
linewidth = 1.5;


if 1 % Figure
  %%  
  xlim = [-20 20];
  zlim = [-5 5];
  S = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[-105:10:110]); % B field lines
  SX = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[0 0]-0.01);  % X line  
  nrows = 3;
  ncols = 1;
  h = gobjects([nrows*ncols,1]);
  ip = 0;
  for irows = 1:nrows;
    for icols = 1:ncols
      ip = ip + 1;
      h(ip) = subplot(nrows,ncols,ip);
    end
  end

  isub = 1;
  
  if 1 % xyz
    hca = h(1); isub = isub + 1;  
    %plot3(hca,SX(1).X,SX(1).X*0,SX(1).Y,'--k',SX(2).X,SX(2).X*0,SX(2).Y,'--k');
    hold(hca,'on')

    %for il = 1:numel(S)
    %  plot3(hca,S(il).X,S(il).X*0,S(il).Y,'color',[0.5 0.5 0.5])
    %end
    for ip = 1:numel(p)
      if ip == 1; hold(hca,'on'); end
      plot3(hca,p(ip).x(1),p(ip).y(1),p(ip).z(1),'Marker','o','color',colors(ip,:))
      plot3(hca,p(ip).x(end),p(ip).y(end),p(ip).z(end),'Marker','x','color',colors(ip,:))    
      plot3(hca,p(ip).x,p(ip).y,p(ip).z,'color',colors(ip,:))
    end
    hold(hca,'off')
    %view(hca,[0 -1 0])
    %view(hca,[0 0 1])
    %axis(hca,'equal')
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'y';
    hca.ZLabel.String = 'z';
  end
  
  if 1 % xz plane
    hca = h(isub); isub = isub + 1;
    plot(hca,SX(1).X,SX(1).Y,'--k',SX(2).X,SX(2).Y,'--k');
    hold(hca,'on')

    for il = 1:numel(S)
      plot(hca,S(il).X,S(il).Y,'color',[0.5 0.5 0.5])
    end
    for ip = 1:numel(p)
      plot(hca,p(ip).x(1),p(ip).z(1),'Marker','o','color',colors(ip,:),'MarkerFaceColor',colors(ip,:))
      plot(hca,p(ip).x(end),p(ip).z(end),'Marker','x','color',colors(ip,:))
      plot(hca,p(ip).x,p(ip).z,'color',colors(ip,:),'linewidth',linewidth)
    end
    hold(hca,'off')  
    hca.XLim = xlim;
    hca.YLim = zlim;
    %axis(hca,'equal')
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'z';
  end
  
  if 1 % xy plane
    hca = h(isub); isub = isub + 1;  

    plot(hca,0,0)
    hold(hca,'on');
    for ip = 1:numel(p)    
      plot(hca,p(ip).x(1),p(ip).y(1),'Marker','o','color',colors(ip,:),'MarkerFaceColor',colors(ip,:))
      plot(hca,p(ip).x(end),p(ip).y(end),'Marker','x','color',colors(ip,:))
      plot(hca,p(ip).x,p(ip).y,'color',colors(ip,:),'linewidth',linewidth)
    end
    plot(hca,[0 0],hca.YLim,'--k')    
    
    hca.XLabel.String = 'x';
    hca.YLabel.String = 'y';
  
    hold(hca,'off')  
    
    hca.XLim = xlim;  
    hca.XGrid = 'on';
    hca.YGrid = 'on';

    irf_legend(hca,'X line',[0.5 1.01],'color',[0.5 0.5 0.5],'fontsize',14,'horizontalalignment','center');
  end

end
