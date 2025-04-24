%% xy
units = irf_units;

% Define magnetic field
b = 20e3; % m

zvec = b*linspace(-4,4,100);

%[X,Y,Z] = ndgrid(xvec,yvec,zvec);


%AY0 = Ay(X,Y,Z);


B0 = 10e-9;
E0 = 20e-3;
% Integrate orbits
m = units.me;
q = -units.e;

Bx = @(x,y,z) B0*tanh(z/b);
By = @(x,y,z) x*0;
Bz = @(x,y,z) x*0;
Ex = @(x,y,z) x*0;
Ey = @(x,y,z) z*0 + 1*E0;
Ez = @(x,y,z) -6*E0*(z/b).*exp(-z.^2/b.^2);
Ay = @(x,y,z,t) B0*b*log(sech(z/b)) - t.*Ey(z,y,z);

options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,xvec([1 end])),...
                 'AbsTol',1e-6);
options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,1,[-40 40]),...
                 'AbsTol',1e-12);
options = odeset('AbsTol',1e-12);
%options = odeset();
EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 0.2;
% Good for y, but do not cross at the same z.


x_init_all = [0 0 30e3 0 0 1000e3];
x_init_all = [0 0 100e3 0 0 -100e3];

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
  p(ip).Ay = Ay(x_sol(:,1),x_sol(:,2),x_sol(:,3),t);
  p(ip).py = m*p(ip).vy + q*p(ip).Ay;
  p(ip).Ey = Ey(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  p(ip).Bx = Bx(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  p(ip).Ez = Ez(x_sol(:,1),x_sol(:,2),x_sol(:,3));
end

colors = pic_colors('matlab');
%colors = [colors(2)];
colors(3,:) = colors(1,:);
%colors(3,:) = colors(5,:).^0.5;
colors(1,:) = colors(1,:).^0.2;
linewidth = 1.5;

fontsize = 16;

if 1 % Figure 1
  nRows = 5; nCols = 2;
  h = setup_subplots(nRows,nCols,'vertical');
  
  %ip = 0; h = gobjects(0);
  %for iRow = 1:nRows
  %  for iCol = 1:nCols
  %    ip = ip + 1;
  %    h(ip) = subplot(nRows,nCols,ip);
  %  end
  %end

  isub = 1;
  
  hca = h(isub); isub = isub + 1;
  plot(hca, ...
    zvec*1e-3,Bx(0,0,zvec)*1e9,...
    zvec*1e-3,Ey(0,0,zvec)*1e3, ...
    zvec*1e-3,Ez(0,0,zvec)*1e3)
  hca.XLabel.String = 'z (km)';
  irf_legend(hca,{'B_x','E_y','E_z'},[0.02 0.98])

  hca = h(isub); isub = isub + 1;
  plot(hca,zvec*1e-3,Ay(0,0,zvec,0))
  hca.XLabel.String = 'z (km)';

  ip = 1;

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).y*1e-3,p(ip).z*1e-3)
  hca.XLabel.String = 'y (km)';
  hca.YLabel.String = 'z (km)';

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).y*1e-3 ,p(ip).t,p(ip).z*1e-3)
  hca.YLabel.String = 'r (km)';
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'y','z'},[0.02 0.1], 'fontsize',fontsize)


  %hca = h(isub); isub = isub + 1;
  %plot(hca,p(ip).y,p(ip).z)

  %hca = h(isub); isub = isub + 1;
  %plot(hca,p(ip).t,p(ip).y)

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vy,p(ip).t,p(ip).Ay*q/m,p(ip).t,p(ip).py)
  hca.XLabel.String = 't (s)';

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vz*1e-3,p(ip).t,-p(ip).Ey.*p(ip).Bx./p(ip).Bx.^2*1e-3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'v_z','E_yB_x/B^2'},[0.02 0.1], 'fontsize',fontsize)

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vy*1e-3,p(ip).t,p(ip).Ez.*p(ip).Bx./p(ip).Bx.^2*1e-3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'v_y','E_zB_x/B^2'},[0.02 0.1], 'fontsize',fontsize)


  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,(p(ip).Ey)*1e3, ...
    p(ip).t,(p(ip).vz.*p(ip).Bx)*1e3, ...
    p(ip).t,(p(ip).Ey + p(ip).vz.*p(ip).Bx)*1e3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'E_y','v_zB_x','E_y+v_zB_x'},[0.02 0.1], 'fontsize',fontsize)

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,(p(ip).Ez)*1e3, ...
    p(ip).t,( - p(ip).vy.*p(ip).Bx)*1e3,...
    p(ip).t,(p(ip).Ez - p(ip).vy.*p(ip).Bx)*1e3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'E_z','-v_yB_x','E_z-v_yB_x'},[0.02 0.1], 'fontsize',fontsize)


  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,(p(ip).Ez)*1e3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'E_z'},[0.02 0.1], 'fontsize',fontsize)
  
  hl = findobj(gcf,'type','line');
  c_eval('hl(?).LineWidth = 2;',1:numel(hl))

  c_eval('h(?).FontSize = 18;',1:numel(h))
  c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

  hlink = linkprop(h(4:end),{'XLim'});
end



if 0 % Figure 1
  S = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[-105:10:110]);
  SX = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[0 0]); 
  %[C,H] = contour3()
  h1(1) = subplot(3,2,1);
  h1(2) = subplot(3,2,3);
  h1(3) = subplot(3,2,5);
  h2 = subplot(3,2,[2 4 6]);


  for ip = 1:numel(p)
    hca = h1(1);
    plot(hca,p(ip).t,p(ip).x)
    hold(hca,'on')

    hca = h1(2);
    plot(hca,p(ip).t,p(ip).y)
    hold(hca,'on')

    hca = h1(3);
    plot(hca,p(ip).t,p(ip).z)
    hold(hca,'on')
  end

  hold(h1(1),'off')
  hold(h1(2),'off')
  hold(h1(3),'off')

  hca = h2(1);
  %plot3(hca,1,1,1);
  plot3(hca,SX(1).X,SX(1).X*0,SX(1).Y,'--k',SX(2).X,SX(2).X*0,SX(2).Y,'--k');
  hold(hca,'on')

  for il = 1:numel(S)
    plot3(hca,S(il).X,S(il).X*0,S(il).Y,'color',[0.5 0.5 0.5])
  end
  for ip = 1:numel(p)
    plot3(hca,p(ip).x,p(ip).y,p(ip).z)
  end
  hold(hca,'off')
  view(hca,[0 -1 0])
  %view(hca,[0 0 1])
  %axis(hca,'equal')
end

if 0 % Figure 2
  xmark = -4.58;
  xlim = [-20 20];
  zlim = [-5 5];
  S = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[-105:10:110]);
  SX = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[0 0]-0.01); 
  %[C,H] = contour3()
  nrows = 2;
  ncols = 1;
  clear h
  h(1) = subplot(nrows,ncols,1);
  h(2) = subplot(nrows,ncols,2);
  %h(3) = subplot(nrows,ncols,3);
  

  isub = 1;
  
  if 0
    hca = h(1); isub = isub + 1;  
    plot3(hca,SX(1).X,SX(1).X*0,SX(1).Y,'--k',SX(2).X,SX(2).X*0,SX(2).Y,'--k');
    hold(hca,'on')

    for il = 1:numel(S)
      plot3(hca,S(il).X,S(il).X*0,S(il).Y,'color',[0.5 0.5 0.5])
    end
    for ip = 1:numel(p)

      plot3(hca,p(ip).x(1),p(ip).y(1),p(ip).z(1),'Marker','o','color',colors(ip,:))
      plot3(hca,p(ip).x(end),p(ip).y(end),p(ip).z(end),'Marker','x','color',colors(ip,:))    
      plot3(hca,p(ip).x,p(ip).y,p(ip).z,'color',colors(ip,:))
    end
    hold(hca,'off')
    %view(hca,[0 -1 0])
    %view(hca,[0 0 1])
    %axis(hca,'equal')
  end
  
  hca = h(isub); isub = isub + 1;
  plot(hca,SX(1).X,SX(1).Y,'--k',SX(2).X,SX(2).Y,'--k');
  hold(hca,'on')

  for il = 1:numel(S)
    plot(hca,S(il).X,S(il).Y,'color',[0.5 0.5 0.5])
  end
  for ip = 1:numel(p)
    plot(hca,p(ip).x(1),p(ip).z(1),'Marker','o','color',colors(ip,:))
    plot(hca,p(ip).x(end),p(ip).z(end),'Marker','x','color',colors(ip,:))
    plot(hca,p(ip).x,p(ip).z,'color',colors(ip,:),'linewidth',linewidth)
  end
  hold(hca,'off')  
  hca.XLim = xlim;
  hca.YLim = zlim;
  %axis(hca,'equal')
  
  hca = h(isub); isub = isub + 1;  
  holdon = 0;
  plot(hca,0,0)
  hold(hca,'on');
  for ip = 1:numel(p)    
    plot(hca,p(ip).x(1),p(ip).y(1),'Marker','o','color',colors(ip,:))
    plot(hca,p(ip).x(end),p(ip).y(end),'Marker','x','color',colors(ip,:))
    plot(hca,p(ip).x,p(ip).y,'color',colors(ip,:),'linewidth',linewidth)
  end
  plot(hca,[0 0],hca.YLim,'--k')
  for ip = 1:numel(p)
    xq = xmark;
    yq = -0.7854;
    [ix] = find(p(ip).x>-10);
    [~,iy] = min(abs(p(ip).y(ix)-yq));
    %ind = intersect(ix,iy);
    ind = ix(iy);
    vs = 10;
    quiver(hca,p(ip).x(ind),p(ip).y(ind),vs*p(ip).vx(ind),vs*p(ip).vy(ind),0,'color',colors(ip,:),'linewidth',1.5*linewidth)
  end
  
  hold(hca,'off')  
  %axis(hca,'equal')
  hca.XLim = xlim;  
  
  %hl = findobj(gcf,'type','line');
  for ip = 1:numel(h)
    hold(h(ip),'on')
    plot(h(ip),xmark*[1 1],hca.YLim,':k')
    hold(h(ip),'off')
  end
  %c_eval('hl(?).LineWidth = 1;',1:numel(hl))
  
  drawnow
  grid on
  if 0  
  axis(h(1),'equal')
  axis(h(2),'equal')
  drawnow
  h(1).YLim = zlim;
  end
  %ylim = h(2).YLim;
  %axis(h(2),'equal')
  h(1).XLim = xlim;
  h(2).XLim = xlim;
  h(2).YLim = ylim;
  ylim = [-23 8];
  h(2).YLim = ylim;
  %h(1).XLim = xlim;
  %h(2).Position(4) = 0.7;
  %h(2).Position(4) = h(2).Position(3)*diff(ylim)/diff(xlim);
  axis(h(2),'equal')
  %h(2).Position(2) = 0.15;
  %compact_panels(0.01)
  h(1).Position(2) = h(2).Position(2)+h(2).Position(4)+0.01;
  h(1).Position(4) = 0.2;
  
  c_eval('h(?).XTickLabels = []; h(?).YTickLabels = [];',1:numel(h))
  c_eval('h(?).XTick = -40:5:40;',1:numel(h))
  h(2).XLabel.String = 'x';
  h(2).YLabel.String = 'y';
  h(1).YLabel.String = 'z';
  c_eval('h(?).FontSize = 14;',1:numel(h))
  c_eval('h(?).LineWidth = 1;',1:numel(h))
end

if 0 % Figure 3
  %%
  xmark = -4.58;
  xlim = [-20 20];
  zlim = [-5 5];
  S = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[-105:10:110]);
  SX = contourcs(xvec,zvec,squeeze(Ay(X,Y,Z))',[0 0]-0.01); 
  %[C,H] = contour3()
  nrows = 3;
  ncols = 1;
  clear h
  h(1) = subplot(nrows,ncols,1);
  h(2) = subplot(nrows,ncols,2);
  h(3) = subplot(nrows,ncols,3);
  

  isub = 1;
  
  if 1
    hca = h(1); isub = isub + 1;  
    plot3(hca,SX(1).X,SX(1).X*0,SX(1).Y,'--k',SX(2).X,SX(2).X*0,SX(2).Y,'--k');
    hold(hca,'on')

    for il = 1:numel(S)
      plot3(hca,S(il).X,S(il).X*0,S(il).Y,'color',[0.5 0.5 0.5])
    end
    for ip = 1:numel(p)

      plot3(hca,p(ip).x(1),p(ip).y(1),p(ip).z(1),'Marker','o','color',colors(ip,:))
      plot3(hca,p(ip).x(end),p(ip).y(end),p(ip).z(end),'Marker','x','color',colors(ip,:))    
      plot3(hca,p(ip).x,p(ip).y,p(ip).z,'color',colors(ip,:))
    end
    hold(hca,'off')
    %view(hca,[0 -1 0])
    %view(hca,[0 0 1])
    %axis(hca,'equal')
  end
  
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
  
  hca = h(isub); isub = isub + 1;  
  
  plot(hca,0,0)
  hold(hca,'on');
  for ip = 1:numel(p)    
    plot(hca,p(ip).x(1),p(ip).y(1),'Marker','o','color',colors(ip,:),'MarkerFaceColor',colors(ip,:))
    plot(hca,p(ip).x(end),p(ip).y(end),'Marker','x','color',colors(ip,:))
    plot(hca,p(ip).x,p(ip).y,'color',colors(ip,:),'linewidth',linewidth)
  end
  plot(hca,[0 0],hca.YLim,'--k')
  for ip = 1:numel(p)
    xq = xmark;
    yq = -0.7854;
    [ix] = find(p(ip).x>-10);
    [~,iy] = min(abs(p(ip).y(ix)-yq));
    %ind = intersect(ix,iy);
    ind = ix(iy);
    vs = 10;
    quiver(hca,p(ip).x(ind),p(ip).y(ind),vs*p(ip).vx(ind),vs*p(ip).vy(ind),0,'color',colors(ip,:),'linewidth',1.5*linewidth)
  end
  
  hold(hca,'off')  
  %axis(hca,'equal')
  hca.XLim = xlim;  
  
  irf_legend(hca,'X line',[0.65 0.9],'color',[0.5 0.5 0.5],'fontsize',14);
  
 
  
  
  hca = h(isub); isub = isub + 1;  
 
  
  plot(hca,[0 0],0.90*[-2 2],'--k')
 
  hold(hca,'on')
  
  for ip = 1:numel(p)
    %xq = xmark;
    yq = -0.7854;
    [ix] = find(p(ip).x>-10);
    [~,iy] = min(abs(p(ip).y(ix)-yq));
    %ind = intersect(ix,iy);
    ind = ix(iy);
    ind = ind + 0;
    %ind = length(p(ip).x);
    
    tBx = Bx(p(ip).x,p(ip).y,p(ip).z);
    vxBz = p(ip).vx.*Bz(p(ip).x,p(ip).y,p(ip).z);
    vzBx = -p(ip).vz.*Bx(p(ip).x,p(ip).y,p(ip).z);
    
    vyBx = -p(ip).vy.*Bx(p(ip).x,p(ip).y,p(ip).z);
    tEy = Ey(p(ip).x,p(ip).y,p(ip).z);
    
    toplot = 1*vxBz+1*vzBx;
    hh=plot(hca,p(ip).x(1),toplot(1),'Marker','o','color',colors(ip,:),'linewidth',linewidth,'MarkerFaceColor',colors(ip,:));  
    plot(hca,p(ip).x(1:ind),toplot(1:ind),'color',colors(ip,:),'linewidth',linewidth,'linestyle','-');
        
    toplot = 1*vxBz+0*vzBx;
    %plot(hca,p(ip).x(1),toplot(1),'Marker','o','color',colors(ip,:),'linewidth',linewidth);   
    plot(hca,p(ip).x(1:ind),toplot(1:ind),'color',colors(ip,:),'linewidth',linewidth,'linestyle','-');
    
    
    plot(hca,xlim,-tEy(1)*[1 1],'color',[0.5 0.5 0.5],'linewidth',linewidth);
    %plot(hca,p(ip).t,vyBx,'color',colors(ip,:).^0.8,'linewidth',linewidth);   
    grid(hca,'on')
    
  end
  hca.XLabel.String = 'x';
  hca.YLabel.String = '-(v_zB_x-v_xB_z)';
  set(hca,'ColorOrder',[0.5 0.5 0.5])
  irf_legend(hca,'-E_y',[0.29 0.72],'color',[0.5 0.5 0.5],'fontsize',14);
  hold(hca,'off')
  
  %hl = findobj(gcf,'type','line');
  for ip = 1:numel(h)
    hold(h(ip),'on')
    plot(h(ip),xmark*[1 1],h(ip).YLim,':k')
    hold(h(ip),'off')
  end
  %c_eval('hl(?).LineWidth = 1;',1:numel(hl))
  
  
  %c_eval('h(?).XTickLabels = []; h(?).YTickLabels = [];',1:numel(h))
  c_eval('h(?).XTick = -40:5:40;',1:numel(h))
  h(2).XLabel.String = 'x';
  h(2).YLabel.String = 'y';
  h(1).YLabel.String = 'z';
  c_eval('h(?).FontSize = 14;',1:numel(h))
  c_eval('h(?).LineWidth = 1;',1:numel(h))
  
  h(2).XGrid = 'on';
  h(2).YGrid = 'on';
  ylim = [-22.5 8];
  %ylim = [-17 8];
  h(2).YLim = ylim;
  h(2).YTick = h(2).XTick;
  
  h(1).Position(2) = 0.72;
  h(2).Position(4) = 0.3;
  h(2).Position(1) = h(1).Position(1);
  h(2).Position(3) = h(1).Position(3);
  
  h(3).Position(2) = 0.2;
  h(3).Position(4) = 0.2;
  h(3).Position(1) = h(1).Position(1);
  h(3).Position(3) = h(1).Position(3);
  
  
  drawnow
  axis(h(2),'equal')
  h(2).YLim = ylim;
  
  c_eval('h(?).XTick = []; h(?).YTick = [];',1:2)
  c_eval('h(?).XTick = []; h(?).YTick = 0;',3)
h(2).Children = h(2).Children(end:-1:1);
  %%
  drawnow
  %ylim = h(2).YLim;
  axis(h(2),'equal')
  drawnow
  h(1).XLim = xlim;
  h(2).XLim = xlim;
  h(2).YLim = ylim;
  drawnow
  %h(1).XLim = xlim;
  %h(2).Position(4) = 0.7;
  %h(2).Position(4) = h(2).Position(3)*diff(ylim)/diff(xlim);
  axis(h(2),'equal')
  h(2).XLim = xlim;
  h(2).Position(3) = h(1).Position(3);
  %h(2).Position(2) = 0.15;
  %compact_panels(0.01)
  %h(1).Position(2) = h(2).Position(2)+h(2).Position(4)+0.01;
  %h(1).Position(4) = 0.2;
  %h(3).Position(2) = 0.2;
  
end



%% xz, Hall fields
units = irf_units;

% Define magnetic field
b = 200e3; % m

zvec = b*linspace(-4,4,100);

%[X,Y,Z] = ndgrid(xvec,yvec,zvec);


%AY0 = Ay(X,Y,Z);


B0 = 10e-9;
E0 = 11e-3;
% Integrate orbits
m = units.mp;
q = units.e;

Bx = @(x,y,z) 0;
By = @(x,y,z) -B0*(z/b).*exp(-z.^2/b.^2);
Bz = @(x,y,z) z*0;
Ex = @(x,y,z) z*0;
Ey = @(x,y,z) z*0;
Ez = @(x,y,z) -E0*(z/b).*exp(-z.^2/b.^2);
Ax = @(x,y,z,t) B0*(b/2).*exp(-z.^2/b.^2);
phiz = @(x,y,z) E0*(b/2).*exp(-z.^2/b.^2);

options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,xvec([1 end])),...
                 'AbsTol',1e-6);
options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,1,[-40 40]),...
                 'AbsTol',1e-12);
options = odeset('AbsTol',1e-12);
%options = odeset();
EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 20;
% Good for y, but do not cross at the same z.


x_init_all = [0 0 30e3 0 0 1000e3];
x_init_all = [0 0 500e3 0 0 00e3];

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
  p(ip).Ax = Ax(x_sol(:,1),x_sol(:,2),x_sol(:,3),t);
  p(ip).px = m*p(ip).vx + q*p(ip).Ax;
  p(ip).Ey = Ey(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  p(ip).Bx = Bx(x_sol(:,1),x_sol(:,2),x_sol(:,3));
  p(ip).Ez = Ez(x_sol(:,1),x_sol(:,2),x_sol(:,3));
end

colors = pic_colors('matlab');
%colors = [colors(2)];
colors(3,:) = colors(1,:);
%colors(3,:) = colors(5,:).^0.5;
colors(1,:) = colors(1,:).^0.2;
linewidth = 1.5;

fontsize = 16;

if 1 % Figure 1
  nRows = 5; nCols = 2;
  h = setup_subplots(nRows,nCols,'vertical');
  
  %ip = 0; h = gobjects(0);
  %for iRow = 1:nRows
  %  for iCol = 1:nCols
  %    ip = ip + 1;
  %    h(ip) = subplot(nRows,nCols,ip);
  %  end
  %end

  isub = 1;
  
  hca = h(isub); isub = isub + 1;
  plot(hca, ...
    zvec*1e-3,By(0,0,zvec)*1e9,...
    zvec*1e-3,Ez(0,0,zvec)*1e3,'--')
  hca.XLabel.String = 'z (km)';
  irf_legend(hca,{'B_y','E_z'},[0.02 0.98], 'fontsize',fontsize)

  if 0
    hca = h(isub); isub = isub + 1;
    plot(hca, zvec*1e-3, -Ez(0,0,zvec)./By(0,0,zvec)*1e-3)
    hca.XLabel.String = 'z (km)';
    hca.XLabel.String = 'v (km/z)';
    irf_legend(hca,{'-E_z/B_y'},[0.02 0.98], 'fontsize',fontsize)
  end

  hca = h(isub); isub = isub + 1;
  plot(hca,zvec*1e-3,Ax(0,0,zvec,0))
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'A_x (...)';

  hca = h(isub); isub = isub + 1;

  phi_from_ez = -gradient(phiz(0,0,zvec),zvec);
  phi_from_ez = cumtrapz(zvec,Ez(0,0,zvec));
  plot(hca,zvec*1e-3,phiz(0,0,zvec), zvec*1e-3,phi_from_ez,'--')
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = '\phi_z (V)';


  ip = 1;

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).x*1e-3,p(ip).z*1e-3)
  hca.XLabel.String = 'x (km)';
  hca.YLabel.String = 'z (km)';
  %axis(hca,'equal')

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).x*1e-3 ,p(ip).t,p(ip).z*1e-3)
  hca.YLabel.String = 'r (km)';
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'x','z'},[0.02 0.1], 'fontsize',fontsize)

  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vx*1e-3 ,p(ip).t,p(ip).vz*1e-3)
  hca.YLabel.String = 'v (km/s)';
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'v_x','v_z'},[0.02 0.1], 'fontsize',fontsize)

  hca = h(isub); isub = isub + 1;
  vzvx = p(ip).vz./p(ip).vx;
  vzvx(vzvx>100) = NaN;
  plot(hca,p(ip).t,vzvx)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'v_z/v_x';

  hca = h(isub); isub = isub + 1;
  vabs = sqrt(p(ip).vx.^2 + p(ip).vz.^2);
  angle = acosd(p(ip).vx./vabs);
  angle2 = atand(p(ip).vx./p(ip).vz);
  plot(hca,p(ip).t,angle2)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'atan(v_x/v_z)';

  hca = h(isub); isub = isub + 1;
  mv22 = m*(p(ip).vx.^2 + p(ip).vz.^2)/2;  
  plot(hca,p(ip).t,mv22/units.eV)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'mv^2/2 (eV)';

  %hca = h(isub); isub = isub + 1;
  %plot(hca,p(ip).y,p(ip).z)

  %hca = h(isub); isub = isub + 1;
  %plot(hca,p(ip).t,p(ip).y)

  hca = h(isub); isub = isub + 1;
  mv22 = m*(p(ip).vx.^2 + p(ip).vz.^2)/2;  
  vabs = sqrt(p(ip).vx.^2 + p(ip).vz.^2);
  angle2 = atand(p(ip).vx./p(ip).vz);
  
  plot(hca,mv22/units.eV,angle2)
  hca.XLabel.String = 'U (eV)';
  hca.YLabel.String = 'angle';


  if 0 % px vx aX
  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vx*1e-3, p(ip).t,p(ip).Ax*q/m*1e-3, p(ip).t,p(ip).px/m*1e-3)
  hca.XLabel.String = 't (s)';
  hca.YLabel.String = 'v_x (km/s)';
  irf_legend(hca,{'v_x','(q/m)A_x','p_x/m'},[0.02 0.1], 'fontsize',fontsize)
  end

  if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vz*1e-3,p(ip).t,-p(ip).Ey.*p(ip).Bx./p(ip).Bx.^2*1e-3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'v_z','E_yB_x/B^2'},[0.02 0.1], 'fontsize',fontsize)
  end
  if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,p(ip).vy*1e-3,p(ip).t,p(ip).Ez.*p(ip).Bx./p(ip).Bx.^2*1e-3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'v_y','E_zB_x/B^2'},[0.02 0.1], 'fontsize',fontsize)
  end
  if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,(p(ip).Ey)*1e3, ...
    p(ip).t,(p(ip).vz.*p(ip).Bx)*1e3, ...
    p(ip).t,(p(ip).Ey + p(ip).vz.*p(ip).Bx)*1e3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'E_y','v_zB_x','E_y+v_zB_x'},[0.02 0.1], 'fontsize',fontsize)
  end
  if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,(p(ip).Ez)*1e3, ...
    p(ip).t,( - p(ip).vy.*p(ip).Bx)*1e3,...
    p(ip).t,(p(ip).Ez - p(ip).vy.*p(ip).Bx)*1e3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'E_z','-v_yB_x','E_z-v_yB_x'},[0.02 0.1], 'fontsize',fontsize)
  end
  if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,p(ip).t,(p(ip).Ez)*1e3)
  hca.XLabel.String = 't (s)';
  irf_legend(hca,{'E_z'},[0.02 0.1], 'fontsize',fontsize)
  end

  hl = findobj(gcf,'type','line');
  c_eval('hl(?).LineWidth = 2;',1:numel(hl))

  c_eval('h(?).FontSize = 18;',1:numel(h))
  c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

  %hlink = linkprop(h(4:end),{'XLim'});
end


%% Analytical opening angle with symbolic expressions, higher level, e.g. Ax and phi not B0 E0 b
syms phi Ax b B0 E0 z q m kB eb U Uz Ux



% Ax = B0*(b/2)*exp(-z^2/b^2);
% phi = E0*(b/2)*exp(-z^2/b^2); % V

% Ax = B0*F0;
% phi = E0*F0;
% phi = (B0/E0)*Ax

phi = U/q;
Ax = phi/eb;
%phi = Ax/eb;
vx = -(q/m)*Ax;

Ux = m*vx^2/2;

%U = Ux + Uz;

Uz = U - Ux;

vz = (2*Uz/m)^0.5;

vxvz = vx/vz;

v = sqrt(2*U/m);

f_phi = matlabFunction(phi);
f_Ax = matlabFunction(Ax);
f_U = matlabFunction(U);
f_Ux = matlabFunction(Ux);
f_Uz = matlabFunction(Uz);
f_vx = matlabFunction(vx);
f_vz = matlabFunction(vz);
f_vxvz = matlabFunction(vxvz);
f_v = matlabFunction(v);



units = irf_units;
m_ = units.mp;
q_ = units.e;
U_ = linspace(0,7000,1050)*units.eV; % eV -> J
U_ = logspace(-1,4,1050)*units.eV; % eV -> J
eb_ =  linspace(0e3,1000e3,249); % E0 = Ez = 10 mV/m, B0 = By = 10 nT -> E0/B0 = 1000 km/s
%eb_ =  logspace(5,8,19);

[UU,EB] = ndgrid(U_,eb_);

VX = f_vx(UU,EB,m_);
VZ = real(f_vz(UU,EB,m_));
VXVZ = real(f_vxvz(UU,EB,m_));
%VXVZ(not(isreal(f_vz(UU,EB,m_)))) = NaN;
%UU = f_U(U);
UX = f_Ux(UU,eb_,m_);
UZ = f_Uz(UU,eb_,m_);

VX(UX > UU) = NaN;
UX(UX > UU) = NaN;
UZ(UZ < 0) = NaN;
VZ(VZ == 0) = NaN;
VXVZ(VXVZ == 0) = NaN;

x_label = 'U (eV)';
y_label = '\phi/A_x (km/s)'; % 'E_0/B_0 (km/s)'

if 1
  %v_ = linspace(0,2000,250); % km/s
  %[VV,EB] = ndgrid(v_,eb_);
  UU = sqrt(2*UU/m_)*1e-3*units.eV;
  %VV*units.eV;
  %v = sqrt(2*U/m);
  x_label = 'v (km/s)';
end
fontsize = 16;

nRows = 1; nCols = 4;
h = setup_subplots(nRows,nCols,'horizontal');

isub = 1;

if 0 % U_x
hca = h(isub); isub = isub + 1;
pcolor(hca, UU/units.eV, EB*1e-3, (UX/units.eV))
shading(hca,'flat')
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;
hcb = colorbar(hca);
hcb.YLabel.String = 'U_x';
end
if 0 % Uz
hca = h(isub); isub = isub + 1;
pcolor(hca, UU/units.eV, EB*1e-3, (UZ/units.eV))
shading(hca,'flat')
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;
hcb = colorbar(hca);
hcb.YLabel.String = 'U_z';
end

hca = h(isub); isub = isub + 1;
pcolor(hca, UU/units.eV, EB*1e-3, VX*1e-3)
shading(hca,'flat')
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;
hcb = colorbar(hca);
hcb.YLabel.String = 'v_x (km/s)';

hca = h(isub); isub = isub + 1;
pcolor(hca, UU/units.eV, EB*1e-3, VZ*1e-3)
shading(hca,'flat')
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;
hcb = colorbar(hca);
hcb.YLabel.String = 'v_z (km/s)';

hca = h(isub); isub = isub + 1;
pcolor(hca, UU/units.eV, EB*1e-3, abs(VXVZ))
shading(hca,'flat')
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;
hcb = colorbar(hca);
hcb.YLabel.String = '|v_x/v_z|';
hca.CLim = [0 prctile(hca.Children.CData(:),99)];

hca = h(isub); isub = isub + 1;
pcolor(hca, UU/units.eV, EB*1e-3, atand(VXVZ))
shading(hca,'flat')
hca.XLabel.String = x_label;
hca.YLabel.String = y_label;
hcb = colorbar(hca);
hcb.YLabel.String = '\theta = tan^{-1}(v_x/v_z)';
%hca.CLim = [0 prctile(hca.Children.CData(:),99)];

%colormap(pic_colors('candy4'))
colormap(irf_colormap('waterfall'))
compact_panels(0.01,0.01)
c_eval('h(?).Color = [0.9 0.9 0.9];',1:numel(h))
c_eval('h(?).Box = ''on''; h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = 16;',1:numel(h))
%c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

%c_eval('h(?).XScale = ''log''; h(?).YScale = ''log'';',1:numel(h))


hb = findobj(gcf,'type','colorbar');
%c_eval('h(?).Position(3) = h(?).Position(3)*1.1;',1:numel(hb))
c_eval('h(?).Position(2) = 0.15;',1:numel(hb))
%c_eval('hb(?).Position(1) = hb(?).Position(1)-0.0;',1:numel(hb))
c_eval('hb(?).Location = ''south''; hb(?).Position(3) = hb(?).Position(3)/2; hb(?).Position(1) = hb(?).Position(1) + hb(?).Position(3);',1:numel(hb))

delete(h(isub:end))
h(isub:end) = [];

for ip = 1:isub-1
  XData = h(ip).Children.XData;
  YData = h(ip).Children.YData;
  CData = h(ip).Children.CData;
  CData(CData>prctile(CData(:),97)) = NaN; 
  clim = hca.CLim;
  hold(h(ip),'on')
  [cc,hh] = contour(h(ip),XData,YData,CData,'k');
  clabel(cc,hh)
  hold(h(ip),'off')
  hca.CLim = clim;
end
%%
hca = h(isub); isub = isub + 1;
pcolor(hca, AX, EB, VZ)
shading(hca,'flat')
hca.XLabel.String = 'A_x';
hca.YLabel.String = 'E_0/B_0';

hca = h(isub); isub = isub + 1;
pcolor(hca, AX, EB, VXVZ)
shading(hca,'flat')
hca.XLabel.String = 'A_x';
hca.YLabel.String = 'E_0/B_0';

%% old plots
hca = h(isub); isub = isub + 1;
plot(hca, ...
  zvec*1e-3, f_U(E0_,b_,q_,zvec)/units.eV, ...
  zvec*1e-3, f_Ux(B0_,b_,m_, q_,zvec)/units.eV, ...
  zvec*1e-3, f_Uz(B0_,E0_,b_,m_,q_,zvec)/units.eV)
hca.XLabel.String = 'z (km)';
hca.YLabel.String = 'U (eV)';
irf_legend(hca,{'U = q\phi','U_x from v_x(A_x)','U_z = U - U_x'}',[0.02 0.98], 'fontsize',fontsize)
%irf_legend(hca,{'U','U_x','U_z'},[0.02 0.1], 'fontsize',fontsize)

if 1 % |vx|, |vz|
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, abs(f_vx(B0_,b_,m_,q_,zvec))*1e-3, zvec*1e-3, abs(f_vz(B0_,E0_,b_,m_,q_,zvec))*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = '|v| (km/s)';
  irf_legend(hca,{'v_x','v_z'},[0.02 0.1], 'fontsize',fontsize)
end
if 0 % vx
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, f_vx(B0_,b_,m_,q_,zvec)*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'v_x (km/s)';
end
if 0 % vz
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, f_vz(B0_,E0_,b_,m_,q_,zvec)*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'v_z (km/s)';
end
if 0 % vx/vz
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, f_vz(B0_,E0_,b_,m_,q_,zvec)./f_vx(B0_,b_,m_,q_,zvec))
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'v_z/v_x';
end

hca = h(isub); isub = isub + 1;
plot(hca, zvec*1e-3, abs(f_vx(B0_,b_,m_,q_,zvec)./f_vz(B0_,E0_,b_,m_,q_,zvec)))
hca.XLabel.String = 'z (km)';
hca.YLabel.String = '|v_x/v_z|';

hca = h(isub); isub = isub + 1;
vx_ = f_vx(B0_,b_,m_,q_,zvec);
vz_ = f_vz(B0_,E0_,b_,m_,q_,zvec);
x_ = cumsum(vx_);
z_ = cumsum(vz_);

plot(hca, x_, z_);
hca.XLabel.String = 'x (...)';
hca.YLabel.String = 'z (...)';

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))

c_eval('h(?).FontSize = 18;',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

%% Analytical opening angle with symbolic expressions
syms phi Ax b B0 E0 z q m kB

Ax = B0*(b/2)*exp(-z^2/b^2);
phi = E0*(b/2)*exp(-z^2/b^2); % V

vx = -(q/m)*Ax;

U = q*phi;

Ux = m*vx^2/2;

Uz = U - Ux;
vz = (2*Uz/m)^0.5;


vzvx = vz/vx;

f_phi = matlabFunction(phi);
f_U = matlabFunction(U);
f_Ux = matlabFunction(Ux);
f_Uz = matlabFunction(Uz);
f_vx = matlabFunction(vx);
f_vz = matlabFunction(vz);




units = irf_units;
kB_ = units.kB;
b_ = 400e3;
B0_ = 32.e-9;
E0_ = 10e-3;
m_ = units.mp;
q_ = units.e;

zvec = 3*b_*linspace(-1,1,100);

fontsize = 16;

nRows = 5; nCols = 1;
h = setup_subplots(nRows,nCols,'vertical');

isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca, zvec*1e-3, f_phi(E0_,b_,zvec))
hca.XLabel.String = 'z (km)';
hca.YLabel.String = '\phi (V)';

hca = h(isub); isub = isub + 1;
plot(hca, ...
  zvec*1e-3, f_U(E0_,b_,q_,zvec)/units.eV, ...
  zvec*1e-3, f_Ux(B0_,b_,m_, q_,zvec)/units.eV, ...
  zvec*1e-3, f_Uz(B0_,E0_,b_,m_,q_,zvec)/units.eV)
hca.XLabel.String = 'z (km)';
hca.YLabel.String = 'U (eV)';
irf_legend(hca,{'U = q\phi','U_x from v_x(A_x)','U_z = U - U_x'}',[0.02 0.98], 'fontsize',fontsize)
%irf_legend(hca,{'U','U_x','U_z'},[0.02 0.1], 'fontsize',fontsize)

if 1 % |vx|, |vz|
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, abs(f_vx(B0_,b_,m_,q_,zvec))*1e-3, zvec*1e-3, abs(f_vz(B0_,E0_,b_,m_,q_,zvec))*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = '|v| (km/s)';
  irf_legend(hca,{'v_x','v_z'},[0.02 0.1], 'fontsize',fontsize)
end
if 0 % vx
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, f_vx(B0_,b_,m_,q_,zvec)*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'v_x (km/s)';
end
if 0 % vz
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, f_vz(B0_,E0_,b_,m_,q_,zvec)*1e-3)
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'v_z (km/s)';
end
if 0 % vx/vz
  hca = h(isub); isub = isub + 1;
  plot(hca, zvec*1e-3, f_vz(B0_,E0_,b_,m_,q_,zvec)./f_vx(B0_,b_,m_,q_,zvec))
  hca.XLabel.String = 'z (km)';
  hca.YLabel.String = 'v_z/v_x';
end

hca = h(isub); isub = isub + 1;
plot(hca, zvec*1e-3, abs(f_vx(B0_,b_,m_,q_,zvec)./f_vz(B0_,E0_,b_,m_,q_,zvec)))
hca.XLabel.String = 'z (km)';
hca.YLabel.String = '|v_x/v_z|';

hca = h(isub); isub = isub + 1;
vx_ = f_vx(B0_,b_,m_,q_,zvec);
vz_ = f_vz(B0_,E0_,b_,m_,q_,zvec);
x_ = cumsum(vx_);
z_ = cumsum(vz_);

plot(hca, x_, z_);
hca.XLabel.String = 'x (...)';
hca.YLabel.String = 'z (...)';

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))

c_eval('h(?).FontSize = 18;',1:numel(h))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
