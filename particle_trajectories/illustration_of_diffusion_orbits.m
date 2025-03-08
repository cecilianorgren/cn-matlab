% illustration_diffusion_orbit

%% xz
% Define magnetic field
a = 5*1e3;
b = 1*1e3;
xvec = 0;
zvec = b*linspace(-10,10,100);
yvec = 0;

[X,Y,Z] = ndgrid(xvec,yvec,zvec);
%dx = x(2) - x(1);
%dy = y(2) - y(1);
%dz = z(2) - z(1);
%x_xline = x;
%y_xline = x*b/a;

Ay = @(x,y,z) (x/a).^2 - (z/b).^2;
AY0 = Ay(X,Y,Z);


% Integrate orbits
units = irf_units;
m = units.me;
q = -units.e;
B0 = 10e-9; % T
E0 = 3e-3; % V/m
lz = 30e3; % m
Bx = @(x,y,z) B0*tanh(z/lz);
By = @(x,y,z) x*0;
Bz = @(x,y,z) 0*x/a^2;
Ex = @(x,y,z) x*0;
Ey = @(x,y,z) x*0 + z*0 + E0; % V/m
Ez = @(x,y,z) x*0;
Ez = @(x,y,z) 10*E0*-1*(z/1/lz).*exp(-(z/(lz)).^2);

options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,xvec([1 end])),...
                 'AbsTol',1e-6);
options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,2,[-400*1e3 400*1e4]),...
                 'AbsTol',1e-12);
%options = odeset();
EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 1;

nP = 10;

vz0 = E0/B0; % m/s
Te0 = 80; % eV
Te = Te0 + 0.3*Te0*randn([nP 1]); % eV
Te(Te<0) = 0;
vt = sqrt(2*Te*units.eV/units.me); % m/s
wce = units.e*B0/units.me;
rhoe = vt/wce;
z0 = 3*lz;
ph = rand(nP,1)*360;
vy =  vt.*cosd(ph);
vz =  vt.*sind(ph);

x_init_all = [repmat([0 0 z0],nP,1) vy*0 vy vz];


clear p
ylim = [0 0];
zlim = [0 0];
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
  ylim(1) = min([ylim(1); p(ip).y]);
  ylim(2) = max([ylim(2); p(ip).y]);
  zlim(1) = min([zlim(1); p(ip).z]);
  zlim(2) = max([zlim(2); p(ip).z]);
end
nP = numel(p);

colors = pic_colors('matlab');
%colors = [colors(2)];
colors(3,:) = colors(1,:);
%colors(3,:) = colors(5,:).^0.5;
colors(1,:) = colors(1,:).^0.2;
linewidth = 1.5;



if 1 % Figure
  nRows = 3;
  nCols = 1;
  h = setup_subplots(nRows,nCols);
  isub = 1;

  if 1 % Bx, Ez, Ey
    hca = h(isub); isub = isub + 1;

    zz = linspace(-3*lz,3*lz,1000);
    hca.ColorOrder = pic_colors('matlab');
    plot(hca, zz*1e-3, Bx(0,0,zz)*1e9, zz*1e-3, Ey(0,0,zz)*1e3, zz*1e-3, Ez(0,0,zz)*1e3)
    irf_legend(hca,{'B_x','E_y','E_z'},[0.02 0.98])
    hca.XLabel.String = 'z (km)';
    hca.YLabel.String = 'E (mV/m), B (nTkm)';
  end    

  iPs = 1:nP; 
  if 1 % (y,z)
    hca = h(isub); isub = isub + 1;

    zz = linspace(-max(zlim),max(zlim),200);
    %zz = linspace(zlim(1),zlim(2),200);
    yy = ylim;
    [Y,Z] = ndgrid(yy,zz);
    hp = pcolor(hca,Y*1e-3,Z*1e-3,Bx(Y*0,Y,Z));
    hp.FaceAlpha = 1;
    shading(hca,'flat')
    colormap(hca,pic_colors('blue_red'))    

  
    hold(hca,'on')
    for ip = iPs
      plot(hca,p(ip).y*1e-3,p(ip).z*1e-3,'LineWidth',2)
    end
    hold(hca,'off')
    
    if 1 % E quivers
      hold(hca,'on')      
      zz = zlim(1):10e3:zlim(2);
      yy = ylim(1):15e3:ylim(2);
      [Y,Z] = ndgrid(yy,zz);
      ey = Ey(Y*0,Y,Z);
      ez = Ez(Y*0,Y,Z);
      hq = quiver(hca,Y*1e-3,Z*1e-3,ey*1e3,ez*1e3,0,'color','k','marker','.','MaxHeadSize',0);
      hold(hca,'off'),
    end

    axis(hca,'equal')
    %hca.YLim = (zz([1 end])*1e-3+ [-20 0]);
    hca.YLim = zlim*1e-3;
    hca.XLim = ylim*1e-3;
    hca.XLabel.String = 'y (km)';
    hca.YLabel.String = 'z (km)';
  end    
  if 0 % E + vzBx
    hca = h(isub); isub = isub + 1;
    
    
    for ip = iPs
      ey = Ey(p(ip).x,p(ip).y,p(ip).z);
      vzbx = p(ip).vz.*Bx(p(ip).x,p(ip).y,p(ip).z);
      plot(hca,p(ip).t,ey+vzbx)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 'time (s)';
    hca.YLabel.String = 'E_y + v_zB_x (mV/m)';
  end    
  if 0 % dU/dt vs t
    hca = h(isub); isub = isub + 1;
    
    
    for ip = iPs
      t = p(ip).t;
      z = p(ip).z;
      U = units.me*( p(ip).vx.^2 + p(ip).vy.^2 + p(ip).vz.^2)/2;      
      dUdt = gradient(U,t);
      plot(hca,t,dUdt)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 'time (s)';
    hca.YLabel.String = 'dU/dt (eV/s)';
  end    
  if 0 % dU/dt vs z
    hca = h(isub); isub = isub + 1;
    
    
    for ip = iPs
      t = p(ip).t;
      z = p(ip).z;
      U = units.me*( p(ip).vx.^2 + p(ip).vy.^2 + p(ip).vz.^2)/2;      
      dUdt = gradient(U,t);
      plot(hca,dUdt,z*1e-3)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.YLabel.String = 'z (km)';
    hca.XLabel.String = 'dU/dt (eV/s)';
  end  
  if 1 % energy
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs      
      U = units.me*( p(ip).vx.^2 + p(ip).vy.^2 + p(ip).vz.^2)/2;      
      plot(hca,p(ip).t,U/units.eV)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 't (s)';
    hca.YLabel.String = 'U (eV)';
  end
  if 0 % (vy,vz)
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs      
      vy = p(ip).vy;
      vz = p(ip).vz;
      plot(hca,vy*1e-3,vz*1e-3)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    axis(hca,'equal')
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = 'v_z (km/s)';
  end
  if 0 % forces
    hca = h(isub); isub = isub + 1;

    
    for ip = iPs
      if ip == 1, hold(hca,'on'); end
      plot(hca,p(ip).t,Ey(p(ip).x,p(ip).y,p(ip).z))
    end
    hold(hca,'off')
  end
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
  %%
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
  h(1).Position(4) = 0.3;
  
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