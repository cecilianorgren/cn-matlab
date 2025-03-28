% illustration_diffusion_orbit
%% Gyromotion, fake

vy = @(z) 0.5*sin(z);
Ay = @(z) -0.5*sin(z);

z = linspace(-pi/2,pi/2,100);

h = setup_subplots(1,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,z,vy(z),z,Ay(z),z,vy(z)+Ay(z))
hca.Box = 'on';
%hca.Visible = 'off';
hca.XLim = z([1 end]);

irf_legend(hca,{'m_iv_y','eA_y','p_y'}',[0.02 0.7],'fontsize',20)

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))
c_eval('h(?).FontSize = 20;',1:numel(h))


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
Ay = @(x,y,z,t) -(B0*lz)*log(cosh((z/lz))) - 1*E0*t;
Bx = @(x,y,z) B0*tanh(z/lz);
By = @(x,y,z) x*0;
Bz = @(x,y,z) 0*x/a^2;
Ex = @(x,y,z) x*0;
Ey = @(x,y,z) x*0 + z*0 + E0; % V/m
Ez = @(x,y,z) x*0;
Ez = @(x,y,z) 10*E0*-1*(z/1/lz).*exp(-(z/(lz)).^2)*1; % finite Ez makes H change with time, and thus with y (since the electron moves along y)
Phi = @(x,y,z) -y.*Ey(x,y,z) - z.*Ez(x,y,z);
H_p = @(x,y,z,px,py,pz,t) 0*(1/2/m)*(px).^2 + (1/2/m)*(px-q*Ay(x,y,z,t)).^2 + 0*(1/2/m)*(pz).^2 + 1*q*Phi(x,y,z);
H_v = @(x,y,z,vx,vy,vz) (m/2)*(vx).^2 + (m/2)*(vy).^2 + (m/2)*(vz).^2 + 1*q*Phi(x,y,z);


options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,xvec([1 end])),...
                 'AbsTol',1e-6);
options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,2,[-400*1e3 400*1e4]),...
                 'AbsTol',1e-22,'RelTol',1e-22);
%options = odeset();
EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 1;

nP = 1;

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
  p(ip).px = units.me*p(ip).vx;
  p(ip).py = units.me*p(ip).vy - units.e*Ay(p(ip).x,p(ip).y,p(ip).z,p(ip).t);
  p(ip).pz = units.me*p(ip).vz;
  p(ip).py_mvy = units.me*p(ip).vy;
  p(ip).py_qAy = - units.e*Ay(p(ip).x,p(ip).y,p(ip).z,p(ip).t);
  p(ip).H_v = H_v(p(ip).x,p(ip).y,p(ip).z,p(ip).vx,p(ip).vy,p(ip).vz);
  p(ip).H_p = H_p(p(ip).x,p(ip).y,p(ip).z,p(ip).px,p(ip).py,p(ip).pz,p(ip).t);
  p(ip).Phi = Phi(p(ip).x,p(ip).y,p(ip).z);
  p(ip).Phi = Phi(p(ip).x-p(ip).x(1),p(ip).y-p(ip).y(1),p(ip).z-p(ip).z(1));
  p(ip).U = (m/2)*(p(ip).vx).^2 + (m/2)*(p(ip).vy).^2 + (m/2)*(p(ip).vz).^2;
end
nP = numel(p);

colors = pic_colors('matlab');
%colors = [colors(2)];
colors(3,:) = colors(1,:);
%colors(3,:) = colors(5,:).^0.5;
colors(1,:) = colors(1,:).^0.2;
linewidth = 1.5;



if 1 % Figure
  %%
  nRows = 5;
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
  
  if 1 % py
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs                
      plot(hca,p(ip).t,p(ip).py)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 't (s)';
    hca.YLabel.String = 'p_y (...)';
  end
  if 1 % mvy
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs                
      plot(hca,p(ip).t,p(ip).py_mvy)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 't (s)';
    hca.YLabel.String = 'm_ev_y (...)';
  end
  if 1 % -qAy
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs                
      plot(hca,p(ip).t,p(ip).py_qAy)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 't (s)';
    hca.YLabel.String = '-eA_y (...)';
  end
  if 0 % H
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs                
      plot(hca,p(ip).t,p(ip).H_v,p(ip).t,p(ip).U,p(ip).t,-p(ip).Phi*units.e)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 't (s)';
    hca.YLabel.String = 'H (...)';
  end
  if 0 % H(y)
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs                
      plot(hca,p(ip).y*1e-3,p(ip).H_v,p(ip).y*1e-3,p(ip).U,p(ip).y*1e-3,-p(ip).Phi*units.e)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 'y (km)';
    hca.YLabel.String = 'H (...)';
  end
  if 0 % H_v, H_p vs t
    hca = h(isub); isub = isub + 1;
    
    for ip = iPs                
      plotyy(hca,p(ip).t,p(ip).H_v,p(ip).t,p(ip).H_p)
      if ip == 1, hold(hca,'on'); end
    end
    hold(hca,'off')
    hca.XLabel.String = 'y (km)';
    hca.YLabel.String = 'H (...)';
  end
  if 0 % energy
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
