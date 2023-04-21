%% Coloumb scattering
e = 1.6022e-19;
eps0 = 8.8542e-12;
kB = 1.3806e-23;
me = 9.1094e-31;
mp = 1836*me;

m = mp;
q =   e;
qtarget = e;

pot_r = @(r,q) q./(4*pi*eps0*r);
pot_xy = @(x,y,q) q./(4*pi*eps0*sqrt(x.^2+y.^2));

U_xy = @(X,Y,q,qtarget) q*pot_xy(X,Y,qtarget);

H_xy = @(x,y,vx,vy,q,m) pot_xy(x,y,q) + (vx.^2 + vy.^2)*m/2;


Ex = @(x,y,q) x.*q./(4*pi*eps0*(x.^2+y.^2).^1.5);
Ey = @(x,y,q) y.*q./(4*pi*eps0*(x.^2+y.^2).^1.5);


nx = 500;
ny = 501;
xvec = 1e0*linspace(-1,1,nx); % m
yvec = 1e0*linspace(-1,1,ny); % m



U_eV = 5e0; % eV
U_J = U_eV*kB/e;
vabs = sqrt(U_eV*e*2/m); % m/s

vx0 = [1e2 2e2 5e2 1e3 2e3 5e3 1e4 2e4 5e4 1e5 2e5 5e5]; % m/s

Ekin0 = m*vx0.^2/2/e; 

x0 = xvec(1)*ones(size(vx0));
y0 = yvec(fix(nx/2)+3)*ones(size(vx0));
y0 = 0.0001*ones(size(vx0));
vy0 = 0*ones(size(vx0)); % m/s

[X,Y] = ndgrid(xvec,yvec);
R = sqrt(X.^2 + Y.^2);

Pxy = pot_xy(X,Y,q);
Pr = pot_r(R,q);
Uxy = U_xy(X,Y,q,qtarget);
EX = Ex(X,Y,q);
EY = Ey(X,Y,q);
%Hxy0 = H_xy(x0,y0,vx0,vy0,q,m);

% integrate particle
clear p
for ip = 1:numel(x0)
  x_init = [x0(ip),y0(ip),vx0(ip),vy0(ip)];
  tstart = 0;
  tstop = 1.5*(xvec(end)-xvec(1))/vx0(ip);
  
  options = odeset('AbsTol',1e-16,'RelTol',1e-16);
  EoM = @(ttt,xxx) eom_coloumb(ttt,xxx,m,q,qtarget); 
  
  
  tic;
  [t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);%,options); % 
  toc
  x_sol_tmp(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  
  %x_sol = [x_sol; x_sol_tmp];
  
  x = x_sol_tmp(:,1);
  y = x_sol_tmp(:,2);
  vx = x_sol_tmp(:,3);
  vy = x_sol_tmp(:,4);

  v = sqrt(vx.^2 + vy.^2);
  deflection_angle = atan2d(vy(end),vx(end));

  p(ip).x0 = x0(ip);
  p(ip).y0 = y0(ip);
  p(ip).vx0 = vx0(ip);
  p(ip).vy0 = vy0(ip);
  p(ip).U0_eV = Ekin0(ip);
  p(ip).t = t;
  p(ip).x = x;
  p(ip).y = y;
  p(ip).vx = vx;
  p(ip).vy = vy;
  p(ip).deflection_angle = deflection_angle;
end

nrows = 1; ncols = 1; ipanel = 0;
h = gobjects([nrows,ncols]);
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;



if 0
  hca = h(isub); isub = isub + 1;  
  [C,hc] = contour(hca,X,Y,log10(abs(Uxy)),20,'linewidth',1,'color',0.5+[0 0 0]);
  %clabel(C,hc,'LabelSpacing',72,'Color',0.5+[0 0 0],'FontWeight','bold');
  if 0 % E
    hold(hca,'on') 
    plEX = EX; plEX(abs(plEX)>prctile(abs(plEX),90)) = NaN;
    plEY = EY; plEY(abs(plEY)>prctile(abs(plEY),90)) = NaN;
    quiver(hca,X,Y,sign(EX).*abs(EX),sign(EY).*abs(EY))
    hold(hca,'off') 
  end
  if 1 % particles
    hold(hca,'on')
    clear hl
    for ip = 1:numel(x0)
      hl(ip) = plot(hca,p(ip).x,p(ip).y,'linewidth',1,'linestyle','-');
    end
    %hleg = legend(hl,arrayfun(@(s) sprintf('v_{x0} = %g m/s',s),cat(1,p.vx0),'UniformOutput',false),'Orientation','vertical');
  end
  hold(hca,'off')
  hca.XLim = xvec([1 end]);
  hca.YLim = yvec([1 end]);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'y (m)';
  %axis(hca,'equal')
  axis(hca,'square')
end

if 1 % deflection angle as a function of initial speed
  hca = h(isub); isub = isub + 1;  
  plot(hca,cat(1,p.vx0)*1e-3,cat(1,p.deflection_angle),'-','color',[0.5 0.5 0.5],'LineWidth',1)
  hold(hca,'on')
  clear hl
  for ip = 1:numel(x0)
    hl(ip) = plot(hca,p(ip).vx0*1e-3,cat(1,p(ip).deflection_angle),'o','markersize',10,'LineWidth',3);
  end
  %hleg = legend(hl,arrayfun(@(s) sprintf('v_{x0} = %g m/s',s.vx0),p.vx0,'UniformOutput',false),'Orientation','vertical');
  hleg = legend(hl,arrayfun(@(s) sprintf('v_{x0} = %g km/s, U_{kin0} = %g eV',s.vx0*1e-3,s.U0_eV),p,'UniformOutput',false),'Orientation','vertical');
  hold(hca,'off')
  hca.XLabel.String = 'v_{x0} (km/s)';
  hca.YLabel.String = 'Deflection angle, \theta (deg)';
  %axis(hca,'square')
  hca.YLim = [0 180];
  hca.YScale = 'log';
  irf_legend(hca,sprintf('y_0 = %1.4f m',y0(1)),[0.02 0.98],'color','k','fontsize',16)
end

c_eval('h(?).FontSize = 16;',1:numel(h))
hl = findobj(h(1),'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))


%% Runaway electrons
e = 1.6022e-19;
eps0 = 8.8542e-12;
kB = 1.3806e-23;
me = 9.1094e-31;
mp = 1836*me;

m = mp;
q =   -e;
qtarget = e;

Ex = @(x,y,q) x.*q./(4*pi*eps0.*(x.^2+y.^2).^1.5);
Ey = @(x,y,q) y.*q./(4*pi*eps0.*(x.^2+y.^2).^1.5);


L = 1e-12;
x1 = -L;
x2 = L;

y1 = -L;
y2 = L;

x = linspace(x1,x2,120);
y = linspace(y1,y2,121);

[X,Y] = ndgrid(x,y);

% Initialize targets
nT = 1000;
xT = x1 + (x2-x1)*rand(nT,1);
yT = y1 + (y2-y1)*rand(nT,1);


EX = zeros(size(X));
EY = zeros(size(X));
for it = 1:nT
  EX = EX + Ex(X-xT(it),Y-yT(it),qtarget);
  EY = EY + Ey(X-xT(it),Y-yT(it),qtarget);
end


imagesc(EX)

%%


U_eV = 5e0; % eV
U_J = U_eV*kB/e;
vabs = sqrt(U_eV*e*2/m); % m/s

vx0 = [1e2 2e2 5e2 1e3 2e3 5e3 1e4 2e4 5e4 1e5 2e5 5e5]; % m/s

Ekin0 = m*vx0.^2/2/e; 

x0 = xvec(1)*ones(size(vx0));
y0 = yvec(fix(nx/2)+3)*ones(size(vx0));
y0 = 0.0001*ones(size(vx0));
vy0 = 0*ones(size(vx0)); % m/s

[X,Y] = ndgrid(xvec,yvec);
R = sqrt(X.^2 + Y.^2);

Pxy = pot_xy(X,Y,q);
Pr = pot_r(R,q);
Uxy = U_xy(X,Y,q,qtarget);
EX = Ex(X,Y,q);
EY = Ey(X,Y,q);
%Hxy0 = H_xy(x0,y0,vx0,vy0,q,m);

% integrate particle
clear p
for ip = 1:numel(x0)
  x_init = [x0(ip),y0(ip),vx0(ip),vy0(ip)];
  tstart = 0;
  tstop = 1.5*(xvec(end)-xvec(1))/vx0(ip);
  
  options = odeset('AbsTol',1e-16,'RelTol',1e-16);
  EoM = @(ttt,xxx) eom_coloumb(ttt,xxx,m,q,qtarget); 
  
  
  tic;
  [t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);%,options); % 
  toc
  x_sol_tmp(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)
  
  %x_sol = [x_sol; x_sol_tmp];
  
  x = x_sol_tmp(:,1);
  y = x_sol_tmp(:,2);
  vx = x_sol_tmp(:,3);
  vy = x_sol_tmp(:,4);

  v = sqrt(vx.^2 + vy.^2);
  deflection_angle = atan2d(vy(end),vx(end));

  p(ip).x0 = x0(ip);
  p(ip).y0 = y0(ip);
  p(ip).vx0 = vx0(ip);
  p(ip).vy0 = vy0(ip);
  p(ip).U0_eV = Ekin0(ip);
  p(ip).t = t;
  p(ip).x = x;
  p(ip).y = y;
  p(ip).vx = vx;
  p(ip).vy = vy;
  p(ip).deflection_angle = deflection_angle;
end

nrows = 1; ncols = 1; ipanel = 0;
h = gobjects([nrows,ncols]);
for irow = 1:nrows, for icol = 1:ncols, ipanel = ipanel + 1; h(irow,icol) = subplot(nrows,ncols,ipanel); end, end
isub = 1;



if 0
  hca = h(isub); isub = isub + 1;  
  [C,hc] = contour(hca,X,Y,log10(abs(Uxy)),20,'linewidth',1,'color',0.5+[0 0 0]);
  %clabel(C,hc,'LabelSpacing',72,'Color',0.5+[0 0 0],'FontWeight','bold');
  if 0 % E
    hold(hca,'on') 
    plEX = EX; plEX(abs(plEX)>prctile(abs(plEX),90)) = NaN;
    plEY = EY; plEY(abs(plEY)>prctile(abs(plEY),90)) = NaN;
    quiver(hca,X,Y,sign(EX).*abs(EX),sign(EY).*abs(EY))
    hold(hca,'off') 
  end
  if 1 % particles
    hold(hca,'on')
    clear hl
    for ip = 1:numel(x0)
      hl(ip) = plot(hca,p(ip).x,p(ip).y,'linewidth',1,'linestyle','-');
    end
    %hleg = legend(hl,arrayfun(@(s) sprintf('v_{x0} = %g m/s',s),cat(1,p.vx0),'UniformOutput',false),'Orientation','vertical');
  end
  hold(hca,'off')
  hca.XLim = xvec([1 end]);
  hca.YLim = yvec([1 end]);
  hca.XLabel.String = 'x (m)';
  hca.YLabel.String = 'y (m)';
  %axis(hca,'equal')
  axis(hca,'square')
end

if 1 % deflection angle as a function of initial speed
  hca = h(isub); isub = isub + 1;  
  plot(hca,cat(1,p.vx0)*1e-3,cat(1,p.deflection_angle),'-','color',[0.5 0.5 0.5],'LineWidth',1)
  hold(hca,'on')
  clear hl
  for ip = 1:numel(x0)
    hl(ip) = plot(hca,p(ip).vx0*1e-3,cat(1,p(ip).deflection_angle),'o','markersize',10,'LineWidth',3);
  end
  %hleg = legend(hl,arrayfun(@(s) sprintf('v_{x0} = %g m/s',s.vx0),p.vx0,'UniformOutput',false),'Orientation','vertical');
  hleg = legend(hl,arrayfun(@(s) sprintf('v_{x0} = %g km/s, U_{kin0} = %g eV',s.vx0*1e-3,s.U0_eV),p,'UniformOutput',false),'Orientation','vertical');
  hold(hca,'off')
  hca.XLabel.String = 'v_{x0} (km/s)';
  hca.YLabel.String = 'Deflection angle, \theta (deg)';
  %axis(hca,'square')
  hca.YLim = [0 180];
  hca.YScale = 'log';
  irf_legend(hca,sprintf('y_0 = %1.4f m',y0(1)),[0.02 0.98],'color','k','fontsize',16)
end

c_eval('h(?).FontSize = 16;',1:numel(h))
hl = findobj(h(1),'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))




%% Help functions
function  x_res = eom_coloumb(t,x_vect,m,q,qtarget)

  eps0 = 8.8542e-12;

  x = x_vect(1);
  y = x_vect(2);
  vx = x_vect(3);
  vy = x_vect(4);
   
  Ex = x.*qtarget./(4*pi*eps0*(x.^2+y.^2).^1.5);
  Ey = y.*qtarget./(4*pi*eps0*(x.^2+y.^2).^1.5);


  % Equations to be solved
  x_res = zeros(4,1);
  x_res(1) = vx; % dx/dt = vx;
  x_res(2) = vy; % dy/dt = vy;
  x_res(3) = (q/m)*(Ex); % dvx/dt = ax;
  x_res(4) = (q/m)*(Ey); % dvy/dt = ay;

end   