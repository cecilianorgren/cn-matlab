


% Integrate trajectory
options = odeset();
options = odeset('AbsTol',1e-1,'RelTol',1e-1);
%options = odeset('Events',@exitBox,varargin{:});
%options = odeset('AbsTol',1e-7,'AbsTol',1e-9,'Events',@exitBox);
%options = odeset('RelTol',1e-6);

RE = 6371200; % m
me = 9.1094e-31;
mp = 1.6726e-27;
e = 1.6022e-19;
kB = 1.3806e-23;

m = me;
q = -e;

EoM = @(ttt,xxx) eom_dipole(ttt,xxx,m,q); 

U_eV = 1e3; % eV
U_J = U_eV*kB/e;
vabs = sqrt(U_eV*e*2/m); % m/s
%res=sqrt(data*Units.eV*2/Units.me)/1000;

pa = 90;
vpar = vabs*cosd(pa);
vperp = vabs*sind(pa);

x_init = [0, 5*RE, 1*RE, vperp,0,vpar];
tstart = 0;
tstop = 0.00001;

[Bx0,By0,Bz0] = phys251_fun_magnetic_dipole_field(x_init(1),x_init(2),x_init(3));
wce = e*sqrt(Bx0^2 + By0^2 + Bz0^2)/me;
wcp = e*sqrt(Bx0^2 + By0^2 + Bz0^2)/mp;
Tce = 2*pi/wce;
Tcp = 2*pi/wcp;

tstart = 0;
stop = Tce;
%disp(sprintf('tstart = %5.2f, tstop = %5.2f, [x0,y0,z0] = [%5.1f, %5.1f, %5.1f], [vx0,vy0,vz0] = [%5.2f, %5.2f, %5.2f]'))
tic;
[t,x_sol_tmp] = ode45(EoM,[tstart tstop],x_init,options);%,options); % 
toc
x_sol_tmp(:,7) = t; % x_sol = (x,y,z,vx,vy,vz,t)

%x_sol = [x_sol; x_sol_tmp];

x = x_sol_tmp(:,1);
y = x_sol_tmp(:,2);
z = x_sol_tmp(:,3);
vz = x_sol_tmp(:,4);
vy = x_sol_tmp(:,5);
vz = x_sol_tmp(:,6);

%%
[SX,SY,SZ] = sphere(20);

hca = subplot(1,1,1);
surf(hca,SX,SY,SZ,'facealpha',0.5)
colormap(hca,[1 0 1])

hold(hca,'on')
plot3(hca,x/RE,y/RE,z/RE)
plot3(hca,x(1)/RE,y(1)/RE,z(1)/RE,'go','MarkerSize',10,'Linewidth',3)
hold(hca,'off')
axis(hca,'equal')



%% Help functions
function  x_res = eom_dipole(t,x_vect,m,q)
  x = x_vect(1);
  y = x_vect(2);
  z = x_vect(3);
  vx = x_vect(4);
  vy = x_vect(5);
  vz = x_vect(6);

  %if isnan(x)
  %  1;
  %end
  %disp(sprintf('t = %g, x = %g, y = %g, z = %g',t,x,y,z))

  %method = 'spline';        

  [Bx,By,Bz] = phys251_fun_magnetic_dipole_field(x,y,z);
  Ex = 0;
  Ey = 0;
  Ez = 0;
  %disp(sprintf('B = [%g, %g, %g]',Bx,By,Bz))
  %plot3(x,y,z,'*g')
  drawnow
  
  %disp(sprintf('%.3f %.3f %.3f, %.3f %.3f %.3f, %.3f %.3f %.3f, %.3f, %.3f %.3f',Ex,Ey,Ez,Bx,By,Bz,x_vect(1),x_vect(2),x_vect(3),x_vect(4),x_vect(5),x_vect(6)))        
  %it = size(x_sol_all,1);
  %x_sol_all(it+1,1) = Ex;
  %x_sol_all(it+1,2) = Ey;
  %x_sol_all(it+1,3) = Ez;
  %x_sol_all(it+1,4) = Bx;
  %x_sol_all(it+1,5) = By;
  %x_sol_all(it+1,6) = Bz;
  %x_sol_all(it+1,7) = t;

  % Equations to be solved
  x_res = zeros(6,1);
  x_res(1) = vx; % dx/dt = vx;
  x_res(2) = vy; % dy/dt = vy;
  x_res(3) = vz; % dz/dt = vz;
  x_res(4) = (q/m)*(Ex + vy*Bz - vz*By);
  x_res(5) = (q/m)*(Ey + vz*Bx - vx*Bz);
  x_res(6) = (q/m)*(Ez + vx*By - vy*Bx);                                              

end   