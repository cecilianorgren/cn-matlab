% fields models
Bx = @(x,y,z) x*0 + y*0 + z*0 + 10*1e-9*tanh(z/10e3); % L
By = @(x,y,z) x*0 + y*0 + z*0 + 5e-9; % M
Bz = @(x,y,z) x*0 + y*0 + z*0; % N
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0 + 1e-3;
Ez = @(x,y,z) x*0 + y*0 + z*0;

% time_reversal
x0 = 0;
y0 = 0;
z0 = 0;
vx0 = 0;
vy0 = 0;
vz0 = 10000;

x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s

% Integrate trajectory forward
T = 0.01; % s
stopfunction = @(t,y) eom.lim(t,y,limN);
options = odeset('RelTol',1e-11);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

EoM = @(ttt,xxx) eom.general(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
%EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
[t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options

x = x_sol(:,1);
y = x_sol(:,2);
z = x_sol(:,3);
vx = x_sol(:,4);
vy = x_sol(:,5);
vz = x_sol(:,6); 
  
% Integrate trajectory backward
figure(45)
plot3(x*1e-3,y*1e-3,z*1e-3)  