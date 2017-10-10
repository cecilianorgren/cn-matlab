global B0 lr lz phi0 L rlim zlim
nParticles = 1;
tint = [0 0.01];  % s

% Initialize starting position
x0 = -rlim + 2*rlim*rand(nParticles,1); 
y0 = -rlim + 2*rlim*rand(nParticles,1); 
z0 = -zlim + 2*zlim*rand(nParticles,1); 

% Initialize starting velocity
vtpar = cn_eV2v(5000,'ev'); % eV -> km/s (parallel thermal velocity)
vtper = cn_eV2v(5000,'ev'); % eV -> km/s (perpendicular thermal velocity)
vx0 = vtper*randn(nParticles,1)/sqrt(2);
vy0 = vtper*randn(nParticles,1)/sqrt(2);
vz0 = vtpar*randn(nParticles,1)/sqrt(2); % shift parallel velocity by electron hole velocity

% Terminate integration when the particle passes outside the box in
% z-direction
options = odeset('Events',@ExB.events);%,'OutputSel',1,'Refine',refine);


% Integrate
tic
[t,x_sol] = ode23(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0]*1e3,options);
toc
tic 
[t2,x_sol2] = ode45(@ExB.EquationOfMotion,tint,[x0;y0;z0;vx0;vy0;vz0]*1e3);
toc

x = x_sol(:,1); 
y = x_sol(:,2);
z = x_sol(:,3);
vx = x_sol(:,4);
vy = x_sol(:,5);
vz = x_sol(:,6);

x2 = x_sol2(:,1); 
y2 = x_sol2(:,2);
z2 = x_sol2(:,3);
vx2 = x_sol2(:,4);
vy2 = x_sol2(:,5);
vz2 = x_sol2(:,6);


plot3(x*1e-3,y*1e-3,z*1e-3,x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go');
%plot3(x2*1e-3,y2*1e-3,z2*1e-3,x2(1)*1e-3,y2(1)*1e-3,z2(1)*1e-3,'go');
legend('23','34')
t(end)
