% Forward and backward integration
% Model parameters
ic = 1;
mms_2015Nov12.Bmodel;
limN = 30e3;
T = 0.1;
options = odeset('AbsTol',1e-3,'RelTol',1e-3);
options = odeset('AbsTol',1e-18);
options = odeset('RelTol',1e-6); % this makes results correct

% Forward integration
x0 = 0; % m
y0 = 0;
z0 = 20e3; 

vx0 = 0000e3; % m/s
vy0 = 0000e3;
vz0 = 5000e3; 

x_init = [x0 y0 z0 vx0 vy0 vz0];

%stopfunction = @(t,y) eom.lim(t,y,limN);
%options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine)

EoM = @(ttt,xxx) eom.general(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
[t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
x_sol_forward = x_sol;
rf = x_sol_forward(:,1:3);
vf = x_sol_forward(:,4:6);
tf = t;
        
% Backward integration
x_init = x_sol_forward(end,:);

%stopfunction = @(t,y) eom.lim(t,y,limN);
%options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine)

EoM = @(ttt,xxx) eom.general_backward(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
[t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
x_sol_backward = x_sol;
rb = x_sol_backward(:,1:3);
vb = x_sol_backward(:,4:6);
tb = t;

% Plot different orbits to compare integration
colors = mms_colors('matlab');
%colors = colors([1 2 4:end],:);

subplot(2,2,1)
plot3(rf(:,1)*1e-3,rf(:,2)*1e-3,rf(:,3)*1e-3,'color',colors(1,:)); hold on
plot3(rf(1,1)*1e-3,rf(1,2)*1e-3,rf(1,3)*1e-3,'go'); 
plot3(rf(end,1)*1e-3,rf(end,2)*1e-3,rf(end,3)*1e-3,'rx'); 
xlabel('L'); ylabel('M'); zlabel('N'); view([0 1 0])
hold off;

subplot(2,2,2)
plot3(rb(:,1)*1e-3,rb(:,2)*1e-3,rb(:,3)*1e-3,'color',colors(3,:)); hold on
plot3(rb(1,1)*1e-3,rb(1,2)*1e-3,rb(1,3)*1e-3,'go'); 
plot3(rb(end,1)*1e-3,rb(end,2)*1e-3,rb(end,3)*1e-3,'rx'); 
xlabel('L'); ylabel('M'); zlabel('N'); view([0 1 0])
hold off;

subplot(2,2,3)
lines = plot(tf,vf,tb(end)-tb,vb);
c_eval('lines(?).Color = colors(?,:);',1:3)
c_eval('lines(?).Color = colors(?-3,:).^3;',4:6)

subplot(2,2,4)
lines = plot(tf,rf,tb(end)-tb,rb);
c_eval('lines(?).Color = colors(?,:).^2;',1:3)
c_eval('lines(?).Color = colors(?-3,:).^(1/2); lines(?).LineStyle = ''-'';',4:6)
