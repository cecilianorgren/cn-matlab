% integrate electron orbits for a given magnetic field

% Define time to run the particle
T = 0.04;
units = irf_units;
electron_energy = 1000; % eV
vt = sqrt(electron_energy*units.eV*2/units.me)/1000; 

velocity_angle = 90;
% Initial positions and velocitites
x0 = 0; 
y0 = 0; % km
z0 = 0;
vx0 = 0; % km/s
vy0 = vt*cosd(velocity_angle);
vz0 = -vt*sind(velocity_angle);

% Parameters
a = 1e-3; % a = E0, eg 1 mV/m
B0 = 10e-9; % b = B0, eg 20 nT
d = 15e3; % d = thickness of current sheet, eg 800 km.
eps = 0; % small number
E0 = 1*2*1e-3;

Bx = @(z) B0*z./z;
By = @(z) z*0;
Bz = @(z) z*0;
%Bz = @(x) b*x/d/15;
Ex = @(z) z*0;
Ey = @(z) z*0;
Ez = @(z) E0*exp(-z.^2*2/(d^2));

hca = subplot(2,1,1);
zz = linspace(-d,d,20);
plotyy(hca,zz*1e-3,Ez(zz)*1e3,zz*1e-3,Bx(zz)*1e9)
hca.XLabel.String = 'z';

stopfunction = @(t,y) ExB.events(t,y,L);
options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
    
% Initial positions and velocities                                   
x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s

% Integrate trajectory
% Matlab ode solver :             
% If the integration terminates beforehand, due to the passing
% of the particle outside of the box, the integration stops, this is
% defined in options and ExB.events.

%EoM = @(ttt,xxx) mr.eom(ttt,xxx,Ex,Ey,Ez,Bx,By,Bz);
EoM = @(ttt,xxx) mms_2015Nov12.eom_separatrix(ttt,xxx,a,b,d,E0,eps);
[t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
x = x_sol(:,1);
y = x_sol(:,2);
z = x_sol(:,3);
vx = x_sol(:,4);
vy = x_sol(:,5);
vz = x_sol(:,6);      


% plotting

%h =
hca = subplot(2,1,2);

xlim = max(abs(x))*1e-3;
ylim = max(abs(y))*1e-3;
zlim = max(abs(z))*1e-3;
xpl = xlim*[-1 -1 1 1];
ypl = ylim*[1 -1 -1 1];
zpl = 0*[-1 -1 1 1];

if 0
% to print with opacity, figure renderer must be set to OpenGl
patch(xpl,ypl,zpl);
hp = findobj(gcf,'type','patch');
set(hp,'facealpha',0.05)
hold on
end

plot3(hca,x*1e-3,y*1e-3,z*1e-3,'k',...
      x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go',...
      x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'rx') % plot in km's
xlabel(hca,'x (GSM)')
ylabel(hca,'y (GSM)')
zlabel(hca,'z (GSM)')
hold off
view([1 0 0])
set(hca,'ylim',[-10 100])
set(hca,'zlim',d*[-1 1]*1e-3)

%subplot(1,3,2)
%plot(y,vy)
%subplot(1,3,3)
%plot(z,vz)


