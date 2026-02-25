%% integrate electron orbits for a given magnetic field
units = irf_units;
% Define time to run the particle
T = 60;


% Initial positions and velocitites
x0 = 00*1e3;
y0 = 00*1e3; % m
z0 = 1000*1e3;
vx0 = 100*1e3; % m/s
vy0 = 0*1e3;
vz0 = 100*1e3;

% Parameters
a = 5e-3; % a = E0, eg 1 mV/m
b = 20e-9; % b = B0, eg 20 nT
d = 2000e3; % d = thickness of current sheet, eg 800 km.
eps = 1e-1; % small number

Bx = @(x,y,z) b*z/d;
By = @(x,y,z) 0;
Bz = @(x,y,z) 1*-b*eps;
%Bz = @(x,y,z) 0*b*x/d/15;
Ex = @(x,y,z) 0;
Ey = @(x,y,z) a;
Ez = @(x,y,z) 0;
   
% Initial positions and velocities                                   
x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s

% Integrate trajectory
EoM = @(ttt,xxx) eom.general_proton(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
[t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
  
x = x_sol(:,1);
y = x_sol(:,2);
z = x_sol(:,3);
vx = x_sol(:,4);
vy = x_sol(:,5);
vz = x_sol(:,6);      

fontsize = 14;

hca = subplot(1,3,1);
plot3(hca,x*1e-3,y*1e-3,z*1e-3);
hold(hca,'on')
  plot3(hca,x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go');
  plot3(hca,x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'r+');
hold(hca,'off')
hca.XLabel.String = 'X (km)';
hca.YLabel.String = 'Y (km)';
hca.ZLabel.String = 'Z (km)';
hca.Box = 'on';
hca.FontSize = fontsize;
%axis(hca,'equal')

hca = subplot(1,3,2);
plot(hca,t,[x y z]*1e-3);
hca.XLabel.String = 't (s)';
hca.YLabel.String = 'distance (km)';
legend('X','Y','Z','location','northwest')
hca.FontSize = fontsize;

hca = subplot(1,3,3);
plot(hca,t,[vx vy vz]*1e-3);
hca.YLabel.String = 'speed (km/s)';
hca.XLabel.String = 't (s)';
legend('v_X','v_Y','v_Z','location','northwest')
hca.FontSize = fontsize;

%% henrik's analytical solutions, first model
C1 = b*units.e/units.mp/d;
C2 = units.e*a/units.mp;
C = C1*C2;
%k1 = (z0-vz0*airy(3,0)/airy(2,0))/(airy(1,0)-airy(0,0)*airy(3,0)/airy(2,0));
%k2 = (vz0-z0*airy(1,0)/airy(0,0))/(airy(3,0)-airy(2,0)*airy(1,0)/airy(0,0));
x_an = @(t) x0 + vx0*t;
y_an = @(t) y0 + (vy0 - C2*z0^2/2)*t + (C1/2)*t.^2;
z_an = @(t) k1*airy(0,C^(1/3)*t)+k2*airy(2,C^(1/3)*t);
z_an = @(t) k1*airy(0,C^(1/3)*t)+k2*airy(2,C^(1/3)*t);

t_an = linspace(0,T,500);
% Ai: airy(0,inp)
% Bi: airy(2,inp)

%% Snalytical solution, first model
% Get the k-constants
colors = mms_colors('matlab');
fontsize = 14;

q = units.e; 
m = units.mp;
C = sqrt((q/m)^2*a*b/d);
zfun1 = @(t) (1/3)*sqrt(t).*besselj(-1/3,(2/3)*C*t.^(3/2));
zfun2 = @(t) (1/3)*sqrt(t).*besselj(1/3,(2/3)*C*t.^(3/2));

t_an = linspace(0,T,500);
dt = t_an(2) - t_an(1);
diff_zfun1 = diff(zfun1(t))/dt;
diff_zfun2 = diff(zfun2(t))/dt;

A1p = diff_zfun1(2);
A2p = diff_zfun2(2);
A1 = zfun1(1e-10);
A2 = zfun2(1e-10);

k1 = (z0-vz0*A2p/A2)/(A1p-A1*A2p/A2);
k2 = (vz0-z0*A1p/A1)/(A2p-A2*A1p/A1);
k1 = 2.4735*z0;
k2 = 4.3995*vz0;
C1 = (b/d)*units.e/units.mp;
C2 = (a)*units.e/units.mp;
C = C1*C2;
%k1 = (z0-vz0*airy(3,0)/airy(2,0))/(airy(1,0)-airy(0,0)*airy(3,0)/airy(2,0));
%k2 = (vz0-z0*airy(1,0)/airy(0,0))/(airy(3,0)-airy(2,0)*airy(1,0)/airy(0,0));
x_an = @(t) x0 + vx0*t;
y_an = @(t) y0 + (vy0 - C1*z0^2/2)*t + (C2/2)*t.^2;
z_an = @(t) (1/3)*sqrt(t).*(k1.*besselj(-1/3,(2/3)*C*t.^(3/2)) + k2.*besselj(1/3,(2/3)*C*t.^(3/2)));


% Ai: airy(0,inp)
% Bi: airy(2,inp)

%% Analytical soluation, second model (with small normal magnetic field)
colors = mms_colors('matlab');
fontsize = 14;

q = units.e; 
m = units.mp;
C = sqrt((q/m)^2*a*b/d);
zfun1 = @(t) (1/3)*sqrt(t).*besselj(-1/3,(2/3)*C*t.^(3/2));
zfun2 = @(t) (1/3)*sqrt(t).*besselj(1/3,(2/3)*C*t.^(3/2));

t_an = linspace(0,T,500);
dt = t_an(2) - t_an(1);
diff_zfun1 = diff(zfun1(t))/dt;
diff_zfun2 = diff(zfun2(t))/dt;

A1p = diff_zfun1(2);
A2p = diff_zfun2(2);
A1 = zfun1(1e-10);
A2 = zfun2(1e-10);

k1 = (z0-vz0*A2p/A2)/(A1p-A1*A2p/A2);
k2 = (vz0-z0*A1p/A1)/(A2p-A2*A1p/A1);
k1 = 2.4735*z0;
k2 = 4.3995*vz0;
C1 = (b/d)*q/m;
C2 = (a)*b/m;
C3 = q*a/m;
C4 = vy0 - 0.5*C1*z0^2 - C2*eps*x0;
C5 = C2*C3/2;
C6 = C2*C4;
C8 = C1*z0^2/2 + C2*eps*x0 + vy0;

C = C1*C2;
%k1 = (z0-vz0*airy(3,0)/airy(2,0))/(airy(1,0)-airy(0,0)*airy(3,0)/airy(2,0));
%k2 = (vz0-z0*airy(1,0)/airy(0,0))/(airy(3,0)-airy(2,0)*airy(1,0)/airy(0,0));
x_an = @(t) x0 + vx0*t + C6*eps*t.^2/2-C5*eps*t.^3/3;
y_an = @(t) y0 + (C4-C2*eps*x0)*t + (C3/2-C2*eps*vx0/2)*t.^2;
z_an = @(t) (1/3)*sqrt(t).*(k1.*besselj(-1/3,(2/3)*C*t.^(3/2)) + k2.*besselj(1/3,(2/3)*C*t.^(3/2)));


% Ai: airy(0,inp)
% Bi: airy(2,inp)

%% henrik's analytical solutions, second model
C1 = units.e*d/units.mp/units.e
C2 = -units.e*a/units.mp;
C3 = 0;
C5 = 0;
C6 = 0;
C10 = 0;
k1 = z0;
k2 = 1;
x_an = @(t) x0+vx0*t-C5*eps*t.^3/3+C6*eps*t.^2/2;
y_an = @(t) y0-C3*t-C2*eps*x0*-C2*eps*vx0*t;
z_an = @(t) k1*airy(0,(-C10).^(1/3)*t)+k2*airy(2,(-C10).^(1/3)*t);

t_an = linspace(0,T/10,500);
% Ai: airy(0,inp)
% Bi: airy(2,inp)
%% plotton numerical and analytical together
nrows = 2;
ncols = 3;
isub = 1;
for ip = 1:nrows*ncols;
  h(ip) = subplot(nrows,ncols,ip);
end

% numerical
hca = h(isub); isub = isub + 1;
plot3(hca,x*1e-3,y*1e-3,z*1e-3);
hold(hca,'on')
  plot3(hca,x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go');
  plot3(hca,x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'r+');
hold(hca,'off')
hca.XLabel.String = 'X (km)';
hca.YLabel.String = 'Y (km)';
hca.ZLabel.String = 'Z (km)';
hca.Box = 'on';
%axis(hca,'equal')

hca = h(isub); isub = isub + 1;
plot(hca,t,[x y z]*1e-3);
hca.XLabel.String = 't (s)';
hca.YLabel.String = 'distance (km)';
legend(hca,'X','Y','Z','location','northeast')

hca = h(isub); isub = isub + 1;
plot(hca,t,[vx vy vz]*1e-3);
hca.YLabel.String = 'speed (km/s)';
hca.XLabel.String = 't (s)';
legend(hca,'v_X','v_Y','v_Z','location','northeast')

% analytical
hca = h(isub); isub = isub + 1;
plot3(hca,x_an(t_an)*1e-3,y_an(t_an)*1e-3,z_an(t_an)*1e-3);
hold(hca,'on')
  plot3(hca,x_an(t_an(1))*1e-3,y_an(t_an(1))*1e-3,z_an(t_an(1))*1e-3,'go');
  plot3(hca,x_an(t_an(end))*1e-3,y_an(t_an(end))*1e-3,z_an(t_an(end))*1e-3,'r+');
hold(hca,'off')
hca.XLabel.String = 'X (km)';
hca.YLabel.String = 'Y (km)';
hca.ZLabel.String = 'Z (km)';
hca.Box = 'on';
%axis(hca,'equal')

hca = h(isub); isub = isub + 1;
plot(hca,t_an,[x_an(t_an)' y_an(t_an)' z_an(t_an)']*1e-3);
hca.XLabel.String = 't (s)';
hca.YLabel.String = 'distance (km)';
legend(hca,'X','Y','Z','location','northeast')
hca.YScale = 'lin';


hca = h(isub); isub = isub + 1;
plot(hca,t_an(1:end-1),diff([x_an(t_an)' y_an(t_an)' z_an(t_an)'])/(t_an(2)-t_an(1))*1e-3);
hca.YLabel.String = 'speed (km/s)';
hca.XLabel.String = 't (s)';
legend(hca,'v_X','v_Y','v_Z','location','northeast')

h(2).YLim = [-0.5 1.3]*1e4;
h(5).YLim = h(2).YLim;
h(3).YLim = [-2000 5000]*1e0;
h(6).YLim = h(3).YLim;

h(1).View = [-70 40];
h(4).View = h(1).View;

for ip = 1:3
  h(ip).FontSize = fontsize;
end

for ip = [1 4]
  h(ip).XLim = [0 3000];
  h(ip).YLim = [-1000 1.7e5];
  h(ip).ZLim = [-1000 1200];
  h(ip).Position(3) = 0.2;
  h(ip).Position(1) = h(ip).Position(1)*1.1;
end
for ip = [2 3 5 6]
  h(ip).XLim = [0 T];  
  h(ip).Position(3) = 0.2;
  h(ip).Position(1) = h(ip).Position(1)*1.1;
end
%% plotting

%h =
%h(1)=subplot(1,3,1);
xlim = max(abs(x))*1e-3;
ylim = max(abs(y))*1e-3;
zlim = max(abs(z))*1e-3;
xpl = xlim*[-1 -1 1 1];
ypl = ylim*[1 -1 -1 1];
zpl = 0*[-1 -1 1 1];

if 1
% to print with opacity, figure renderer must be set to OpenGl
patch(xpl,ypl,zpl);
hp = findobj(gcf,'type','patch');
set(hp,'facealpha',0.05)
hold on
end

plot3(x*1e-3,y*1e-3,z*1e-3,'k',...
      x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go',...
      x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'rx') % plot in km's
hca=gca;
xlabel(hca,'x (GSM)')
ylabel(hca,'y (GSM)')
zlabel(hca,'z (GSM)')
hold off
%set(hca,'ylim',[])

%subplot(1,3,2)
%plot(y,vy)
%subplot(1,3,3)
%plot(z,vz)

%% Plot of bessel functions
x=0.01:0.01:15;
colors = mms_colors('matlab');
plot(x,besselj(1/3,x),x,besselj(-1/3,x),x,x*0,'k')
fontsize = 20;
ht1 = text(2,besselj(1/3,2),'J_{1/3}','HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(1,:));
ht2 = text(6,besselj(-1/3,6),'J_{-1/3}','HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(2,:));
hca = gca;
hca.FontSize = fontsize;
hca.YLim = [-0.5 2];

%% Plot of sqrt(t)*bessel functions
x=0:0.001:15;
colors = mms_colors('matlab');
plot(x,sqrt(x).*besselj(1/3,x.^(3/2)),x,sqrt(x).*besselj(-1/3,x.^(3/2)),x,x*0,'k')
fontsize = 20;
ht1 = text(2,sqrt(2).*besselj(1/3,2),'t^{1/2}J_{1/3}(t)','HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(1,:));
%ht2 = text(9,sqrt(9).*besselj(-1/3,9),'t^{1/2}J_{-1/3}(t)','HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(2,:));
ht2 = text(8,0.6,'t^{1/2}J_{-1/3}(t)','HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(2,:));
hca = gca;
hca.FontSize = fontsize;
hca.YLim = [-1 1];
hca.XLabel.String = 't';


%%
x=0:0.001:15;
colors = mms_colors('matlab');
plot(x,sqrt(x).*besselj(1/3,x.^(3/2)),x,sqrt(x).*besselj(-1/3,x.^(3/2)),x,x*0,'k')
fontsize = 20;
ht1 = text(2,sqrt(2).*besselj(1/3,2),'t^{1/2}J_{1/3}(t)','HorizontalAlignment','left','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(1,:));
%ht2 = text(9,sqrt(9).*besselj(-1/3,9),'t^{1/2}J_{-1/3}(t)','HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(2,:));
ht2 = text(8,0.6,'t^{1/2}J_{-1/3}(t)','HorizontalAlignment','right','VerticalAlignment','bottom','fontsize',fontsize,'color',colors(2,:));
hca = gca;
hca.FontSize = fontsize;
hca.YLim = [-1 1];
hca.XLabel.String = 't';

