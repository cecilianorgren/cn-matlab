% integrate electron orbits for a given magnetic field

% Define time to run the particle
T = 10;

% Initial positions and velocitites
x0 = 0; % km
y0 = 0; 
z0 = 0;
vx0 = 1; % km/s
vy0 = 0;
vz0 = 0;

% Parameters
d = 10; % spatial scale of sine field
E0 = 100; % sin electric field

Ez = @(z) E0*sin(pi/d*z);

hca = subplot(2,1,1);
zz = linspace(-d,d,20);
plot(hca,zz,Ez(zz))
hca.XLabel.String = 'z';
   
% Initial positions and velocities                                   
x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s

% Integrate trajectory
% Matlab ode solver :             
% If the integration terminates beforehand, due to the passing
% of the particle outside of the box, the integration stops, this is
% defined in options and ExB.events.
options = odeset('RelTol',1e-10,'AbsTol',1e-12);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

EoM = @(ttt,xxx) thesis.spatial_electric_field(ttt,xxx,a,E0);
[t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
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

plot(hca,x*1e-3,vx*1e-3)
%%
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


%% Plot simultaneously the electron orbits of several particles
nRows = 2;
nCols = 3;
units = irf_units;

tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
CS_normal_velocity = 70; % km/s

ic = 1;

c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp.resample(obsB);'...
],ic)

if 0
obsB = gseB1.tlim(tintModel);
obsE = gseE1.tlim(tintModel);
obsEpar = gseE1par.tlim(tintModel);
obsEperp = gseE1perp.tlim(tintModel);
obsVepar = gseVe1par.tlim(tintModel);
obsVeperp = gseVe1perp.tlim(tintModel);
end

zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
hca = subplot(nRows,nCols,[1 2]);
colors = hca.ColorOrder;
linesObs = plot(hca,zObs,[obsB.data],'.');


% Model parameters
a = 1e-3; % a = E0, eg 1 mV/m
B0 = 10e-9; % b = B0, eg 20 nT
d = 28e3; % d = thickness of current sheet, eg 800 km.
eps = 0.01; % small number
E0 = 2*1e-3;

%Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(z) -B0*tanh(z*pi/d).*(1-exp(-z.^2*5/(d^2)))-1*B0/6;
By = @(z) 0.5*B0-3*0.5*B0*sin(2/3*pi/d*z).*exp(-z.^2*2/(d^2)).*(1-exp(-z.^2*5/(d^2)));
Bz = @(z) -B0*eps+z*0;
Ex = 0;
Ey = 0;
Ez = @(z) -E0*sin(pi/d*(z));

% First plot data and model fit to data
hca.ColorOrder = colors;
hold(hca,'on')
hca = subplot(nRows,nCols,[1 2]);
zMod = linspace(-d*1.5,d*1.5,20);
%plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
lineMod = plot(hca,zMod*1e-3,[Bx(zMod)]*1e9,'--','color',colors(1,:));
plot(hca,zMod*1e-3,[By(zMod)]*1e9,'--','color',colors(2,:))
plot(hca,zMod*1e-3,[Bz(zMod)]*1e9,'--','color',colors(3,:))
hold(hca,'off')
legend([linesObs(1) lineMod],{'Observed data','Model fit'})

hca.Title.String = 'Magnetic field';
hca.YGrid = 'on';
hca.YLabel.String = 'B (nT)';
hca.XLabel.String = 'N (km)';


% Integration
T = 0.15;
electron_energy = 200; % eV
vt = sqrt(electron_energy*units.eV*2/units.me)/1000; 

velocity_angle = -40;
% Initial positions and velocitites
x0 = 0; 
y0 = 0; % km
z0 = 20;
vx0 = 0; % km/s
vy0 = vt*cosd(velocity_angle);
vz0 = -vt*sind(velocity_angle);

% Plot the electron distribution in the center of the current sheet, and
% where the test particle lies in that distribution.
%time = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt');
time = tintObs(1) + 0.5*(tintObs.stop-tintObs.start)+z0/CS_normal_velocity;
hca = subplot(nRows,nCols,[1 2]);
hold(hca,'on')
h_estart = plot(z0*[1 1],hca.YLim,'k:');
hold(hca,'off')

vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0.5 4.5];  
%projclim = [4 7];  
vlabels = {'v_M','v_N','v_L'};
hca = subplot(nRows,nCols,3);
%mms.plot_projection(hca,ePDist1.convertto('1/(cm^2 s sr keV)'),'tint',time,'xyz',[M;N;L],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scPot1,'vlabel',vlabels);
mms.plot_projection(hca,ePDist1.convertto('s^3/km^6'),'tint',time,'xyz',[M;N;L],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scPot1,'vlabel',vlabels);
colormap('jet')
hold(hca,'on')
v0 = [vx0,vy0,vz0];
v0LMN = v0*1e-3;
h_v0 = plot3(hca,v0LMN(2),v0LMN(3),v0LMN(1),'sk'); 
h_v0.MarkerSize = 10;
h_v0.MarkerFaceColor = 0+[1 1 1];
hold(hca,'off')
hca.Title.String = time.utc;

vlabels = {'v_L','v_N','-v_M'};
hca = subplot(nRows,nCols,6);
%mms.plot_projection(hca,ePDist1.convertto('1/(cm^2 s sr keV)'),'tint',time,'xyz',[M;N;L],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scPot1,'vlabel',vlabels);
c_eval('dist = ePDist?;',ic)
mms.plot_projection(hca,dist.convertto('s^3/km^6'),'tint',time,'xyz',[L;N;-M],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scPot1,'vlabel',vlabels);
colormap('jet')
hold(hca,'on')
v0 = [vx0,vy0,vz0];
v0LMN = v0*1e-3;
h_v0 = plot3(hca,v0LMN(1),v0LMN(3),-v0LMN(1),'sk'); 
h_v0.MarkerSize = 10;
h_v0.MarkerFaceColor = 0+[1 1 1];
hold(hca,'off')
hca.Title.String = time.utc;


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
EoM = @(ttt,xxx) mms_2015Nov12.eom(ttt,xxx,a,B0,d,E0,eps);
[t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
x = x_sol(:,1);
y = x_sol(:,2);
z = x_sol(:,3);
vx = x_sol(:,4);
vy = x_sol(:,5);
vz = x_sol(:,6);      


% plotting

hca = subplot(nRows,nCols,[4 5]);

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
axis(hca,'equal')
hca.ZLim = [-20 20];
xlabel(hca,'L (km)')
ylabel(hca,'M (km)')
zlabel(hca,'N (km)')
hold off
view([1 0 0])
hca.XGrid = 'on';
hca.YGrid = 'on';
hca.ZGrid = 'on';
%set(hca,'ylim',[-10 100])
%set(hca,'zlim',d*[-1 1]*1e-3)
hca.Title.String = ['Electron test particle trajectory, E_{e0} = ' num2str(electron_energy) ' eV'];

