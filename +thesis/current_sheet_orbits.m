% integrate electron orbits for a given magnetic field
% Set up plot
nCols = 1;
nRows = 4;

h(1) = subplot(nRows,nCols,1);
h(2) = subplot(nRows,nCols,[2 3]);
h(3) = subplot(nRows,nCols,4);


% Harris sheet without guide field

% Define time to run the particle
T = 0.1;
units = irf_units;
electron_energy = 200; % eV
vt = sqrt(electron_energy*units.eV*2/units.me)/1000; 
velocity_angle = 90;
m = units.me;
q = -units.e;

% Initial positions and velocitites
x0 = 0*1e3; % m
y0 = 0*1e3; 
z0 = 15*1e3;
vx0 = 0*1e3; % m/s
vy0 = vt*cosd(velocity_angle)*1e3;
vz0 = -vt*sind(velocity_angle)*1e3;

% Parameters
Er = 0e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 5*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 0;

% Harris current sheet
Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z) x*0 + y*0 + z*0;

% Generalized potential
% must rescale, matlab can handle the small numbers
Ay = @(z) B0*d*log(cosh(z/d));
Az = @(z) -Bg*z;
pz = m*vz0 + q*Az(z0);
py = m*vy0 + q*Ay(z0);
Gamma = @(z,py,pz,phi0,d) 0.5/m*((py-q*Ay(z)).^2+(pz-q*Az(z)).^2);% + q*phi0;
Gamma0 = Gamma(0,py,pz,phi0,d);

zz = linspace(-5*d,5*d,100);
% Plot magnetic field
hca = h(1);
plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9)
hca.XLabel.String = 'z';

% Plot generalized potential
hca = h(3);
plotyy(hca,zz*1e-3,[Gamma(zz,py,pz,phi0,d); zz*0+pz]',...
           zz*1e-3,[py-q*Ay(zz)]')
%hca.YLim = pz*[-2 2]; 
hca.XLabel.String = 'z';

% Integrate particle orbit
x_init = [x0;y0;z0;vx0;vy0;vz0]; % *m, m/s
EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
[t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6);      

hca = h(1,2);
xlim = max(abs(x))*1e-3;
ylim = max(abs(y))*1e-3;
zlim = max(abs(z))*1e-3;
xpl = xlim*[-1 -1 1 1];
ypl = ylim*[1 -1 -1 1];
zpl = 0*[-1 -1 1 1];

if 0 % 3D
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
else % 2D
  plot(hca,z*1e-3,y*1e-3,'k',...
        z(1)*1e-3,y(1)*1e-3,'go',...
        z(end)*1e-3,y(end)*1e-3,'rx') % plot in km's
  xlabel(hca,'y')
  ylabel(hca,'z') 
end

zlim = 30;
for iPanel = 1:3  
  h(1,iPanel).XLim = zlim*[-1 1];
end

%% integrate electron orbits for a given magnetic field: several particles
% Set up p  lot
fig = figure(33);
fig.Position(4) = 1e3*0.5;
fig.Position(3) = 1e3*0.3;
nCols = 1;
nRows = 4;

clear h
h(1) = subplot(nRows,nCols,1);
h(2) = subplot(nRows,nCols,[2:3]);
h(3) = subplot(nRows,nCols,4);

units = irf_units;

% Parameters
Er = 0e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 2*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 0;

% Harris current sheet
Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z) x*0 + y*0 + z*0;
Ay = @(z) B0*d*log(cosh(z/d));
Az = @(z) -Bg*z;

Gamma = @(z,py,pz,phi0,d) 0.5/m*((py-q*Ay(z)).^2+(pz-q*Az(z)).^2);

% Plot magnetic field
hca = h(1);
zz = linspace(-5*d,5*d,100);
plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9,'k');
hca.YTick = 0;
hca.XLabel.String = 'z';
if Bg == 0
  hca.Title.String = 'No guide field';
else
  hca.Title.String = 'Finite guide field';
end
hca.Title.Position(2) = B0*1e9 + 2;
  
%irf_legend(hca,'B',[0.01 1.05],'fontsize',14)
%irf_legend(hca,'B',[0.01 1.05],'fontsize',14)

% Initiate particles
particles = {};
iParticle = 0;
if 1 % no drift
  iParticle = iParticle + 1;
  particle.T = 0.02;
  particle.energy = 50;
  particle.z0 = -15; 
  particle.y0 = 0;
  particle.velocity_angle = 00;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end
if 1 % meandering
  iParticle = iParticle + 1;
  particle.T = 0.01;
  particle.energy = 150;
  particle.z0 = 5;  
  particle.y0 = 10;
  particle.velocity_angle = 70;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end
if 1 % meandering 8's
  iParticle = iParticle + 1;
  particle.T = 0.02;
  particle.energy = 70;
  particle.z0 = -1;  
  particle.y0 = -10;
  particle.velocity_angle = 160;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end
if 1 % grad B
  iParticle = iParticle + 1;
  particle.T = 0.02;
  particle.energy = 80;
  particle.y0 = 0;
  particle.z0 = 8;  
  particle.velocity_angle = -70;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end

for iParticle = 1:numel(particles)
  % Load particle data
  particle = particles{iParticle};
  electron_energy = particle.energy;
  velocity_angle = particle.velocity_angle;
  q = particle.q;
  m = particle.m;
  T = particle.T;
  
  vt = sqrt(electron_energy*units.eV*2/m)/1000;
  % Initial positions and velocitites
  x0 = 0*1e3; % m
  y0 = particle.y0*1e3; 
  z0 = particle.z0*1e3;
  vx0 = 0*1e3; % m/s
  vy0 = vt*cosd(velocity_angle)*1e3;
  vz0 = -vt*sind(velocity_angle)*1e3;

  % Generalized potential
  pz = m*vz0 + q*Az(z0);
  py = m*vy0 + q*Ay(z0);
  Gamma0 = Gamma(0,py,pz,phi0,d);
  
  % Integrate particle orbit
  x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
  EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
  [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
  x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
  vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6);      

  hca = h(2);
  hold(hca,'on')
  horb = plot(hca,z*1e-3,y*1e-3,...
                  z(1)*1e-3,y(1)*1e-3,'go',...
                  z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
  xlabel(hca,'z')
  ylabel(hca,'y') 
  hold(hca,'off')
  
  % Plot generalized potential
  if 1
    hca = h(3);
    hold(hca,'on')
    vecGamma = Gamma(zz,py,pz,phi0,d);
    hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    limGamma = vecGamma; 
    zzGamma = zz;
    zzGamma(limGamma>electron_energy*units.eV) = NaN;
    limGamma(limGamma>electron_energy*units.eV) = NaN;
    hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    %hpot.Color = ;
    hca.XLabel.String = 'z';
    hold(hca,'off')
  end
end
hold(h(2),'off')
%axis(h(2),'equal')
zlim = 16;
for iPanel = 1:numel(h) 
  h(iPanel).XLim = zlim*[-1 1];
end
for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  %h(ii).Box = 'off';
  h(ii).FontSize = 14;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
hold(h(2),'on')
plot(h(2),[0 0],h(2).YLim,'k')
h(2).XLim = zlim*[-1 1];
hold(h(2),'off')
h(3).XLim = zlim*[-1 1];
%h(3).YLim = [0 2]*1e-16;
h(3).YLim = [0 500];
h(3).YScale = 'lin';
%% Plot one electron distribution and associated orbit
nRows = 2;
nCols = 3;
units = irf_units;

tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
tintObs = tintObs;
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
Bx = @(z) -B0*tanh(z*pi/d).*(1-exp(-z.^2*5/(d^2)))-0*B0/6;
By = @(z) 1*0.5*B0-3*0.5*B0*sin(4/3*pi/d*z).*exp(-z.^2*2/(d^2)).*(1-exp(-z.^2*5/(d^2)));
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
T = 0.03;
electron_energy = 150; % eV
vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

velocity_angle = -30;
% Initial positions and velocitites
x0 = 0; 
y0 = 0; % km
z0 = 15;
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
hca.ZLim = [-25 25];
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

%% Plot simultaneously the electron orbits of several particles
% Set up plot
nRows = 3;
nCols = 3;
units = irf_units;

for ii = 1:nRows*nCols
  h(ii) = subplot(nRows,nCols,ii); 
end

xfactor = 0.9; 
h_field = subplot(nRows,nCols,[1 2]); 
h_NM = subplot(nRows,nCols,[4 5]); h_NM.Position(3) = h_NM.Position(3)*xfactor;
h_NL = subplot(nRows,nCols,[7 8]); h_NL.Position(3) = h_NL.Position(3)*xfactor;
h_psd(1) = subplot(nRows,nCols,3);
h_psd(2) = subplot(nRows,nCols,6);
h_psd(3) = subplot(nRows,nCols,9);

tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
tintObs = tintObs;
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


% Colors
colors = mms_colors('xyz');

zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
hca = h_field;%subplot(nRows,nCols,[1 2]);
set(hca,'colororder',colors)
hca.ColorOrder = colors;
linesObs = plot(hca,zObs,[obsB.data],'.');
linesObs(1).Color = colors(1,:);
linesObs(2).Color = colors(2,:);
linesObs(3).Color = colors(3,:);
set(hca,'colororder',colors)
irf_legend(hca,{'B_L','B_M','B_N'},[0.95 0.95],'fontsize',14)

% Model parameters
a = 1e-3; % a = E0, eg 1 mV/m
B0 = 10e-9; % b = B0, eg 20 nT
d = 9e3; % d = thickness of current sheet, eg 800 km.
b = 12e3;

eps = 0.05; % Normal magnetic field ratio Bn = B0*eps;
E0 = 2*1e-3;

%Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(z) -B0*tanh(z/d).*(1-exp(-z.^2/(b^2)))-0*B0/15;
By = @(z) 0.5*B0-2*B0*sin(4/3*z/d).*exp(-z.^2/(b^2)).*(1-exp(-z.^2/(b^2)));
By = @(z) 0.5*B0-2*B0*sin(4/3*z/d).*(exp(-z.^2/(b^2))-exp(-2*z.^2/(b^2)));
Bz = @(z) B0*eps+z*0;
Ex = 0;
Ey = 0;
Ez = @(z) -E0*sin(pi/d*(z));

T = 0.1; % integration time

% Electron test particles
if 0
  electron_energy = 300*rand(1,3); % eV
  vt = -sqrt(electron_energy*units.eV*2/units.me)/1000;

  velocity_angle = 90*(2*rand(1,3) - 1);
  velocity_angle_L= 90+10*(2*rand(1,3) - 1);
  % Initial positions and velocitites
  x0 = [0 0 0]; 
  y0 = [0 0 0]; % km
  z0 = 15*(2*rand(1,3) - 1);
  vx0 = vt.*cosd(velocity_angle_L); % km/s
  vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
  vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
else
  electron_energy = [150 150 150]+50; % eV
  vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

  velocity_angle = -[-45 0 45];
  % Initial positions and velocitites
  x0 = [0 0 0]; 
  y0 = [0 0 0]; % km
  z0 = -[-5 -10 -15];
  z0 = -[-15 -15 -15];
  vx0 = [0 0 0]; % km/s
  vy0 = vt.*cosd(velocity_angle);
  vz0 = -vt.*sind(velocity_angle);
end


% First plot data and model fit to data
hold(hca,'on')
hca = h_field;
zMod = linspace(-d*1.5,d*1.5,20)*3;
%plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
lineMod = plot(hca,zMod*1e-3,[Bx(zMod)]*1e9,'--','color',colors(1,:));
plot(hca,zMod*1e-3,[By(zMod)]*1e9,'--','color',colors(2,:))
plot(hca,zMod*1e-3,[Bz(zMod)]*1e9,'--','color',colors(3,:))
hold(hca,'off')
%legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
hca.Title.String = 'Magnetic field';
hca.YGrid = 'on';
hca.YLabel.String = 'B (nT)';
hca.XLabel.String = 'N (km)';

colors = [0.9290    0.6940    0.1250;...
    0.4940    0.1840    0.5560;
    0.3010    0.7450    0.9330;...
    0.6350    0.0780    0.1840]; 
  %%
nParticles = 3;




% Integration
x_sol_all = [];

for iParticle= 1:nParticles 
  %pause
  
  % Plot the electron distribution in the center of the current sheet, and
  % where the test particle lies in that distribution.
  %time = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt');
  time = tintObs(1) + 0.5*(tintObs.stop-tintObs.start)+z0(iParticle)/CS_normal_velocity;  
  hold(h_field,'on')
  h_estart(iParticle) = plot(h_field,z0(iParticle)*[1 1],h_field.YLim,'-.','color',colors(iParticle,:));
  hold(h_field,'off')

  vlim = 12*1e3;
  elevlim = 15;
  strCMap = 'jet';
  projclim = [0.5 4.5];  
  %projclim = [4 7];  
  vlabels = {'v_M','v_N','v_L'};
  hca = h_psd(iParticle);
  %mms.plot_projection(hca,ePDist1.convertto('1/(cm^2 s sr keV)'),'tint',time,'xyz',[M;N;L],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scPot1,'vlabel',vlabels);
  mms.plot_projection(hca,ePDist1.convertto('s^3/km^6'),'tint',time,'xyz',[M;N;L],'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scPot1,'vlabel',vlabels);
  colormap('jet')
  hold(hca,'on')
  v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)];
  v0LMN = v0*1e-3;
  h_v0 = plot3(hca,v0LMN(2),v0LMN(3),v0LMN(1),'sk'); 
  h_v0.MarkerSize = 10;
  h_v0.MarkerFaceColor = colors(iParticle,:);0+[1 1 1];
  hold(hca,'off')
  timeUTC = time.utc;
  hca.Title.String = timeUTC(12:23);


  stopfunction = @(t,y) ExB.events(t,y,L);
  options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

  % Initial positions and velocities                                   
  x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s
  x_init = x_init(:,iParticle);
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

  x_sol_all = [x_sol_all;x_sol(:,1:3)];
  
  % plotting
  hca = h_NM; hold(hca,'on');  

  if 0
  % to print with opacity, figure renderer must be set to OpenGl
  patch(xpl,ypl,zpl);
  hp = findobj(gcf,'type','patch');
  set(hp,'facealpha',0.05)
  hold on
  end

  plot(hca,y*1e-3,z*1e-3,'color',colors(iParticle,:))
  hold(hca,'on')
  plot(hca,y(1)*1e-3,z(1)*1e-3,'go',...
           y(end)*1e-3,z(end)*1e-3,'rx') % plot in km's
  hold(hca,'off')
  axis(hca,'equal')
  hca.YLim = [-25 25];
  ylabel(hca,'N (km)')
  xlabel(hca,'M (km)')  
  hold off
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %set(hca,'ylim',[-10 100])
  %set(hca,'zlim',d*[-1 1]*1e-3)
  hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV'];
  hold(hca,'on');  
  
  
  hca = h_NL; hold(hca,'on');

  if 0
  % to print with opacity, figure renderer must be set to OpenGl
  patch(xpl,ypl,zpl);
  hp = findobj(gcf,'type','patch');
  set(hp,'facealpha',0.05)
  hold on
  end

  plot(hca,x*1e-3,z*1e-3,'color',colors(iParticle,:))
  hold(hca,'on')
  plot(hca,x(1)*1e-3,z(1)*1e-3,'go',...
           x(end)*1e-3,z(end)*1e-3,'rx') % plot in km's
  hold(hca,'off')
  %axis(hca,'equal')
  hca.YLim = [-25 25];
  ylabel(hca,'N (km)')
  xlabel(hca,'L (km)')  
  hold off
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  %set(hca,'ylim',[-10 100])
  %set(hca,'zlim',d*[-1 1]*1e-3
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV'];
 
  
  
  hold(hca,'on');
end

Llim = [min(x_sol_all(:,1)) max(x_sol_all(:,1))]*1e-3;
Mlim = [min(x_sol_all(:,2)) max(x_sol_all(:,2))]*1e-3;
Nlim = [min(x_sol_all(:,3)) max(x_sol_all(:,3))]*1e-3;
if Nlim(1)>-10, NLim(1)=-10; end

limfactor = 1.2;
h_NM.XLim = Mlim*limfactor;
h_NM.YLim = Nlim*limfactor;
h_NL.XLim = Llim*limfactor;
h_NL.YLim = Nlim*limfactor;

h_NL.Box = 'on';
h_NM.Box = 'on';
h_field.Position(3) = h_field.Position(3)*xfactor;
h_field.XTick = [-50:10:50];
h_field.YLim= [-13 13];

h_field.XLim = [-40 40];
h_axes = findobj(gcf,'type','axes');
for iax = 1:numel(h_axes)
  h_axes(iax).FontSize = 14;
end
h_colorbar = findobj(gcf,'type','ColorBar');
for icb = 1:numel(h_colorbar)
  h_colorbar(icb).Position(1) = h_colorbar(icb).Position(1)+0.07;
end
legend(h_field,[linesObs(1) lineMod h_estart],{'Observed data','Model fit','Particle starting position'},'location','best','fontsize',12);
set(h_NL,'colororder',[0 0 0; colors(1:3,:); 0 0 0]);
irf_legend(h_NL,{'E = [',sprintf('%.0f ',electron_energy(1)),sprintf('%.0f ',electron_energy(2)),sprintf('%.0f',electron_energy(3)),'] eV'},[0.05 1.1],'fontsize',14)