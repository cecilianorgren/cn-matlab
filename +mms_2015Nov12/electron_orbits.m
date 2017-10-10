% integrate electron orbits for a given magnetic field
if 0
% Define time to run the particle
T = 0.04;
units = irf_units;
electron_energy = 200; % eV
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
E0 = 2*1e-3;

Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(z) -B0*tanh(z*pi/d).*(1-exp(-z.^2*5/(d^2)));
By = @(z) 0.5*B0-3*0.5*B0*sin(2/3*pi/d*z).*exp(-z.^2*2/(d^2)).*(1-exp(-z.^2*5/(d^2)));
Bz = -B0*eps;
%Bz = @(x) b*x/d/15;
Ex = 0;
Ey = 0;
Ez = @(z) -E0*sin(pi/d*(z));

hca = subplot(2,1,1);
zz = linspace(-d,d,20);
plotyy(hca,zz*1e-3,Ez(zz)*1e3,zz*1e-3,[Bx(zz); By(zz); sqrt(Bx(zz).^2 + By(zz).^2)]*1e9)
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
EoM = @(ttt,xxx) mms_2015Nov12.eom(ttt,xxx,a,B0,d,E0,eps);
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

end
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
h_v0 = plot3(hca,v0LMN(1),v0LMN(3),1+0*v0LMN(1),'sk'); 
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
limN = 30e3;
Er = 0*-1e-3; % reconnection electric field, V/m
E0 = 2e-3;
B0 = 10e-9; % asymptotical magnetic field, T
d = 12e3; % thickness of current sheet, m
b = 8e3; % bifurcation length scale, m
Bg = 5*1e-9; % guide field, T
Bn = 1.0*1e-9; % normal field, T
BH = 5e-9;

%Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(x,y,z) x*0 + y*0 + z*0 - B0*tanh(z/d).*(1-exp(-z.^2/(2*b^2)));
By = @(x,y,z) x*0 + y*0 + z*0 + Bg-5*BH*sin(4/3*z/d).*exp(-z.^2/(2*b^2)).*(1-exp(-z.^2/(2*b^2)));
%By = @(x,y,z) x*0 + y*0 + z*0 + Bg-7*BH*tanh(z/d).*exp(-z.^2/d^2).*(1-exp(-z.^2/(2*b^2))); 
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z) x*0 + y*0 + z*0 - E0*sin(z/d);

T = 0.1; % integration time

% Electron test particles
particle_set = 4;
switch particle_set 
  case 1 % random
    electron_energy = 200*rand(1,3); % eV
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
  case 2
    electron_energy = [100 50 120]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = [34 30 -25];
    velocity_angle_L = [35 90 90];
    % Initial positions and velocitites
    x0 = [0 0 0]; 
    y0 = [0 0 0]; % km
    z0 = -[-0 -10 -15];
    %z0 = -[-15 -15 -15];
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 3
    electron_energy = [100 150 200]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = -[-45 30 25];
    % Initial positions and velocitites
    x0 = [0 0 0]; 
    y0 = [0 0 0]; % km
    z0 = -[-0 -10 -15];
    %z0 = -[-15 -15 -15];
    vx0 = [0 0 0]; % km/s
    vy0 = vt.*cosd(velocity_angle);
    vz0 = -vt.*sind(velocity_angle);
   case 4
    limN = 40e3;
    T = 0.3;
    electron_energy = [100 200 300]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = [25 30 35];
    % Initial positions and velocitites
    x0 = [0 0 0]; 
    y0 = [0 0 0]; % km
    z0 = [30 30 30];
    %z0 = -[-15 -15 -15];
    
    velocity_angle_L = -[30 60 30];
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
end


% First plot data and model fit to data
hold(hca,'on')
hca = h_field;
zMod = linspace(-d*1.5,d*1.5,20)*3;
%plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',colors(1,:));
plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',colors(2,:))
plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',colors(3,:))
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
  
nParticles = 3;

% Integration
x_sol_all = [];
%saveParticle = cell(1,nParticles);
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
  h_v0 = plot3(hca,v0LMN(2),v0LMN(3),1+0*v0LMN(1),'sk'); 
  h_v0.MarkerSize = 10;
  h_v0.MarkerFaceColor = colors(iParticle,:);0+[1 1 1];
  hold(hca,'off')
  timeUTC = time.utc;
  hca.Title.String = timeUTC(12:23);

  % Initial positions and velocities                                   
  x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s
  x_init = x_init(:,iParticle); 
  % Integrate trajectory
  stopfunction = @(t,y) eom.lim(t,y,limN);
  options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

  EoM = @(ttt,xxx) eom.general(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
  [t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
  x = x_sol(:,1);
  y = x_sol(:,2);
  z = x_sol(:,3);
  vx = x_sol(:,4);
  vy = x_sol(:,5);
  vz = x_sol(:,6); 
  
  %saveParticle(iParticle) = x_sol;
  x_sol_all = [x_sol_all;x_sol(:,1:3)];
  
  % plotting
  hca = h_NM; hold(hca,'on');  


  plot(hca,y*1e-3,z*1e-3,'color',colors(iParticle,:))
  hold(hca,'on')
  plot(hca,y(1)*1e-3,z(1)*1e-3,'go',...
           y(end)*1e-3,z(end)*1e-3,'rx') % plot in km's
  hold(hca,'off')
  %axis(hca,'equal')
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
set(h_NM,'colororder',[0 0 0; colors(1:3,:); 0 0 0]);
irf_legend(h_NM,{'E = [',sprintf('%.0f ',electron_energy(1)),sprintf('%.0f ',electron_energy(2)),sprintf('%.0f',electron_energy(3)),'] eV'},[0.99 0.99],'fontsize',14)
irf_legend(h_NM,{'\theta_L = [',sprintf('%.0f ',velocity_angle_L(1)),sprintf('%.0f ',velocity_angle_L(2)),sprintf('%.0f',velocity_angle_L(3)),'] ^o'},[0.99 0.84],'fontsize',14)
%set(h_NL,'colororder',[0 0 0; colors(1:3,:); 0 0 0]);
%irf_legend(h_NL,{'E = [',sprintf('%.0f ',electron_energy(1)),sprintf('%.0f ',electron_energy(2)),sprintf('%.0f',electron_energy(3)),'] eV'},[0.05 0.6],'fontsize',14)
%irf_legend(h_NL,{'\theta_L = [',sprintf('%.0f ',velocity_angle_L(1)),sprintf('%.0f ',velocity_angle_L(2)),sprintf('%.0f',velocity_angle_L(3)),'] ^o'},[0.05 0.4],'fontsize',14)


%% Make figure for pitchangle and stuff
%x_sol = saveParticles{1};
x = x_sol(:,1);
y = x_sol(:,2);
z = x_sol(:,3);
vx = x_sol(:,4);
vy = x_sol(:,5);
vz = x_sol(:,6); 
  
Bxyz = [Bx(x,y,z),By(x,y,z),Bz(x,y,z)]; normBxyz = irf_norm(Bxyz);
Vxyz = [vx vy vz]; normVxyz = irf_norm(Vxyz);  
pitchangle = acosd(normBxyz(:,1).*normVxyz(:,1) + ...
                   normBxyz(:,2).*normVxyz(:,2) + ...
                   normBxyz(:,3).*normVxyz(:,3));
                 
fig = figure(39);

%% Electron orbits of several particles, not slice distributions but pitchangle plot
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
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;

% Colors
colors = mms_colors('xyz');

% zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
% hca = h_field;%subplot(nRows,nCols,[1 2]);
% set(hca,'colororder',colors)
% hca.ColorOrder = colors;
% linesObs = plot(hca,zObs,[obsB.data],'.');
% linesObs(1).Color = colors(1,:);
% linesObs(2).Color = colors(2,:);
% linesObs(3).Color = colors(3,:);
% set(hca,'colororder',colors)
% irf_legend(hca,{'B_L','B_M','B_N'},[0.95 0.95],'fontsize',14)

% Model parameters
limN = 30e3;
Er = 0*-1e-3; % reconnection electric field, V/m
E0 = 2e-3;
B0 = 10e-9; % asymptotical magnetic field, T
d = 12e3; % thickness of current sheet, m
b = 8e3; % bifurcation length scale, m
Bg = 5*1e-9; % guide field, T
Bn = 1.0*1e-9; % normal field, T
BH = 5e-9;

%Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(x,y,z) x*0 + y*0 + z*0 - B0*tanh(z/d).*(1-exp(-z.^2/(2*b^2)));
By = @(x,y,z) x*0 + y*0 + z*0 + Bg-5*BH*sin(4/3*z/d).*exp(-z.^2/(2*b^2)).*(1-exp(-z.^2/(2*b^2)));
%By = @(x,y,z) x*0 + y*0 + z*0 + Bg-7*BH*tanh(z/d).*exp(-z.^2/d^2).*(1-exp(-z.^2/(2*b^2))); 
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z) x*0 + y*0 + z*0 - E0*sin(z/d);

T = 0.2; % integration time

% Electron test particles
particle_set = 5;
switch particle_set 
  case 1 % random
    electron_energy = 200*rand(1,3); % eV
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
  case 2
    electron_energy = [100 50 120]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = [34 30 -25];
    velocity_angle_L = [35 90 90];
    % Initial positions and velocitites
    x0 = [0 0 0]; 
    y0 = [0 0 0]; % km
    z0 = -[-0 -10 -15];
    %z0 = -[-15 -15 -15];
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 3
    electron_energy = [100 150 200]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = -[-45 30 25];
    % Initial positions and velocitites
    x0 = [0 0 0]; 
    y0 = [0 0 0]; % km
    z0 = -[-0 -10 -15];
    %z0 = -[-15 -15 -15];
    vx0 = [0 0 0]; % km/s
    vy0 = vt.*cosd(velocity_angle);
    vz0 = -vt.*sind(velocity_angle);
  case 4
    limN = 40e3;
    T = 0.5;
    electron_energy = [100 200 300]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = [25 30 35];
    % Initial positions and velocitites
    x0 = [0 0 0]; 
    y0 = [0 0 0]; % km
    z0 = [30 30 30];
    %z0 = -[-15 -15 -15];
    
    velocity_angle_L = -[30 60 30];
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 5
    limN = 40e3;
    T = 0.5;
    electron_energy = 50+[10 40 10 40]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = 30+[15 20 25 10];
    % Initial positions and velocitites
    x0 = [0 0 0 0]; 
    y0 = [0 0 0 0]; % km
    z0 = [30 30 -30 -30];
    %z0 = -[-15 -15 -15];
    
    velocity_angle_L = [-20 -20 20 20];
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 6 % energy and pitch angle
    limN = 40e3;
    T = 0.5;
    electron_energy = [100 100 300 100]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;
    
    % Initial positions and velocitites
    x0 = [0 0 0 0]; 
    y0 = [0 0 0 0]; % km
    z0 = [30 30 30 -30];
    
    pa0 = [10 10 10 10];
    B0 = [Bx(x0,y0,z0); By(x0,y0,z0); Bz(x0,y0,z0)];
    b0 = irf_norm(B0')';
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
end
nParticles = numel(electron_energy);

% Integration
x_sol_all = [];
saveParticle = cell(1,nParticles);
for iParticle= 1:nParticles 

  % Initial positions and velocities                                   
  x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s
  x_init = x_init(:,iParticle); 
  % Integrate trajectory
  stopfunction = @(t,y) eom.lim(t,y,limN);
  options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);

  EoM = @(ttt,xxx) eom.general(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
  [t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
  x = x_sol(:,1);
  y = x_sol(:,2);
  z = x_sol(:,3);
  vx = x_sol(:,4);
  vy = x_sol(:,5);
  vz = x_sol(:,6); 
  
  Bxyz = [Bx(x,y,z),By(x,y,z),Bz(x,y,z)]; normBxyz = irf_norm(Bxyz);
  Vxyz = [vx vy vz]; normVxyz = irf_norm(Vxyz);  
  pitchangle = acosd(normBxyz(:,1).*normVxyz(:,1) + ...
                     normBxyz(:,2).*normVxyz(:,2) + ...
                     normBxyz(:,3).*normVxyz(:,3));
      
  saveParticle{iParticle}.r = x_sol(:,1:3);
  saveParticle{iParticle}.r0 = [x0(iParticle),y0(iParticle),z0(iParticle)];
  saveParticle{iParticle}.v = x_sol(:,4:6);
  saveParticle{iParticle}.v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)];
  saveParticle{iParticle}.B = Bxyz;
  saveParticle{iParticle}.pa = pitchangle;
end

% Plot, 4 panels, incl B
colors = mms_colors('1234');
% Set up plot
nRows = 5;
nCols = 1;
units = irf_units;

h(1) = subplot(nRows,nCols,1); 
h(2) = subplot(nRows,nCols,2); 
h(3) = subplot(nRows,nCols,3); 
h(4) = subplot(nRows,nCols,[4 5]); 

isub = 1;
if 1 % magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,20)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
end
if 1 % Time: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'x',tintObs)
end

if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [10 400];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.CLim = h(2).CLim;  
  xlabel(hca,'N (km)')  
end
for iP = 1:numel(saveParticle) % Electron pitchangles
  % Pick out the data
  z = saveParticle{iP}.r(:,3); % N
  pa = saveParticle{iP}.pa; % pitch angles    
  
  step = 2;
  hold(hca,'on');
  plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
  hold(hca,'off') 
end

for iP = 1:numel(saveParticle) % Electron trajectories
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  
  hca = h(isub);
  
  if 1
    hl(iP) = plot(hca,x*1e-3,z*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,x(1)*1e-3,z(1)*1e-3,'go',...
              x(end)*1e-3,z(end)*1e-3,'rx'); % plot in km's
  ylabel(hca,'N (km)'); hca.YLim = [-30 30];
  xlabel(hca,'L (km)')  
  else
    hl(iP) =  plot(hca,z*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,z(1)*1e-3,x(1)*1e-3,'go',...
             z(end)*1e-3,x(end)*1e-3,'rx') % plot in km's
  ylabel(hca,'L (km)')
  xlabel(hca,'N (km)'); hca.XLim = [-30 30];
  hca.XLim = h(2).XLim;
  end
  if iP == numel(saveParticle); hold(hca,'off'); end  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
legend(hl,labels{:},'location','west')

h(1).Position(2) = h(1).Position(2)-h(1).Position(4);
h(1).Position(4) = h(1).Position(4)*2;
h(2).Position(3) = h(1).Position(3);
h(2).Position(2) = h(1).Position(2)-h(1).Position(4)*0.5;
h(3).Position(2) = h(2).Position(2)-h(2).Position(4);
h(3).Position(3) = h(1).Position(3);

for ii = 1:4;
  h(ii).Position(3) = 0.7;
end

h(2).XTick = [];
h(1).XGrid = 'on';
%h(1).XAxisLocation= 'top';
h(1).XLim = h(3).XLim;
h(1).YLim = [-13 13];
h(1).Title.String = '';
h(4).XLim(1) = 0;
h(4).YGrid = 'on';
h(4).Position(2) = 0.18;

hylab = ylabel(h(2),{'Pitchangle (\circ)'},'interpreter','tex');

hcb = findobj(fig.Children,'type','colorbar');
hcb.Position(1) = h(2).Position(1) + h(2).Position(3)+0.01;
hcb.Position(2) = h(3).Position(2);
hcb.Position(4) = h(2).Position(4)*2;
hcb.FontSize = 10;
hcb.YLabel.FontSize = 10;

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)'};
legshift = 0; % the two sc configuration plots

for ii = 1:4
  irf_legend(h(ii),legends{ii},[0.01 0.9],'color',[0 0 0],'fontsize',14);
  h(ii).FontSize = 14;  
  h(ii).YLabel.FontSize = 14;  
end

set(gcf, 'InvertHardCopy', 'off');
set(gcf,'paperpositionmode','auto');
set(gcf,'color','white');

%% Plot, 3 panels, no B
colors = mms_colors('1234');
% Set up plot
nRows = 4;
nCols = 1;
units = irf_units;

h(1) = subplot(nRows,nCols,1); 
h(2) = subplot(nRows,nCols,2); 
h(3) = subplot(nRows,nCols,[3 4]); 

isub = 1;
if 1 % Time: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  % find magnetic mirror angle
  tref4 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength4 = [-0.5 0];
  tref3 = irf_time('2015-11-12T07:19:21.40Z','utc>epochtt'); tlength3 = [-0.5 0];
  tref2 = irf_time('2015-11-12T07:19:21.50Z','utc>epochtt'); tlength2 = [-0.60 0];
  %tref2 = irf_time('2015-11-12T07:19:21.20Z','utc>epochtt'); tlength2 = [-0.30 0.3];
  tref1 = irf_time('2015-11-12T07:19:20.97Z','utc>epochtt'); tlength1 = [0 0.40];
  tref1 = irf_time('2015-11-12T07:19:21.5Z','utc>epochtt'); tlength1 = [-0.60 0.0];
  c_eval('tref = tref?;',ic)
  c_eval('tintB = tref+tlength?;',ic)
  thetaref = -00;
  c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,thetaref+asind((mvaB?.tlim(tintB).abs.data/mvaB?.abs.resample(tref).data).^0.5));')
  %c_eval('alphaB = irf.ts_scalar(mvaB?.tlim(tintB).time,asind((mvaB?.tlim(tintB).abs.data.^-1*mvaB?.abs.resample(tref).data).^0.5));')
  elim = [10 400];
  try
    c_eval('irf_spectrogram(hca,ePitch?.tlim(tint+[-5 5]).elim(elim).deflux.specrec(''pa''),''log'');',ic)
  catch
    %irf_spectrogram(hca,ePDist4.tlim(tint).pitchangles(dmpaB4,20).elim(elim).deflux.specrec('pa'),'log');
  end
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  irf_plot(hca,alphaB,'k--');
  irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  irf_zoom(hca,'x',tintObs)
end

if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [10 400];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hold(hca,'on')
  set(hca,'ColorOrder',mms_colors('11'))
  %irf_plot(hca,alphaB,'k--');
  %irf_plot(hca,180-alphaB,'k--');  
  hold(hca,'off')
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  hca.YLabel.String = {'\theta_{PA,e}','(\circ)'};
  ylabel(hca,{'\theta_{PA,e}','(\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  hca.CLim = h(1).CLim;  
  xlabel(hca,'N (km)')  
end
for iP = 1:numel(saveParticle) % Electron pitchangles
  % Pick out the data
  z = saveParticle{iP}.r(:,3); % N
  pa = saveParticle{iP}.pa; % pitch angles    
  
  step = 5;
  hold(hca,'on');
  plot(hca,z(1:step:end)*1e-3,pa(1:step:end),'.','color',colors(iP,:))
  hold(hca,'off') 
end

for iP = 1:numel(saveParticle) % Electron trajectories
  % Pick out the data
  x = saveParticle{iP}.r(:,1); % L
  y = saveParticle{iP}.r(:,2); % M  
  z = saveParticle{iP}.r(:,3); % N  
  pa = saveParticle{iP}.pa; % pitch angles  
  
  hca = h(isub);
  
  if 1
    hl(iP) = plot(hca,x*1e-3,z*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,x(1)*1e-3,z(1)*1e-3,'go',...
              x(end)*1e-3,z(end)*1e-3,'rx'); % plot in km's
  ylabel(hca,'N (km)'); hca.YLim = [-30 30];
  xlabel(hca,'L (km)')  
  else
    hl(iP) =  plot(hca,z*1e-3,x*1e-3,'color',colors(iP,:));
    if iP == 1; hold(hca,'on'); end
    plot(hca,z(1)*1e-3,x(1)*1e-3,'go',...
             z(end)*1e-3,x(end)*1e-3,'rx') % plot in km's
  ylabel(hca,'L (km)')
  xlabel(hca,'N (km)'); hca.XLim = [-30 30];
  hca.XLim = h(2).XLim;
  end
  if iP == numel(saveParticle); hold(hca,'off'); end  
  %hca.Title.String = ['Electron test particle trajectory'];%, E_{e0} = ' num2str(electron_energy) ' eV']; 
end
labels = arrayfun(@(x) {['E = ' num2str(x) ' eV']},electron_energy);
legend(hl,labels{:},'location','west')

h(2).Position(3) = h(1).Position(3);
h(2).Position(2) = h(1).Position(2)-h(1).Position(4);

h(3).Position(3) = h(1).Position(3);

%%


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
set(h_NM,'colororder',[0 0 0; colors(1:3,:); 0 0 0]);
irf_legend(h_NM,{'E = [',sprintf('%.0f ',electron_energy(1)),sprintf('%.0f ',electron_energy(2)),sprintf('%.0f',electron_energy(3)),'] eV'},[0.99 0.99],'fontsize',14)
irf_legend(h_NM,{'\theta_L = [',sprintf('%.0f ',velocity_angle_L(1)),sprintf('%.0f ',velocity_angle_L(2)),sprintf('%.0f',velocity_angle_L(3)),'] ^o'},[0.99 0.84],'fontsize',14)
%set(h_NL,'colororder',[0 0 0; colors(1:3,:); 0 0 0]);
%irf_legend(h_NL,{'E = [',sprintf('%.0f ',electron_energy(1)),sprintf('%.0f ',electron_energy(2)),sprintf('%.0f',electron_energy(3)),'] eV'},[0.05 0.6],'fontsize',14)
%irf_legend(h_NL,{'\theta_L = [',sprintf('%.0f ',velocity_angle_L(1)),sprintf('%.0f ',velocity_angle_L(2)),sprintf('%.0f',velocity_angle_L(3)),'] ^o'},[0.05 0.4],'fontsize',14)
