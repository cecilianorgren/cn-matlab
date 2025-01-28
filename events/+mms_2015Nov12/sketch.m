% Make LN plane sketch of electron velocities
ic = 1;
CS_normal_velocity = 70; % km/s
units = irf_units;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
colors = mms_colors('xyz');

% Set up plot
nRows = 1;
nCols = 1;
for ii = 1:nRows*nCols
  h(ii) = subplot(nRows,nCols,ii); 
end
h(1).Position = h(1).Position.*[1.0 1.6 0.74 0.8];
isub = 1;

% Select and prepare data
if 1
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs);'...
'obsVe = mvaVe?.tlim(tintObs);'...
'obsVeperp = mvaVe?perp.tlim(tintObs);'...
'obsBslow = mvaB?.tlim(tintObs).resample(obsVe);'...
],ic)
else
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaJ?par.tlim(tintObs);'...
'obsVe = mvaJ?.tlim(tintObs);'...
'obsVeperp = mvaJ?perp.tlim(tintObs);'...
'obsBslow = mvaB?.tlim(tintObs).resample(obsVe);'...
],ic)
end
obsVepar_ = obsVe-obsVeperp;%par;

normB = obsBslow/obsBslow.abs;
normVe = obsVe/obsVe.abs;
pitchangle = acosd(normB.dot(normVe).data);

% length = t*v_cs
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsVe = (obsVe.time.epochUnix-mean(obsVe.time.epochUnix))*CS_normal_velocity;
nL = 10; limL = 15; obsL = linspace(-limL,limL,nL); % length in L-direction - should scale with B
[obsN,obsL] = meshgrid(zObs,obsL);
obsBM = repmat(obsB.y.data,1,nL)';

if 0 % Plot data in spatial scale
  hca = h(isub); isub = isub + 1;
  set(hca,'colororder',colors)
  hca.ColorOrder = colors;
  linesObs = plot(hca,zObs,[obsB.data],'.');
  linesObs(1).Color = colors(1,:);
  linesObs(2).Color = colors(2,:);
  linesObs(3).Color = colors(3,:);
  set(hca,'colororder',colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.95 0.95],'fontsize',14)
end

colors = mms_colors('matlab');

% Plot BM in LN plane
hca = h(isub); isub = isub + 1;
hp = pcolor(obsN,obsL,obsBM);
shading(hca,'flat')
hcb = colorbar('peer',hca); hcb.YLim = [-2 12];
hcb.YLabel.String = 'B_M (nT)';
hcb.FontSize = 16;
colormap(hca,flipdim(cn.cmap('bluered3'),1)); hca.CLim = 25*[-1 1]; colorBL = mms_colors('matlab'); colorBL = colorBL(1,:)*0.8;
%colormap(hca,cn.cmap('bluered3')); hca.CLim = 25*[-1 1]; colorBL = mms_colors('matlab'); colorBL = colorBL(4,:)*0.8;
%colormap(hca,tocolumn(1:-0.01:0.5)*[1 1 1]); hca.CLim = [0 15]; colorBL = [0.5 0.5 0.5];
%colormap(hca,'hot'); hca.CLim = [-2 12];

% Plot quivers
hold(hca,'on')
hqB = plot_quivers_tmp(hca,[obsB.data(:,3) obsB.data(:,1)],[zObs zObs*0],0,colorBL);
hqBabs = plot(hca,zObs,obsB.abs.data,'linewidth',2,'color',0.8*colors(1,:));
colorVe = mms_colors('matlab'); colorVe = colorVe(2,:); [0.7 0.1 0.5];
colorVe_ = mms_colors('matlab'); colorVe_ = 0*colorVe_(7,:);
vescale = 0.015;
%plot_quivers_tmp(hca,[obsVe.data(:,3) obsVe.data(:,1)]*vescale,[zObsVe zObsVe*0],0,colorVe)
%plot_xo(hca,obsVe.data(:,2),[zObsVe zObsVe*0],colorVe);
hqVpar = plot_quivers_tmp(hca,[obsVepar_.data(:,3) obsVepar_.data(:,1)]*vescale,[zObsVe zObsVe*0],0,colorVe_);
hqVparXO = plot_xo(hca,obsVepar_.data(:,2),[zObsVe zObsVe*0],colorVe_);
hqVperp = plot_quivers_tmp(hca,[obsVeperp.data(:,3) obsVeperp.data(:,1)]*vescale,[zObsVe zObsVe*0],0,colorVe);
hqVperpXO = plot_xo(hca,obsVeperp.data(:,2),[zObsVe zObsVe*0],colorVe);
hold(hca,'off')

% BL axes
hca.YLabel.String = 'B_L (nT)';
hca.YAxisLocation = 'right';
hca.XLabel.String = 'N (km)';

% Make new axis for Ve
ax2 = axes('Position',hca.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
ax2.YLabel.String = 'V_{e} (km/s)';
ax2.YLim = hca.YLim/vescale;
ax2.XTick = [];

hca.FontSize = 18;
ax2.FontSize = 18;
hcb.FontSize = 18;

hleg = legend([hqBabs hqB(1) hqVperp(1) hqVperpXO(23) hqVpar(1) hqVparXO(15)],{'|B|','B_{LN}','V_{e,\perp,LN}','V_{e,\perp,M}','V_{e,||,LN}','V_{e,||,M}'},'location','southwest');

hleg.Position(1) = 0.9;
hleg.Position(2) = 0.4;
hcb.Position(1) = 0.8;
%hca.Position(3) = 0.67;

   set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'paperpositionmode','auto');
    set(gcf,'color','white');
%% Make LN plane sketch of current sheet
ic = 1;
CS_normal_velocity = 70; % km/s
units = irf_units;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
colors = mms_colors('xyz');

% Set up plot
nRows = 1;
nCols = 1;
for ii = 1:nRows*nCols
  h(ii) = subplot(nRows,nCols,ii); 
end
h(1).Position = h(1).Position.*[1.0 1.6 0.74 0.8];
isub = 1;

% Select and prepare data
if 0
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs);'...
'obsVe = mvaVe?.tlim(tintObs);'...
'obsVeperp = mvaVe?perp.tlim(tintObs);'...
'obsBslow = mvaB?.tlim(tintObs).resample(obsVe);'...
],ic)
else
c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp.resample(obsB);'...
'obsVepar = mvaJ?par.tlim(tintObs);'...
'obsVe = mvaJ?.tlim(tintObs);'...
'obsVeperp = mvaJ?perp.tlim(tintObs);'...
'obsBslow = mvaB?.tlim(tintObs).resample(obsVe);'...
],ic)
end
obsVepar_ = obsVe-obsVeperp;%par;

normB = obsBslow/obsBslow.abs;
normVe = obsVe/obsVe.abs;
pitchangle = acosd(normB.dot(normVe).data);

% length = t*v_cs
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsVe = (obsVe.time.epochUnix-mean(obsVe.time.epochUnix))*CS_normal_velocity;
nL = 10; limL = 15; obsL = linspace(-limL,limL,nL); % length in L-direction - should scale with B
[obsN,obsL] = meshgrid(zObs,obsL);
obsBM = repmat(obsB.y.data,1,nL)';

if 0 % Plot data in spatial scale
  hca = h(isub); isub = isub + 1;
  set(hca,'colororder',colors)
  hca.ColorOrder = colors;
  linesObs = plot(hca,zObs,[obsB.data],'.');
  linesObs(1).Color = colors(1,:);
  linesObs(2).Color = colors(2,:);
  linesObs(3).Color = colors(3,:);
  set(hca,'colororder',colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.95 0.95],'fontsize',14)
end

colors = mms_colors('matlab');

% Plot BM in LN plane
hca = h(isub); isub = isub + 1;
hp = pcolor(obsN,obsL,obsBM);
shading(hca,'flat')
hcb = colorbar('peer',hca); hcb.YLim = [-2 12];
hcb.YLabel.String = 'B_M (nT)';
hcb.FontSize = 16;
colormap(hca,flipdim(cn.cmap('bluered3'),1)); hca.CLim = 25*[-1 1]; colorBL = mms_colors('matlab'); colorBL = colorBL(1,:)*0.8;
%colormap(hca,cn.cmap('bluered3')); hca.CLim = 25*[-1 1]; colorBL = mms_colors('matlab'); colorBL = colorBL(4,:)*0.8;
%colormap(hca,tocolumn(1:-0.01:0.5)*[1 1 1]); hca.CLim = [0 15]; colorBL = [0.5 0.5 0.5];
%colormap(hca,'hot'); hca.CLim = [-2 12];
%hcb.Position(1) = 0.9;
%hca.Position(3) = 0.67;

% Plot quivers
hold(hca,'on')
hqB = plot_quivers_tmp(hca,[obsB.data(:,3) obsB.data(:,1)],[zObs zObs*0],0,colorBL);
hqBabs = plot(hca,zObs,obsB.abs.data,'linewidth',2,'color',0.8*colors(1,:));
colorVe = mms_colors('matlab'); colorVe = colorVe(2,:); [0.7 0.1 0.5];
colorVe_ = mms_colors('matlab'); colorVe_ = 0*colorVe_(7,:);
vescale = 0.015;
%plot_quivers_tmp(hca,[obsVe.data(:,3) obsVe.data(:,1)]*vescale,[zObsVe zObsVe*0],0,colorVe)
%plot_xo(hca,obsVe.data(:,2),[zObsVe zObsVe*0],colorVe);
hqVpar = plot_quivers_tmp(hca,[obsVepar_.data(:,3) obsVepar_.data(:,1)]*vescale,[zObsVe zObsVe*0],0,colorVe_);
hqVparXO = plot_xo(hca,obsVepar_.data(:,2),[zObsVe zObsVe*0],colorVe_);
hqVperp = plot_quivers_tmp(hca,[obsVeperp.data(:,3) obsVeperp.data(:,1)]*vescale,[zObsVe zObsVe*0],0,colorVe);
hqVperpXO = plot_xo(hca,obsVeperp.data(:,2),[zObsVe zObsVe*0],colorVe);
hold(hca,'off')

% BL axes
hca.YLabel.String = 'B_L (nT)';
hca.YAxisLocation = 'right';
hca.XLabel.String = 'N (km)';

% Make new axis for Ve
ax2 = axes('Position',hca.Position,...
    'XAxisLocation','top',...
    'YAxisLocation','left',...
    'Color','none');
ax2.YLabel.String = 'J (nA/m^2)';
ax2.YLim = hca.YLim/vescale;
ax2.XTick = [];

hca.FontSize = 18;
ax2.FontSize = 18;
hcb.FontSize = 18;

hleg = legend([hqBabs hqB(1) hqVperp(1) hqVperpXO(23) hqVpar(1) hqVparXO(15)],{'|B|','B_{LN}','J_{\perp,LN}','J_{\perp,M}','J_{||,LN}','J_{||,M}'},'location','southwest');

hleg.Position(1) = 0.9;
hleg.Position(2) = 0.4;
hcb.Position(1) = 0.8;
%hca.Position(3) = 0.67;

   set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'paperpositionmode','auto');
    set(gcf,'color','white');
    
%% Plot simultaneously the electron orbits of several particles
% Set up plot
nRows = 2;
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
d = 28e3; % d = thickness of current sheet, eg 800 km.
eps = 0.05; % Normal magnetic field ratio Bn = B0*eps;
E0 = 2*1e-3;

%Bx = @(z) -B0*tanh(z*pi/d);
Bx = @(z) -B0*tanh(z*pi/d).*(1-exp(-z.^2*5/(d^2)))-0*B0/6;
By = @(z) 1*0.5*B0-3*0.5*B0*sin(4/3*pi/d*z).*exp(-z.^2*2/(d^2)).*(1-exp(-z.^2*5/(d^2)));
Bz = @(z) B0*eps+z*0;
Ex = 0;
Ey = 0;
Ez = @(z) -E0*sin(pi/d*(z));

T = 0.3; % integration time

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
zMod = linspace(-d*1.5,d*1.5,20);
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

   set(gcf, 'InvertHardCopy', 'off');
    set(gcf,'paperpositionmode','auto');
    set(gcf,'color','white');