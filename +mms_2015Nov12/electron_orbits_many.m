%% Electron orbits of several particles, not slice distributions but pitchangle plot
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
tintObs = tintObs;
CS_normal_velocity = 70; % km/s

ic = 1;

c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar = obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp = obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar = obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp = obsVeperp.resample(obsB);'...
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
% limN = 30e3;
% Er = 0*-1e-3; % reconnection electric field, V/m
% E0 = -2e-3;
% B0 = 10e-9; % asymptotical magnetic field, T
% d = 12e3; % thickness of current sheet, m
% b = 8e3; % bifurcation length scale, m
% Bg = 5*1e-9; % guide field, T
% Bn = 1.0*1e-9; % normal field, T
% BH = 5e-9;
% 
% %Bx = @(z) -B0*tanh(z*pi/d);
% Bx = @(x,y,z) x*0 + y*0 - 1*abs(z)/d*0.05*B0 - B0*tanh(z/d).*(1-exp(-z.^2/(2*b^2)));
% By = @(x,y,z) x*0 + y*0 + z*0 + Bg-5*BH*sin(4/3*z/d).*exp(-z.^2/(2*b^2)).*(1-exp(-z.^2/(2*b^2)));
% %By = @(x,y,z) x*0 + y*0 + z*0 + Bg-7*BH*tanh(z/d).*exp(-z.^2/d^2).*(1-exp(-z.^2/(2*b^2))); 
% Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
% Ex = @(x,y,z) x*0 + y*0 + z*0;
% Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
% Ez = @(x,y,z) x*0 + y*0 + z*0 - E0*sin(z/d).*exp(-z.^2/(2*b^2)).*(1-exp(-z.^2/(2*b^2)));

mms_2015Nov12.Bmodel;

T = 0.02; % integration time

% Electron test particles
%particle_set = 62;
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
    electron_energy = [50 120 50 120]; % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = [-10 20 -10 20];
    % Initial positions and velocitites
    x0 = [0 0 0 0]; 
    y0 = [50 50 -50 -50]; % km
    z0 = [30 30 -30 -30];
    %z0 = -[-15 -15 -15];
    
    velocity_angle_L = [-20 -10 20 10];
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 52 % starting from +N
    limN = 40e3;
    T = 0.5;
    nP = 50;
    electron_energy = 50+20*randn(1,nP); % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = 360*rand(1,nP); % NM-plane angle
    % Initial positions and velocitites
    x0 = zeros(1,nP);
    y0 = zeros(1,nP); % km
    z0 = zeros(1,nP)+30;
    
    velocity_angle_L = 00 + 90*randm(1,nP);
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 53 % starting from -N
    %%
    limN = 40e3;
    T = 0.5;
    nP = 5;
    electron_energy = 50+20*randn(1,nP); % eV
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000;

    velocity_angle = 360*rand(1,nP); % NM-plane angle
    % Initial positions and velocitites
    x0 = zeros(1,nP);
    y0 = zeros(1,nP); % km
    z0 = zeros(1,nP)-30;
    
    velocity_angle_L = 00 - 90*rand(1,nP);
    vx0 = vt.*cosd(velocity_angle_L); % km/s
    vy0 = vt.*cosd(velocity_angle).*sind(velocity_angle_L);
    vz0 = -vt.*sind(velocity_angle).*sind(velocity_angle_L);
  case 61 % energy and pitch angle, % starting from -N
    %%
    nP = 40;
    limN = 40e3;
    T = 0.5;
    electron_energy = 40 + 50*randn(1,nP); % eV
    electron_energy(electron_energy<0) = [];
    nP = numel(electron_energy);    
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000; % km/s
    
    % Initial positions and velocities
    x0 = zeros(1,nP); % km
    y0 = zeros(1,nP); 
    z0 = zeros(1,nP)-30;
    
    pa0 = 0 + 20*randn(1,nP);
    az0 = 360*rand(1,nP);
    B0 = [Bx(x0*1e3,y0*1e3,z0*1e3); By(x0*1e3,y0*1e3,z0*1e3); Bz(x0*1e3,y0*1e3,z0*1e3)]; B0 = B0(:,1)'; % same for all particles since they have same startin position
    b0 = irf_norm(B0);
    b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
    b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);
    thetaBL = acosd(b0*[1 0 0]');
    
    % B coordinate system
    vx0_pa = vt.*cosd(pa0); % km/s
    vy0_pa = vt.*cosd(az0).*sind(pa0);
    vz0_pa = vt.*sind(az0).*sind(pa0);
    
    v0_pa = [vx0_pa' vy0_pa' vz0_pa']';
    
    rI = [1 0 0; 0 1 0; 0 0 1];
    rA = [b0; b0_perp1; b0_perp2];
    rB = [1 0 0; 0 1 0; 0 0 1];
    
    % A>B, v' = B(A^-1)v
    
    % LMN (xyz) coordinate system
    v0 = rB*rA^-1*v0_pa;
    vx0 = v0(1,:);
    vy0 = v0(2,:);
    vz0 = v0(3,:);
  case 62 % energy and pitch angle, % starting from +N
    nP = 40;
    limN = 40e3;
    T = 0.5;
    electron_energy = 100 + 30*randn(1,nP); % eV
    electron_energy(electron_energy<0) = [];
    nP = numel(electron_energy);    
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000; % km/s
    
    % Initial positions and velocities
    x0 = zeros(1,nP); % km
    y0 = zeros(1,nP); 
    z0 = zeros(1,nP)+30;
    
    pa0 = 180 + 20*randn(1,nP);
    az0 = 360*rand(1,nP);
    B0 = [Bx(x0*1e3,y0*1e3,z0*1e3); By(x0*1e3,y0*1e3,z0*1e3); Bz(x0*1e3,y0*1e3,z0*1e3)]; B0 = B0(:,1)'; % same for all particles since they have same startin position
    b0 = irf_norm(B0);
    b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
    b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);
    thetaBL = acosd(b0*[1 0 0]');
    
    % B coordinate system
    vx0_pa = vt.*cosd(pa0); % km/s
    vy0_pa = vt.*cosd(az0).*sind(pa0);
    vz0_pa = vt.*sind(az0).*sind(pa0);
    
    v0_pa = [vx0_pa' vy0_pa' vz0_pa']';
    
    rI = [1 0 0; 0 1 0; 0 0 1];
    rA = [b0; b0_perp1; b0_perp2];
    rB = [1 0 0; 0 1 0; 0 0 1];
    
    % A>B, v' = B(A^-1)v
    
    % LMN (xyz) coordinate system
    v0 = rB*rA^-1*v0_pa;
    vx0 = v0(1,:);
    vy0 = v0(2,:);
    vz0 = v0(3,:);
  case 7 % crescents, % starting from +N
    %%
    nP = 20;
    limN = 40e3;
    T = 0.5;
    electron_energy = 200 + 50*randn(1,nP); % eV
    electron_energy(electron_energy<0) = [];
    nP = numel(electron_energy);    
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000; % km/s
    
    % Initial positions and velocities
    x0 = zeros(1,nP); % km
    y0 = zeros(1,nP); 
    z0 = 10+2*randn(1,nP);
    
    pa0 = 90 + 10*randn(1,nP);
    az0 = 0+90*rand(1,nP);
    B0 = [Bx(x0*1e3,y0*1e3,z0*1e3); By(x0*1e3,y0*1e3,z0*1e3); Bz(x0*1e3,y0*1e3,z0*1e3)]; B0 = B0(:,1)'; % same for all particles since they have same startin position
    b0 = irf_norm(B0);
    b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
    b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);
    thetaBL = acosd(b0*[1 0 0]');
    
    % B coordinate system
    vx0_pa = vt.*cosd(pa0); % km/s
    vy0_pa = vt.*cosd(az0).*sind(pa0);
    vz0_pa = vt.*sind(az0).*sind(pa0);
    
    v0_pa = [vx0_pa' vy0_pa' vz0_pa']';
    
    rI = [1 0 0; 0 1 0; 0 0 1];
    rA = [b0; b0_perp1; b0_perp2];
    rB = [1 0 0; 0 1 0; 0 0 1];
    
    % A>B, v' = B(A^-1)v
    
    % LMN (xyz) coordinate system
    v0 = rB*rA^-1*v0_pa;
    vx0 = v0(1,:);
    vy0 = v0(2,:);
    vz0 = v0(3,:);  
end
limN = 32e3;
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
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
  [t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
  if 1 % interpolate to even time steps    
    tlim = 1e-8;
    dtcut = 10*(t(end)-t(end-1));
    [x_sol, t] = resample(x_sol, t, 2/mean(diff(t))); 
    x_sol = x_sol(t<(t(end)-dtcut),:); t = t(t<(t(end)-dtcut));
    x_sol = x_sol(t>(t(1)+dtcut),:); t = t(t>(t(1)+dtcut));
  end
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
      
  saveParticle{iParticle}.t = t;
  saveParticle{iParticle}.r = x_sol(:,1:3);
  saveParticle{iParticle}.r0 = [x0(iParticle),y0(iParticle),z0(iParticle)];
  saveParticle{iParticle}.v = x_sol(:,4:6);
  saveParticle{iParticle}.v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)];
  saveParticle{iParticle}.B = Bxyz;
  saveParticle{iParticle}.pa = pitchangle;
end

% Bin particles to get 'probability density' of pitchangles
all_z = [];
all_pa = [];
all_pa0 = [];
all_E = [];
all_T = [];
for iP = 1:numel(saveParticle) % Electron pitchangles
  % Pick out the data
  ind = 1:numel(saveParticle{iP}.pa);500; 
  z = saveParticle{iP}.r(ind,3); % N
  pa = saveParticle{iP}.pa(ind); % pitch angle
  v = saveParticle{iP}.v(ind,:);
  E = v(:,1).^2 + v(:,2).^2 + v(:,3).^2;
  all_z = [all_z; z];
  all_pa = [all_pa; pa];  
  all_pa0 = [all_pa0 pa(1)];
  all_E = [all_E; E];
  all_T = [all_T; saveParticle{iP}.t(end)-saveParticle{iP}.t(1)];
end
edges_pa = 0:10:180;
edges_z = [-30:1:30];
[EDGES_PA,EDGES_Z] = meshgrid(edges_pa,edges_z);
[nt,edges,mid,loc] = histcn([all_z all_pa],edges_z*1e3,edges_pa);

if 0 % plot
  %%
surf(edges_z,edges_pa,EDGES_PA',(nt'))
colorbar
shading('flat')
view([0 0 1])
%set(gca,'clim',[0 2000])
elseif 0
  %%
  hca = subplot(2,1,1);
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt'))
  colorbar('peer',hca)
  shading(hca,'flat')
  view([0 0 1])
  
  hca = subplot(2,1,2);
  for ii = 1:numel(saveParticle)
    if ii == 1; hold(hca,'on'); end
    plot(hca,saveParticle{ii}.r(:,1)*1e-3,saveParticle{ii}.r(:,3)*1e-3)
  end
  hold(hca,'off')
elseif 0
  %%
  hca = subplot(2,1,1);
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt'))
  colorbar('peer',hca)
  shading(hca,'flat')
  view([0 0 1])
  
  hca = subplot(2,1,2);
  for ii = 1:numel(saveParticle)
    if ii == 1; hold(hca,'on'); end
    plot3(hca,saveParticle{ii}.r(:,1)*1e-3,saveParticle{ii}.r(:,2)*1e-3,saveParticle{ii}.r(:,3)*1e-3)
  end
  hold(hca,'off')  
elseif 0
  %%
  nPanels = 4;
  hca = subplot(nPanels,1,1);
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt'))
  colorbar('peer',hca)
  shading(hca,'flat')
  view([0 0 1])
  
  hca = subplot(nPanels,1,2);
  for ii = 1:numel(saveParticle)
    if ii == 1; hold(hca,'on'); end
    plot(hca,saveParticle{ii}.r(:,1)*1e-3,saveParticle{ii}.r(:,3)*1e-3)
  end
  hold(hca,'off')  
  hca = subplot(nPanels,1,3);
  for ii = 1:numel(saveParticle)
    if ii == 1; hold(hca,'on'); end
    plot(hca,saveParticle{ii}.r(:,1)*1e-3,saveParticle{ii}.r(:,2)*1e-3)
  end
  hold(hca,'off')  
  
  hca = subplot(nPanels,1,4);
  zz = linspace(-3*d,3*d,100);
  plot(zz,Bx(0,0,zz)*1e9,zz,By(0,0,zz)*1e9,zz,Bz(0,0,zz)*1e9,zz,sqrt(Bx(0,0,zz).^2+By(0,0,zz).^2+Bz(0,0,zz).^2)*1e9)
  %hca.YTick = -10:5:10;
  hca.YGrid = 'on';
end