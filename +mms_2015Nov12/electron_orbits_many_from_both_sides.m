%% Electron orbits of several particles, not slice distributions but pitchangle plot
ic = 1;
toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
c_eval('tintObs = tintObs+-toffset?;',ic)
CS_normal_velocity = 70; % km/s
tintObs = tintObs + 3*[-1 1];


c_eval([...
'obsB = mvaB?.tlim(tintObs);'...
'obsCurvB = mvaCurvB.resample(obsB).tlim(tintObs);'...
'obsE = mvaE?.tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsE = mvaE?_new.tlim(tintObs); obsE = obsE.resample(obsB);'... %'obsE = mvaE?_new.tlim(tintObs) - mvaEht.resample(mvaE?_new).tlim(tintObs); obsE = obsE.resample(obsB);'...
'obsEpar = mvaE?par.tlim(tintObs); obsEpar = obsEpar.resample(obsB);'...
'obsEperp = mvaE?perp.tlim(tintObs); obsEperp = obsEperp.resample(obsB);'...
'obsVepar = mvaVe?par.tlim(tintObs); obsVepar = obsVepar.resample(obsB);'...
'obsVeperp = mvaVe?perp.tlim(tintObs); obsVeperp = obsVeperp.resample(obsB);'...
'obsPDist = ePDist?.tlim(tintObs);'...
'obsPitch = ePitch?.tlim(tintObs);'...
],ic)
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;
zObsPitch = (obsPitch.time.epochUnix-mean(obsPitch.time.epochUnix))*CS_normal_velocity;

% Colors
colors = mms_colors('xyz');


% Model parameters
mms_2015Nov12.Bmodel;

%plot(zObs,obsB.data,'.'); hold on; plot(zMod*1e-3,[Bx(0,0,zMod)' By(0,0,zMod)' Bz(0,0,zMod)']*1e9,'--'); hold off;
%plot(zObs,obsE.data,'.'); hold on; plot(zMod*1e-3,[Ex(0,0,zMod)' Ey(0,0,zMod)' Ez(0,0,zMod)']*1e3,'--'); hold off;

% See what part of the model field that is par and perp
B_ = [Bx(0,0,zObs*1e3) By(0,0,zObs*1e3) Bz(0,0,zObs*1e3)]; 
Bnorm = irf_norm(B_);
Evec = [Ex(0,0,zObs*1e3) Ey(0,0,zObs*1e3) Ez(0,0,zObs*1e3)];
Epar =  sum(Evec.*Bnorm,2);
Eperp = Evec - Bnorm.*repmat(Epar,1,3);

% model curvature of B
modCurvB = zeros(numel(zObs),3);
dz = zObs(2) - zObs(1);
modCurvB(2:end,:) = repmat(Bnorm(2:end,3),1,3).*diff(Bnorm)/dz; %nT/km

modCurvB_2 = zeros(numel(zObs),3);
dz = zObs(2) - zObs(1);
modCurvB_2(2:end,:) = repmat(B_(2:end,3),1,3).*diff(B_)/dz; %nT/km
if 0
h(1) = subplot(4,1,1); plot(zObs,[B_ sqrt(sum(B_.^2,2))]*1e9);% hold on; plot(zObs,obsB.data,'.'); hold off;
h(2) = subplot(4,1,2); plot(zObs,Bnorm)
h(3) = subplot(4,1,3); plot(zObs,modCurvB)
h(4) = subplot(4,1,4); plot(zObs,modCurvB_2)
linkaxes(h,'x')
c_eval('h(?).YGrid = ''on'';',1:4)
c_eval('h(?).XLim = [-30 30];',1:4)
end

% Electron test particles
particle_sets = [61 62]; % the particle sets to initialize
nP = 100;
limN = 40e3;
T = 0.4; % integration time

edges_pa = 0:10:180;
edges_z = [-limN-5:1:limN+5];
edges_z = linspace(-limN-5*1e3,limN+5*1e3,100)*1e-3;
[EDGES_PA,EDGES_Z] = meshgrid(edges_pa,edges_z);

tic
for iDist = 1:2;
particle_set = particle_sets(iDist);
switch particle_set % Particle initialization
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

    velocity_angle = [15 20 15 20];
    % Initial positions and velocitites
    x0 = [0 0 0 0]; 
    y0 = [0 0 0 0]; % km
    z0 = [30 30 -30 -30];
    %z0 = -[-15 -15 -15];
    
    velocity_angle_L = [-20 -20 20 20];
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
    E0 = 10; ET = 50;
    electron_energy = E0 + ET*randn(1,nP); % eV
    electron_energy(electron_energy<0) = [];
    while numel(electron_energy)< nP      
      electron_energy = [electron_energy E0 + 50*randn(1,(nP-numel(electron_energy)))];
      electron_energy(electron_energy<0) = [];
    end
      
    nP = numel(electron_energy);    
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000; % km/s
    
    % Initial positions and velocities
    x0 = zeros(1,nP); % km
    y0 = zeros(1,nP); 
    z0 = zeros(1,nP)-32;
    
    pa0 = 0 + 30*randn(1,nP);
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
    E0 = 10; ET = 50;
    electron_energy = E0 + ET*randn(1,nP); % eV
    electron_energy(electron_energy<0) = [];
   while numel(electron_energy)< nP      
      electron_energy = [electron_energy E0 + 50*randn(1,(nP-numel(electron_energy)))];
      electron_energy(electron_energy<0) = [];
   end
          
    nP = numel(electron_energy);    
    vt = sqrt(electron_energy*units.eV*2/units.me)/1000; % km/s
    
    % Initial positions and velocities
    x0 = zeros(1,nP); % km
    y0 = zeros(1,nP); 
    z0 = zeros(1,nP)+32;
    
    pa0 = 180 + 30*randn(1,nP);
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
c_eval('electron_energy? = electron_energy;',iDist)
nParticles = numel(electron_energy);
disp(sprintf('nParticles: %.0f, z0: %g km',nParticles, z0(1,1)));

% Integration, do this separately for two populations
x_sol_all = [];
saveParticle = cell(1,nParticles);
for iParticle= 1:nParticles 

  % Initial positions and velocities                                   
  x_init = [x0;y0;z0;vx0;vy0;vz0]*1e3; % m, m/s
  x_init = x_init(:,iParticle); 
  % Integrate trajectory
  stopfunction = @(t,y) eom.lim(t,y,limN);
  options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine)

  EoM = @(ttt,xxx) eom.general(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
  %EoM = @(ttt,xxx) eom.interp_data(ttt,xxx,0,0,zObs,obsB.x.data,obsB.y.data,obsB.z.data,obsE.x.data,obsE.y.data,obsE.z.data);
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
      
  saveParticle{iParticle}.t = t;
  saveParticle{iParticle}.T = t(end);
  saveParticle{iParticle}.r = x_sol(:,1:3);
  saveParticle{iParticle}.r0 = [x0(iParticle),y0(iParticle),z0(iParticle)];
  saveParticle{iParticle}.v = x_sol(:,4:6);
  saveParticle{iParticle}.v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)];
  saveParticle{iParticle}.B = Bxyz;
  saveParticle{iParticle}.pa = pitchangle;  
end

% Bin particles to get 'probability density' of pitchangles
all_T = [];
all_z = [];
all_pa = [];
all_pa0 = [];
all_energy = [];
if 1
  for iP = 1:numel(saveParticle) % Electron pitchangles
    % Pick out the data
    ind = 1:numel(saveParticle{iP}.pa);
    z = saveParticle{iP}.r(ind,3); % N
    pa = saveParticle{iP}.pa(ind); % pitch angle
    v = saveParticle{iP}.v(ind,:);    
    energy = v(:,1).^2 + v(:,2).^2 + v(:,3).^2;
    all_T = [all_T; saveParticle{iP}.T];
    all_z = [all_z; z];
    all_pa = [all_pa; pa];  
    all_pa0 = [all_pa0 pa(1)];
    all_energy = [all_energy; energy];
  end
else
  %[all_r,all_r0,all_v,all_v0,all_B,all_pa] = deal(saveParticle);
end
c_eval('[nt?,edges?,mid?,loc?] = histcn([all_z all_pa],edges_z*1e3,edges_pa);',iDist)
c_eval('nt?(nt?==0) = NaN;',iDist);
c_eval('all_T? = all_T;',iDist)
if 0

c_eval('nt = nt?; loc = loc?;',iDist)
for iz = 1:numel(unique(loc(:,1)));
  for ipa = 1:numel(unique(loc(:,2)));
    %iz 
    %ipa
    rows = intersect(loc1,[iz ipa],'rows');
    if ~isempty(rows)
      nt_E(iz,ipa) = sum(all_energy(rows,:));
    else
      nt_E(iz,ipa) = NaN;
    end
  end
end
c_eval('nt?E = nt_E;',iDist)
end
allParticleSets{iDist} = saveParticle;
toc
end

%% Plot
figure(31)
nPanels = 9;
clear h;
c_eval('h(?) = subplot(nPanels,1,?);',1:nPanels)
linkaxes(h,'x')
clim = 500;
cmap = cn.cmap('bluepink');
isub = 1;
if 1 % Magnetic field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data],'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)

  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 0 % Magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
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
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 1 % Electric field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data],'.');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)

  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3 3];
end  
if 0 % Electric field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data obsE.abs.data],'.');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)

  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Ex(0,0,zMod).^2+Ey(0,0,zMod).^2+Ez(0,0,zMod).^2)*1e3,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3 3];
end  
if 0 % Electric field, perp par
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1b');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsEperp.data obsEpar.data],'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'E_{\perp,L}','E_{\perp,M}','E_{\perp,N}','E_{||}'},[0.01 0.2],'fontsize',14)

  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  plot(hca,zObs,Eperp(:,1)*1e3,'--','color',B_colors(1,:));
  plot(hca,zObs,Eperp(:,2)*1e3,'--','color',B_colors(2,:))
  plot(hca,zObs,Eperp(:,3)*1e3,'--','color',B_colors(3,:))
  plot(hca,zObs,Epar*1e3,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3 3];
end  
if 0 % plot model "phase space", both from left and right at both times, separately
  hca = h(isub); isub = isub + 1;
  hs1 = surf(hca,edges_z,edges_pa,EDGES_PA',-nt1'); hs1.FaceAlpha = 0.8;
  hold(hca,'on');
  hs2 = surf(hca,edges_z,edges_pa,EDGES_PA',nt2'); hs2.FaceAlpha = 0.8;
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];  
  colormap(hca,cmap(end:-1:1,:));
  hca.CLim = clim*[-1 1];
  hcb.YLim = hca.CLim;
end
if 1 % plot model "phase space", both from left and right at both times, mized/added up
  hca = h(isub); isub = isub + 1;
  nt1_ = nt1; nt1_(isnan(nt1_))=0;
  nt2_ = nt2; nt2_(isnan(nt2_))=0;
  nt12 = nt1_+nt2_;
  hs1 = surf(hca,edges_z,edges_pa,EDGES_PA',nt12'); hs1.FaceAlpha = 0.8;  
  hcb = colorbar('peer',hca);
  colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];  
  colormap(hca,cmap(end:-1:1,:));
  hca.CLim = 2*clim*[-1 1];
  hcb.YLim = [0 hca.CLim(2)];
end
if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [10 100];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPitch,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')

  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(4).CLim;  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [100 200];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPitch,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')

  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(4).CLim;  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [200 400];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPitch,plotPitch.f,log10(plotPitch.p'))
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
  %hca.CLim = h(4).CLim;  
  xlabel(hca,'N (km)')  
end
if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [250 1000];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPitch,plotPitch.f,log10(plotPitch.p'))
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
  %hca.CLim = h(4).CLim;  
  xlabel(hca,'N (km)')  
end
if 1 % plot model "phase space"
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt1'))
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];
  hca.CLim = clim*[-1 1];
  colormap(hca,cmap);
  
  hcb.YLim = [0 hca.CLim(2)];
end
if 1 % plot model "phase space"
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt2'))
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180]; 
  colormap(hca,cmap(end:-1:1,:));
  hca.CLim = clim*[-1 1];
  hcb.YLim = [0 hca.CLim(2)];
end
if 0 % plot model "phase space", deflux, multiply nt with v^2
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt1E')*units.e/units.me*units.eV)
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];
  %hca.CLim = clim*[-1 1];
  colormap(hca,cmap);  
  %hcb.YLim = [0 hca.CLim(2)];
end
if 0 % plot model "phase space", deflux, multiply nt with v^2
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt2E'))
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180]; 
  colormap(hca,cmap(end:-1:1,:));
  %hca.CLim = clim*[-1 1];
  %hcb.YLim = [0 hca.CLim(2)];
end

c_eval('h(?).XLim = [-30 30];',1:3)
irf_plot_axis_align


c_eval('h(?).Position(4) = h(?).Position(4)*1.5;',1:numel(h))

%% Plot, paper
figure(33)
nPanels = 9;
clear h;
c_eval('h(?) = subplot(nPanels,1,?);',1:nPanels)
c_eval('h(?).Position(4) = h(?).Position(4)*1.3;',1:numel(h))
%h(2:2:2*(nPanels-2)) = [];
%h(7) = subplot(nPanels,2,nPanels*2-1);
%h(8) = subplot(nPanels,2,nPanels*2);
linkaxes(h(1:nPanels-1),'x')
clim = 499;
cmap = cn.cmap('bluepink');
isub = 1;
if 1 % Magnetic field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,obsB.data,'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 0 % Curvature of magnetic field, obs+mod
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,obsCurvB.data,'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  modfac = 0.1;
  lineMod = plot(hca,zObs,modCurvB(:,1)*modfac,'--','color',B_colors(1,:));
  plot(hca,zObs,modCurvB(:,2)*modfac,'--','color',B_colors(2,:))
  plot(hca,zObs,modCurvB(:,3)*modfac,'--','color',B_colors(3,:))  
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'curb B (1/km)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30]; 
end
if 1 % Curvature of magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1'); 
  set(hca,'colororder',B_colors)
  lines = plot(hca,zObs,obsCurvB.data,'-');
  c_eval('lines(?).Color = B_colors(?,:);',1:3)
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'L','M','N'},[0.01 0.2],'fontsize',14)
  hca.YGrid = 'on';
  hca.YLabel.String = {'curv B','(1/km)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  %hca.YLim = [-13 13];
end
if 0 % Magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
 
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
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'B (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-13 13];
end
if 1 % Electric field, no absolute value
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data],'-');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  set(hca,'colororder',B_colors)
  %irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)
  irf_legend(hca,{'L','M','N'},[0.01 0.95],'fontsize',14)
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'E','(mV/m)'};
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3.9 3.9];
end  
if 0 % Electric field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsE.data obsE.abs.data],'.');
  linesObs(1).Color = B_colors(1,:); 
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'E_L','E_M','E_N'},[0.01 0.2],'fontsize',14)

  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Ex(0,0,zMod)]*1e3,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[Ey(0,0,zMod)]*1e3,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Ez(0,0,zMod)]*1e3,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Ex(0,0,zMod).^2+Ey(0,0,zMod).^2+Ez(0,0,zMod).^2)*1e3,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3 3];
end  
if 0 % Electric field, perp par
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1b');
  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsEperp.data obsEpar.data],'.');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'E_{\perp,L}','E_{\perp,M}','E_{\perp,N}','E_{||}'},[0.01 0.2],'fontsize',14)

  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,100)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  plot(hca,zObs,Eperp(:,1)*1e3,'--','color',B_colors(1,:));
  plot(hca,zObs,Eperp(:,2)*1e3,'--','color',B_colors(2,:))
  plot(hca,zObs,Eperp(:,3)*1e3,'--','color',B_colors(3,:))
  plot(hca,zObs,Epar*1e3,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = 'E (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [-30 30];
  hca.YLim = [-3 3];
end  
if 1 % eDist omni 64 perp
  hca = h(isub); isub = isub + 1;
  pas = [75 105];
  c_eval('ePitch?lim = ePDist?.pitchangles(dmpaB?,pas);',ic)
  c_eval('obsPitchPerp = ePitch?lim.tlim(tintObs).deflux.specrec;',ic)  
  if 1
    plot_z1 = zObsPitch;
    plot_E1 = obsPitchPerp.f(1,:);
    plot_p1 = log10(obsPitchPerp.p);
    plot_p1(2:2:end,:) = NaN;
    pcolor(hca,plot_z1,plot_E1,plot_p1')
    hold(hca,'on')
    plot_z2 = zObsPitch;
    plot_E2 = obsPitchPerp.f(2,:);
    plot_p2 = log10(obsPitchPerp.p);
    plot_p2(1:2:end,:) = NaN;
    pcolor(hca,plot_z2,plot_E2,plot_p2')
    hold(hca,'off')
  else
    pcolor(hca,plot_z1,obsPitchPerp.f(1,:),log10(obsPitchPerp.p'))
  end
  set(hca,'yscale','log');
  set(hca,'ytick',[1e1 1e2 1e3 1e4]);
  shading(hca,'flat');
  hca.YLabel.String = {'E_e','(eV)'}; 
  colormap(hca,'jet')   
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = obsPitchPerp.p_label;  
  irf_legend(hca,{[num2str(pas(1),'%.0f') '<\theta<' num2str(pas(2),'%.0f')]},[0.02 0.1],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 2e3];
  hca.CLim = [5.5 8.2];
  hcb.YLabel.String = ' ';
end
if 1 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [10 1000];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPitch,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hcb = colorbar('peer',hca); hca.CLim = [7.4 8.1];
  hcb.YLabel.String = plotPitch.p_label;  
  hcbref = hcb;
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  hca.YLabel.String = {'\theta_{PA} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.02 0.7],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(4).CLim;  
  xlabel(hca,'N (km)')  
  
end
if 0 % Distance: ePDist pa low energies
  hca = h(isub); isub = isub + 1;
  elim = [400 1000];
  plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
  pcolor(hca,zObsPitch,plotPitch.f,log10(plotPitch.p'))
  shading(hca,'flat')
  hcb = colorbar('peer',hca); %hca.CLim = [6.4 7.5];  
  hcb.YLabel.String = plotPitch.p_label;  
  %irf_pl_mark(hca,tref,'k')
  hca.XGrid = 'off';
  hca.YGrid = 'off';    
  %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
  %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
  hca.YTick = [45 90 135];   
  colormap(hca,'jet')
  %irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.98 0.15],'fontsize',12,'color',[0 0 0]);
  %hca.CLim = h(4).CLim;  
  xlabel(hca,'N (km)')  
end
if 1 % plot model "phase space"
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt1'))
  hcb = colorbar('peer',hca);
  hcb_ = [];
  hcb_ = [hcb_ hcb];
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];
  hca.CLim = clim*[-1 1];
  colormap(hca,cmap);
  
  hcb.YLim = [0 hca.CLim(2)];
  hca.Box = 'on';
  hca.ZTick = [45 90 135]; 
  hca.ZLabel.String = {'\theta_{PA} (\circ)'};
  hcb.YLabel.String = 'counts';  
  hcb.Position([1 3 4]) = hcbref.Position([1 3 4]); 
end
if 1 % plot model "phase space"
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt2'))
  hcb = colorbar('peer',hca);
  hcb_ = [hcb_ hcb];
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180]; 
  colormap(hca,cmap(end:-1:1,:));
  hca.CLim = clim*[-1 1];
  hcb.YLim = [0 hca.CLim(2)];
  hca.Box = 'on';
  hca.ZTick = [45 90 135];  
  hca.ZLabel.String = {'\theta_{PA} (\circ)'};
  hcb.YLabel.String = 'counts';  
  hcb.Position([1 3 4]) = hcbref.Position([1 3 4]); 
end
if 0 % plot model "phase space", both from left and right at both times, separately
  hca = h(isub); isub = isub + 1;
  hs1 = surf(hca,edges_z,edges_pa,EDGES_PA',-nt1'); hs1.FaceAlpha = 0.8;
  hold(hca,'on');
  hs2 = surf(hca,edges_z,edges_pa,EDGES_PA',nt2'); hs2.FaceAlpha = 0.8;
  hold(hca,'off');
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];  
  colormap(hca,cmap(end:-1:1,:));
  hca.CLim = clim*[-1 1];
  hcb.YLim = hca.CLim;
end
if 1 % plot model "phase space", both from left and right at both times, mized/added up
  hca = h(isub); isub = isub + 1;
  nt1_ = nt1; nt1_(isnan(nt1_))=0;
  nt2_ = nt2; nt2_(isnan(nt2_))=0;
  nt12 = nt1_+nt2_;
  hs1 = surf(hca,edges_z,edges_pa,EDGES_PA',nt12'); hs1.FaceAlpha = 0.8;  
  hcb = colorbar('peer',hca);
  hcb_ = [hcb_ hcb];
  cmap_ = cn.cmap('purplepink');
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];  
  colormap(hca,cmap_(end:-1:1,:));
  hca.CLim = 1*clim*[-1 1];
  hcb.YLim = [0 hca.CLim(2)];
  hca.Box = 'on';
  hca.ZTick = [45 90 135];  
  hca.ZLabel.String = {'\theta_{PA} (\circ)'};
  hcb.YLabel.String = 'counts';  
  hcb.Position([1 3 4]) = hcbref.Position([1 3 4]); 
end
drawnow
if 0 % plot model "phase space", deflux, multiply nt with v^2
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt1E')*units.e/units.me*units.eV)
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180];
  %hca.CLim = clim*[-1 1];
  colormap(hca,cmap);  
  %hcb.YLim = [0 hca.CLim(2)];
  hca.Box = 'on';
end
if 0 % plot model "phase space", deflux, multiply nt with v^2
  hca = h(isub); isub = isub + 1;
  surf(hca,edges_z,edges_pa,EDGES_PA',(nt2E'))
  hcb = colorbar('peer',hca);
  %colormap(hca,'jet')
  shading(hca,'flat')
  hca.View = [0 0];
  hca.ZLim = [0 180]; 
  colormap(hca,cmap(end:-1:1,:));
  %hca.CLim = clim*[-1 1];
  %hcb.YLim = [0 hca.CLim(2)];
  hca.Box = 'on';
end
irf_plot_axis_align(h(1:end-1))
hca.XLabel.String = 'N (km)';
c_eval('h(?).XTick = [];',1:isub-2)
if 1 % Electrons injected, histogram
  hca = h(isub); isub = isub + 1;
  edges_e = linspace(0,2*ET+E0,21);
  c_eval('[nt?_e,loc?_e] = histc(electron_energy?,edges_e);',1:2)
  isub = 1;  
  if 0
    hh1 = histogram(hca,electron_energy1,edges_e);  hh1.FaceColor = [1 0 0.5]; hh1.FaceAlpha = 0.5; % set(get(hh1,'Children'),'FaceAlpha',0.3);
    hold(hca,'on')
    hh2 = histogram(hca,electron_energy2,edges_e); hh2.FaceColor = [0 0 1]; hh2.FaceAlpha = 0.3;    
    hold(hca,'off')
    
    hca.XLim = edges_e([1 end]);
  else
    hca.Position(3) = hca.Position(3)*0.45;    
    hca1 = hca;
    hca2 = subplot(nPanels,2,nPanels*2);
    hca2.Position(2) = hca1.Position(2);
    hca2.Position(4) = hca1.Position(4);
    hca2.Position(1) = h(1).Position(1) + h(1).Position(3) - hca1.Position(3)+ 0.02;
    
    hh1 = bar(hca1,edges_e,nt1_e); hh1.FaceColor = [1.0    0.5    0.7]; %hh1.FaceAlpha = 0.5; % set(get(hh1,'Children'),'FaceAlpha',0.3);
    hold(hca,'on')    
    hh2 = bar(hca2,edges_e,nt2_e); hh2.FaceColor = [0.2 0.5 0.9];[0.5625    0.7581    0.8868];[0.7 0.7 1];
    hca1.XLim = edges_e([1 end]);
    hca2.XLim = edges_e([1 end]);    
    hca2.YTickLabel = [];
    c_eval('hca?.XLabel.String = ''Energy (eV)'';',1:2)
    hca1.YLabel.String = '# electrons';
    c_eval('hca?.Position(2) = hca?.Position(2)-0.05;',1:2);
    %c_eval('hca?.YLim = [0 20];',1:2)
    hca1.YLim = [0 max([hca1.YLim(2) hca2.YLim(2)])];
    hca2.YLim = hca1.YLim;
  end
  
end
c_eval('h(?).XLim = [-30 30];',1:7)
c_eval('hcb_(?).Position([1 3 4]) = hcbref.Position([1 3 4]);',1:3)
h(2).YLim = [-0.019 0.019];

irf_legend(h(6),'e^- coming from left (-N)',[0.02 0.98]);
irf_legend(h(7),'e^- coming from right (+N)',[0.02 0.98]);
irf_legend(h(8),'all e^-',[0.02 0.98]);
c_eval('h(?).YLabel.String = '' '';',5)
c_eval('h(?).ZLabel.String = '' '';',6:8)
h(5).YLabel.String = 'Pitchangle (\circ)';
h(5).YLabel.Position(2) = -1;
c_eval('h(?).FontSize = 13;',1:9)
hcbref.YLabel.Position(2) = hcbref.YLabel.Position(2)+0.4;

%%
figure(32)
nPanels = 2;
c_eval('h(?) = subplot(nPanels,1,?);',1:nPanels)
edges_e = linspace(0,100,21);
c_eval('[nt?_e,loc?_e] = histc(electron_energy?,edges_e);',1:2)
isub = 1;
hca = h(isub); isub = isub + 1;
hh = bar(hca,edges_e,nt1_e);
hca = h(isub); isub = isub + 1;
hh = bar(hca,edges_e,nt2_e);









