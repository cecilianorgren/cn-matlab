%% Observed data and model
units = irf_units;
ic = 1;

% Time interval of observed data
t_center = irf_time('2015-11-12T07:19:21.175000012Z','utc>EpochTT'); 
tintObs = t_center + 1*[-1 1];
toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
c_eval('tintObs = tintObs+-toffset?;',ic)
CS_normal_velocity = 70; % km/s

% Observed data
c_eval(['obsPDist = ePDist?.tlim(tintObs);' ...
        'obsB = mvaB?.tlim(tintObs);' ...
        'obsE = mvaE?.tlim(tintObs);'],ic)

zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;
zObsB = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsE = (obsE.time.epochUnix-mean(obsE.time.epochUnix))*CS_normal_velocity;
zObs = zObsB;

% Model parameters
mms_2015Nov12.Bmodel;
% TSeries of model data
modE = irf.ts_vec_xyz(obsE.time,[Ex(zObsE*0,zObsE*0,zObsE*1e3) Ey(zObsE*0,zObsE*0,zObsE*1e3) Ez(zObsE*0,zObsE*0,zObsE*1e3)]*1e3);
modB = irf.ts_vec_xyz(obsB.time,[Bx(zObsB*0,zObsB*0,zObsB*1e3) By(zObsB*0,zObsB*0,zObsB*1e3) Bz(zObsB*0,zObsB*0,zObsB*1e3)]*1e9);
if 1 % enforce Ez (EN) to be perpendicular to B, only par E can then come from Ey (EM)  
  Bnorm = modB.norm;  
  Eperp = cross(Bnorm.resample(modE),cross(modE,Bnorm.resample(modE)));  
  modEenf = Eperp;
  modEenf.data(:,2) = modEenf.data(:,2) + modE.data(:,2);    
end

[modEenfpar,modEenfperp] = irf_dec_parperp(modB,modEenf); modEenfpar.name = 'mod Eenf par'; modEenfperp.name = 'mod Eenf perp';
[modEpar,modEperp] = irf_dec_parperp(modB,modE); modEpar.name = 'mod E par'; modEperp.name = 'mod E perp';
[obsEpar,obsEperp] = irf_dec_parperp(obsB,obsE); obsEpar.name = 'obs E par'; obsEperp.name = 'obs E perp';
%irf_plot({obsE,modE,modEperp,modEpar,modEenf,modEenfperp,modEenfpar})
%irf_plot({obsE,modE,modEperp,modEenf,modEenfperp},'comp')
figure(101); 
h=irf_plot(6);
irf_plot(irf_panel('BL'),{obsB.x,modB.x},'comp')
irf_plot(irf_panel('BM'),{obsB.y,modB.y},'comp')
irf_plot(irf_panel('BN'),{obsB.z,modB.z},'comp')
irf_plot(irf_panel('EL'),{obsE.x,modE.x},'comp')
irf_plot(irf_panel('EM'),{obsE.y,modE.y},'comp')
irf_plot(irf_panel('EN'),{obsE.z,modE.z},'comp')


%% Initialize particles
solidangle = 0;
do64 = 1;
% Define pitch (pol) and azimuthal angle for population stating on each
% side, total nP is 2*nE*nAz*nPol
nPol = 12;
nAz = 24;
%nE = 16; 
iE = 3:34; nE = numel(iE);
nP = 2*nE*nAz*nPol;

% Initial position
x0 = zeros(nP,1); % km
y0 = zeros(nP,1);
z0 = zeros(nP,1);
z0(1:end/2) = -30; % km, half starting on left
z0(end/2+1:end) = 30; % km, half starting on right


% Initial energy and velocity 
if do64
  energyElectron64 = ePDist1.tlim(tintObs).einterp; % eV
  energyElectron = energyElectron64.depend{1}(1,iE); % eV
else
  energyElectron = obsPDist.depend{1}(1,iE); % eV
end

edgesAngleAzimuthal = linspace(0,360,nAz+1); % perp-plane angle
dAngleAzimuthal = diff(edgesAngleAzimuthal(1:2));
angleAzimuthal = edgesAngleAzimuthal(2:end) - 0.5*dAngleAzimuthal;

edgesAnglePolar = [linspace(0,90,nPol+1) linspace(90,180,nPol+1)]; % angle to B
dAnglePolar = diff(edgesAnglePolar(1:2));
anglePolar = edgesAnglePolar(2:end) - 0.5*dAnglePolar;
anglePolar(nPol+1) = [];

[eENERGY,aAZIM,aPOLAR] = meshgrid(energyElectron,angleAzimuthal,anglePolar);
eVEL = sqrt(eENERGY*units.eV*2/units.me)/1000; % km/s    
  
% Initial velocitites
vx0_pa = eVEL.*cosd(aPOLAR); % km/s
vy0_pa = eVEL.*cosd(aAZIM).*sind(aPOLAR);
vz0_pa = -eVEL.*sind(aAZIM).*sind(aPOLAR);
v0_pa = [vx0_pa(:) vy0_pa(:) vz0_pa(:)];

v0 = v0_pa*0;

% Turn from field-aligned coordinate system (above), to LMN system
B0 = [Bx(x0*1e3,y0*1e3,z0*1e3), By(x0*1e3,y0*1e3,z0*1e3), Bz(x0*1e3,y0*1e3,z0*1e3)]; % magnetic field at each particles starting location
absB0 = sqrt(sum(B0.^2,2)); % T
% get gyroradius for each particle, used for adaptive stopping condition
vperp0 = ((vy0_pa(:)*1e3).^2+(vz0_pa(:)*1e3).^2).^0.5; % m/s
rho0 = units.me*vperp0(:)/units.e./absB0;

for iP = 1:nP
  b0 = irf_norm(B0(iP,:));
  b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
  b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);
  % A>B, v' = B(A^-1)v
  rI = [1 0 0; 0 1 0; 0 0 1];
  rA = [b0; b0_perp1; b0_perp2];
  rB = [1 0 0; 0 1 0; 0 0 1];
  % LMN (xyz) coordinate system  
  v0_pa(iP,:);
  v0_ = torow(rB*rA^-1*v0_pa(iP,:)');
  v0(iP,:) = v0_';
end

vx0 = v0(:,1);
vy0 = v0(:,2);
vz0 = v0(:,3);

x_init_all = double([x0, y0, z0, vx0, vy0, vz0])*1e3; % m, m/s


% Get f from obsPDist
f0 = nan(nP,1); 

% Take the direction of the bins to be same for all times
% f_center = obsPDist(find(abs(obsPDist.time-t_center)==min(abs(obsPDist.time-t_center))));
% [VX_obs,VY_obs,VZ_obs] = f_center.v(lmn,'squeeze'); % obs, this gives the direction of the 
%                                               % bin, which is opposite to the 
%                                               % particles that entered the bin  
% vx0_obs = -VX_obs; % change to direction in which the particles are going
% vy0_obs = -VY_obs;
% vz0_obs = -VZ_obs;

t_left = t_center + min(z0)/CS_normal_velocity;
tind_left  = find(abs(obsPDist.time-t_left)==min(abs(obsPDist.time-t_left)));
%f_left = obsPDist(find(abs(obsPDist.time-t_left)==min(abs(obsPDist.time-t_left))));
if do64
  f_left = obsPDist(tind_left+[-1 0]).e64;
else
  f_left = obsPDist(tind_left);
end
[VX_obs_left,VY_obs_left,VZ_obs_left] = f_left.v(lmn,'squeeze');
vx0_obs_left = VX_obs_left; % change to direction in which the particles are going, THIS WAS WRONG, SHOULD BE NO MINUS
vy0_obs_left = VY_obs_left;
vz0_obs_left = VZ_obs_left;

t_right = t_center + max(z0)/CS_normal_velocity;
tind_right = find(abs(obsPDist.time-t_right)==min(abs(obsPDist.time-t_right)));
%f_right = obsPDist(find(abs(obsPDist.time-t_right)==min(abs(obsPDist.time-t_right))));
if do64
  f_right = obsPDist(tind_right+[0 +1]).e64;
else
  f_right = obsPDist(tind_right);
end
[VX_obs_right,VY_obs_right,VZ_obs_right] = f_right.v(lmn,'squeeze');
vx0_obs_right = VX_obs_right; % change to direction in which the particles are going, THIS WAS WRONG, SHOULD BE NO MINUS
vy0_obs_right = VY_obs_right;
vz0_obs_right = VZ_obs_right;

disp('Initializing particles')
fprintf('iP = %5.0f/%5.0f\n',0,nP)

% t0_left = t_center + -30/CS_normal_velocity;
% t0_right = t_center + 30/CS_normal_velocity;
% tind_left  = find(abs(obsPDist.time-t0_left)==min(abs(obsPDist.time-t0_left)));
% tind_right = find(abs(obsPDist.time-t0_right)==min(abs(obsPDist.time-t0_right)));
%f0_left32 = obsPDist(tind_left);
%f0_left64 = obsPDist(tind_left+[0 1]).einterp; f0_left64 = f0_left64(1);
%f0_right32 = obsPDist(tind_right);
%f0_right64 = obsPDist(tind_right+[0 1]).einterp; f0_right64 = f0_right64(1);

%f0_left64 = obsPDist(tind_left+[-1 0]).e64;
%f0_right64 = obsPDist(tind_right+[0 1]).e64;

f0_left_scale = 1;
f0_right_scale = 1.1;
for iP = 1:nP     
  %if mod(iP,100) == 0, disp(sprintf('iParticle = %g/%g',iP,nP)); end % display progress
  if mod(iP,100) == 0, fprintf([repmat('\b', 1, 12) '%5.0f/%5.0f\n'],iP,nP); end % display progress  
  % time and distribution at particle starting location
  %t0 = t_center + z0(iP)/CS_normal_velocity;
%  f0_tmp = obsPDist(find(abs(obsPDist.time-t0)==min(abs(obsPDist.time-t0))));
  %disp(f0_tmp.time.utc)
  
  % find closest bin
  if z0(iP)<0
    vx0_obs = vx0_obs_left;
    vy0_obs = vy0_obs_left;
    vz0_obs = vz0_obs_left;
    f0_tmp = f_left;
  else
    vx0_obs = vx0_obs_right;
    vy0_obs = vy0_obs_right;
    vz0_obs = vz0_obs_right;  
    f0_tmp = f_right;
  end
  vx0_diff = vx0_obs - vx0(iP);
  vy0_diff = vy0_obs - vy0(iP);
  vz0_diff = vz0_obs - vz0(iP);
  v0_diff = sqrt(vx0_diff.^2 + vy0_diff.^2 + vz0_diff.^2);
  min_ind = find(v0_diff == min(v0_diff(:)));
  %[min_en,min_az,min_pol] = ind2sub(size(squeeze(f_right.data(1,:,:,:))),min_ind);      
  [min_en,min_az,min_pol] = ind2sub(size(squeeze(f0_tmp.data(1,:,:,:))),min_ind);      
  
  %quiver3([0 0],[0 0],[0 0],[vx0_obs(min_en,min_az,min_pol) vx0(iP)],[vy0_obs(min_en,min_az,min_pol) vy0(iP)],[vz0_obs(min_en,min_az,min_pol) vz0(iP)])
  %axis equal
  
  % assign f
  f0(iP) = f0_tmp.data(1,min_en,min_az,min_pol);
end
pause(0.1); beep

%% Integration
T = 1; % Integration time
limN = 40e3; % outer N limit for stopping integration
% unfortunately some of the higher energy electrons which have larger
% gyroradius exit immediately, might need to update exit condition, like
% minimum integration time
nParticles = numel(f0);

tic

x_sol_all = [];
saveParticle = cell(1,nParticles);
disp('Integrating particles')
fprintf('iP = %5.0f/%5.0f\n',0,nP)
for iParticle = 1:nParticles
  if 0%f0(iParticle)==0
    continue;
  end
    
  if mod(iParticle,100) == 0, fprintf([repmat('\b', 1, 12) '%5.0f/%5.0f\n'],iParticle,nParticles); end % display progress  
  %if mod(iParticle,100) == 0, disp(sprintf('iParticle = %g/%g',iParticle,nParticles)); end
  % Initial positions and velocities                                   
  x_init = x_init_all(iParticle,:); % m, m/s
  
  
  % Integrate trajectory
  if 0 % fix limN
    stopfunction = @(t,y) eom.lim(t,y,limN);
  else % adaptive limN based on rho0    
    limNrho = double(30e3+1.5*rho0(iParticle));
    if limNrho < limN, limN_ = limN;
    else limN_ = limNrho;
    end
    stopfunction = @(t,y) eom.lim(t,y,limN_); % adaptive
  end
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
      
  if 0
    figure(30)
    hca = subplot(2,1,1);
    plot3(hca,xyz(:,1),xyz(:,2),xyz(:,3))
    hold(hca,'on')
    plot3(hca,x(1),y(1),z(1),'go')
    plot3(hca,x(end),y(end),z(end),'rx')
    plot3(hca,x,y,z)
    hold(hca,'off')
    hca = subplot(2,1,2);
    plot(hca,z,pitchangle)
    pause
  end
  saveParticle{iParticle}.t = t;
  saveParticle{iParticle}.T = t(end);
  saveParticle{iParticle}.r = x_sol(:,1:3);
  saveParticle{iParticle}.r0 = [x0(iParticle),y0(iParticle),z0(iParticle)];
  saveParticle{iParticle}.v = x_sol(:,4:6); % m/s
  saveParticle{iParticle}.v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)]; % km/s
  saveParticle{iParticle}.B = Bxyz;
  saveParticle{iParticle}.pa = pitchangle;
  saveParticle{iParticle}.energy = units.me*sum(x_sol(:,4:6).^2,2)/2/units.eV; % eV
  saveParticle{iParticle}.f = f0(iParticle);
  %disp(sprintf('energy = %g, energy = %g ',saveParticle{iParticle}.energy(1),eENERGY(iParticle)))
  %saveParticle{iParticle}.dE = dEE_col(iParticle);
  %saveParticle{iParticle}.dv = dVV_col(iParticle);
  %eVEL = sqrt(eENERGY*units.eV*2/units.me)/1000;
end
toc

saveParticle_many = saveParticle;
beep

%% Binning
saveParticle = saveParticle_many;
doleft = 1;
doright = 1;
if all([doleft doright]); doall = 1; end

limN = 40e3;

% Decide binning, match with real observed bins if possible.
times = obsPDist.tlim(tintObs).time;
posN = (times-times(1) +- (times(end)-times(1))/2)*CS_normal_velocity;
posN(abs(posN)>1.1*limN*1e-3) = [];
dN = posN(2)-posN(1);
edgesN = [posN(1)-dN/2; posN+dN/2];
%edgesE = [0 f_center.depend{1}];
edgesE = [0 min(ePDist1.depend{1}(1:2,:))]; % bin according to lowest energy levels since most particles are there

edgesPolarAngle = 0:11.25:180;
edgesAzimuthalAngle = 0:11.25:360;

f_binned_all = zeros(numel(edgesN)-1,numel(edgesE)-1,numel(edgesAzimuthalAngle)-1,numel(edgesPolarAngle)-1);
f_occupied_all = zeros(numel(edgesN)-1,numel(edgesE)-1,numel(edgesAzimuthalAngle)-1,numel(edgesPolarAngle)-1);

allPa0 = nan(numel(saveParticle),1);
allF = nan(numel(saveParticle),1);
allE0 = zeros(numel(saveParticle),1);
allEend = zeros(numel(saveParticle),1);
allv0 = zeros(numel(saveParticle),1);
allTend = zeros(numel(saveParticle),1);
allzend = zeros(numel(saveParticle),1);
allzmin = zeros(numel(saveParticle),1);
allz0 = zeros(numel(saveParticle),1);
alla = zeros(numel(saveParticle),1);

disp('Binning particles')
fprintf('iP = %5.0f/%5.0f\n',0,nP)
tic;
for iP = 1:nP%(nP/2+1):nP %:numel(saveParticle) % step through particles
  if mod(iP,100) == 0, fprintf([repmat('\b', 1, 12) '%5.0f/%5.0f\n'],iP,nP); end % display progress  
  %if mod(iP,100) == 0, disp(sprintf('iParticle = %g/%g',iP,nParticles)); end

  if f0(iP)==0
    %continue
  end
  if saveParticle{iP}.t(end)==1
    %disp(sprintf('skipping: %g',iP))
    continue
  end
  if 0%saveParticle{iP}.energy(1)<35
    continue;
  end
  if 0%saveParticle{iP}.energy(1)>250
    continue;
  end
  % Pick out the data
  ind = 1:numel(saveParticle{iP}.t);
  x = saveParticle{iP}.r(ind,1); % L
  y = saveParticle{iP}.r(ind,2); % M  
  z = saveParticle{iP}.r(ind,3); % N  
  allzend(iP) = z(end);
  allzmin(iP) = sign(z(1))*min(abs(z));
  allTend(iP) = saveParticle{iP}.t(end);
  pa = saveParticle{iP}.pa(ind); % pitch angle  
  
  allz0(iP) = z(1);
  allPa0(iP) = pa(1);
  v = saveParticle{iP}.v(ind,:); % m/s
  energy = saveParticle{iP}.energy;
  allE0(iP) = energy(1);
  allEend(iP) = energy(end);
  
  if 0%z(1)<0 
    continue
    end
  
  if z(1)>0 && min(pa) > 90
    %break
  end


  % make selections on which particles to include
  if 0%energy(1)<100
    continue;
  end
  if doleft && doright
  elseif doleft
    if z(1)>0, continue; end
  elseif doright
    if z(1)<0, continue; end
  else 
    continue
  end
  % [count edges mid loc] = histcn(X, edge1, edge2, ..., edgeN)
  zbin = [25 30]*1e3;
  pabin = [75 105];
  [a b c d] = histcn([z, pa],zbin,pabin);
  if any(a>0)
    continue;
  end
  if 0%z(1)>0 && z(end)>0
    %%
    nrows = 4;
    ncols = 1;
    npanels = nrows*ncols;    
    for ipanel = 1:npanels
      h(ipanel) = subplot(nrows,ncols,ipanel);
    end
    difft =  diff(saveParticle{iP}.t);
    diffe =  diff(saveParticle{iP}.energy);
    isub = 1;
    if 1 % particle orbit, LMN
      hca = h(isub); isub = isub + 1;
      hs = plot3(hca,x*1e-3,y*1e-3,z*1e-3);
      hca.Title.String = sprintf('iP = %g',iP);    
      %hca.XLim = [-40 40];
      %hca.YLim = [0 180];
      hca.XLabel.String = 'L (km)';
      hca.YLabel.String = 'M (km)';
      hca.ZLabel.String = 'N (km)';      
    end
    if 1 % particle orbit, NL
      hca = h(isub); isub = isub + 1;
      hs = plot(hca,z*1e-3,x*1e-3);
      hca.Title.String = sprintf('iP = %g',iP);    
      %hca.XLim = [-40 40];
      %hca.YLim = [0 180];
      hca.XLabel.String = 'N (km)';
      hca.YLabel.String = 'L (km)';      
    end
    if 1
      hca = h(isub); isub = isub + 1;
      hs = scatter(hca,z(2:end)*1e-3,pa(2:end),z(2:end)*0+5,diffe./difft);
      hca.Title.String = sprintf('iP = %g',iP);    
      hca.XLim = [-40 40];
      hca.YLim = [0 180];
      hca.XLabel.String = 'km';
      hold(hca,'on')
      hp = patch(hca,[zbin zbin(2) zbin(1)]*1e-3,[pabin(1) pabin(1) pabin(2) pabin(2)],mms_colors('x'));
      hp.FaceAlpha = 0.5;
      hold(hca,'off')
      hcb = colorbar('peer',hca);
      colormap(hca,'parula');
      hcb.YLabel.String = 'diff(energy)/diff(t)';
    end
    if 1
      hca = h(isub); isub = isub + 1;
      scatter(hca,saveParticle{iP}.t(2:end),saveParticle{iP}.energy(2:end),saveParticle{iP}.t(2:end)*0+2,diffe) 
      hca.Title.String = sprintf('iP = %g',iP);    
      hca.XLabel.String = 't';
    end
    %pause
  end
  %alla(iP) = a;

  vtot = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2); % m/s
  allv0(iP) = vtot(1);
  normv = [v(:,1)./vtot v(:,2)./vtot v(:,3)./vtot];
  
  lmnxyz = saveParticle{iP}.r;
  xyzB = [Bx(lmnxyz(:,1),lmnxyz(:,2),lmnxyz(:,3)) By(lmnxyz(:,1),lmnxyz(:,2),lmnxyz(:,3)) Bz(lmnxyz(:,1),lmnxyz(:,2),lmnxyz(:,3))];
  xyzBnorm = xyzB./repmat(sqrt(sum(xyzB.^2,2)),1,3);
  oldxyz = saveParticle{iP}.r*inv(lmn'); % from LMN to dsl/gse xyz
  oldx = oldxyz(:,1);
  oldy = oldxyz(:,2);
  oldz = oldxyz(:,3);
  oldvxyz = normv*inv(lmn');  % from LMN to dsl/gse
  oldvx = oldvxyz(:,1);
  oldvy = oldvxyz(:,2);
  oldvz = oldvxyz(:,3); 
    
  polarAngle = acosd(-oldvz); % angle from z, oldvz is normalized and ind gse/dsl coordinates
  azimuthalAngle = atan2d(-oldvy,-oldvx); % angle is 0 if vx>0 (vy=0)
  % atan gives angles
  negangles = find(azimuthalAngle<0);
  azimuthalAngle(negangles) = 360+azimuthalAngle(negangles);
  
  % Sort into N and E bins
  [bins_occupied,~,mid,loc] = histcn([z energy azimuthalAngle polarAngle],edgesN*1e3,edgesE,edgesAzimuthalAngle,edgesPolarAngle);
  
  if 1%z(1)>0
  f_occupied = bins_occupied(1:(numel(edgesN)-1),1:32,1:32,1:16);
  f_occupied(f_occupied~=0) = 1; % only count one particle once, no matter how many time steps it stays in/passes through there
  f_occupied_all = f_occupied_all + f_occupied; % total number of particles passing through
  f_binned = bins_occupied; 
  this_f = saveParticle{iP}.f;
  if z(1)>0
    this_f = this_f*1.0;  
  end
  f_binned(f_binned~=0) = this_f; % put value of occupied bins to f of each particle
  f_binned_all = f_binned_all + f_binned(1:(numel(edgesN)-1),1:32,1:32,1:16); % add all f's
  
  allF(iP) = saveParticle{iP}.f;
  end
  if energy(1) < 200 && energy(end) > 300
    
  end
  if 0%energy(1) < 200 && energy(end) > 300%0 % diagnostic plotting
    pitchangle = acosd(sum(xyzBnorm.*normv,2));
    
    figure(40)
    nrows = 6;
    ncols = 1;
    npanels = nrows*ncols;
    for ip = 1:npanels
      h(ip) = subplot(nrows,ncols,ip);
    end
    isub = 1;

    colors = mms_colors('matlab');
    if 1
      hca = h(isub); isub = isub + 1;  
      scatter(hca,lmnxyz(:,3)*1e-3,pitchangle)
      hca.XLabel.String = 'N';
      hca.YLabel.String = 'PA';      
      hca.YLim = [0 180];
    end    
    if 1
      hca = h(isub); isub = isub + 1;  
      ax = plotyy(hca,z,polarAngle,z,azimuthalAngle);      
      hca.XLabel.String = 'N';
      hca.YLabel.String = 'Angle';      
      legend('polar','azimuthal')
      ax(1).YLim = [0 180];
      ax(2).YLim = [0 360];
    end   
    if 0
      hca = h(isub); isub = isub + 1;  
      plot3(hca,lmnxyz(:,1)*1e-3,lmnxyz(:,2)*1e-3,lmnxyz(:,3)*1e-3)
      hold(hca,'on')
      plot3(hca,lmnxyz(1,1)*1e-3,lmnxyz(1,2)*1e-3,lmnxyz(1,3)*1e-3,'go')
      plot3(hca,lmnxyz(end,1)*1e-3,lmnxyz(end,2)*1e-3,lmnxyz(end,3)*1e-3,'rx')
      hold(hca,'off')      
      hca.XLabel.String = 'L';
      hca.YLabel.String = 'M';
      hca.ZLabel.String = 'N';  
      view(hca,[0 1 0]); 
      axis(hca,'equal')
    end  
    if 0
      hca = h(isub); isub = isub + 1;  
      plot3(hca,oldx*1e-3,oldy*1e-3,oldz*1e-3)
      hold(hca,'on')
      plot3(hca,oldx(1)*1e-3,oldy(1)*1e-3,oldz(1)*1e-3,'go')
      plot3(hca,oldx(end)*1e-3,oldy(end)*1e-3,oldz(end)*1e-3,'rx')
      hold(hca,'off')      
      hca.XLabel.String = 'X';
      hca.YLabel.String = 'Y';
      hca.ZLabel.String = 'Z';   
      view(hca,double(lmn(2,:)));
      axis(hca,'equal')
    end
    if 1 % pitch angle from tsPDist, psd
      midN_ = mid{1};
      midTime_ = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
      timesMap_ = midTime_ + midN_*1e-3/CS_normal_velocity;
      timesMap_ = timesMap_([1 end]);%+0.015*[-1 1];
      tsFmap_ = obsPDist.tlim(timesMap_([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
      tsFmap_.data = f_binned;
      
      hca = h(isub); isub = isub + 1;
      elim = [0 1000];
      irf_spectrogram(hca,tsFmap_.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
    end
    if 1 % pitch angle from tsPDist, deflux
      midN_ = mid{1};
      midTime_ = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
      timesMap_ = midTime_ + midN_*1e-3/CS_normal_velocity;
      timesMap_ = timesMap_([1 end]);%+0.015*[-1 1];
      tsFmap_ = obsPDist.tlim(timesMap_([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
      tsFmap_.data = f_binned;
      
      hca = h(isub); isub = isub + 1;
      elim = [0 1000];
      irf_spectrogram(hca,tsFmap_.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
    end
    if 1      
      hca = h(isub); isub = isub + 1;      
      edgesPA = linspace(0,180,16);
      [bins_occupied,~,mid,loc] = histcn([z energy pitchangle],edgesN*1e3,edgesE,edgesPA);
      xx = mid{1};
      yy = mid{3};
      cc = squeeze(sum(bins_occupied,2))*saveParticle{iP}.f;
      pcolor(hca,xx,yy,log10(cc)');      
      hb = colorbar('peer',hca);
    end
  end
end
f_binned_all = f_binned_all./f_occupied_all;
toc
%f_binned_all(f_binned_all==0) = NaN;

% Making PDist from test particles
midN = mid{1};
midTime = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
timesMap = midTime + midN*1e-3/CS_normal_velocity;
timesMap_ = timesMap([1 end])+0.015*[-1 1];

tsFmap = obsPDist.tlim(timesMap([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
tsFmap.data = f_binned_all;
tsFmap = obsPDist.clone(timesMap,f_binned_all);
tsFmap.depend{1} = repmat(mid{2},tsFmap.length,1); % energy
tsFmap.depend{2} = repmat(mid{3},tsFmap.length,1); % azimuthal angle
tsFmap.depend{3} = mid{4};                         % polar angle
tsFmap.ancillary.energy0 = mid{2};
tsFmap.ancillary.energy1 = mid{2};
tsFmap.ancillary.energy = tsFmap.depend{1};
tsFmap.ancillary.esteptable = ones(tsFmap.length,1);

if doleft && doright
  tsFmap_all = tsFmap;
elseif doleft
  tsFmap_left = tsFmap;
elseif doright
  tsFmap_right = tsFmap;  
end
pause(0.1); beep

%% Plot
npanels = 10;
%[h1,h2] = initialize_combined_plot(npanels,1,1,0.4,'vertical');
h1 = irf_plot(npanels,'newfigure');
elim = [20 1000];
ylim = [10 1000];
colors = mms_colors('matlab');

isub = 1;

if 1 % ts B
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',[colors(1:3,:); 0 0 0])
  irf_plot(hca,{obsB.x,obsB.y,obsB.z,obsB.abs},'comp')
  hold(hca,'on')
  set(hca,'ColorOrder',[colors(1:3,:); 0 0 0])
  irf_plot(hca,{modB.x,modB.y,modB.z,modB.abs},'comp','--')
  hold(hca,'off')  
  hca.YLabel.String = 'B (nT)';
end
if 1 % ts E
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  irf_plot(hca,{obsE.x,obsE.y,obsE.z},'comp')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  irf_plot(hca,{modE.x,modE.y,modE.z},'comp','--')
  hold(hca,'off')  
  hca.YLabel.String = 'E (mV/m)';    
end

if 1 % deflux obs  omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.omni.tlim(tintObs).elim(elim).deflux.specrec('energy'),'log'); hca.YScale = 'log';
  hca.YLim = [10 1000];
end
if 1 % deflux map  omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.omni.tlim(tintObs).elim(elim).deflux.specrec('energy'),'log'); hca.YScale = 'log';
  hca.YLim = [10 1000];
end
if 0 % deflux obs  par
  hca = h1(isub); isub = isub + 1;
  palim = [0 30];
  irf_spectrogram(hca,obsPDist.pitchangles(gseB1,palim).tlim(tintObs).deflux.specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % deflux map  par
  hca = h1(isub); isub = isub + 1;
  palim = [0 30];
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,palim).tlim(tintObs).deflux.specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % deflux obs  perp
  hca = h1(isub); isub = isub + 1;
  palim = [75 105];
  irf_spectrogram(hca,obsPDist.pitchangles(gseB1,palim).tlim(tintObs).elim(elim).deflux.specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % deflux map  perp
  hca = h1(isub); isub = isub + 1;
  palim = [75 105];
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,palim).tlim(tintObs).elim(elim).deflux.specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % psd obs  perp
  hca = h1(isub); isub = isub + 1;
  palim = [75 105];
  irf_spectrogram(hca,obsPDist.pitchangles(gseB1,palim).tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % psd map  perp
  hca = h1(isub); isub = isub + 1;
  palim = [75 105];
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,palim).tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % deflux obs  apar
  hca = h1(isub); isub = isub + 1;
  palim = [150 180];
  irf_spectrogram(hca,obsPDist.pitchangles(gseB1,palim).tlim(tintObs).deflux.specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % deflux map  apar
  hca = h1(isub); isub = isub + 1;
  palim = [150 180];
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,palim).tlim(tintObs).deflux.specrec('energy'),'log'); hca.YScale = 'log';  
  hca.YLim = ylim;  
  irf_legend(hca,{[num2str(palim(1),'%.0f') '<\theta_B<' num2str(palim(2),'%.0f')]},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % deflux obs pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log');
  irf_legend(hca,{sprintf('%g<E<%g eV',elim)},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % delux map pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
  irf_legend(hca,{sprintf('%g<E<%g eV',elim)},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 0 % psd obs omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.omni.tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';
  irf_legend(hca,{sprintf('%g<E<%g eV',elim)},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 1000];
end
if 0 % psd obs omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.omni.tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';
  irf_legend(hca,{sprintf('%g<E<%g eV',elim)},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
  hca.YLim = [10 1000];
end
if 1 % psd obs pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.einterp.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa'),'log');
  irf_legend(hca,{sprintf('%g<E<%g eV',elim)},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end
if 1 % psd map pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
  irf_legend(hca,{sprintf('%g<E<%g eV',elim)},[0.98 0.10],'fontsize',12,'color',[0 0 0]);
end

irf_zoom(h1,'x',tsFmap.time) % tintObs
irf_plot_axis_align

for ip = 5:6
  %h1(ip).CLim = [7.2 8.2];
end
% for ip = 5:10
%   h1(ip).CLim = [7.0 8.1];
% end

for ip = 7:8
  %h1(ip).CLim = [-26.7 -25.9];
  %h1(ip).CLim = [-28.7 -27.2];
end 

h1(10).CLim = h1(9).CLim;
h1(8).CLim = h1(7).CLim;
