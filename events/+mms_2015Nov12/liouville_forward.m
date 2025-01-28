% liouville_forward

units = irf_units;
% Observed data
ic = 1;
toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
c_eval('tintObs = tintObs+-toffset?;',ic) % tinObs(2) correspond
CS_normal_velocity = 70; % km/s

tintObs = tintObs + 1*[-1 1];

c_eval(['obsPDist = ePDist?.tlim(tintObs);' ...
        'obsB = mvaB?.tlim(tintObs);' ...
        'obsE = mvaE?.tlim(tintObs);'],ic)

zf0 = 0; % starting/ending point of where to get the f
nPSD = numel(zf0); % number of different Liouville mappings to do;
time = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
% c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0(?)/CS_normal_velocity;',1:nPSD)
tind = find(abs(obsPDist.time-time)==min(abs(obsPDist.time-time)));
  
c_eval('gseBref = mean(gseB?.tlim(time+[-0.005 0.005]).data,1);',ic)
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;
zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsB = zObs;
zObsE = (obsE.time.epochUnix-mean(obsE.time.epochUnix))*CS_normal_velocity;

% Model parameters
mms_2015Nov12.Bmodel;

%xyzB = [Bx(zz*0,zz*0,zz*1) By(zz*0,zz*0,zz*1) Bz(zz*0,zz*0,zz*1)];
%xyzE = [Ex(zz*0,zz*0,zz*1) Ey(zz*0,zz*0,zz*1) Ez(zz*0,zz*0,zz*1)];

tsEmod = irf.ts_vec_xyz(obsE.time,[Ex(zObsE*0,zObsE*0,zObsE*1e3) Ey(zObsE*0,zObsE*0,zObsE*1e3) Ez(zObsE*0,zObsE*0,zObsE*1e3)]*1e3);
tsBmod = irf.ts_vec_xyz(obsB.time,[Bx(zObsB*0,zObsB*0,zObsB*1e3) By(zObsB*0,zObsB*0,zObsB*1e3) Bz(zObsB*0,zObsB*0,zObsB*1e3)]*1e9);
[modEpar,modEperp] = irf_dec_parperp(tsBmod,tsEmod); modEpar.name = 'mod E par'; modEperp.name = 'mod E perp';
[obsEpar,obsEperp] = irf_dec_parperp(obsB,obsE); obsEpar.name = 'obs E par'; obsEperp.name = 'obs E perp';
%irf_plot({obsE,tsEmod,modEperp,modEpar})
%irf_plot({obsB,tsBmod},'comp')

%% Initialize particles
solidangle = 0;
if 1 % Initialize particles, from scratch
  zf0 = 30; % km
  time = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
  tind = find(abs(obsPDist.time-time)==min(abs(obsPDist.time-time)));
  f0 = obsPDist(tind);
  
  Pol = 8;
  nAz = 16;
  nE = 32; iE = 1:32; nE = numel(iE);
  nP = nE*nAz*nPol;

  edgesAnglePolar = linspace(150,180,nPol+1); % angle to B
  %edgesAnglePolar = linspace(0,30,nPol+1); % angle to B
  dAnglePolar = diff(edgesAnglePolar(1:2));
  anglePolar = edgesAnglePolar(2:end) - 0.5*dAnglePolar;
  %anglePolar = ePDist1.depend{1}(1,14:16); nPol = numel(anglePolar);
  edgesAngleAzimuthal = linspace(0,360,nAz+1); % perp-plane angle
  dAngleAzimuthal = diff(edgesAngleAzimuthal(1:2));
  angleAzimuthal = edgesAngleAzimuthal(2:end) - 0.5*dAngleAzimuthal;
  
  energyElectron = ePDist1.depend{1}(1,iE); % eV

  [eENERGY,aAZIM,aPOLAR] = meshgrid(energyElectron,angleAzimuthal,anglePolar);
  eVEL = sqrt(eENERGY*units.eV*2/units.me)/1000; % km/s
  
  % Solid angle for the different particles
  dAngleAzimuthal = diff(angleAzimuthal(1:2)); % constant 
  solidAngle = nan(nPol,1);
  
  if solidangle
  for ipol = 1:nPol
    pol_tmp = linspace(edgesAnglePolar(ipol),edgesAnglePolar(ipol+1),100);
    %dSolidAngle = sind(pol_tmp);
    solidAngle(ipol) = trapz(pol_tmp,sind(pol_tmp))*dAngleAzimuthal*pi/180; % radians
    %solidAngle(ipol) = (cosd(edgesAnglePolar(ipol+1))-cosd(edgesAnglePolar(ipol)))*dAngleAzimuthal*pi/180; % radians
    %solidAngle_(ipol) = dAnglePolar*pi/180*dAngleAzimuthal*pi/180;
    %plot(pol_tmp,dSolidAngle);
    %title(sprintf('solid angle = %g',solidAngle(ipol)))
    %pause
  end
  [~,~,SOLIDANGLE] = meshgrid(energyElectron, angleAzimuthal,solidAngle);
  end
  
  % Initial positions and velocitites
  x0 = zeros(1,nP); % km
  y0 = zeros(1,nP);
  z0 = zeros(1,nP) + zf0;

  vx0_pa = eVEL.*cosd(aPOLAR); % km/s
  vy0_pa = eVEL.*cosd(aAZIM).*sind(aPOLAR);
  vz0_pa = -eVEL.*sind(aAZIM).*sind(aPOLAR);

  % Turn from field-aligned coordinate system (above), to LMN system
  v0_pa = [vx0_pa(:) vy0_pa(:) vz0_pa(:)]';

  B0 = [Bx(x0*1e3,y0*1e3,z0*1e3); By(x0*1e3,y0*1e3,z0*1e3); Bz(x0*1e3,y0*1e3,z0*1e3)]; B0 = B0(:,1)'; % same for all particles since they have same startin position
  b0 = irf_norm(B0);
  b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
  b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);


  rI = [1 0 0; 0 1 0; 0 0 1];
  rA = [b0; b0_perp1; b0_perp2];
  rB = [1 0 0; 0 1 0; 0 0 1];

  % A>B, v' = B(A^-1)v

  % LMN (xyz) coordinate system
  v0 = rB*rA^-1*v0_pa;
  vx0 = v0(1,:); % km/s
  vy0 = v0(2,:);
  vz0 = v0(3,:);

  VX0 = reshape(vx0,nE,nAz,nPol); % km/s
  VY0 = reshape(vy0,nE,nAz,nPol);
  VZ0 = reshape(vz0,nE,nAz,nPol);

  X0 = reshape(x0,nE,nAz,nPol);
  Y0 = reshape(y0,nE,nAz,nPol);
  Z0 = reshape(z0,nE,nAz,nPol);

  x_init_all = double([x0;y0;z0;vx0;vy0;vz0])'*1e3; % m, m/s
  
  if 0
    figure(18); 
    hca = subplot(1,1,1); 
    quiver3(hca,vx0*0,vy0*0,vz0*0,vx0,vy0,vz0); 
    hold(hca,'on')
    hq = quiver3(hca,0,0,0,b0(1)*max(max(v0)),b0(2)*max(max(v0)),b0(3)*max(max(v0)),'r'); 
    hq.LineWidth = 2;
    hold(hca,'off')
    hca.XLabel.String = 'L';
    hca.YLabel.String = 'M';
    hca.ZLabel.String = 'N';
    axis(hca,'equal');
  end
  
  F0obs = squeeze(f0.data(1,:,:,:));
  [VX_obs,VY_obs,VZ_obs] = f0.v(lmn,'squeeze'); % obs, this gives the direction of the 
                                      % bin, which is opposite to the 
                                      % particles that entered the bin  
  vx0_obs = -VX_obs; % change to direction in which the particles are going
  vy0_obs = -VY_obs;
  vz0_obs = -VZ_obs;
  
  nParticles = nP;
  f0_all = nan(nP,1); 
  all_solidAngle = nan(nParticles,1);
  all_polarAngle = nan(nParticles,1);  
  solidangle_ratio = nan(nParticles,1);
  % need to normalize f0 using the solid angle of the particle and the
  % solida angle of the bin from which f_init is taken
  for ip = 1:nParticles    
%     energysize = size(obj.depend{1});
%     theta = obj.depend{3};
%     dangle = pi/16;
%     lengthphi = 32;
% 
%     z2 = ones(lengthphi,1)*sind(theta);
%     solida = dangle*dangle*z2;      
%     allsolida = repmat(solida,1,1,length(dist.time), energysize(2));
%     allsolida = squeeze(permute(allsolida,[3 4 1 2]));
%     dists = dist.data.*allsolida;
%     omni = squeeze(irf.nanmean(irf.nanmean(dists,3),4))/(mean(mean(solida)));
      
    vx0_diff = vx0_obs - vx0(ip);
    vy0_diff = vy0_obs - vy0(ip);
    vz0_diff = vz0_obs - vz0(ip);
    v0_diff = sqrt(vx0_diff.^2 + vy0_diff.^2 + vz0_diff.^2);
    min_ind = find(v0_diff == min(v0_diff(:)));
    [min_en,min_az,min_pol] = ind2sub(size(F0obs),min_ind);    
    if solidangle
      min_pol_tmp = linspace(f0.depend{3}(min_pol)-0.5*diff(f0.depend{3}(1:2)),f0.depend{3}(min_pol)+0.5*diff(f0.depend{3}(1:2)),100);            
      this_solidAngle =    trapz(min_pol_tmp,sind(min_pol_tmp))*diff(f0.depend{2}(1:2))*pi/180; % radians   
      % solidAngle(ipol) = trapz(pol_tmp,sind(pol_tmp))*dAngleAzimuthal*pi/180; % radians    
      all_solidAngle(ip) = this_solidAngle;
      solidangle_ratio(ip) = SOLIDANGLE(ip)/this_solidAngle;
    end
    all_polarAngle(ip) = f0.depend{3}(min_pol);
    f0_all(ip) = F0obs(min_ind);%/this_solidAngle*SOLIDANGLE(ip);
    
  end
  
%   iDist = 1;
% 
%   % Assign phase space density to each particle
%   z0s = unique(z0);
%   c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + z0s(?)/CS_normal_velocity;',1:numel(z0s))
%   c_eval('t?ind = find(abs(obsPitch.time-t?)==min(abs(obsPitch.time-t?)));',1:numel(z0s))
%   c_eval('f0? = obsPitch(t?ind);',1:numel(z0s))
% 
%   c_eval('obsPAbin! = f0?.depend{2}; obsPAbin! = [0 obsPAbin! + obsPAbin!(1)]; ?;',ic,1:numel(z0s))
%   c_eval('obsEbin! = f0?.depend{1}; obsEbin! = [0 obsEbin! + obsEbin!(1)]; ? ;',ic,1:numel(z0s))
%   
%   for iAz = 1:nAz
%     for iPol = 1:nPol
%       [~,zind]=intersect(z0s,Z0(1,iAz,iPol));
%       c_eval('obsEbin = obsEbin?; obsPAbin = obsPAbin?; f0 = f0?;',zind)    
%       for iE = 1:nE    
%         iBinE = iE;
%         iBinPa = find(PA0(iE,iAz,iPol)>obsPAbin,1,'last'  );
%         F0(iE,iPa) = f0.data(1,iBinE,iBinPa);
%       end
%     end
%   end
% 
%   f0mod = reshape(F0,nE*nPa,1);
end
if 0 % Initialize particles, from FPI bins, easier to assign f to a particle like this 
  zf0 = 30; % from where to send in particles
  time = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
  tind = find(abs(obsPDist.time-time)==min(abs(obsPDist.time-time)));
  f0 = obsPDist(tind);

  [VX,VY,VZ] = f0.v(lmn,'squeeze'); % obs, this gives the direction of the 
                                      % bin, which is opposite to the 
                                      % particles that entered the bin  
  vx0 = -VX; % change to direction in which the particles are going
  vy0 = -VY;
  vz0 = -VZ;

  % Transform into LMN coordinate system
  xX = reshape(vx0,size(vx0,1)*size(vx0,2)*size(vx0,3),1);
  yY = reshape(vy0,size(vy0,1)*size(vy0,2)*size(vy0,3),1);
  zZ = reshape(vz0,size(vz0,1)*size(vz0,2)*size(vz0,3),1);

  newTmpX = [xX yY zZ]*lmn(1,:)';
  newTmpY = [xX yY zZ]*lmn(2,:)';
  newTmpZ = [xX yY zZ]*lmn(3,:)';

  vx0 = reshape(newTmpX,size(vx0,1),size(vx0,2),size(vx0,3));
  vy0 = reshape(newTmpY,size(vx0,1),size(vx0,2),size(vx0,3));
  vz0 = reshape(newTmpZ,size(vx0,1),size(vx0,2),size(vx0,3));

  nP = numel(vx0);
  x0 = zeros(1,nP); % km
  y0 = zeros(1,nP); 
  z0 = zeros(1,nP) + zf0;

  normB = gseBref./sqrt(sum(gseBref.^2));
  vtot = sqrt(vx0.^2 + vy0.^2 + vz0.^2);
  normvx0 = vx0./vtot; normvy0 = vy0./vtot; normvz0 = vz0./vtot;

  thetab = acosd(normvx0.*normB(1) + normvy0.*normB(2) + normvz0.*normB(3));

  % only run for electrons that are going in towards the current sheet.
  ingoing = find(thetab>160);find(thetab>130);
  % remove higher energy particles that have f = 0
  empty_bin = find(f0.data == 0);
  keep_ind = setdiff(ingoing(:),empty_bin);

  f0_all = f0.data(keep_ind);
  x_init_all = [double(tocolumn(x0(keep_ind))) double(tocolumn(y0(keep_ind))) double(tocolumn(z0(keep_ind))) ....
                double(vx0(keep_ind))          double(vy0(keep_ind))          double(vz0(keep_ind))]*1e3; % m, m/s
end
  
%% Integration
T = 2; % Integration time
limN = 40e3; % outer N limit for stopping integration
nParticles = numel(f0_all);

tic

x_sol_all = [];
saveParticle = cell(1,nParticles);
%disp(sprintf('nParticles = %g',nParticles))
for iParticle = 1:nParticles;
  if mod(iParticle,100) == 0
    disp(sprintf('iParticle = %g/%g',iParticle,nParticles))
  end
  % Initial positions and velocities                                   
  x_init = x_init_all(iParticle,:); % m, m/s
  
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
  saveParticle{iParticle}.f = f0_all(iParticle);
  %disp(sprintf('energy = %g, energy = %g ',saveParticle{iParticle}.energy(1),eENERGY(iParticle)))
  %saveParticle{iParticle}.dE = dEE_col(iParticle);
  %saveParticle{iParticle}.dv = dVV_col(iParticle);
  %eVEL = sqrt(eENERGY*units.eV*2/units.me)/1000;
end
toc

%% Binning
% Decide binning, match with real observed bins if possible.
times = obsPDist.tlim(tintObs).time;
posN = (times-times(1) +- (times(end)-times(1))/2)*CS_normal_velocity;
posN(abs(posN)>1.1*limN*1e-3) = [];
dN = posN(2)-posN(1);
edgesN = [posN(1)-dN/2; posN+dN/2];
edgesE = [0 f0.depend{1}];
edgesPolarAngle = 0:11.25:180;
edgesAzimuthalAngle = 0:11.25:360;

f_binned_all = zeros(numel(edgesN)-1,numel(edgesE)-1,numel(edgesAzimuthalAngle)-1,numel(edgesPolarAngle)-1);
f_occupied_all = zeros(numel(edgesN)-1,numel(edgesE)-1,numel(edgesAzimuthalAngle)-1,numel(edgesPolarAngle)-1);

allPa = nan(numel(saveParticle),1);
allF = nan(numel(saveParticle),1);
allE0 = zeros(numel(saveParticle),1);
allEend = zeros(numel(saveParticle),1);
allv0 = zeros(numel(saveParticle),1);
%tsFlmap = 
  
%times = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
% c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0(?)/CS_normal_velocity;',1:nPSD)
%tind = find(abs(obsPDist.time-time) == min(abs(obsPDist.time-time)));
nBins = numel(edgesN)-1;

iPlotParticles = 1:numel(saveParticle);
iPincluded = [];
plotIncluded = 0;

tic;
for iP = iPlotParticles; %:numel(saveParticle) % step through particles
  if mod(iP,100) == 0
    disp(sprintf('iParticle = %g/%g',iP,nParticles))
  end
  % Pick out the data
  ind = 1:numel(saveParticle{iP}.t);
  z = saveParticle{iP}.r(ind,3); % N  
  %pa = saveParticle{iP}.pa(ind); % pitch angle  
  
  %allPa(iP) = pa(1);
  %if pa(1)< (130), pa(1); continue; end  
    
%   if any(abs(z)<abs(z(1))*0.6)
%     iPincluded = [iPincluded; iP];
%     if plotIncluded
%       colors = mms_colors('matlab'); color = colors(1,:);
% 
%       thisParticle = saveParticle{iP};
%       thisX = thisParticle.r(:,1)*1e-3; % L
%       thisY = thisParticle.r(:,2)*1e-3; % M
%       thisZ = thisParticle.r(:,3)*1e-3; % N
%       hca = subplot(1,1,1);
%       plot3(hca,thisX,thisY,thisZ,'color',color)  
%       hca.XLabel.String = 'L';
%       hca.YLabel.String = 'M';
%       hca.ZLabel.String = 'N';
%       hca.Title.String = sprintf('iP=%g, Included',iP);
%       view(hca,[0 1 0])
%       pause
%     end
%   else
%     if plotIncluded
%       colors = mms_colors('matlab'); color = colors(2,:);    
%       thisParticle = saveParticle{iP};
%       thisX = thisParticle.r(:,1)*1e-3; % L
%       thisY = thisParticle.r(:,2)*1e-3; % M
%       thisZ = thisParticle.r(:,3)*1e-3; % N
%       hca = subplot(1,1,1);
%       plot3(hca,thisX,thisY,thisZ,'color',color)  
%       hca.XLabel.String = 'L';
%       hca.YLabel.String = 'M';
%       hca.ZLabel.String = 'N';
%       hca.Title.String = sprintf('iP=%g, NOT Included',iP);
%       view(hca,[0 1 0])
%       pause
%     end
%     continue
%   end
%   
  %iP
  v = saveParticle{iP}.v(ind,:); % m/s
  energy = saveParticle{iP}.energy;
  allE0(iP) = energy(1);
  allEend(iP) = energy(end);
  if energy>2000;
    1;
  end
  
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
  

  % plot3(oldx,oldy,oldz,saveParticle{iP}.r(:,1),saveParticle{iP}.r(:,2),saveParticle{iP}.r(:,3))
  
  polarAngle = acosd(-oldvz); % angle from z, oldvz is normalized and ind gse/dsl coordinates
  azimuthalAngle = atan2d(-oldvy,-oldvx); % angle is 0 if vx>0 (vy=0)
  % atan gives angles
  negangles = find(azimuthalAngle<0);
  azimuthalAngle(negangles) = 360+azimuthalAngle(negangles);
  
  % Sort out if angles are wrong, the particle should be put into the box
  % it enters into, not the direction its heading. This could invert the
  % order.
  
  %plot(saveParticle{iP}.t,[polarAngle,azimuthalAngle])
  %plot(saveParticle{iP}.t,[polarAngle/180,oldvz])
  %plot(saveParticle{iP}.t,[azimuthalAngle/360,oldvx]); legend('azimuthal angle/360','vx DSL')
 
 
%   all_T = [all_T; saveParticle{iP}.T];
%   all_z = [all_z; z];
%   all_pa = [all_pa; pa];  
%   all_pa0 = [all_pa0 pa(1)];
%   all_energy = [all_energy; energy];
%   all_zstop = [all_zstop; z(end)];
%   all_zstart = [all_zstart; z(1)];

  % Sort into N and E bins
  [bins_occupied,~,mid,loc] = histcn([z energy azimuthalAngle polarAngle],edgesN*1e3,edgesE,edgesAzimuthalAngle,edgesPolarAngle);
  
  f_occupied = bins_occupied(1:(numel(edgesN)-1),1:32,1:32,1:16);
  f_occupied(f_occupied~=0) = 1;
  f_occupied_all = f_occupied_all + f_occupied;
  f_binned = bins_occupied; f_binned(f_binned~=0) = saveParticle{iP}.f;
  f_binned_all = f_binned_all + f_binned(1:42,1:32,1:32,1:16);
  
  allF(iP) = saveParticle{iP}.f;
  
  if 0 % diagnostic plotting
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
  
%   hca = subplot(3,1,1);
%   pcolor(hca,mid{1}*1e-3,mid{2},bins_occupied')
%   hca.YScale = 'log';
%   hca = subplot(3,1,2);
%   pcolor(hca,mid{1}*1e-3,mid{2},f_binned')
%   hca.YScale = 'log';
%   
%   % Step through bins
%   % loc: one value for each time step with the number of the bin that that
%   % step corresponds to.
%   for iBin = 1:nBins
%     
%   end
%   
%   thisParticle = saveParticle{iP};
%   thisX = thisParticle.r(:,1)*1e-3; % L
%   thisY = thisParticle.r(:,2)*1e-3; % M
%   thisZ = thisParticle.r(:,3)*1e-3; % N
%   
%   hca = subplot(3,1,3);
%   plot3(hca,thisX,thisY,thisZ)    
  
  %bins_occupied(bins_occupied>0) = 1;
  %all_npart = all_npart + bins_occupied;
  %this_vol = sum(v(1,:),2)^2;
  %all_vvol_previous = all_vvol;
  %all_vvol = all_vvol + bins_occupied*this_vol;
  %all_psd_2 = (all_psd_2.*all_vvol_previous + bins_occupied*saveParticle{iP}.f*this_vol)./all_vvol;

  % Add phase space density to N/PA grid
  %all_psd = all_psd + bins_occupied*saveParticle{iP}.f;%/saveParticle{iP}.dv;   
  %all_psd = all_psd.*all_npart + bins_occupied*saveParticle{iP}.f;%/saveParticle{iP}.dv;   
  %all_psd = all_psd./(all_npart+bins_occupied);
  % Add DEF to N/PA grid, more complicated since energy changes    
  %all_def = all_def + bins_occupied*energy(1).^2;%/saveParticle{iP}.dE;
  if numel(iPincluded) > 99
    break;
  end
end
f_binned_all = f_binned_all./f_occupied_all;
toc
%f_binned_all(f_binned_all==0) = NaN;

%% Making PDist from test particles
midN = mid{1};
midTime = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
timesMap = midTime + midN*1e-3/CS_normal_velocity;
timesMap_ = timesMap([1 end])+0.015*[-1 1];
%tsFmap = obsPDist.tlim(timesMap([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
tsFmap.data = f_binned_all;
tsFmap = obsPDist.clone(timesMap,f_binned_all);
tsFmap.depend{1} = repmat(mid{2},tsFmap.length,1); % energy
tsFmap.depend{2} = repmat(mid{3},tsFmap.length,1); % azimuthal angle
tsFmap.depend{3} = mid{4};                         % polar angle
tsFmap.ancillary.energy0 = mid{2};
tsFmap.ancillary.energy1 = mid{2};
tsFmap.ancillary.energy = tsFmap.depend{1};
tsFmap.ancillary.esteptable = ones(tsFmap.length,1);

% Make TSeries for the model E and B


%[VXdb,VYdb,VZdb] = tsFmap.v(lmn,'squeeze');
%a = find(tsFmap.data~=0);
%scatter3(VXdb(a),VYdb(a),VZdb(a));
%hist([thetab(keep_ind) allPa]); legend('initializing','post integration')

%% Plot
npanels = 10;
%[h1,h2] = initialize_combined_plot(npanels,1,1,0.4,'vertical');
h1 = irf_plot(npanels);
elim = [000 1000];
colors = mms_colors('matlab');

zz = zObs*1e3;%tocolumn(linspace(-30,30,200)*1e3);
xyzB = [Bx(zz*0,zz*0,zz*1) By(zz*0,zz*0,zz*1) Bz(zz*0,zz*0,zz*1)];
xyzE = [Ex(zz*0,zz*0,zz*1) Ey(zz*0,zz*0,zz*1) Ez(zz*0,zz*0,zz*1)];
isub = 1;

if 0 % xyzB vs z
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,zz*1e-3,xyzB*1e9,'--')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  plot(hca,zObs,obsB.data)
  hold(hca,'off')
  hca.YLabel.String = 'B_{@B} (nT)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(zz) max(zz)]*1e-3;
end
if 0 % xyzE vs z
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  plot(hca,zz*1e-3,xyzE*1e3,'--')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  plot(hca,zObsE,obsE.data)
  hold(hca,'off')
  hca.YLabel.String = 'E_{@B} (mV/m)';
  hca.XLabel.String = 'N (km)';
  hca.XLim = [min(zz) max(zz)]*1e-3;
end


if 1 % ts B
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  irf_plot(hca,{obsB.x,obsB.y,obsB.z},'comp')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  irf_plot(hca,{tsBmod.x,tsBmod.y,tsBmod.z},'comp','--')
  hold(hca,'off')  
  hca.YLabel.String = 'B (nT)';
end
if 1 % ts E
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',colors(1:3,:))
  irf_plot(hca,{obsE.x,obsE.y,obsE.z},'comp')
  hold(hca,'on')
  set(hca,'ColorOrder',colors(1:3,:))  
  irf_plot(hca,{tsEmod.x,tsEmod.y,tsEmod.z},'comp','--')
  hold(hca,'off')  
  hca.YLabel.String = 'E (mV/m)';    
end


if 1 % deflux obs  omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.omni.tlim(tintObs).deflux.specrec('energy'),'log'); hca.YScale = 'log';
  hca.YLim = [10 1000];
end
if 1 % deflux map  omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.omni.tlim(tintObs).deflux.specrec('energy'),'log'); hca.YScale = 'log';
  hca.YLim = [10 1000];
end
if 1 % deflux obs pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log');
end
if 1 % delux map pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,15).tlim(tintObs).deflux.elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
end
if 1 % psd obs omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.omni.tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';
  hca.YLim = [10 1000];
end
if 1 % psd obs omni
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.omni.tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';
  hca.YLim = [10 1000];
end
if 1 % psd map pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.einterp.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa'),'log');
end
if 1 % psd map pitchangles
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,tsFmap.pitchangles(gseB1,15).tlim(tintObs).elim(elim).specrec('pa'),'log'); %hca.CLim = h1(isub-2).CLim;
end
irf_zoom(h1,'x',tintObs)
irf_plot_axis_align
for ip = 3:6
  h1(ip).CLim = [7 8.3];
end
for ip = 9:10
  h1(ip).CLim = [-26.5 -26.0];
end

%isub = 1;
%hca = h2;
%holdon = 0;
%%
if 0
nPlotP = numel(iPincluded);
indPlot = 1:10;
for iP_ = 1:1:nPlotP
  iP = iPincluded(iP_);
  thisParticle = saveParticle{iP};
  indPlot = 1:numel(thisParticle.r(:,1));
  thisX = thisParticle.r(indPlot,1)*1e-3; % L
  thisY = thisParticle.r(indPlot,2)*1e-3; % M
  thisZ = thisParticle.r(indPlot,3)*1e-3; % N
   
  if 1 % colormap according to f
    cmap = colormap(h1(1)); % parula
    ncmap = size(cmap,1); % 64
    clim = h1(4).CLim;
    colorind = ceil(ncmap*(log10(thisParticle.f)-clim(1))/diff(clim));
    if colorind == 0
      plotcolor = [1 1 1];
    elseif colorind > ncmap
      colorind = 64;
      plotcolor = cmap(fix(colorind),:);
    else
      plotcolor = cmap(fix(colorind),:);
    end
  else % colormap according to ending position
    if thisZ(end)<0;
      plotcolor = [0.1000    0.4000    1.0000];
    else
      plotcolor = [0.9500    0.7000         0];
    end
  end
      
  plot3(hca,thisX,thisY,thisZ,'color',plotcolor)
  hca.XLabel.String = 'L';
  hca.YLabel.String = 'M';
  hca.ZLabel.String = 'N';
  if holdon == 0
    holdon = 1;
    hold(hca,'on')
  end
  pause(0.1)
end
end
