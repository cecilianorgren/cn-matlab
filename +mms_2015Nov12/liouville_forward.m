% liouville_forward

units = irf_units;
% Observed data
ic = 1;
toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
c_eval('tintObs = tintObs+-toffset?;',ic) % tinObs(2) correspond
CS_normal_velocity = 70; % km/s

tintObs = tintObs + 1*[-1 1];

c_eval(['obsPDist = ePDist?.tlim(tintObs);'],ic)

zf0 = 0; % starting/ending point of where to get the f
nPSD = numel(zf0); % number of different Liouville mappings to do;
time = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
% c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0(?)/CS_normal_velocity;',1:nPSD)
tind = find(abs(obsPDist.time-time)==min(abs(obsPDist.time-time)));
  
c_eval('gseBref = mean(gseB?.tlim(time+[-0.005 0.005]+-1).data,1);',ic)
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;

% Model parameters
mms_2015Nov12.Bmodel;

%% % Initialize particles, from scratch
% zf0 = 30; % km
% 
% nPol = 8;
% nAz = 32;
% nE = 32;
% nP = nE*nAz*nPol;
% 
% anglePolar = linspace(0,90,nPol); % angle to B
% angleAzimuthal = linspace(0,360,nAz); % perp-plane angle
% energyElectron = ePDist1.depend{1}(1,:); % eV
% 
% [eENERGY,aAZIM,aPOLAR] = meshgrid(energyElectron,angleAzimuthal,anglePolar);
% eVEL = sqrt(eENERGY*units.eV*2/units.me)/1000;
% 
% Initial positions and velocitites
% x0 = zeros(1,nP); % km
% y0 = zeros(1,nP);
% z0 = zeros(1,nP) + zf0;
% 
% vx0_pa = eVEL.*cosd(aPOLAR); % km/s
% vy0_pa = eVEL.*cosd(aAZIM).*sind(aPOLAR);
% vz0_pa = -eVEL.*sind(aAZIM).*sind(aPOLAR);
% 
% Turn from field-aligned coordinate system (above), to LMN system
% v0_pa = [vx0_pa(:) vy0_pa(:) vz0_pa(:)]';
% 
% B0 = [Bx(x0*1e3,y0*1e3,z0*1e3); By(x0*1e3,y0*1e3,z0*1e3); Bz(x0*1e3,y0*1e3,z0*1e3)]; B0 = B0(:,1)'; % same for all particles since they have same startin position
% b0 = irf_norm(B0);
% b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
% b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);
%     
%     
% rI = [1 0 0; 0 1 0; 0 0 1];
% rA = [b0; b0_perp1; b0_perp2];
% rB = [1 0 0; 0 1 0; 0 0 1];
% 
% A>B, v' = B(A^-1)v
% 
% LMN (xyz) coordinate system
% v0 = rB*rA^-1*v0_pa;
% vx0 = v0(1,:);
% vy0 = v0(2,:);
% vz0 = v0(3,:);
% 
% VX0 = reshape(vx0,nE,nAz,nPol);
% VY0 = reshape(vy0,nE,nAz,nPol);
% VZ0 = reshape(vz0,nE,nAz,nPol);
% 
% X0 = reshape(x0,nE,nAz,nPol);
% Y0 = reshape(y0,nE,nAz,nPol);
% Z0 = reshape(z0,nE,nAz,nPol);
% 
% x_init_all = [x0;y0;z0;vx0;vy0;vz0]'*1e3; % m, m/s
% 
% 
% nParticles = nP;
% iDist = 1;
% 
% Assign phase space density to each particle
% z0s = unique(z0);
% c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + z0s(?)/CS_normal_velocity;',1:numel(z0s))
% c_eval('t?ind = find(abs(obsPitch.time-t?)==min(abs(obsPitch.time-t?)));',1:numel(z0s))
% c_eval('f0? = obsPitch(t?ind);',1:numel(z0s))
% 
% c_eval('obsPAbin! = f0?.depend{2}; obsPAbin! = [0 obsPAbin! + obsPAbin!(1)]; ?;',ic,1:numel(z0s))
% c_eval('obsEbin! = f0?.depend{1}; obsEbin! = [0 obsEbin! + obsEbin!(1)]; ? ;',ic,1:numel(z0s))
% F0 = nan(nE,nAz,nPol);
% 
% for iAz = 1:nAz
%   for iPol = 1:nPol
%     [~,zind]=intersect(z0s,Z0(1,iAz,iPol));
%     c_eval('obsEbin = obsEbin?; obsPAbin = obsPAbin?; f0 = f0?;',zind)    
%     for iE = 1:nE    
%       iBinE = iE;
%       iBinPa = find(PA0(iE,iAz,iPol)>obsPAbin,1,'last'  );
%       F0(iE,iPa) = f0.data(1,iBinE,iBinPa);
%     end
%   end
% end
% 
% f0mod = reshape(F0,nE*nPa,1);
% 
% quiver3(vx0*0,vy0*0,vz0*0,vx0,vy0,vz0); axis equal;

%% Initialize particles, from FPI bins, easier to assign f to a particle like this 
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
ingoing = find(thetab<90);
% remove higher energy particles that have f = 0
empty_bin = find(f0.data == 0);
keep_ind = setdiff(ingoing,empty_bin);

f0_all = f0.data(keep_ind);
x_init_all = [double(tocolumn(x0(keep_ind))) double(tocolumn(y0(keep_ind))) double(tocolumn(z0(keep_ind))) double(vx0(keep_ind)) double(vy0(keep_ind)) double(vz0(keep_ind))]*1e3; % m, m/s

%% Integration
T = 2; % Integration time
limN = 40e3; % outer N limit for stopping integration
nParticles = numel(f0_all);

tic

x_sol_all = [];
saveParticle = cell(1,nParticles);
disp(sprintf('nParticles = %g',nParticles))
for iParticle = 1:1:2000;nParticles;
  if mod(iParticle,10) == 0
    disp(sprintf('iParticle = %g',iParticle))
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
      
  saveParticle{iParticle}.t = t;
  saveParticle{iParticle}.T = t(end);
  saveParticle{iParticle}.r = x_sol(:,1:3);
  saveParticle{iParticle}.r0 = [x0(iParticle),y0(iParticle),z0(iParticle)];
  saveParticle{iParticle}.v = x_sol(:,4:6);
  saveParticle{iParticle}.v0 = [vx0(iParticle),vy0(iParticle),vz0(iParticle)];
  saveParticle{iParticle}.B = Bxyz;
  saveParticle{iParticle}.pa = pitchangle;
  saveParticle{iParticle}.energy = units.me*sum(x_sol(:,4:6).^2,2)/2/units.eV; % eV
  saveParticle{iParticle}.f = f0_all(iParticle);
  %saveParticle{iParticle}.dE = dEE_col(iParticle);
  %saveParticle{iParticle}.dv = dVV_col(iParticle);
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

%tsFlmap = 
  
%times = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
% c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0(?)/CS_normal_velocity;',1:nPSD)
tind = find(abs(obsPDist.time-time)==min(abs(obsPDist.time-time)));
nBins = numel(edgesN)-1;

iPlotParticles = 1:2000;%numel(saveParticle);
iPincluded = [];
plotIncluded = 0;

for iP = iPlotParticles %:numel(saveParticle) % step through particles
  % Pick out the data
  ind = 1:numel(saveParticle{iP}.t);
  z = saveParticle{iP}.r(ind,3); % N
 
  if any(abs(z)<abs(z(1))*0.6)
    iPincluded = [iPincluded; iP];
    if plotIncluded
      colors = mms_colors('matlab'); color = colors(1,:);

      thisParticle = saveParticle{iP};
      thisX = thisParticle.r(:,1)*1e-3; % L
      thisY = thisParticle.r(:,2)*1e-3; % M
      thisZ = thisParticle.r(:,3)*1e-3; % N
      hca = subplot(1,1,1);
      plot3(hca,thisX,thisY,thisZ,'color',color)  
      hca.XLabel.String = 'L';
      hca.YLabel.String = 'M';
      hca.ZLabel.String = 'N';
      hca.Title.String = sprintf('iP=%g, Included',iP);
      view(hca,[0 1 0])
      pause
    end
  else
    if plotIncluded
      colors = mms_colors('matlab'); color = colors(2,:);    
      thisParticle = saveParticle{iP};
      thisX = thisParticle.r(:,1)*1e-3; % L
      thisY = thisParticle.r(:,2)*1e-3; % M
      thisZ = thisParticle.r(:,3)*1e-3; % N
      hca = subplot(1,1,1);
      plot3(hca,thisX,thisY,thisZ,'color',color)  
      hca.XLabel.String = 'L';
      hca.YLabel.String = 'M';
      hca.ZLabel.String = 'N';
      hca.Title.String = sprintf('iP=%g, NOT Included',iP);
      view(hca,[0 1 0])
      pause
    end
    continue
  end
  
  iP
  pa = saveParticle{iP}.pa(ind); % pitch angle
  v = saveParticle{iP}.v(ind,:);    
  energy = saveParticle{iP}.energy;
  
  vtot = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);
  normv = [v(:,1)./vtot v(:,2)./vtot v(:,3)./vtot];
  
  oldxyz = saveParticle{iP}.r*inv(lmn');
  oldx = oldxyz(:,1);
  oldy = oldxyz(:,2);
  oldz = oldxyz(:,3);
  oldvxyz = normv*inv(lmn');
  oldvx = oldvxyz(:,1);
  oldvy = oldvxyz(:,2);
  oldvz = oldvxyz(:,3); 

  % plot3(oldx,oldy,oldz,saveParticle{iP}.r(:,1),saveParticle{iP}.r(:,2),saveParticle{iP}.r(:,3))
  
  polarAngle = acosd(oldvz); % angle from z, oldvz is normalized
  azimuthalAngle = atan2d(oldvy,oldvx); % angle is 0 if vx>0 (vy=0)
  
  
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
  
  f_binned = bins_occupied; f_binned(f_binned~=0) = saveParticle{iP}.f;
  f_binned_all = f_binned_all + f_binned(1:42,1:32,1:32,1:16);
  
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
end
%
midN = mid{1};
midTime = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5;
timesMap = midTime + midN*1e-3/CS_normal_velocity;
timesMap_ = timesMap([1 end])+0.015*[-1 1];
tsFmap = obsPDist.tlim(timesMap([1 end])+0.015*[-1 1]); % need to add one sampling period to get the right number
tsFmap.data = f_binned_all;

%% Plot
npanels = 2;
[h1,h2] = initialize_combined_plot(npanels,1,1,0.4,'vertical');

isub = 1;
hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,tsFmap.pitchangles(gseB1,15).tlim(tintObs).deflux.specrec('pa'),'log');
hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.pitchangles(gseB1,15).tlim(tintObs).deflux.specrec('pa'),'log');

isub = 1;
hca = h2;
holdon = 0;

nPlotP = numel(iPincluded);
for iP_ = 1:nPlotP;
  iP = iPincluded(iP_);
  thisParticle = saveParticle{iP};
  thisX = thisParticle.r(:,1)*1e-3; % L
  thisY = thisParticle.r(:,2)*1e-3; % M
  thisZ = thisParticle.r(:,3)*1e-3; % N
   
  plot3(hca,thisX,thisY,thisZ)  
  hca.XLabel.String = 'L';
  hca.YLabel.String = 'M';
  hca.ZLabel.String = 'N';
  if holdon == 0
    holdon = 1;
    hold(hca,'on')
  end
  pause(0.1)
end

