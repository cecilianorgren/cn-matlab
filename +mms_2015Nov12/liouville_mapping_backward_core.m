% Liouville mapping backward
%% Observed fields, for comparison
units = irf_units;
% Observed data
ic = 1;
toffset1 = 0; toffset2 = -0.03; toffset3 = 0.1; toffset4 = 0.1;
tintObs = irf.tint('2015-11-12T07:19:20.65Z/2015-11-12T07:19:21.70Z');
c_eval('tintObs = tintObs+-toffset?;',ic) % tinObs(2) correspond
CS_normal_velocity = 70; % km/s
%CS_normal_velocity = 70; % km/s

tintObs = tintObs + 1*[-1 1];

c_eval([...
'obsPDist = ePDist?.tlim(tintObs);'...
],ic)
%zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;
zObsPDist = (obsPDist.time.epochUnix-mean(obsPDist.time.epochUnix))*CS_normal_velocity;
%zObsPitch = (obsPitch.time.epochUnix-mean(obsPitch.time.epochUnix))*CS_normal_velocity;

% Model parameters
mms_2015Nov12.Bmodel;

%% Initialize 'final' distributions
% We will start from these and integrate each point in the mms data 
% (32x32x16) back in time until they reach either + or + 30 km. We will
% then use the f they arrive at to assign f at the starting point.

limN = 25e3;

% Outer distributions to map from
t_plus  = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 +  limN/CS_normal_velocity;
t_minus = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 +- limN/CS_normal_velocity;
ind_t_plus  = find(abs(obsPDist.time-t_plus)  == min(abs(obsPDist.time-t_plus)));
ind_t_minus = find(abs(obsPDist.time-t_minus) == min(abs(obsPDist.time-t_minus)));
f0_plus = obsPDist(ind_t_plus); %  f0_plus.ancillary.esteptable = f0_plus.ancillary.esteptable(ind_t_plus,:);
f0_minus = obsPDist(ind_t_minus); % f0_minus.ancillary.esteptable = f0_minus.ancillary.esteptable(ind_t_plus,:);
[VX_plus,VY_plus,VZ_plus] = f0_plus.v(lmn,'squeeze');
[VX_minus,VY_minus,VZ_minus] = f0_minus.v(lmn,'squeeze');

if 0 % get for one z
  zf0 = [10 15]; % starting/ending point of where to get the f
  nPSD = numel(zf0); % number of different Liouville mappings to do;
  times = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0/CS_normal_velocity;
  % c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + zf0(?)/CS_normal_velocity;',1:nPSD)
  c_eval('t?ind = find(abs(obsPDist.time-times(?))==min(abs(obsPDist.time-times(?))));',1:nPSD)
  tind = []; c_eval('tind = [tind t?ind];',1:nPSD)    
  %c_eval('f0? = obsPDist(t?ind); f0?.ancillary.esteptable = f0?.ancillary.esteptable(t?ind,:);',1:nPSD)
  %c_eval('f0? = obsPDist(t?ind);',1:nPSD)
else % many, defined by the indices
  tind = 44:1:59;  
  nPSD = numel(tind);
  %c_eval('t?ind = tind(?);',1:nPSD)
  %c_eval('t? = obsPDist(t?ind);',1:nPSD)
  %c_eval('f0? = obsPDist(t?ind);',1:nPSD) 
  zf0 = zObsPDist(tind);  
end
F = obsPDist(tind);


%% Integration
iE_ = [2:19];
iAz_  = [1:1:32];
iPol_  = [1:1:16];

nE = 32;numel(iE_);
nAz = 32;numel(iAz_);
nPol = 16;numel(iPol_);

nP = nE*nAz*nPol;

allParticles = cell(1,nPSD); % will contain saveParticle
allF = cell(1,nPSD);
allf = zeros(nPSD,32,32,16);

for iPSD = 1:nPSD % do many trajectory integrations for each 'final' distribution  
  %% Initialization of phase space grid
  % measured distribution
  f0 = obsPDist(tind(iPSD));    
  
  % initialize psd  
  f = zeros(32,32,16);
  pa_start = nan(32,32,16);
  pa_stop = nan(32,32,16);
  z_stop = nan(32,32,16);  
  
  % initialze starting position and velocities
  x0 = zeros(nE,nAz,nPol);
  y0 = zeros(nE,nAz,nPol);
  z0 = zeros(nE,nAz,nPol) + zf0(iPSD);
  
  [VX,VY,VZ] = f0.v(lmn,'squeeze');
  
  % Choose one energy level, initially time saving when building script  
  % VX = VX(iE_,iAz_,iPol_);
  % VY = VY(iE_,iAz_,iPol_);
  % VZ = VZ(iE_,iAz_,iPol_);
    
  vx0 = -VX;
  vy0 = -VY;
  vz0 = -VZ;
  
  %% Trajectory integration
  tic
  T = 2; % s, integration time
  
  saveParticle = cell(nE,nAz,nPol);
  disp('iE = ')
  for iE = iE_
    disp(sprintf('%.0f ',iE))
    for iAz = iAz_
      for iPol = iPol_
        % Initial positions and velocities                                   
        x_init = [x0(iE,iAz,iPol),...
                  y0(iE,iAz,iPol),...
                  z0(iE,iAz,iPol),...
                  vx0(iE,iAz,iPol),...
                  vy0(iE,iAz,iPol),...
                  vz0(iE,iAz,iPol)]*1e3; % m, m/s        
        x_init = double(x_init);
        
        % Integrate trajectory
        stopfunction = @(t,y) eom.lim(t,y,limN);        
        options = odeset('Events',stopfunction,'RelTol',1e-6);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine)

        EoM = @(ttt,xxx) eom.general_backward(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);        
        [t,x_sol] = ode45(EoM,[0 T],x_init,options); % ,options
        x = x_sol(:,1);
        y = x_sol(:,2);
        z = x_sol(:,3);
        vx = x_sol(:,4);
        vy = x_sol(:,5);
        vz = x_sol(:,6); 
                
        
        if 0; %t(end) > 0.98*T%(abs(z(end))*1e-3 - 30) < 1 % Diagnostics
          disp(sprintf('t(end) = %.5f',t(end)))
          plot3(x*1e-3,y*1e-3,z*1e-3); hold on
          plot3(x(1)*1e-3,y(1)*1e-3,z(1)*1e-3,'go'); hold on
          plot3(x(end)*1e-3,y(end)*1e-3,z(end)*1e-3,'rx'); hold on
          xlabel('L'); ylabel('M'); zlabel('N'); view([0 1 0])
          pause
        end % // diagnostics

        energy = units.me*sum(x_sol(:,4:6).^2,2)/2/units.eV; % eV
        Bxyz = [Bx(x,y,z),By(x,y,z),Bz(x,y,z)]; normBxyz = irf_norm(Bxyz);
        Vxyz = [vx vy vz]; normVxyz = irf_norm(Vxyz);  
        pitchangle = acosd(normBxyz(:,1).*normVxyz(:,1) + ...
                           normBxyz(:,2).*normVxyz(:,2) + ...
                           normBxyz(:,3).*normVxyz(:,3));

        saveParticle{iE,iAz,iPol}.t = t;
        saveParticle{iE,iAz,iPol}.r = x_sol(:,1:3);
        saveParticle{iE,iAz,iPol}.v = x_sol(:,4:6);
        saveParticle{iE,iAz,iPol}.B = Bxyz;
        saveParticle{iE,iAz,iPol}.pa = pitchangle;
        saveParticle{iE,iAz,iPol}.energy = energy; % eV
        %saveParticle{iE,iAz,iPol}.f = f0(iParticle);
        
        % Assign phase space density
        if t(end) > 0.99*T % particle stays a long time inside current sheet
          f_assign = 0;
          disp(sprintf('iE = %.0f, iAz = %.0f, iPol = %.0f',iE,iAz,iPol))
        else           
          if x_sol(:,6) > 0.99*limN; % it ended up to the left
            f0 = f0_plus;
            vx_ = VX_plus*1e-3; % 10^3 km/s
            vy_ = VY_plus*1e-3;
            vz_ = VZ_plus*1e-3;
          else
            f0 = f0_minus;
            vx_ = VX_minus*1e-3; % 10^3 km/s
            vy_ = VY_minus*1e-3;
            vz_ = VZ_minus*1e-3;
          end

          v_end = -x_sol(end,4:6)*1e-6; % integration end velocity
          v_mag = sqrt(vx_.^2 + vy_.^2 + vz_.^2);

          v_diff = sqrt((vx_-v_end(1)).^2 + (vy_-v_end(2)).^2 + (vz_-v_end(3)).^2);
          [v_diff_min,i1] = min(v_diff(:));
          [i1x,i1y,i1z] = ind2sub(size(v_diff),i1);

          if 0 % Diagnositcs
            %%
            colors = mms_colors('matlab');
            hca = subplot(1,1,1);
            quiver3(hca,0,0,0,v_end(1),v_end(2),v_end(3),'color',colors(2,:),'linewidth',2)
            hold(hca,'on')
            quiver3(hca,vx_(i1)*0,vy_(i1)*0,vz_(i1)*0,vx_(i1),vy_(i1),vz_(i1),1,'color',colors(1,:),'linewidth',1)


            %quiver3(hca,vx_(i12)*0,vy_(i12)*0,vz_(i12)*0,vx_(i12),vy_(i12),vz_(i12),1,'color',colors(5,:),'linewidth',2)
            %quiver3(hca,vx_(:)*0,vy_(:)*0,vz_(:)*0,vx_(:),vy_(:),vz_(:),'color',colors(3,:))                    
            %pause
          end
          f_assign = f0.data(i1);
        end        
        saveParticle{iE,iAz,iPol}.f = f_assign;
        f(iE,iAz,iPol) = f_assign;
        pa_start(iE,iAz,iPol) = pitchangle(1);
        pa_stop(iE,iAz,iPol) = pitchangle(end);
        z_stop(iE,iAz,iPol) = x_sol(end,3);
      end
    end
    toc
  end
  
  allf(iPSD,:,:,:) = f;
  %c_eval('tsFmap = f0.clone(t?,reshape(f,[1 size(f)]));',iPSD)
  %c_eval('tsFmap? = tsFmap;',iPSD)
%   if iPDS == 1;
%     ePDist_LM = tsFmap;
%   else
%     ePDist_LM = combine(tsFmap);
%   end 
  allParticles{iPSD}.particles = saveParticle;
  allParticles{iPSD}.f = f;
  %allParticles{iPSD}.tsf = tsFmap;
  %toc
  
end
tsFmap = F;
tsFmap.data = allf;
save(['/home/cecilia/Data/' datestr(now,'yyyy-mm-dd_HH_MM_ss') '.mat'])
