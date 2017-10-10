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

if 1 % get for one z
  zf0 = [0]; % starting/ending point of where to get the f
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
iE_ = [6 9 12];[5:13];
iAz_  = [1:1:32];
iPol_  = [1:1:16];

nE = 32;numel(iE_);
nAz = 32;numel(iAz_);
nPol = 16;numel(iPol_);

nP = nE*nAz*nPol;

allParticles = cell(1,nPSD); % will contain saveParticle
allF = cell(1,nPSD);
allf = zeros(nPSD,32,32,16);
allzstop = zeros(nPSD,32,32,16);

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
    
  vx0 = VX;
  vy0 = VY;
  vz0 = VZ;
  
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
  allzstop(iPSD,:,:,:) = z_stop;
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

tsZstop = F;
tsZstop.data = allzstop;

%% Bin phase space density of particles to get 'probability density' of pitchangles (Liouville's theorem)
mm = units.me/units.mp;
all_T = [];
all_z = [];
all_pa = [];
all_pa0 = [];
all_energy = [];
all_zstop = [];
all_zstart = [];
edges_z = -(limN*1e-3-1):1:(limN*1e-3+1);
edges_pa = -0:10:180;
all_psd = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_psd_2 = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_def = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_npart = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_vvol = zeros(numel(edges_z)-1,numel(edges_pa)-1);
dvt_all = 0;

for iiEE = iE_
  for iiAz = iAz_
    for iiPol = iPol_
      Ps = saveParticle{iiEE,iiAz,iiPol};
      % Pick out the data
      ind = 1:numel(Ps.pa);
      z = Ps.r(ind,3); % N
      pa = Ps.pa(ind); % pitch angle
      v = Ps.v(ind,:);    
      energy = Ps.energy;  
      all_z = [all_z; z];
      all_pa = [all_pa; pa];  
      all_pa0 = [all_pa0 pa(1)];
      all_energy = [all_energy; energy];
      all_zstop = [all_zstop; z(end)];
      all_zstart = [all_zstart; z(1)];

      [bins_occupied,~,mid,loc] = histcn([z pa],edges_z*1e3,edges_pa);
      bins_occupied(bins_occupied>0) = 1;
      all_npart = all_npart + bins_occupied;
      this_vol = sum(v(1,:),2)^2;
      all_vvol_previous = all_vvol;
      all_vvol = all_vvol + bins_occupied*this_vol;  

      % Add phase space density to N/PA grid
      all_psd = all_psd + bins_occupied*Ps.f;%/saveParticle{iP}.dv;     
      % Add DEF to N/PA grid, more complicated since energy changes    
      %all_def = all_def + bins_occupied*energy(1).^2;%/saveParticle{iP}.dE;
    end
  end
end

return


%iP=1;plot(saveParticle{iP}.r(:,3),saveParticle{iP}.energy); hold on; for iP=1:10:nParticles, plot(saveParticle{iP}.r(:,3),saveParticle{iP}.energy); end

%% Plot comparison between observed and mapped distribution
iPSD = 1;

nrows = 4;
ncols = 3;
for ipl = 1:nrows*ncols
  h(ipl) = subplot(nrows,ncols,ipl);
end

  
c_eval('scpot = scPot?;',ic)
c_eval('dmpaB?slow = dmpaB?.resample(gseVe?);',ic)
c_eval('dslE?slow = dslE?.resample(gseVe?);',ic)
c_eval('ePitch = obsPitch;',ic)
 
% Plot format input
vlim = 12*1e3;
elevlim = 15;
strCMap = 'jet';
projclim = [0 5];  
  
  c_eval('hatE = dslE?slow.resample(t!).data/dslE?slow.resample(t!).abs.data;',ic,iPSD)
  c_eval('hatB = dmpaB?slow.resample(t!).data/dmpaB?slow.resample(t!).abs.data;',ic,iPSD)
  c_eval('hatExB = cross(hatE,hatB);',ic)
  par = hatB;
  perp1 = hatExB;
  perp2 = cross(par,perp1);  
  
  zz0 = zf0(iPSD);
  c_eval('time = t?;',iPSD)
  c_eval('timeUTC = t?.utc;',iPSD)
  isub = 1;

  %vectors = {hatExB,'ExB'; hatE,'E'; hatB,'B'};
  vectors = {[1 0 0],'x'; [0 1 0],'y'; [0 0 1],'z'};
  
  c_eval('psd1 = f0?;',iPSD)
  c_eval('psd2 = tsFmap?;',iPSD)  
  c_eval('ePitch = obsPitch(t?ind);',iPSD); it = 1;
  %psd3 = psd2; psd3.data(psd3.data == 0) = NaN; psd3.data = abs(psd3.data-psd1.data);
  c_eval('psd3 = f0_plus;',iPSD)
  c_eval('psd4 = f0_minus;',iPSD)
  
  for idist = 1:4
    c_eval('dist = psd?.convertto(''s^3/km^6'');',idist)
    % Perpendicular plane    
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;perp2;par]; vlabels = {'v_{ExB}','v_{perp2}','v_{||}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);                
    hca.Title.String = {sprintf('N = %.1f km',zz0),timeUTC(1:23)};

    % B plane 1
    hca = h(isub); isub = isub + 1; 
    xyz = [perp1;par;-perp2]; vlabels = {'v_{ExB}','v_{||}','-v_{perp2}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';

    % B plane 2
    hca = h(isub); isub = isub + 1;
    xyz = [perp2;par;perp1]; vlabels = {'v_{perp2}','v_{||}','v_{ExB}'};
    mms.plot_projection(hca,dist,'tint',time,'xyz',xyz,'elevationlim',elevlim,'vlim',vlim,'clim',projclim,'scpot',scpot,'vlabel',vlabels,'vectors',vectors);        
    hca.Title.String = '';

    if 0 % scatter plot, error
      hca = h(isub); isub = isub + 1;
      [vxx,vyy,vzz] = dist.v(lmn,'squeeze');
      [xxx,yyy,zzz] = dist.xyz(lmn,'squeeze');
      EEx = sign(vxx).*units.me.*(vxx.^2)/2/units.eV;
      EEy = sign(vyy).*units.me.*(vyy.^2)/2/units.eV;
      EEz = sign(vzz).*units.me.*(vzz.^2)/2/units.eV;
      EE = units.me*(vxx.^2 + vyy.^2 + vzz.^2)/2/units.eV;
      step = 2;
      scatter3(hca,EEx(1:step:end),EEy(1:step:end),EEz(1:step:end),EE(1:step:end)*0+10,EE(1:step:end));
      hca.XLabel.String = 'L';
      hca.YLabel.String = 'M';
      hca.ZLabel.String = 'N';
      axis(hca,'equal')
    end
    if 0 % Pitchangle distribution
      hca = h(isub); isub = isub + 1;
      plot(hca,ePitch.depend{1}(it,:),squeeze(ePitch.data(it,:,[1 ceil(numel(ePitch.depend{2})/2) numel(ePitch.depend{2})])));
      hca.YScale = 'log'; hca.XScale = 'log';
      hca.YLabel.String = ['f_e (' ePitch.units ')'];
      hca.XLabel.String = 'E (eV)';
      hca.XLim = [10 1000];
      legend(hca,{'0','90','180'})
      hca.YTick = 10.^[0:5];
      hca.YLim = [1e0 1e5];
    end
  end


  for ii = [1:nrows*ncols]
    colormap(h(ii),strCMap)
  end
%  cn.print(['e_proj_fix_mms' num2str(ic) '_' timeUTC '_opengl'],'opengl','path',eventPath)

%% Plot various diagnostics
iiPSD = 1;
iiE = 11;
iiAz = 2;
iiPol = 10;

ileft = find(z_stop(:)<0);
iright = find(z_stop(:)>0);

nrows = 3;
ncols = 3;
for ipl = 1:nrows*ncols
  h(ipl) = subplot(nrows,ncols,ipl);
end

isub = 1;

if 1 % Position at limN
  hca = h(isub); isub = isub + 1;
  hist(hca,z_stop(:)*1e-3)
  hca.Title.String = 'Origin of e';
  hca.XLabel.String = 'N (km)';
  hca.YLabel.String = 'Number';
end
if 1 % Position at limN vs. pa inside current sheet
  hca = h(isub); isub = isub + 1;
  scatter(hca,z_stop(ileft)*1e-3,pa_start(ileft)); hold(hca,'on')
  scatter(hca,z_stop(iright)*1e-3,pa_start(iright)); hold(hca,'off')  
  hca.YLabel.String = 'Pitchangle inside (deg.)';
  hca.XLabel.String = 'N (km)';
  hca.YLim = [0 180];  
  hca.YTick = [0 45 90 135 180];
end
if 1 % PA vs PA
  hca = h(isub); isub = isub + 1;
  scatter(hca,pa_stop(ileft),pa_start(ileft)); hold(hca,'on')
  scatter(hca,pa_stop(iright),pa_start(iright)); hold(hca,'off')
  hca.YLabel.String = 'Pitchangle inside (deg.)';
  hca.XLabel.String = 'Pitchangle outside (deg.)';
  hca.XLim = [0 180];
  hca.YLim = [0 180];
  hca.XTick = [0 45 90 135 180];
  hca.YTick = [0 45 90 135 180];
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},(all_npart)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  hca.XLabel.String = 'N (km)';
  hca.YLabel.String = 'Pitchangle (deg.)';
end
if 1
  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},log10(all_psd)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  hca.XLabel.String = 'N (km)';
  hca.YLabel.String = 'Pitchangle (deg.)';
end



%% Plot different orbits to diagnostisize the method
iiPSD = 1;
iiE = 11;
iiAz = 2;
iiPol = 10;

Ps = allParticles{iPSD}.particles{iiE,iiAz,iiPol}; % particle structure


plot3(Ps.r(:,1)*1e-3,Ps.r(:,2)*1e-3,Ps.r(:,3)*1e-3); hold on
plot3(Ps.r(1,1)*1e-3,Ps.r(1,2)*1e-3,Ps.r(1,3)*1e-3,'go'); 
plot3(Ps.r(end,1)*1e-3,Ps.r(end,2)*1e-3,Ps.r(end,3)*1e-3,'rx'); 

qs = 4e-4;
quiver3(Ps.r(1,1)*1e-3,Ps.r(1,2)*1e-3,Ps.r(1,3)*1e-3,...
        Ps.v(1,1)*1e-3,Ps.v(1,2)*1e-3,Ps.v(1,3)*1e-3,...
        qs,'g'); 
quiver3(Ps.r(end,1)*1e-3,Ps.r(end,2)*1e-3,Ps.r(end,3)*1e-3,...
        Ps.v(end,1)*1e-3,Ps.v(end,2)*1e-3,Ps.v(end,3)*1e-3,...
        qs,'r'); 
hold off;
axis equal
title(sprintf('pa_{start} = %.0f, pa_{stop} = %.0f, z_{stop} = %.1f',pa_start(iiE,iiAz,iiPol),pa_stop(iiE,iiAz,iiPol),z_stop(iiE,iiAz,iiPol)*1e-3))
xlabel('L'); ylabel('M'); zlabel('N'); view([0 1 0])

%% 
ipl = find(z_stop(:)>0);
ipl = find(pa_start(:) < 100 & pa_start(:) > 90);

hca = subplot(1,1,1); hold(hca,'on')
axis(hca,'equal')
%title(sprintf('pa_{start} = %.0f, pa_{stop} = %.0f, z_{stop} = %.1f',pa_start(iiE,iiAz,iiPol),pa_stop(iiE,iiAz,iiPol),z_stop(iiE,iiAz,iiPol)*1e-3))
xlabel('L'); ylabel('M'); zlabel('N'); view([0 1 0])
for ii = 1:numel(ipl)
  Ps = allParticles{iPSD}.particles{ipl(ii)}; % particle structure
  plot3(hca,Ps.r(:,1)*1e-3,Ps.r(:,2)*1e-3,Ps.r(:,3)*1e-3); 
  %plot3(Ps.r(1,1)*1e-3,Ps.r(1,2)*1e-3,Ps.r(1,3)*1e-3,'go'); 
  %plot3(Ps.r(end,1)*1e-3,Ps.r(end,2)*1e-3,Ps.r(end,3)*1e-3,'rx'); 
pause(0.1)
%  qs = 4e-4;
%   quiver3(Ps.r(1,1)*1e-3,Ps.r(1,2)*1e-3,Ps.r(1,3)*1e-3,...
%           Ps.v(1,1)*1e-3,Ps.v(1,2)*1e-3,Ps.v(1,3)*1e-3,...
%           qs,'g'); 
%   quiver3(Ps.r(end,1)*1e-3,Ps.r(end,2)*1e-3,Ps.r(end,3)*1e-3,...
%           Ps.v(end,1)*1e-3,Ps.v(end,2)*1e-3,Ps.v(end,3)*1e-3,...
%           qs,'r'); 
end
hold(hca,'off');
axis equal
%title(sprintf('pa_{start} = %.0f, pa_{stop} = %.0f, z_{stop} = %.1f',pa_start(iiE,iiAz,iiPol),pa_stop(iiE,iiAz,iiPol),z_stop(iiE,iiAz,iiPol)*1e-3))
xlabel('L'); ylabel('M'); zlabel('N'); view([0 1 0])

%% Plot
  nrows = 9; ncols = 1;
  for isub = 1:nrows*ncols
    h(isub) = subplot(nrows,ncols,isub);
  end
  h = irf_plot(nrows);
  isub = 1;
  
  elim = [40 1000];
  
  if 1 % Magnetic field
  hca = h(isub); isub = isub + 1;
  % Colors
  B_colors = mms_colors('xyz1');

  zObs = (obsB.time.epochUnix-mean(obsB.time.epochUnix))*CS_normal_velocity;  
  set(hca,'colororder',B_colors)
  hca.ColorOrder = B_colors;
  linesObs = plot(hca,zObs,[obsB.data obsB.abs.data],'-');
  linesObs(1).Color = B_colors(1,:);
  linesObs(2).Color = B_colors(2,:);
  linesObs(3).Color = B_colors(3,:);
  linesObs(4).Color = B_colors(4,:);
  set(hca,'colororder',B_colors)
  irf_legend(hca,{'B_L','B_M','B_N'},[0.01 0.2],'fontsize',14)
  
  hold(hca,'on')
  zMod = linspace(-d*1.5,d*1.5,50)*3;
  %plot(hca,zMod*1e-3,[Bx(zMod); By(zMod); sqrt(Bx(zMod).^2 + By(zMod).^2)]*1e9)
  lineMod = plot(hca,zMod*1e-3,[Bx(0,0,zMod)]*1e9,'--','color',B_colors(1,:));
  plot(hca,zMod*1e-3,[By(0,0,zMod)]*1e9,'--','color',B_colors(2,:))
  plot(hca,zMod*1e-3,[Bz(0,0,zMod)]*1e9,'--','color',B_colors(3,:))
  plot(hca,zMod*1e-3,sqrt(Bx(0,0,zMod).^2+By(0,0,zMod).^2+Bz(0,0,zMod).^2)*1e9,'--','color',B_colors(4,:))
  hold(hca,'off')
  %legend(hca,[linesObs(1) lineMod],{'Observed data','Model fit'})
  %hca.Title.String = 'Magnetic field';
  hca.YGrid = 'on';
  hca.YLabel.String = {'B','(nT)'};
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

  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},(all_npart)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  hcb.YLabel.String = {'# unique trajectories','passing this point'};

  
  if 1 % Distance: ePDist PSD pa low energies
    hca = h(isub); isub = isub + 1;    
    plotPitch = obsPitch.elim(elim).specrec('pa');
    pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
    shading(hca,'flat')
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('11'))
    hold(hca,'off')
    %irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';   
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = plotPitch.p_label;

    %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
    %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
    hca.YLabel.String = {'Pitchangle','(\circ)'};
    hca.YTick = [45 90 135];   
    colormap(hca,'jet')
    irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
    %hca.CLim = h(2).CLim;  
    xlabel(hca,'N (km)')  
  end
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},log10(all_psd)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  %hca.CLim = h(isub-2).CLim;
  hcb.YLabel.String = {'PSD','s^3/km^6'};
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},log10(all_psd_2)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  
  if 1 % Distance: ePDist DEF pa low energies
    hca = h(isub); isub = isub + 1;    
    plotPitch = obsPitch.elim(elim).deflux.specrec('pa');
    pcolor(hca,zObsPDist,plotPitch.f,log10(plotPitch.p'))
    shading(hca,'flat')
    hold(hca,'on')
    set(hca,'ColorOrder',mms_colors('11'))
    hold(hca,'off')
    %irf_pl_mark(hca,tref,'k')
    hca.XGrid = 'off';
    hca.YGrid = 'off';   
    hcb = colorbar('peer',hca);
    hcb.YLabel.String = plotPitch.p_label;

    %hca.YLabel.String = {'\theta_{PA,e} (\circ)'};
    %ylabel(hca,{'\theta_{PA,e} (\circ)'},'interpreter','tex')
    hca.YLabel.String = {'Pitchangle','(\circ)'};
    hca.YTick = [45 90 135];   
    colormap(hca,'jet')
    irf_legend(hca,{[num2str(elim(1),'%.0f') '<E_e<' num2str(elim(2),'%.0f') ' eV']},[0.01 0.5],'fontsize',12,'color',[0 0 0]);
    %hca.CLim = h(2).CLim;
    %hca.CLim = [7.1 8];
    xlabel(hca,'N (km)')  
  end
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},log10(all_def)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  %hca.CLim = [7 11];
  
  hca = h(isub); isub = isub + 1;
  pcolor(hca,mid{1}*1e-3,mid{2},log10(all_def_2)')
  hcb = colorbar('peer',hca);
  shading(hca,'flat')
  %hca.CLim = [7 11];
  

  

  
  h1pos = h(1).Position(3);
  for ip = 1:nrows
    h(ip).Position(3) = h1pos*0.95;
    colormap(h(ip),'jet')
    h(ip).XLim = [-30 30];
  end
  %h(4).CLim = h(3).CLim;
  %h(9).CLim = [3 5.5];
  %h(5).CLim = [6 10];

