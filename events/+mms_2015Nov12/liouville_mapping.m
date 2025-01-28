% Liouville
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

%% Initialize electron test particles
particle_set = '+N'; % the particle sets to initialize: 1: +N, 2: -N
nE = 3;
nPa = 1;
nP = nE*nPa;
limN = 40e3;
T = 1; % integration t  ime
pa1 = 10;
pa2 = 10;


switch particle_set % Particle initialization
  case '+N' % energy and pitch angle, % starting from +N or -N, logspaced in energys        
    pa0_ = 180;
    z0_ = 30;
    pa0 = pa0_ - linspace(pa1,pa2,nPa);
    z0 = zeros(1,nP) + z0_;
  case '-N'
    pa0_ = 0;
    z0_ = -30;
    pa0 = pa0_ + linspace(pa1,pa2,nPa);
    z0 = zeros(1,nP) + z0_;
  case '+-N'
    pa0_minus = 0;
    pa0_plus = 180;
    z0_minus = -30;    
    z0_plus = 30;    
    z0_ = max([z0_plus z0_minus]);
    
    pa0_minus_ = pa0_minus + linspace(pa1,pa2,nPa);
    pa0_plus_ =  pa0_plus - linspace(pa1,pa2,nPa);    
    pa0 = [pa0_plus_ pa0_minus_];
    
    z0_minus_ = zeros(1,nP) + z0_minus;
    z0_plus_ = zeros(1,nP) + z0_plus;
    z0 = [z0_plus_ z0_minus_];
    
    nPa = numel(pa0);
end

nP = nE*nPa;

E1 = 10; E2 = 300; % eV
vt1 = sqrt(E1*units.eV*2/units.me)/1000; % km/s
vt2 = sqrt(E2*units.eV*2/units.me)/1000;
if 1 % equi-spaced v's
  vt = linspace(vt1,vt2,nE); % eV   
  dvt = diff(vt);
  vt_edges = [vt-0.5*dvt(1) vt(end)+0.5*dvt(1)];  
  dvt = diff(vt_edges);
  electron_energy = (vt*1000).^2*units.me/units.eV; % eV   
  E_edges = (vt_edges*1000).^2*units.me/units.eV; % eV 
  dE = diff(E_edges);
else % equi-spaced E's
  %electron_energy = linspace(log10(E1),log10(E2),nE); % eV  
  electron_energy = linspace(E1,E2,nE); % eV    
  dE = diff(electron_energy);
  E_edges = [electron_energy-0.5*dE(1) + electron_energy(end)+0.5*dE(1)];
  dE = diff(E_edges);
  vt = sqrt(electron_energy*units.eV*2/units.me)/1000; % km/s
  vt_edges = sqrt(E_edges*units.eV*2/units.me)/1000; % km/s
end  



% Initial positions and velocities
x0 = zeros(1,nP); % km
y0 = zeros(1,nP); 
%z0 = zeros(1,nP) + z0_;

% Start with field aligned set up
%pa0 = pa0_ + 30*randn(1,nPa);
%az0 = 360*rand(1,nPa);
AZ0 = 360*rand(nE,nPa);
%thetaBL = acosd(b0*[1 0 0]');

[PA0,VT0] = meshgrid(pa0,vt);
[PA0,E0] = meshgrid(pa0,electron_energy);
dEE = repmat(tocolumn(dE),1,nPa);
dVV = repmat(tocolumn(dvt),1,nPa);


dEE_col = reshape(dEE,nE*nPa,1);
dVV_col = reshape(dVV,nE*nPa,1);

% B/parperp coordinate system
% vx0_pa = vt.*cosd(pa0); % km/s
% vy0_pa = vt.*cosd(az0).*sind(pa0);
% vz0_pa = vt.*sind(az0).*sind(pa0);    
% v0_pa = [vx0_pa' vy0_pa' vz0_pa']';
VX0_pa = VT0.*cosd(PA0); % km/s
VY0_pa = VT0.*cosd(AZ0).*sind(PA0);
VZ0_pa = VT0.*sind(AZ0).*sind(PA0);

vx0_pa = reshape(VX0_pa,nE*nPa,1);
vy0_pa = reshape(VY0_pa,nE*nPa,1);
vz0_pa = reshape(VZ0_pa,nE*nPa,1);

v0_pa = [vx0_pa vy0_pa vz0_pa]';
v0 = v0_pa*0;

% Transform to LMN system, % A>B, v' = B(A^-1)v  
rI = [1 0 0; 0 1 0; 0 0 1]; % unit matrix
rB = [1 0 0; 0 1 0; 0 0 1]; % LMN (xyz), final coordinate system

z0s = unique(z0);
for iz0 = 1:numel(z0s)
  %[~,zindA,zindB] = intersect(z0s,z0);  
  zind = find(z0 == z0s(iz0));
  B0 = [Bx(x0(1)*1e3,y0(1)*1e3,z0s(iz0)*1e3); By(x0(1)*1e3,y0(1)*1e3,z0s(iz0)*1e3); Bz(x0(1)*1e3,y0(1)*1e3,z0s(iz0)*1e3)]; % B0 = B0(:,1)'; % NOT TRUE: same for all particles since they have same starting position
  b0 = irf_norm(B0');
  b0_perp1 = cross(b0,[1 0 0]); b0_perp1 = b0_perp1/norm(b0_perp1);
  b0_perp2 = cross(b0,b0_perp1); b0_perp2 = b0_perp2/norm(b0_perp2);

  rA = [b0; b0_perp1; b0_perp2]; % FAC, initial coordinate system
  v0(:,zind) = rB*rA^-1*v0_pa(:,zind);
end

vx0 = v0(1,:);
vy0 = v0(2,:);
vz0 = v0(3,:);

VX0 = reshape(vx0,nE,nPa);
VY0 = reshape(vy0,nE,nPa);
VZ0 = reshape(vz0,nE,nPa);

X0 = reshape(x0,nE,nPa);
Y0 = reshape(y0,nE,nPa);
Z0 = reshape(z0,nE,nPa);

x_init_all = [x0;y0;z0;vx0;vy0;vz0]'*1e3; % m, m/s


nParticles = nP;
iDist = 1;

% Assign phase space density to each particle
z0s = unique(z0);
c_eval('t? = tintObs(2) +- (tintObs(2)-tintObs(1))*0.5 + z0s(?)/CS_normal_velocity;',1:numel(z0s))
c_eval('t?ind = find(abs(obsPitch.time-t?)==min(abs(obsPitch.time-t?)));',1:numel(z0s))
c_eval('f0? = obsPitch(t?ind);',1:numel(z0s))

c_eval('obsPAbin! = f0?.depend{2}; obsPAbin! = [0 obsPAbin! + obsPAbin!(1)]; ?;',ic,1:numel(z0s))
c_eval('obsEbin! = f0?.depend{1}; obsEbin! = [0 obsEbin! + obsEbin!(1)]; ? ;',ic,1:numel(z0s))
F0 = nan(nE,nPa);
for iPa = 1:nPa
  [~,zind]=intersect(z0s,Z0(1,iPa));
  c_eval('obsEbin = obsEbin?; obsPAbin = obsPAbin?; f0 = f0?;',zind)    
  for iE = 1:nE    
    iBinE = find(E0(iE,iPa)<obsEbin,1,'first');
    iBinPa = find(PA0(iE,iPa)>obsPAbin,1,'last');
    F0(iE,iPa) = f0.data(1,iBinE,iBinPa);
  end
end

f0mod = reshape(F0,nE*nPa,1);

if 0 % plot initial phase space density and DEF
  %%
  hca = subplot(3,1,1); h(1) = hca;
  pcolor(hca,f0.depend{2},f0.depend{1},log10(squeeze(f0.data)))
  hca.YScale = 'log';
  hca.XLabel.String = 'Pitchangle (deg.)';
  hca.YLabel.String = 'Energy (eV)';
  
  
  hca = subplot(3,1,2); h(2) = hca;
  pcolor(hca,PA0,E0,log10(F0))
  hca.YScale = 'log';
  hca.XLabel.String = 'Pitchangle (deg.)';
  hca.YLabel.String = 'Energy (eV)';
 
  hca = subplot(3,1,3); h(3) = hca;
  pcolor(hca,PA0,E0,log10(E0));
  hca.YScale = 'log';
  hca.XLabel.String = 'Pitchangle (deg.)';
  hca.YLabel.String = 'Energy (eV)';
  
  xlim = [0 180];
  ylim = [10 1000];
  for ip = 1:3;
    h(ip).XLim = xlim;
    h(ip).YLim = ylim;  
    colorbar('peer',h(ip))
  end
  clim = [0 5];
  for ip = 1:2
    h(ip).CLim = clim;     
  end
end

%% Integration
tic

x_sol_all = [];
saveParticle = cell(1,nParticles);
for iParticle= 1:nParticles 

  % Initial positions and velocities                                   
  x_init = x_init_all(iParticle,:); % m, m/s
  %x_init = x_init(:,iParticle); 
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
  saveParticle{iParticle}.f = f0mod(iParticle);
  saveParticle{iParticle}.dE = dEE_col(iParticle);
  saveParticle{iParticle}.dv = dVV_col(iParticle);
end
toc
%iP=1;plot(saveParticle{iP}.r(:,3),saveParticle{iP}.energy); hold on; for iP=1:10:nParticles, plot(saveParticle{iP}.r(:,3),saveParticle{iP}.energy); end

%% Bin phase space density of particles to get 'probability density' of pitchangles (Liouville's theorem)
mm = units.me/units.mp;
all_T = [];
all_z = [];
all_pa = [];
all_pa0 = [];
all_energy = [];
all_zstop = [];
all_zstart = [];
edges_z = -1-abs(z0_):1:abs(z0_)+1;
edges_pa = -0:10:180;
all_psd = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_psd_2 = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_def = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_npart = zeros(numel(edges_z)-1,numel(edges_pa)-1);
all_vvol = zeros(numel(edges_z)-1,numel(edges_pa)-1);
dvt_all = 0;

for iP = 1:numel(saveParticle) % Electron pitchangles
  % Pick out the data
  ind = 1:numel(saveParticle{iP}.pa);
  z = saveParticle{iP}.r(ind,3); % N
  pa = saveParticle{iP}.pa(ind); % pitch angle
  v = saveParticle{iP}.v(ind,:);    
  energy = saveParticle{iP}.energy;
  all_T = [all_T; saveParticle{iP}.T];
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
  all_psd_2 = (all_psd_2.*all_vvol_previous + bins_occupied*saveParticle{iP}.f*this_vol)./all_vvol;

  % Add phase space density to N/PA grid
  all_psd = all_psd + bins_occupied*saveParticle{iP}.f;%/saveParticle{iP}.dv;   
  %all_psd = all_psd.*all_npart + bins_occupied*saveParticle{iP}.f;%/saveParticle{iP}.dv;   
  %all_psd = all_psd./(all_npart+bins_occupied);
  % Add DEF to N/PA grid, more complicated since energy changes    
  all_def = all_def + bins_occupied*energy(1).^2;%/saveParticle{iP}.dE;
end

all_def = all_def;%./all_npart;%/nPa;
all_psd = all_psd;%./all_npart;%/nPa;

all_def_2 = all_def.*all_psd;%./all_npart;%/nPa;
%all_def = all_psd.*all_def/1e6/mm^2/0.53707;
all_def(all_def == 0) = NaN;
all_def_2(all_def_2 == 0) = NaN;


% Plot
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

