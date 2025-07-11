if 0
%%
B0 = 20e-9;

tA = 0.1;
fBx = @(z,L) B0*tanh(z/L);
fAy  = @(z,L) L*B0*log(cosh(z/L)); % + t/tA;
fPhi  = @(z) z*0;

fpx = @(q,v,z,L,m) m*vx;
fpy = @(q,v,z,L,m) m*vy + q*fAy(z,L);
fpz = @(q,v,z,L,m) m*vz;

% H = m*vx^2/2 + m*vy^2/2 + m*vz^2/2 + q*Phi
% vx = px/m
% vy = py/m - qA/m
% vz = pz/m
% H = px^2/2m + m*(py/m - qA/m)^2/2 + pz^2/2m + p*Phi
fPA = @(q,vy,vz,z,L,m) 0.5*m*(fpy(q,vy,z,L,m)/m).^2 + ...
                       0.5*m*(fpz(q,vz,z,L,m)/m - q*fAy(z,L)/m).^2;
                     
                      
fPA = @(q,py,pz,z,L,m) 0.5*m*(pz/m).^2 + ...
                       0.5*m*(py/m - q*fAy(z,L)/m).^2;


% fH = @(q,vx,vz,vy,z,L,m) m*(fpx(q,vx,z,L,m)/m).^2 + ... 
%                          m*fPA(q,vy,vz,z,L,m) + q*fPhi(z);
   
units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

zv = 4*L*linspace(-1,1,200);
dz = zv(2)-zv(1);
z0 = 2*L;
vy0 = 10000e3;
vz0 = 10000e3;
py0 = m*vy0 - q*fAy(z0,L);
pz0 = m*vz0; 

%py0=0;

h = setup_subplots(2,1);
isub = 1;
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,zv/L,fAy(zv,L))
  hca.YLabel.String = 'A_y';
  hca.XLabel.String = 'z/L';
end
if 0
  hca = h(isub); isub = isub + 1;
  plot(hca,zv(2:end)/L,diff(fAy(zv,L))/dz)
  hca.YLabel.String = 'diff A_y/dz';
  hca.XLabel.String = 'z/L';
end
if 1
  hca = h(isub); isub = isub + 1;
  plot(hca,zv/L,fPA(q,py0,pz0,zv,L,m))
  hca.XLabel.String = 'z/L';
end


%% integrate electron orbits for a given magnetic field: several particles
% Set up p  lot
fig = figure(33);
%fig.Position(4) = 1e3*0.5;
%fig.Position(3) = 1e3*0.3;
nCols = 1;
nRows = 4;

clear h
h(1) = subplot(nRows,nCols,1);
h(2) = subplot(nRows,nCols,[2:3]);
h(3) = subplot(nRows,nCols,4);

units = irf_units;

% Parameters
Er = 0e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 1*5*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 11000*units.e;
Lphi = L;

% Harris current sheet
Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z) x*0 + y*0 + z*0;
Ax = @(z) Bg*z;
Ay = @(z) B0*d*log(cosh(z/d));
Az = @(z) 0*z;%Bg*z;
Phi = @(z) phi0*exp(-(z/Lphi).^2);

%Gamma = @(z,py,px,phi0,d) 0.5/m*((py-q*Ay(z)).^2+(pz-q*Az(z)).^2);
Gamma = @(z,py,px,phi0,d) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2);
Gamma = @(z,py,px,phi0,d) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2) + q*Phi(z);

% Plot magnetic field
hca = h(1);
zz = linspace(-5*d,5*d,100);
plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)].^2,1))*1e9,'k');
%hca.YTick = 0;
hca.XLabel.String = 'z';
hca.YLabel.String = 'B';
irf_legend(h(1),{'B_x','B_y','B_z'},[0.98 0.98])
if Bg == 0
  hca.Title.String = 'No guide field';
else
  hca.Title.String = 'Finite guide field';
end
hca.Title.Position(2) = B0*1e9 + 2;
  
%irf_legend(hca,'B',[0.01 1.05],'fontsize',14)
%irf_legend(hca,'B',[0.01 1.05],'fontsize',14)

% Initiate particles
particles = {};
iParticle = 0;
if 1 % no drift
  iParticle = iParticle + 1;
  particle.T = 0.02;
  particle.energy = 50;
  particle.z0 = -15; 
  particle.y0 = 0;
  particle.velocity_angle = 00;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end
if 1 % meandering
  iParticle = iParticle + 1;
  particle.T = 0.01;
  particle.energy = 150;
  particle.z0 = 5;  
  particle.y0 = 10;
  particle.velocity_angle = 70;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end
if 1 % meandering 8's
  iParticle = iParticle + 1;
  particle.T = 0.02;
  particle.energy = 70;
  particle.z0 = -1;  
  particle.y0 = -10;
  particle.velocity_angle = 160;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end
if 1 % grad B
  iParticle = iParticle + 1;
  particle.T = 0.02;
  particle.energy = 80;
  particle.y0 = 0;
  particle.z0 = 8;  
  particle.velocity_angle = -70;
  particle.m = units.me;
  particle.q = -units.e;
  particles{iParticle} = particle;
end

for iParticle = 1:numel(particles)
  % Load particle data
  particle = particles{iParticle};
  electron_energy = particle.energy;
  velocity_angle = particle.velocity_angle;
  q = particle.q;
  m = particle.m;
  T = particle.T;
  
  vt = sqrt(electron_energy*units.eV*2/m)/1000;
  % Initial positions and velocitites
  x0 = 0*1e3; % m
  y0 = particle.y0*1e3; 
  z0 = particle.z0*1e3;
  vx0 = 0*1e3; % m/s
  vy0 = vt*cosd(velocity_angle)*1e3;
  vz0 = -vt*sind(velocity_angle)*1e3;

  % Generalized potential
  px = m*vx0 + q*Ax(z0);
  pz = m*vz0 + q*Az(z0);
  py = m*vy0 + q*Ay(z0);
  Gamma0 = Gamma(0,py,px,phi0,d);
  
  % Integrate particle orbit
  x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
  EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
  [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
  x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
  vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
  particles{iParticle}.x = x;
  particles{iParticle}.y = y;
  particles{iParticle}.z = z;
  particles{iParticle}.vx = vx;
  particles{iParticle}.vy = vy;
  particles{iParticle}.vz = vz;

  hca = h(2);
  hold(hca,'on')
  horb = plot(hca,z*1e-3,y*1e-3,...
                  z(1)*1e-3,y(1)*1e-3,'go',...
                  z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
  xlabel(hca,'z')
  ylabel(hca,'y') 
  hold(hca,'off')
  
  % Plot generalized potential
  if 1 % along orbit
    hca = h(3);
    hold(hca,'on')
        
    vecGamma = Gamma(zz,py,px,phi0,d);
    hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color.^0.3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')
    
    vecGamma = Gamma(particles{iParticle}.z,py,px,phi0,d);
    hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    
    
%     limGamma = vecGamma; 
%     zzGamma = zz;
%     zzGamma(limGamma>electron_energy*units.eV) = NaN;
%     limGamma(limGamma>electron_energy*units.eV) = NaN;
%     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
%     %hpot.Color = ;
    hca.XLabel.String = 'z';
    hold(hca,'off')
  end
    if 1 % patch in orbitrange
    hca = h(3);
    hold(hca,'on')
        
    vecGamma = Gamma(zz,py,px,phi0,d);
    hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color.^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')
    
    vecGamma = Gamma(particles{iParticle}.z,py,px,phi0,d);
    hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    
    zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
    vecGamma = Gamma(zz_patch,py,px,phi0,d);
    
    hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,horb(1).Color);
    hp.FaceAlpha = 0.5;
    
%     limGamma = vecGamma; 
%     zzGamma = zz;
%     zzGamma(limGamma>electron_energy*units.eV) = NaN;
%     limGamma(limGamma>electron_energy*units.eV) = NaN;
%     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
%     %hpot.Color = ;
    hca.XLabel.String = 'z';
    hold(hca,'off')
  end
  if 0
    hca = h(3);
    hold(hca,'on')
    vecGamma = Gamma(zz,py,px,phi0,d);
    hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    limGamma = vecGamma; 
    zzGamma = zz;
    zzGamma(limGamma>electron_energy*units.eV) = NaN;
    limGamma(limGamma>electron_energy*units.eV) = NaN;
    hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    %hpot.Color = ;
    hca.XLabel.String = 'z';
    hold(hca,'off')
  end
end
hold(h(2),'off')
%axis(h(2),'equal')
zlim = 16;
for iPanel = 1:numel(h) 
  h(iPanel).XLim = zlim*[-1 1];
end
for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  %h(ii).Box = 'off';
  h(ii).FontSize = 14;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
hold(h(2),'on')
plot(h(2),[0 0],h(2).YLim,'k')
h(2).XLim = zlim*[-1 1];
hold(h(2),'off')
h(3).XLim = zlim*[-1 1];
%h(3).YLim = [0 2]*1e-16;
h(3).YLim = [0 500];
h(3).YScale = 'lin';
hold(h(3),'off')

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))

%% integrate electron orbits for a given magnetic field: several particles
units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;
colors = [0 0 1; 0.5 0.5 1; 0 0 0.5; 0 0 0.2];
colors = pic_colors('matlab');

% Set up plot
fig = figure(34);
%fig.Position(4) = 1e3*0.5;
%fig.Position(3) = 1e3*0.3;
nCols = 1;
nRows = 4;

clear h
h = setup_subplots(3,3,'horizontal');
isub = 1;

units = irf_units;

if 0 % No guide field
  % Parameters
  Er = 0e-3; % reconnection electric field, V/m
  B0 = 10e-9; % asymptotical magnetic field, T
  d = 5e3; % thickness of current sheet, m
  Bg = 1*5*1e-9; % guide field, T
  Bn = 0.0*1e-9; % normal field, T
  phi0 = 11000*units.e;
  Lphi = L;
  
  % Harris current sheet
  Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
  By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
  Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
  Ex = @(x,y,z) x*0 + y*0 + z*0;
  Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
  Ez = @(x,y,z) x*0 + y*0 + z*0;
  Ax = @(z) Bg*z;
  Ay = @(z) B0*d*log(cosh(z/d));
  Az = @(z) 0*z;%Bg*z;
  Phi = @(z) phi0*exp(-(z/Lphi).^2);
  
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2);
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2) + q*Phi(z);
  
  % Initiate particles
  particles = {};
  iParticle = 0;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = -15; 
    particle.y0 = 0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.01;
    particle.energy = 150;
    particle.z0 = 5;  
    particle.y0 = 10;
    particle.velocity_angle = 70;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 70;
    particle.z0 = -1;  
    particle.y0 = -10;
    particle.velocity_angle = 160;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 80;
    particle.y0 = 0;
    particle.z0 = 8;  
    particle.velocity_angle = -70;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  
  % Integrate particle trajectories
  for iParticle = 1:numel(particles)
    % Load particle data
    particle = particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0);
    Gamma0 = Gamma(0,py,px);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;
      
  end
  
  % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,100);
  plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)].^2,1))*1e9,'k');
  %hca.YTick = 0;
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'B';
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  if Bg == 0
    hca.Title.String = 'No guide field';
  else
    hca.Title.String = 'Finite guide field';
  end
  hca.Title.Position(2) = B0*1e9 + 2;

  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)
  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)

  % Plot generalized potential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 1 % along orbit
      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')

      if iParticle == 1, hold(hca,'on'); end

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      hca.XLabel.String = 'z';      
    end
    if 1 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    if 0      
      hold(hca,'on')
      vecGamma = Gamma(zz,py,px,phi0,d);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      limGamma = vecGamma; 
      zzGamma = zz;
      zzGamma(limGamma>electron_energy*units.eV) = NaN;
      limGamma(limGamma>electron_energy*units.eV) = NaN;
      hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      %hpot.Color = ;
      hca.XLabel.String = 'z';
      hold(hca,'off')
    end
  end
  hold(hca,'off')
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:))
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z')
    ylabel(hca,'y')     
    %if ip == 1, hold(hca,'on'); end
  end
  hold(hca,'off')
  %axis(h(2),'equal')
  zlim = 16;
  for iPanel = 1:numel(h) 
    h(iPanel).XLim = zlim*[-1 1];
  end
  for ii = 1:numel(h) 
    %hca = subplot(nRows,nCols,ii);
    %h(ii).Box = 'off';
    h(ii).FontSize = 14;
    axis(h(ii),'tight');
    %axis(h(ii),'off');
    %h(ii).Title.Position(2) = -0.2;
    %h(ii).Position(2) = 0.2;
    %h(ii).Position(4) = 0.7;
  end
  %hold(h(2),'on')
  %plot(h(2),[0 0],h(2).YLim,'k')
  %h(2).XLim = zlim*[-1 1];
  %hold(h(2),'off')
  %h(3).XLim = zlim*[-1 1];
  %h(3).YLim = [0 2]*1e-16;
  %h(3).YLim = [0 500];
  %h(3).YScale = 'lin';
  %hold(h(3),'off')
end

if 1 % No guide field, same starting position
  str_title = 'No guide field, same starting position and energy but different phase';
  % Parameters
  Er = 0e-3; % reconnection electric field, V/m
  B0 = 10e-9; % asymptotical magnetic field, T
  d = 5e3; % thickness of current sheet, m
  Bg = 0*5*1e-9; % guide field, T
  Bn = 0.0*1e-9; % normal field, T
  phi0 = 11000*units.e;
  Lphi = L;
  
  % Harris current sheet
  Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
  By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
  Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
  Ex = @(x,y,z) x*0 + y*0 + z*0;
  Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
  Ez = @(x,y,z) x*0 + y*0 + z*0;
  Ax = @(z) Bg*z;
  Ay = @(z) B0*d*log(cosh(z/d));
  Az = @(z) 0*z;%Bg*z;
  Phi = @(z) phi0*exp(-(z/Lphi).^2);
  
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2);
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2) + q*Phi(z);
  
  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 300;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 180;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 270;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  
  % Integrate particle trajectories
  for iParticle = 1:numel(particles)
    % Load particle data
    particle = particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0);
    Gamma0 = Gamma(0,py,px);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;
      
  end
  
  % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,100);
  plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)].^2,1))*1e9,'k');
  %hca.YTick = 0;
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'B';
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  %hca.Title.Position(2) = B0*1e9 + 2;

  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)
  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)

  % Plot generalized potential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 1 % along orbit
      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')

      if iParticle == 1, hold(hca,'on'); end

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      hca.XLabel.String = 'z';      
    end
    if 1 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    if 0      
      hold(hca,'on')
      vecGamma = Gamma(zz,py,px,phi0,d);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      limGamma = vecGamma; 
      zzGamma = zz;
      zzGamma(limGamma>electron_energy*units.eV) = NaN;
      limGamma(limGamma>electron_energy*units.eV) = NaN;
      hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      %hpot.Color = ;
      hca.XLabel.String = 'z';
      hold(hca,'off')
    end
  end
  hold(hca,'off')
  
  hca.Title.String = str_title;
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z')
    ylabel(hca,'y')     
    %if ip == 1, hold(hca,'on'); end
  end
  hold(hca,'off')
  %axis(h(2),'equal')
  zlim = 16;
  for iPanel = 1:numel(h) 
    h(iPanel).XLim = zlim*[-1 1];
  end
  for ii = 1:numel(h) 
    %hca = subplot(nRows,nCols,ii);
    %h(ii).Box = 'off';
    h(ii).FontSize = 14;
    axis(h(ii),'tight');
    %axis(h(ii),'off');
    %h(ii).Title.Position(2) = -0.2;
    %h(ii).Position(2) = 0.2;
    %h(ii).Position(4) = 0.7;
  end
  %hold(h(2),'on')
  %plot(h(2),[0 0],h(2).YLim,'k')
  %h(2).XLim = zlim*[-1 1];
  %hold(h(2),'off')
  %h(3).XLim = zlim*[-1 1];
  %h(3).YLim = [0 2]*1e-16;
  %h(3).YLim = [0 500];
  %h(3).YScale = 'lin';
  %hold(h(3),'off')
end

if 1 % No guide field, same starting position
  str_title = 'No guide field, same starting position and phase but different energies';
  % Parameters
  Er = 0e-3; % reconnection electric field, V/m
  B0 = 10e-9; % asymptotical magnetic field, T
  d = 5e3; % thickness of current sheet, m
  Bg = 0*5*1e-9; % guide field, T
  Bn = 0.0*1e-9; % normal field, T
  phi0 = 11000*units.e;
  Lphi = L;
  
  % Harris current sheet
  Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
  By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
  Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
  Ex = @(x,y,z) x*0 + y*0 + z*0;
  Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
  Ez = @(x,y,z) x*0 + y*0 + z*0;
  Ax = @(z) Bg*z;
  Ay = @(z) B0*d*log(cosh(z/d));
  Az = @(z) 0*z;%Bg*z;
  Phi = @(z) phi0*exp(-(z/Lphi).^2);
  
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2);
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2) + q*Phi(z);
  
  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = -90;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  
  % Integrate particle trajectories
  for iParticle = 1:numel(particles)
    % Load particle data
    particle = particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0);
    Gamma0 = Gamma(0,py,px);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;
      
  end
  
  % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,100);
  plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)].^2,1))*1e9,'k');
  %hca.YTick = 0;
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'B';
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  %hca.Title.Position(2) = B0*1e9 + 2;

  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)
  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)

  % Plot generalized potential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 1 % along orbit
      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')

      if iParticle == 1, hold(hca,'on'); end

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      hca.XLabel.String = 'z';      
    end
    if 1 % line in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      %hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    if 0 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    if 0      
      hold(hca,'on')
      vecGamma = Gamma(zz,py,px,phi0,d);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      limGamma = vecGamma; 
      zzGamma = zz;
      zzGamma(limGamma>electron_energy*units.eV) = NaN;
      limGamma(limGamma>electron_energy*units.eV) = NaN;
      hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      %hpot.Color = ;
      hca.XLabel.String = 'z';
      hold(hca,'off')
    end
  end
  hold(hca,'off')
  
  hca.Title.String = str_title;
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z')
    ylabel(hca,'y')     
    %if ip == 1, hold(hca,'on'); end
  end
  hold(hca,'off')
  %axis(h(2),'equal')
  zlim = 16;
  for iPanel = 1:numel(h) 
    h(iPanel).XLim = zlim*[-1 1];
  end
  for ii = 1:numel(h) 
    %hca = subplot(nRows,nCols,ii);
    %h(ii).Box = 'off';
    h(ii).FontSize = 14;
    axis(h(ii),'tight');
    %axis(h(ii),'off');
    %h(ii).Title.Position(2) = -0.2;
    %h(ii).Position(2) = 0.2;
    %h(ii).Position(4) = 0.7;
  end
  %hold(h(2),'on')
  %plot(h(2),[0 0],h(2).YLim,'k')
  %h(2).XLim = zlim*[-1 1];
  %hold(h(2),'off')
  %h(3).XLim = zlim*[-1 1];
  %h(3).YLim = [0 2]*1e-16;
  %h(3).YLim = [0 500];
  %h(3).YScale = 'lin';
  %hold(h(3),'off')
end

if 1 % Finite guide field, same starting position
  str_title = 'No guide field, same starting position and phase but different energies';
  % Parameters
  Er = 0e-3; % reconnection electric field, V/m
  B0 = 10e-9; % asymptotical magnetic field, T
  d = 5e3; % thickness of current sheet, m
  Bg = 1*3*1e-9; % guide field, T
  Bn = 0.0*1e-9; % normal field, T
  phi0 = 11000*units.e;
  Lphi = L;
  
  % Harris current sheet
  Bx = @(x,y,z) x*0 + y*0 - B0*tanh(z/d);
  By = @(x,y,z) x*0 + y*0 + z*0 + Bg;
  Bz = @(x,y,z) x*0 + y*0 + z*0 + Bn;
  Ex = @(x,y,z) x*0 + y*0 + z*0;
  Ey = @(x,y,z) x*0 + y*0 + z*0 + Er;
  Ez = @(x,y,z) x*0 + y*0 + z*0;
  Ax = @(z) Bg*z;
  Ay = @(z) B0*d*log(cosh(z/d)); % f(x) = Bn*x -> Bz
  Az = @(z) 0*z;%Bg*z;
  Phi = @(z) phi0*exp(-(z/Lphi).^2);
  
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2);
  Gamma = @(z,py,px) 0.5/m*((py-q*Ay(z)).^2+(px-q*Ax(z)).^2) + q*Phi(z);
  
  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = -90;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  
  % Integrate particle trajectories
  for iParticle = 1:numel(particles)
    % Load particle data
    particle = particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0);
    Gamma0 = Gamma(0,py,px);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,Bx,By,Bz,Ex,Ey,Ez);
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;
      
  end
  
  % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,100);
  plot(hca,zz*1e-3,[Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz); By(0,0,zz); Bz(0,0,zz)].^2,1))*1e9,'k');
  %hca.YTick = 0;
  hca.XLabel.String = 'z';
  hca.YLabel.String = 'B';
  irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  %hca.Title.Position(2) = B0*1e9 + 2;

  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)
  %irf_legend(hca,'B',[0.01 1.05],'fontsize',14)

  % Plot generalized potential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 1 % along orbit
      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')

      if iParticle == 1, hold(hca,'on'); end

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      hca.XLabel.String = 'z';      
    end
    if 1 % line in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      %hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    if 0 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    if 0      
      hold(hca,'on')
      vecGamma = Gamma(zz,py,px,phi0,d);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      limGamma = vecGamma; 
      zzGamma = zz;
      zzGamma(limGamma>electron_energy*units.eV) = NaN;
      limGamma(limGamma>electron_energy*units.eV) = NaN;
      hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
      %hpot.Color = ;
      hca.XLabel.String = 'z';
      hold(hca,'off')
    end
  end
  hold(hca,'off')
  
  hca.Title.String = str_title;
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z')
    ylabel(hca,'y')     
    %if ip == 1, hold(hca,'on'); end
  end
  hold(hca,'off')
  %axis(h(2),'equal')
  zlim = 16;
  for iPanel = 1:numel(h) 
    h(iPanel).XLim = zlim*[-1 1];
  end
  for ii = 1:numel(h) 
    %hca = subplot(nRows,nCols,ii);
    %h(ii).Box = 'off';
    h(ii).FontSize = 14;
    axis(h(ii),'tight');
    %axis(h(ii),'off');
    %h(ii).Title.Position(2) = -0.2;
    %h(ii).Position(2) = 0.2;
    %h(ii).Position(4) = 0.7;
  end
  %hold(h(2),'on')
  %plot(h(2),[0 0],h(2).YLim,'k')
  %h(2).XLim = zlim*[-1 1];
  %hold(h(2),'off')
  %h(3).XLim = zlim*[-1 1];
  %h(3).YLim = [0 2]*1e-16;
  %h(3).YLim = [0 500];
  %h(3).YScale = 'lin';
  %hold(h(3),'off')
end



c_eval('h(?).XLim = [-20 20];',1:numel(h))

hl = findobj(gcf,'type','line');

c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).Visible = ''off'';',1:numel(h))

%% integrate electron orbits for a given magnetic field: several particles, common plotting with setup before

units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

% General current sheet structure:
% Harris current sheet, incl. parameters in functions to make more
% versatile later on
Bx = @(x,y,z,B0,d) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z,Bg) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z,Bn) x*0 + y*0 + z*0 + Bn;
%Ex = @(x,y,z) x*0 + y*0 + z*0;
%Ey = @(x,y,z,Er) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z,phi0,Lphi) x*0 + y*0 + z*0 -phi0/Lphi;
Ax = @(z,Bg) Bg*z;
Ay = @(z,B0,d) B0*d*log(cosh(z/d));
Az = @(z) 0*z;
Phi = @(z,phi0,Lphi) phi0*exp(-(z/Lphi).^2);
Phi = @(z,phi0,Lphi) phi0*z/Lphi;

% Pseudo-potential
Gamma = @(z,py,px,B0,Bg,d,phi0,Lphi) 0.5/m*((py-q*Ay(z,B0,d)).^2+(px-q*Ax(z,Bg)).^2) + q*Phi(z,phi0,Lphi);
  
% What colors to use for lines
colors = [0 0 1; 0.5 0.5 1; 0 0 0.5; 0 0 0.2];
colors = pic_colors('matlab');
fontsize = 12;

% Define our different current sheets and particles
ics = 0;
cs = struct([]);

if 0 % No guide field, different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 300;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 180;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 270;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % No guide field, slightly different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 50;
  T = 0.02;
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -120;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -70;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -20;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % Finite guide field
  ics = ics + 1;  
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 1.6e-9; % guide field, T
  cs(ics).Bn = 0.0*1e-9; % normal field, T
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  cs(ics).phi0 = 0; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'Weak guide field, same starting position and phase but different energies';

    % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = -90;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end  
  cs(ics).particles = particles;
end
if 1 % No guide field, finite constant En
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m    
  cs(ics).phi0 = 22000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, finite normal electric field';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = 90;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end

% Set up plot
fig = figure(35);
nCols = 3;
nRows = 3;
h = setup_subplots(nCols,nRows,'horizontal');
isub = 1;

% Loop through different settings
for ics = 1:numel(cs)  
  % Extract parameters for ease of writing below
  B0 = cs(ics).B0; % asymptotical magnetic field, T
  Bg = cs(ics).Bg; % guide field, T
  Bn = cs(ics).Bn; % normal field, T
  %Er = cs(ics).Er = 0e-3; % reconnection electric field, V/m
  d = cs(ics).d; % thickness of current sheet, m
  phi0 = cs(ics).phi0; % potential
  Lphi = cs(ics).Lphi; % length scale of potential
  
  % Integrate particle trajectories
  particles = cell(numel(cs(ics).particles),1);
  for iParticle = 1:numel(cs(ics).particles)
    % Load particle data
    particle = cs(ics).particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0,Bg);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0,B0,d);
    Gamma0 = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,@(x,y,z)Bx(x,y,z,B0,d),@(x,y,z)By(x,y,z,Bg),@(x,y,z)Bz(x,y,z,Bn),Ex,Ey,@(x,y,z)Ez(x,y,z,phi0,Lphi));
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;      
  end
  
   % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,200);
  if 1 % Bx, By, Ez
    set(hca,'ColorOrder',[colors(1:2,:); 0 0 0; colors(3,:)])
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k',zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    set(hca,'ColorOrder',[colors(1:2,:); 0 0 0; colors(3,:)])
    %axes(hca)
    %yyaxis right;
    %ax2 = gca;
    %ax2.Position = hca.Position;
    %ax2.
%     ax2 = axes('XAxisLocation','top',...
%          'YAxisLocation','right',...
%          'Color','none',...
%          'XColor','k','YColor','k',...
%          'Position',hca.Position);
    %plot(ax2,zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    %ax2.YLabel.String = 'E_z (mV/m)';
    irf_legend(hca,{'B_x','B_y','|B|','E_z'},[0.02 0.08],'fontsize',fontsize)
    hca.YLabel.String = 'B (nT), E (mV/m)';
  else
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k');
    irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  hca.YLabel.String = 'B (nT)';
  end
  %hca.YTick = 0;
  hca.XLabel.String = 'z (km)';
  

  % Plot generalized spotential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 0 % along orbit
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')
      if iParticle == 1, hold(hca,'on'); end
      %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
      %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')    
    end
    if 1 % line in orbitrange
      if iParticle == 1, hold(hca,'on'); end
      
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
      %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px,B0,Bg,d,phi0,Lphi);

      %hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      %hp.FaceAlpha = 0.5;
      hca.XLabel.String = 'z';      
    end
    if 0 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    
    hca.XLabel.String = 'z (km)';  
    hca.YLabel.String = '\Lambda (eV)';  
  end
  hold(hca,'off')
  
  hca.Title.String = cs(ics).str_title;
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z (km)')
    ylabel(hca,'y (km)')         
  end
  hold(hca,'off')
end


for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  h(ii).Box = 'on';
  h(ii).FontSize = 12;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
c_eval('h(?).XLim = [-13 13];',1:numel(h))
%c_eval('h(?).YLim = [00 600];',2:3:numel(h))
c_eval('h(?).YLim(2) = 180;',2:3:numel(h))

hl = findobj(gcf,'type','line');

c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).Visible = ''off'';',1:numel(h))
compact_panels(0.04,0.06)


end



%% including 
%% integrate electron orbits for a given magnetic field: several particles, common plotting with setup before

units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

Er = 0e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 0*5*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 1000*units.e;
Lphi = L;



% General current sheet structure:
% Harris current sheet, incl. parameters in functions to make more
% versatile later on
Bx = @(x,y,z,B0,d) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z,Bg) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z,Bn) x*0 + y*0 + z*0 + Bn;
%Ex = @(x,y,z) x*0 + y*0 + z*0;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z,Er) x*0 + y*0 + z*0 + Er;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ez = @(x,y,z,phi0,Lphi) x*0 + y*0 + z*0 -phi0/Lphi;
Ax = @(z,Bg) Bg*z;
Ay = @(z,B0,d) B0*d*log(cosh(z/d));
Az = @(z) 0*z;
Phi = @(z,phi0,Lphi) phi0*exp(-(z/Lphi).^2);
Phi = @(z,phi0,Lphi) phi0*z/Lphi;

% Pseudo-potential
Gamma = @(z,py,px,B0,Bg,d,phi0,Lphi) 0.5/m*((py-q*Ay(z,B0,d)).^2+(px-q*Ax(z,Bg)).^2) + q*Phi(z,phi0,Lphi);
  
% What colors to use for lines
colors = [0 0 1; 0.5 0.5 1; 0 0 0.5; 0 0 0.2];
colors_matlab = pic_colors('matlab');
colors = [0 0 0; colors_matlab(1:3,:)];
colors_BE = [colors_matlab([1 5],:); 0 0 0; colors(3,:)];
fontsize = 12;

% Define our different current sheets and particles
ics = 0;
cs = struct([]);

if 0 % No guide field, different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 300;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 180;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 270;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % No guide field, slightly different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 50;
  T = 0.02;
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -120;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -70;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -20;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % Finite guide field
  ics = ics + 1;  
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 3.0e-9; % guide field, T
  cs(ics).Bn = 0.0*1e-9; % normal field, T
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  cs(ics).phi0 = 0; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'Weak guide field, same starting position and phase but different energies';

    % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = -90;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end  
  cs(ics).particles = particles;
end
if 1 % No guide field, finite constant En
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m    
  cs(ics).phi0 = 20000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, finite normal electric field';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = 90;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
% Test the same particles for all configurations
if 1
  T = 0.015;
  energy = 30;
  z0 = -5;
  iParticle = 0;
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -60;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -65;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  for ics = 1:3
    cs(ics).particles = particles;
  end
end

%cs(ics).str_title = 'No guide field, same starting position and energy but different phase';
% Set up plot
fig = figure(35);
nCols = 3;
nRows = 4;
h = setup_subplots(nRows,nCols,'vertical');
isub = 1;

% Loop through different settings
for ics = 1:numel(cs)  
  % Extract parameters for ease of writing below
  B0 = cs(ics).B0; % asymptotical magnetic field, T
  Bg = cs(ics).Bg; % guide field, T
  Bn = cs(ics).Bn; % normal field, T
  %Er = cs(ics).Er = 0e-3; % reconnection electric field, V/m
  d = cs(ics).d; % thickness of current sheet, m
  phi0 = cs(ics).phi0; % potential
  Lphi = cs(ics).Lphi; % length scale of potential
  
  % Integrate particle trajectories
  particles = cell(numel(cs(ics).particles),1);
  for iParticle = 1:numel(cs(ics).particles)
    % Load particle data
    particle = cs(ics).particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0,Bg);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0,B0,d);
    Gamma0 = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,@(x,y,z)Bx(x,y,z,B0,d),@(x,y,z)By(x,y,z,Bg),@(x,y,z)Bz(x,y,z,Bn),Ex,@(x,y,z)Ey(x,y,z,Er),@(x,y,z)Ez(x,y,z,phi0,Lphi));
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.t = t;
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;      
  end
  
   % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,200);
  if 1 % Bx, By, Ez
    set(hca,'ColorOrder',colors_BE)
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,':',zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    set(hca,'ColorOrder',colors_BE)
    %axes(hca)
    %yyaxis right;
    %ax2 = gca;
    %ax2.Position = hca.Position;
    %ax2.
%     ax2 = axes('XAxisLocation','top',...
%          'YAxisLocation','right',...
%          'Color','none',...
%          'XColor','k','YColor','k',...
%          'Position',hca.Position);
    %plot(ax2,zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    %ax2.YLabel.String = 'E_z (mV/m)';
    
    % irf_legend(hca,{'B_x','B_y','|B|','E_z'},[0.02 0.08],'fontsize',fontsize)
    hca.YLabel.String = 'B (nT), E (mV/m)';
  else
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k');
    irf_legend(hca,{'B_L','B_M','B_N'},[0.98 0.98])
  hca.YLabel.String = 'B (nT)';
  end
  %hca.YTick = 0;
  hca.XLabel.String = 'N (km)';
  %hca.Title.String = cs(ics).str_title;
  
  
  % Plot particles
  hca = h(isub); isub = isub + 1;  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'N (km)')
    ylabel(hca,'M (km)')         
  end
  hold(hca,'off')


  % Plot generalized spotential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 0 % along orbit
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')
      if iParticle == 1, hold(hca,'on'); end
      %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
      %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')    
    end
    if 1 % line in orbitrange
      if iParticle == 1, hold(hca,'on'); end
        
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      
      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px,B0,Bg,d,phi0,Lphi);

     % hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:),'facealpha',0.2);
     % hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      %hp.FaceAlpha = 0.5;
      hca.XLabel.String = 'z';      
    end
    if 0 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    
    hca.XLabel.String = 'N (km)';  
    hca.YLabel.String = '\Lambda (eV)';     
  end
  hold(hca,'off')
  
  

  % Plot phase space separatrix  
  if 1
    hca = h(isub); isub = isub + 1;  
    vx = 0*1e3; 
    vy = 0.99*5000*linspace(-1,1,102)*1e3; 
    vz = 0.99*5000*linspace(-1,1,103)*1e3; 
    %vy = particles{1}.vy(1);
    %vz = particles{1}.vz(1);

    [VX,VY,VZ] = ndgrid(vx,vy,vz);
    z0 = -5*1e3; 
    Ay0 = Ay(z0,B0,d);
    Ax0 = Ax(z0,Bg);
    Px0 = m*VX + q*Ax0;
    Py0 = m*VY + q*Ay0;
    Gamma_0 = Gamma(0,Py0,Px0,B0,Bg,d,phi0,Lphi);

    

    Uk = 0.5*m*(VX.^2 + VY.^2 + VZ.^2);
    Up = -Phi(z0,phi0,Lphi)*units.e;

    H = Uk + Up;

    Pxy2 = (Px0.^2 + Py0.^2)/2/m;
    Gamma_0_Uk = Gamma_0 - 1*Uk - Up;
    Gamma_0_Uk = H - Pxy2;
    [C,hc] = contour(hca,squeeze(VY)*1e-3,squeeze(VZ)*1e-3,squeeze(Gamma_0_Uk)/units.eV,[0 0],'k--','linewidth',1);
    %[C,hc] = contourf(hca,squeeze(VY)*1e-3,squeeze(VZ)*1e-3,squeeze(Uk)/units.eV);
    %clabel(C,hc)
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';

    hold(hca,'on')
    for ip = 1:numel(particles)
      vxp = particles{ip}.vx(1);
      vyp = particles{ip}.vy(1);
      vzp = particles{ip}.vz(1);
      z0_ = particles{ip}.z(1);
      hm = plot(hca,vyp*1e-3,vzp*1e-3,'o','color',colors(ip,:),'MarkerFaceColor',colors(ip,:),'MarkerSize',8);

      px = m*vxp + q*Ax(z0,Bg);
      pz = m*vzp + q*Az(z0);
      py = m*vyp + q*Ay(z0,B0,d);
      Gamma0_ = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    end
    hold(hca,'off')
    axis(hca,'equal')
    % 1;
  else
    isub = isub + 1;  
  end
end


for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  h(ii).Box = 'on';
  h(ii).FontSize = 12;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
%c_eval('h(?).XLim = [-13 13];',1:numel(h))
%c_eval('h(?).YLim = [00 600];',2:3:numel(h))
c_eval('h(?).YLim(2) = 180;',3:nRows:numel(h))

hl = findobj(gcf,'type','line');

c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).Visible = ''off'';',1:numel(h))
hx = [1:3 nRows+(1:3) 2*nRows+(1:3)];
for ip = hx
  h(ip).XLim = [-12 12];
end
irf_legend(h(nRows*(nCols-1)+1),{'B_L','B_M','|B|','E_N'}',[1.05 0.98],'fontsize',fontsize)
%compact_panels(h(hx),0.04,0.06)
%h(1).Title.String = 'Only';
%compact_panels(h([1:3 nRows+(1:3) 2*nRows+(1:3)]),0.01,0.01)
%compact_panels(h(nRows*[1:3]),0.01,0.01)

c_eval('h(?).YLim = [-11 11];',[1 nRows+1 2*nRows+1])

%c_eval('h(?).YLabel = [];',[nRows+1:nRows*nCols])
%c_eval('h(?).YTickLabel = [];',[nRows+1:nRows*nCols])

c_eval('h(?).YLim = [-30 100];',3+nRows*[0:2])

c_eval('h(?).Box = "off"; h(?).XAxisLocation = "origin"; h(?).YAxisLocation = "origin";',1:12)
c_eval('h(?).XLabel.VerticalAlignment = "middle"; h(?).XLabel.VerticalAlignment = "bottom"; h(?).YLabel.HorizontalAlignment = "center"; h(?).YLabel.VerticalAlignment = "bottom";',1:12)
%%
legs = {'a1','a2','a3','a4','b1','b2','b3','b4','c1','c2','c3','c3'};
legs = {'a1','b1','c1','d1','a2','b2','c3','d2','a3','b3','c3','d3'};
for ip = 1:numel(h)
  hleg(ip) = irf_legend(h(ip),[legs{ip} ,')'],[-0.01 1.01],'color',[0 0 0],'fontsize',12);
end

%% Including phase space separatrix
units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

Er = 0e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 1*5*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 11000*units.e;
Lphi = L;

% General current sheet structure:
% Harris current sheet, incl. parameters in functions to make more
% versatile later on
Bx = @(x,y,z,B0,d) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z,Bg) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z,Bn) x*0 + y*0 + z*0 + Bn;
%Ex = @(x,y,z) x*0 + y*0 + z*0;
%Ey = @(x,y,z,Er) x*0 + y*0 + z*0 + Er;
Ez = @(x,y,z,phi0,Lphi) x*0 + y*0 + z*0 -phi0/Lphi;
Ax = @(z,Bg) Bg*z;
Ay = @(z,B0,d) B0*d*log(cosh(z/d));
Az = @(z) 0*z;
Phi = @(z,phi0,Lphi) phi0*exp(-(z/Lphi).^2);
Phi = @(z,phi0,Lphi) phi0*z/Lphi;

% Pseudo-potential
Gamma = @(z,py,px,B0,Bg,d,phi0,Lphi) 0.5/m*((py-q*Ay(z,B0,d)).^2+(px-q*Ax(z,Bg)).^2) + q*Phi(z,phi0,Lphi);
  
% What colors to use for lines
colors = [0 0 1; 0.5 0.5 1; 0 0 0.5; 0 0 0.2];
colors = pic_colors('matlab');
fontsize = 12;

% Define our different current sheets and particles
ics = 0;
cs = struct([]);

if 0 % No guide field, different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 300;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 180;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 270;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % No guide field, slightly different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 50;
  T = 0.02;
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -120;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -70;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -20;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % Finite guide field
  ics = ics + 1;  
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 1.6e-9; % guide field, T
  cs(ics).Bn = 0.0*1e-9; % normal field, T
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  cs(ics).phi0 = 0; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'Weak guide field, same starting position and phase but different energies';

    % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = -90;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end  
  cs(ics).particles = particles;
end
if 1 % No guide field, finite constant En
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m    
  cs(ics).phi0 = 22000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, finite normal electric field';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = 90;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end

% Set up plot
fig = figure(36);
nCols = 3;
nRows = 4;
h = setup_subplots(nCols,nRows,'vertical');
isub = 1;

% Loop through different settings
for ics = 1:numel(cs)  
  % Extract parameters for ease of writing below
  B0 = cs(ics).B0; % asymptotical magnetic field, T
  Bg = cs(ics).Bg; % guide field, T
  Bn = cs(ics).Bn; % normal field, T
  %Er = cs(ics).Er = 0e-3; % reconnection electric field, V/m
  d = cs(ics).d; % thickness of current sheet, m
  phi0 = cs(ics).phi0; % potential
  Lphi = cs(ics).Lphi; % length scale of potential
  
  % Integrate particle trajectories
  particles = cell(numel(cs(ics).particles),1);
  for iParticle = 1:numel(cs(ics).particles)
    % Load particle data
    particle = cs(ics).particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0,Bg);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0,B0,d);
    Gamma0 = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,@(x,y,z)Bx(x,y,z,B0,d),@(x,y,z)By(x,y,z,Bg),@(x,y,z)Bz(x,y,z,Bn), Ex,@(x,y,z) Ey, @(x,y,z,Er)Ez(x,y,z,phi0,Lphi));
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;      
  end
  
   % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,200);
  if 1 % Bx, By, Ez
    set(hca,'ColorOrder',[colors(1:2,:); 0 0 0; colors(3,:)])
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k',zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    set(hca,'ColorOrder',[colors(1:2,:); 0 0 0; colors(3,:)])
    %axes(hca)
    %yyaxis right;
    %ax2 = gca;
    %ax2.Position = hca.Position;
    %ax2.
%     ax2 = axes('XAxisLocation','top',...
%          'YAxisLocation','right',...
%          'Color','none',...
%          'XColor','k','YColor','k',...
%          'Position',hca.Position);
    %plot(ax2,zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    %ax2.YLabel.String = 'E_z (mV/m)';
    irf_legend(hca,{'B_x','B_y','|B|','E_z'},[0.02 0.08],'fontsize',fontsize)
    hca.YLabel.String = 'B (nT), E (mV/m)';
  else
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k');
    irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  hca.YLabel.String = 'B (nT)';
  end
  %hca.YTick = 0;
  hca.XLabel.String = 'z (km)';
  

  % Plot generalized spotential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 0 % along orbit
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')
      if iParticle == 1, hold(hca,'on'); end
      %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
      %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')    
    end
    if 1 % line in orbitrange
      if iParticle == 1, hold(hca,'on'); end
      
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
      %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px,B0,Bg,d,phi0,Lphi);

      %hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      %hp.FaceAlpha = 0.5;
      hca.XLabel.String = 'z';      
    end
    if 0 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    
    hca.XLabel.String = 'z (km)';  
    hca.YLabel.String = '\Lambda (eV)';  
  end
  hold(hca,'off')
  
  hca.Title.String = cs(ics).str_title;
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z (km)')
    ylabel(hca,'y (km)')         
  end
  hold(hca,'off')

  % Plot phase space separatrix
  hca = h(isub); isub = isub + 1;  
  vx = linspace(-3000,3000,100); 
  vy = linspace(-3000,3000,100); 
  [VX,VY] = ndgrid(vx,vy);
  z0 = 5; 
  Ay0 = Ay(z0);
  Ax0 = Ax(z0);
  Px0 = m*VX + q*Ax0;
  Py0 = m*VY + q*Ay0;
  Gamma_ = Gamma(z0,Px0,Py0,B0,Bg,d,phi0,Lphi);
  contourf(hca,VX,BY,Gamma_,[0 0]);
  1;
end


for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  h(ii).Box = 'on';
  h(ii).FontSize = 12;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
c_eval('h(?).XLim = [-13 13];',1:numel(h))
%c_eval('h(?).YLim = [00 600];',2:3:numel(h))
c_eval('h(?).YLim(2) = 180;',2:3:numel(h))

hl = findobj(gcf,'type','line');

c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).Visible = ''off'';',1:numel(h))
compact_panels(0.04,0.06)


%% integrate electron orbits for a given magnetic field: several particles, common plotting with setup before

units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

Er = 0.5e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 0*5*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 1000*units.e;
Lphi = L;



% General current sheet structure:
% Harris current sheet, incl. parameters in functions to make more
% versatile later on
Bx = @(x,y,z,B0,d) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z,Bg) x*0 + y*0 + z*0 + Bg;
Bz = @(x,y,z,Bn) x*0 + y*0 + z*0 + Bn;
%Ex = @(x,y,z) x*0 + y*0 + z*0;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z,Er) x*0 + y*0 + z*0 + Er;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ez = @(x,y,z,phi0,Lphi) x*0 + y*0 + z*0 -phi0/Lphi;
Ax = @(z,Bg) Bg*z;
Ay = @(z,B0,d) B0*d*log(cosh(z/d));
Az = @(z) 0*z;
Phi = @(z,phi0,Lphi) phi0*exp(-(z/Lphi).^2);
Phi = @(z,phi0,Lphi) phi0*z/Lphi;

% Pseudo-potential
Gamma = @(z,py,px,B0,Bg,d,phi0,Lphi) 0.5/m*((py-q*Ay(z,B0,d)).^2+(px-q*Ax(z,Bg)).^2) + q*Phi(z,phi0,Lphi);
  
% What colors to use for lines
colors = [0 0 1; 0.5 0.5 1; 0 0 0.5; 0 0 0.2];
colors_matlab = pic_colors('matlab');
colors = [0 0 0; colors_matlab(1:3,:)];
fontsize = 12;

% Define our different current sheets and particles
ics = 0;
cs = struct([]);

if 1 % No guide field
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = -0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m    
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 300;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 180;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 270;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % Finite guide field
  ics = ics + 1;  
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 3.0e-9; % guide field, T
  cs(ics).Bn = 0.0*1e-9; % normal field, T
  cs(ics).Er = -0.2e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  cs(ics).phi0 = 0; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'Weak guide field, same starting position and phase but different energies';
end
if 1 % No guide field, finite constant En
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m    
  cs(ics).phi0 = 20000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, finite normal electric field';
end

cs = cs(1:3);
% Initiate a bunch of particles
if 1
  T = 0.015;
  energy = 50;
  y0 = 0;
  z0 = -3;
  nParticles = 7;
  iParticle = 0;
  phase = linspace(0,360-360/nParticles,nParticles);  
  particles = {};

  for iP = 1:nParticles
    %particle = struct([]);
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = phase(iP);
    particle.m = units.me;
    particle.q = -units.e;
    particles{iP} = particle;
  end
  for ics = 1:numel(cs)
    cs(ics).particles = particles;
  end
end

%cs(ics).str_title = 'No guide field, same starting position and energy but different phase';
% Set up plot
fig = figure(35);
nCols = 3;
nRows = 5;
h = setup_subplots(nRows,nCols,'vertical');
isub = 1;

% Loop through different settings
for ics = 1:numel(cs)  
  % Extract parameters for ease of writing below
  B0 = cs(ics).B0; % asymptotical magnetic field, T
  Bg = cs(ics).Bg; % guide field, T
  Bn = cs(ics).Bn; % normal field, T
  Er = cs(ics).Er; % reconnection electric field, V/m
  d = cs(ics).d; % thickness of current sheet, m
  phi0 = cs(ics).phi0; % potential
  Lphi = cs(ics).Lphi; % length scale of potential
  
  % Integrate particle trajectories
  particles = cell(numel(cs(ics).particles),1);
  for iParticle = 1:numel(cs(ics).particles)
    % Load particle data
    particle = cs(ics).particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0,Bg);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0,B0,d);
    Gamma0 = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,@(x,y,z)Bx(x,y,z,B0,d),@(x,y,z)By(x,y,z,Bg),@(x,y,z)Bz(x,y,z,Bn),Ex,@(x,y,z)Ey(x,y,z,Er),@(x,y,z)Ez(x,y,z,phi0,Lphi));
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;      
  end
  
   % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,200);
  if 1 % Bx, By, Ez
    set(hca,'ColorOrder',[colors(1:2,:); 0 0 0; colors(3,:)])
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k',zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    set(hca,'ColorOrder',[colors(1:2,:); 0 0 0; colors(3,:)])
    %axes(hca)
    %yyaxis right;
    %ax2 = gca;
    %ax2.Position = hca.Position;
    %ax2.
%     ax2 = axes('XAxisLocation','top',...
%          'YAxisLocation','right',...
%          'Color','none',...
%          'XColor','k','YColor','k',...
%          'Position',hca.Position);
    %plot(ax2,zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    %ax2.YLabel.String = 'E_z (mV/m)';
    
    % irf_legend(hca,{'B_x','B_y','|B|','E_z'},[0.02 0.08],'fontsize',fontsize)
    hca.YLabel.String = 'B (nT), E (mV/m)';
  else
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k');
    irf_legend(hca,{'B_x','B_y','B_z'},[0.98 0.98])
  hca.YLabel.String = 'B (nT)';
  end
  %hca.YTick = 0;
  hca.XLabel.String = 'z (km)';
  %hca.Title.String = cs(ics).str_title;
  

  % Plot generalized spotential
  if 0
    hca = h(isub); isub = isub + 1;
    for iParticle = 1:numel(particles)
      px = particles{iParticle}.px;
      py = particles{iParticle}.py;
      if 0 % along orbit
        vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
        hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')
        if iParticle == 1, hold(hca,'on'); end
        %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
        %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')    
      end
      if 1 % line in orbitrange
        if iParticle == 1, hold(hca,'on'); end
          
        vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
        hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')
  
        
        zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
        vecGamma = Gamma(zz_patch,py,px,B0,Bg,d,phi0,Lphi);
  
       % hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
        hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:),'facealpha',0.2);
       % hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
        %hp.FaceAlpha = 0.5;
        hca.XLabel.String = 'z';      
      end
      if 0 % patch in orbitrange
  
        vecGamma = Gamma(zz,py,px);
        hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')
  
        vecGamma = Gamma(particles{iParticle}.z,py,px);
        hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')
  
        zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
        vecGamma = Gamma(zz_patch,py,px);
  
        hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
        hp.FaceAlpha = 0.5;
  
    %     limGamma = vecGamma; 
    %     zzGamma = zz;
    %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
    %     limGamma(limGamma>electron_energy*units.eV) = NaN;
    %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
    %     %hpot.Color = ;
        hca.XLabel.String = 'z';      
      end
      
      hca.XLabel.String = 'z (km)';  
      hca.YLabel.String = '\Lambda (eV)';     
    end
    hold(hca,'off')
  end
  
  
  % Plot particles
  hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3);
    colors(ip,1:3) = horb.Color;
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'z (km)')
    ylabel(hca,'y (km)')         
  end
  hold(hca,'off')


  % Plot phase space separatrix  
  if 1
    hca = h(isub); isub = isub + 1;  
    vx = 0*1e3; 
    vy = 0.99*5000*linspace(-1,1,102)*1e3; 
    vz = 0.99*5000*linspace(-1,1,103)*1e3; 
    %vy = particles{1}.vy(1);
    %vz = particles{1}.vz(1);

    [VX,VY,VZ] = ndgrid(vx,vy,vz);
    z0 = -5*1e3; 
    Ay0 = Ay(z0,B0,d);
    Ax0 = Ax(z0,Bg);
    Px0 = m*VX + q*Ax0;
    Py0 = m*VY + q*Ay0;
    Gamma_0 = Gamma(0,Py0,Px0,B0,Bg,d,phi0,Lphi);

    
    Uk = 0.5*m*(VX.^2 + VY.^2 + VZ.^2);
    Up = -Phi(z0,phi0,Lphi)*units.e;
    Gamma_0_Uk = Gamma_0 - 1*Uk - Up;
    [C,hc] = contour(hca,squeeze(VY)*1e-3,squeeze(VZ)*1e-3,squeeze(Gamma_0_Uk)/units.eV,[0 0],'k--','linewidth',1);
    %[C,hc] = contourf(hca,squeeze(VY)*1e-3,squeeze(VZ)*1e-3,squeeze(Uk)/units.eV);
    %clabel(C,hc)
    hca.XLabel.String = 'v_y (km/s)';
    hca.YLabel.String = 'v_z (km/s)';

    hold(hca,'on')
    for ip = 1:numel(particles)
      vxp = particles{ip}.vx(1);
      vyp = particles{ip}.vy(1);
      vzp = particles{ip}.vz(1);
      z0_ = particles{ip}.z(1);
      hm = plot(hca,vyp*1e-3,vzp*1e-3,'o','color',colors(ip,:),'MarkerFaceColor',colors(ip,:),'MarkerSize',8);

      px = m*vxp + q*Ax(z0,Bg);
      pz = m*vzp + q*Az(z0);
      py = m*vyp + q*Ay(z0,B0,d);
      Gamma0_ = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    end
    hold(hca,'off')
    axis(hca,'equal')
    % 1;
  else
    isub = isub + 1;  
  end

  % Plot phase space
   hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    vx = particles{ip}.vx;
    vy = particles{ip}.vy;
    vz = particles{ip}.vz;
    horb = plot(hca,vy*1e-3+0*ip,vz*1e-3+0*ip,'color',colors(ip,:));    
    
    if ip == 1, hold(hca,'on'); end
    %plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
    %         z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'v_y (km/s)')
    ylabel(hca,'v_z (km/s)')
  end
  hold(hca,'off')
  axis(hca,'equal')

  % Plot phase space
   hca = h(isub); isub = isub + 1;
  
  for ip = 1:numel(particles)
    vx = particles{ip}.vx;
    vy = particles{ip}.vy;
    vz = particles{ip}.vz;
    horb = plot(hca,vx*1e-3+0*ip,vz*1e-3+0*ip,'color',colors(ip,:));    
    
    if ip == 1, hold(hca,'on'); end
    %plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
    %         z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'v_x (km/s)')
    ylabel(hca,'v_z (km/s)')
  end
  hold(hca,'off')
  axis(hca,'equal')

end


for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  h(ii).Box = 'on';
  h(ii).FontSize = 12;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
%c_eval('h(?).XLim = [-13 13];',1:numel(h))
%c_eval('h(?).YLim = [00 600];',2:3:numel(h))
%c_eval('h(?).YLim(2) = 180;',2:nRows:numel(h))

hl = findobj(gcf,'type','line');

c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).Visible = ''off'';',1:numel(h))
%hx = [1:3 nRows+(1:3) 2*nRows+(1:3)];
irf_legend(h(nRows*(nCols-1)+1),{'B_x','B_y','|B|','E_z'}',[1.05 0.98],'fontsize',fontsize)
%compact_panels(h(hx),0.04,0.06)
%h(1).Title.String = 'Only';
%compact_panels(h([1:3 nRows+(1:3) 2*nRows+(1:3)]),0.01,0.01)
%compact_panels(h(nRows*[1:3]),0.01,0.01)

%c_eval('h(?).YLim = [-11 11];',[1 nRows+1 2*nRows+1])

%c_eval('h(?).YLabel = [];',[nRows+1:nRows*nCols])
%c_eval('h(?).YTickLabel = [];',[nRows+1:nRows*nCols])

%c_eval('h(?).YLim = [-30 100];',2+nRows*[0:2])

c_eval('h(?).Box = "on"; h(?).XAxisLocation = "origin"; h(?).YAxisLocation = "origin";',1:numel(h))
c_eval('h(?).XLabel.VerticalAlignment = "middle"; h(?).XLabel.HorizontalAlignment = "left"; h(?).YLabel.HorizontalAlignment = "center"; h(?).YLabel.VerticalAlignment = "bottom";',1:numel(h))




%% Hall fields: integrate electron orbits for a given magnetic field: several particles, common plotting with setup before

units = irf_units;
L = 5000e3;
q = -units.e;
m = units.me;

Er = 0e-3; % reconnection electric field, V/m
B0 = 10e-9; % asymptotical magnetic field, T
d = 5e3; % thickness of current sheet, m
Bg = 0*5*1e-9; % guide field, T
Bn = 0.0*1e-9; % normal field, T
phi0 = 1000*units.e;
Lphi = L;



% General current sheet structure:
% Harris current sheet, incl. parameters in functions to make more
% versatile later on
Bx = @(x,y,z,B0,d) x*0 + y*0 - B0*tanh(z/d);
By = @(x,y,z,Bg) x*0 + y*0 + z*0 + Bg*z/d;
Bz = @(x,y,z,Bn) x*0 + y*0 + z*0 + Bn;
%Ex = @(x,y,z) x*0 + y*0 + z*0;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ey = @(x,y,z,Er) x*0 + y*0 + z*0 + Er;
Ex = @(x,y,z) x*0 + y*0 + z*0;
Ez = @(x,y,z,phi0,Lphi) x*0 + y*0 + z*0 -phi0/Lphi;
Ax = @(z,Bg) Bg*z.^2/2;
Ay = @(z,B0,d) B0*d*log(cosh(z/d));
Az = @(z) 0*z;
Phi = @(z,phi0,Lphi) phi0*exp(-(z/Lphi).^2);
Phi = @(z,phi0,Lphi) phi0*z/Lphi;

% Pseudo-potential
Gamma = @(z,py,px,B0,Bg,d,phi0,Lphi) 0.5/m*((py-q*Ay(z,B0,d)).^2+(px-q*Ax(z,Bg)).^2) + q*Phi(z,phi0,Lphi);
  
% What colors to use for lines
colors = [0 0 1; 0.5 0.5 1; 0 0 0.5; 0 0 0.2];
colors_matlab = pic_colors('matlab');
colors = [0 0 0; colors_matlab(1:3,:)];
colors_BE = [colors_matlab([1 5],:); 0 0 0; colors(3,:)];
fontsize = 12;

% Define our different current sheets and particles
ics = 0;
cs = struct([]);

if 0 % No guide field, different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 300;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 180;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 270;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % No guide field, slightly different phase
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  
  cs(ics).phi0 = 0*11000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, same starting position and energy but different phase';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  energy = 50;
  T = 0.02;
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -120;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -70;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -20;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
if 1 % Finite guide field
  ics = ics + 1;  
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 3.0e-9; % guide field, T
  cs(ics).Bn = 0.0*1e-9; % normal field, T
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m  
  cs(ics).phi0 = 10000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'Weak guide field, same starting position and phase but different energies';

    % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = -90;
  if 1 % no drift
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % meandering 8's
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % grad B
    iParticle = iParticle + 1;
    particle.T = 0.02;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end  
  cs(ics).particles = particles;
end
if 1 % No guide field, finite constant En
  ics = ics + 1;
  % Define current sheet
  cs(ics).B0 = 10e-9; % asymptotical magnetic field, T
  cs(ics).Bg = 0e-9; % guide field
  cs(ics).Bn = 0.0*1e-9; % normal field, TT
  cs(ics).Er = 0e-3; % reconnection electric field, V/m
  cs(ics).d = 5e3; % thickness of current sheet, m    
  cs(ics).phi0 = 20000; % potential
  cs(ics).Lphi = L; % length scale of potential
  cs(ics).str_title = 'No guide field, finite normal electric field';

  % Initiate particles
  particles = {};
  iParticle = 0;
  z0 = -5;
  y0 = 0;
  %energy = 300;
  velocity_angle = 90;
  T = 0.02;
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 20;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 50;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 100;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = 300;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = velocity_angle;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  cs(ics).particles = particles;
end
% Test the same particles for all configurations
if 1
  T = 0.015;
  energy = 30;
  z0 = -5;
  iParticle = 0;
  if 1 % 3
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 2
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -60;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -65;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 1 % 1
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = -00;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  if 0 % 4
    iParticle = iParticle + 1;
    particle.T = T;
    particle.energy = energy;
    particle.z0 = z0; 
    particle.y0 = y0;
    particle.velocity_angle = 90;
    particle.m = units.me;
    particle.q = -units.e;
    particles{iParticle} = particle;
  end
  for ics = 1:3
    cs(ics).particles = particles;
  end
end

%cs(ics).str_title = 'No guide field, same starting position and energy but different phase';
% Set up plot
fig = figure(35);
nCols = 3;
nRows = 4;
h = setup_subplots(nRows,nCols,'vertical');
isub = 1;

% Loop through different settings
for ics = 1:numel(cs)  
  % Extract parameters for ease of writing below
  B0 = cs(ics).B0; % asymptotical magnetic field, T
  Bg = cs(ics).Bg; % guide field, T
  Bn = cs(ics).Bn; % normal field, T
  %Er = cs(ics).Er = 0e-3; % reconnection electric field, V/m
  d = cs(ics).d; % thickness of current sheet, m
  phi0 = cs(ics).phi0; % potential
  Lphi = cs(ics).Lphi; % length scale of potential
  
  % Integrate particle trajectories
  particles = cell(numel(cs(ics).particles),1);
  for iParticle = 1:numel(cs(ics).particles)
    % Load particle data
    particle = cs(ics).particles{iParticle};
    electron_energy = particle.energy;
    velocity_angle = particle.velocity_angle;
    q = particle.q;
    m = particle.m;
    T = particle.T;

    vt = sqrt(electron_energy*units.eV*2/m)/1000;
    % Initial positions and velocitites
    x0 = 0*1e3; % m
    y0 = particle.y0*1e3; 
    z0 = particle.z0*1e3;
    vx0 = 0*1e3; % m/s
    vy0 = vt*cosd(velocity_angle)*1e3;
    vz0 = -vt*sind(velocity_angle)*1e3;

    % Generalized potential
    px = m*vx0 + q*Ax(z0,Bg);
    pz = m*vz0 + q*Az(z0);
    py = m*vy0 + q*Ay(z0,B0,d);
    Gamma0 = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    % Integrate particle orbit
    x_init = [x0;y0;z0;vx0;vy0;vz0]; % m, m/s
    EoM = @(ttt,xxx) thesis.eom_harris(ttt,xxx,@(x,y,z)Bx(x,y,z,B0,d),@(x,y,z)By(x,y,z,Bg),@(x,y,z)Bz(x,y,z,Bn),Ex,@(x,y,z)Ey(x,y,z,Er),@(x,y,z)Ez(x,y,z,phi0,Lphi));
    [t,x_sol] = ode45(EoM,[0 T],x_init); % ,options
    x =  x_sol(:,1);  y = x_sol(:,2);  z = x_sol(:,3);
    vx = x_sol(:,4); vy = x_sol(:,5); vz = x_sol(:,6); 
    particles{iParticle}.x = x;
    particles{iParticle}.y = y;
    particles{iParticle}.z = z;
    particles{iParticle}.vx = vx;
    particles{iParticle}.vy = vy;
    particles{iParticle}.vz = vz;
    particles{iParticle}.px = px;
    particles{iParticle}.py = py;
    particles{iParticle}.pz = pz;
    particles{iParticle}.Gamma0 = Gamma0;      
  end
  
   % Plot magnetic field
  hca = h(isub); isub = isub + 1;
  zz = linspace(-5*d,5*d,200);
  if 1 % Bx, By, Ez
    set(hca,'ColorOrder',colors_BE)
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,':',zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    set(hca,'ColorOrder',colors_BE)
    %axes(hca)
    %yyaxis right;
    %ax2 = gca;
    %ax2.Position = hca.Position;
    %ax2.
%     ax2 = axes('XAxisLocation','top',...
%          'YAxisLocation','right',...
%          'Color','none',...
%          'XColor','k','YColor','k',...
%          'Position',hca.Position);
    %plot(ax2,zz*1e-3,Ez(0,0,zz,phi0,Lphi)*1e3);
    %ax2.YLabel.String = 'E_z (mV/m)';
    
    % irf_legend(hca,{'B_x','B_y','|B|','E_z'},[0.02 0.08],'fontsize',fontsize)
    hca.YLabel.String = 'B (nT), E (mV/m)';
  else
    plot(hca,zz*1e-3,[Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)]*1e9,zz*1e-3,sqrt(sum([Bx(0,0,zz,B0,d); By(0,0,zz,Bg); Bz(0,0,zz,Bn)].^2,1))*1e9,'k');
    irf_legend(hca,{'B_L','B_M','B_N'},[0.98 0.98])
  hca.YLabel.String = 'B (nT)';
  end
  %hca.YTick = 0;
  hca.XLabel.String = 'N (km)';
  %hca.Title.String = cs(ics).str_title;
  
  
  % Plot particles
  hca = h(isub); isub = isub + 1;  
  for ip = 1:numel(particles)
    z = particles{ip}.z;
    y = particles{ip}.y;
    horb = plot(hca,z*1e-3,y*1e-3,'color',colors(ip,:));
    
    if ip == 1, hold(hca,'on'); end
    plot(hca,z(1)*1e-3,y(1)*1e-3,'go',...
             z(end)*1e-3,y(end)*1e-3,'rx'); % plot in km's
    xlabel(hca,'N (km)')
    ylabel(hca,'M (km)')         
  end
  hold(hca,'off')


  % Plot generalized spotential
  hca = h(isub); isub = isub + 1;
  for iParticle = 1:numel(particles)
    px = particles{iParticle}.px;
    py = particles{iParticle}.py;
    if 0 % along orbit
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^3,'Linestyle','--');%,zz*1e-3,[py-q*Ay(zz)]')
      if iParticle == 1, hold(hca,'on'); end
      %vecGamma = Gamma(particles{iParticle}.z,py,px,B0,Bg,d,phi0,Lphi);
      %hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')    
    end
    if 1 % line in orbitrange
      if iParticle == 1, hold(hca,'on'); end
        
      vecGamma = Gamma(zz,py,px,B0,Bg,d,phi0,Lphi);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      
      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px,B0,Bg,d,phi0,Lphi);

     % hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:),'facealpha',0.2);
     % hp = plot(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,'color',colors(iParticle,:));
      %hp.FaceAlpha = 0.5;
      hca.XLabel.String = 'z';      
    end
    if 0 % patch in orbitrange

      vecGamma = Gamma(zz,py,px);
      hpot = plot(hca,zz*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:).^1,'Linestyle','-');%,zz*1e-3,[py-q*Ay(zz)]')

      vecGamma = Gamma(particles{iParticle}.z,py,px);
      hpot = plot(hca,particles{iParticle}.z*1e-3,vecGamma'/units.eV,'color',colors(iParticle,:));%,zz*1e-3,[py-q*Ay(zz)]')

      zz_patch = linspace(min(particles{iParticle}.z),max(particles{iParticle}.z),100);
      vecGamma = Gamma(zz_patch,py,px);

      hp = patch(hca,[zz_patch zz_patch(end:-1:1)]*1e-3,[vecGamma vecGamma*0+max(vecGamma)]/units.eV,colors(iParticle,:));
      hp.FaceAlpha = 0.5;

  %     limGamma = vecGamma; 
  %     zzGamma = zz;
  %     zzGamma(limGamma>electron_energy*units.eV) = NaN;
  %     limGamma(limGamma>electron_energy*units.eV) = NaN;
  %     hpot = plot(hca,zzGamma*1e-3,electron_energy+zzGamma*0,'color',horb(1).Color);%,zz*1e-3,[py-q*Ay(zz)]')
  %     %hpot.Color = ;
      hca.XLabel.String = 'z';      
    end
    
    hca.XLabel.String = 'N (km)';  
    hca.YLabel.String = '\Lambda (eV)';     
  end
  hold(hca,'off')
  
  

  % Plot phase space separatrix  
  if 1
    hca = h(isub); isub = isub + 1;  
    vx = 0*1e3; 
    vy = 0.99*5000*linspace(-1,1,102)*1e3; 
    vz = 0.99*5000*linspace(-1,1,103)*1e3; 
    %vy = particles{1}.vy(1);
    %vz = particles{1}.vz(1);

    [VX,VY,VZ] = ndgrid(vx,vy,vz);
    z0 = -5*1e3; 
    Ay0 = Ay(z0,B0,d);
    Ax0 = Ax(z0,Bg);
    Px0 = m*VX + q*Ax0;
    Py0 = m*VY + q*Ay0;
    Gamma_0 = Gamma(0,Py0,Px0,B0,Bg,d,phi0,Lphi);

    

    Uk = 0.5*m*(VX.^2 + VY.^2 + VZ.^2);
    Up = -Phi(z0,phi0,Lphi)*units.e;

    H = Uk + Up;

    Pxy2 = (Px0.^2 + Py0.^2)/2/m;
    Gamma_0_Uk = Gamma_0 - 1*Uk - Up;
    Gamma_0_Uk = H - Pxy2;
    [C,hc] = contour(hca,squeeze(VY)*1e-3,squeeze(VZ)*1e-3,squeeze(Gamma_0_Uk)/units.eV,[0 0],'k--','linewidth',1);
    %[C,hc] = contourf(hca,squeeze(VY)*1e-3,squeeze(VZ)*1e-3,squeeze(Uk)/units.eV);
    %clabel(C,hc)
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';

    hold(hca,'on')
    for ip = 1:numel(particles)
      vxp = particles{ip}.vx(1);
      vyp = particles{ip}.vy(1);
      vzp = particles{ip}.vz(1);
      z0_ = particles{ip}.z(1);
      hm = plot(hca,vyp*1e-3,vzp*1e-3,'o','color',colors(ip,:),'MarkerFaceColor',colors(ip,:),'MarkerSize',8);

      px = m*vxp + q*Ax(z0,Bg);
      pz = m*vzp + q*Az(z0);
      py = m*vyp + q*Ay(z0,B0,d);
      Gamma0_ = Gamma(0,py,px,B0,Bg,d,phi0,Lphi);

    end
    hold(hca,'off')
    axis(hca,'equal')
    % 1;
  else
    isub = isub + 1;  
  end
end


for ii = 1:numel(h) 
  %hca = subplot(nRows,nCols,ii);
  h(ii).Box = 'on';
  h(ii).FontSize = 12;
  axis(h(ii),'tight');
  %axis(h(ii),'off');
  %h(ii).Title.Position(2) = -0.2;
  %h(ii).Position(2) = 0.2;
  %h(ii).Position(4) = 0.7;
end
%c_eval('h(?).XLim = [-13 13];',1:numel(h))
%c_eval('h(?).YLim = [00 600];',2:3:numel(h))
c_eval('h(?).YLim(2) = 180;',3:nRows:numel(h))

hl = findobj(gcf,'type','line');

c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).Visible = ''off'';',1:numel(h))
hx = [1:3 nRows+(1:3) 2*nRows+(1:3)];
for ip = hx
  h(ip).XLim = [-12 12];
end
irf_legend(h(nRows*(nCols-1)+1),{'B_L','B_M','|B|','E_N'}',[1.05 0.98],'fontsize',fontsize)
%compact_panels(h(hx),0.04,0.06)
%h(1).Title.String = 'Only';
%compact_panels(h([1:3 nRows+(1:3) 2*nRows+(1:3)]),0.01,0.01)
%compact_panels(h(nRows*[1:3]),0.01,0.01)

c_eval('h(?).YLim = [-11 11];',[1 nRows+1 2*nRows+1])

%c_eval('h(?).YLabel = [];',[nRows+1:nRows*nCols])
%c_eval('h(?).YTickLabel = [];',[nRows+1:nRows*nCols])

c_eval('h(?).YLim = [-30 100];',3+nRows*[0:2])

c_eval('h(?).Box = "off"; h(?).XAxisLocation = "origin"; h(?).YAxisLocation = "origin";',1:12)
c_eval('h(?).XLabel.VerticalAlignment = "middle"; h(?).XLabel.VerticalAlignment = "bottom"; h(?).YLabel.HorizontalAlignment = "center"; h(?).YLabel.VerticalAlignment = "bottom";',1:12)
