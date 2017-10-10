function check_orbits(modelnr,gyrocenter0,T)
nParticles = numel(gyrocenter0);
ExB.model;
re = sqrt(10*Tper)/B0;
L = 14*lz;
W = 8*lr + 5*re; % width of box, km
zlim = L/2; 
rlim = W/2; 
eh_x0 = 0; % m
eh_y0 = 0; % m

% Initialize starting velocity and position
vtpar = cn_eV2v(Tpar,'ev'); % eV -> km/s (parallel thermal velocity)
vtper = cn_eV2v(Tper,'ev'); % eV -> km/s (perpendicular thermal velocity)
[x0,y0,vx0,vy0] = new_position(gyrocenter0,vtper);
z0 = x0;
vz0 = vx0*0;
% Integration setup
stopfunction = @(t,y) ExB.events(t,y,L);
options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
EoM = @(ttt,xxx) ExB.EquationOfMotion(ttt,xxx,B0,lr,lz,phi0,L,rlim,zlim,eh_x0,eh_y0);                                  

% Loop over particles
for jj = 1:nParticles
    x_init = [x0(jj);y0(jj);z0(jj);vx0(jj);vy0(jj);vz0(jj)]*1e3; % m, m/s
    
    % Integration
    [t,x_sol] = ode45(EoM,[0 T],x_init,options);
    x_end{jj} = x_sol; 
    x{jj} = x_sol(:,1);
    y{jj} = x_sol(:,2);
    z{jj} = x_sol(:,3);
    vx{jj} = x_sol(:,4);
    vy{jj} = x_sol(:,5);
    vz{jj} = x_sol(:,6);     
end
    
plot_trajectories(x,y,rlim)
%plot_eh(lr,phi0)

% Nested functions

end
% Not nested functions
function [x0,y0,vx0,vy0] = new_position(gyrocenter,vper)
r0 = gyrocenter;
x0 = r0.*cosd(0);
y0 = r0.*sind(0);  
vx0 = vper.*sind(0)*ones(numel(gyrocenter));
vy0 = vper.*cosd(0)*ones(numel(gyrocenter));
end
function plot_trajectories(x,y,rlim)
strPlot = 'plot(x{1}*1e-3,y{1}*1e-3';
for ii = 2:numel(x)
    strPlot = [strPlot ',x{' num2str(ii) '}*1e-3,y{' num2str(ii) '}*1e-3'];
end
strPlot = [strPlot ')'];    
eval([strPlot ';'])
axis equal
set(gca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])
end

    
