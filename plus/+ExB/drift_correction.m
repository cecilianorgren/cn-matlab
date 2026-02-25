function varargout = drift_correction(modelnr,gyrocenter,gyroradius,nGyroOrbits,doPlot)

% gyrocenter position and gyroradius is given in 'lr'
% start with just one
me = 9.10939999999e-31;
e = 1.6022e-19;
% = 100;
nParticles = numel(gyrocenter);
ExB.model;
gyrocenter = gyrocenter*lr;
gyroradius = gyroradius*lr;
fce = e*B0*1e-9/me/2/pi; % Hz, take B = 10 nT as model and get the vperp from that
vperp = 2*pi*gyroradius*fce; % km/s , gyroradius given in km
T = nGyroOrbits/fce;

% Initialize starting velocity and position
[x0,y0,z0,vx0,vy0,vz0] = new_position(gyrocenter,gyroradius,vperp);

% Integration setup
EoM = @(ttt,xxx) eom(ttt,xxx,B0,lr,lz,phi0);                                  

% Loop over particles
for jj = 1:nParticles
    x_init = [x0(jj);y0(jj);z0(jj);vx0(jj);vy0(jj);vz0(jj)]*1e3; % m, m/s
    
    % Integration
    [t,x_sol] = ode45(EoM,[0 T],x_init);
    x_end{jj} = x_sol; 
    x{jj} = x_sol(:,1);
    y{jj} = x_sol(:,2);
    z{jj} = x_sol(:,3);
    vx{jj} = x_sol(:,4);
    vy{jj} = x_sol(:,5);
    vz{jj} = x_sol(:,6);  
    % Get gyrocenter motion    
    [xc,yc,vc,vazc] = find_max(x{jj},y{jj},t);
    %disp([num2str(gyroradius) ' ' num2str(gyrocenter)])
    Ex = gyrocenter(jj)./(lr.^2).*phi0.*exp(-0.5*(gyrocenter(jj)/lr).^2-0.5*(0/lr).^2-0.5*(0/lz).^2); % mV/m    
    v_ExB(jj) = -1e-3*Ex*1e-3/(B0*1e-9); % km/s
    
    %disp(['v_ExB = ' num2str(mean(v_ExB)) ', vaz = ' num2str(mean(vazc)) ', vc = ' num2str(mean(vc))])
end
    
if doPlot
    plot_trajectories(x,y,4*lr); hold on;
    plot(xc*1e-3,yc*1e-3,'r*'); hold off
    titleStr{1} = ['v_c = ' num2str(mean(vc)*1e-3,'%.0f') ' km/s,  v_{ExB} = ' num2str(v_ExB,'%.0f') ' km/s,  v_{az} = ' num2str(mean(vazc)*1e-3,'%.0f') ' km/s,'];
    titleStr{2} = ['r_c/l_r = ' num2str(gyrocenter/lr) ', r_L/l_r = ' num2str(gyroradius/lr) ', lr = ' num2str(lr) ' km'];
    title(titleStr)
end
%plot_eh(lr,phi0)

if nargout == 0; return;
else varargout = {mean(vc)*1e-3,v_ExB,mean(vazc)*1e-3};
end

% Nested functions
end
% Not nested functions
function [x0,y0,z0,vx0,vy0,vz0] = new_position(gyrocenter,gyroradius,vperp)
    me = 9.10939999999e-31;
    e = 1.6022e-19;
    %vperp = gyroradius*B0*1e-9*e/me;    
    x0 = gyrocenter.*cosd(0)+gyroradius; % x = r_c + r_L
    y0 = gyrocenter.*sind(0);  % y = 0
    z0 = x0*0;
    vx0 = vperp.*sind(0)*ones(numel(gyrocenter)); % vx = 0
    vy0 = vperp.*cosd(0)*ones(numel(gyrocenter)); % vy = vperp
    vz0 = vx0*0;    
end
function [xc,yc,vc,vaztmp] = find_max(x,y,t)
    r = sqrt(x.^2+y.^2); 
    rmax = findpeaks(r);
    nOrbits = numel(rmax)-1;
    indOrbits = cell(nOrbits,1);
    tOrbits = zeros(nOrbits,1);
    for kk = 1:nOrbits
        indRmaxStart = find(r==rmax(kk));
        indRmaxStop = find(r==rmax(kk+1));
        indOrbits{kk} = indRmaxStart:indRmaxStop;
        tOrbits(kk) = mean(t(indRmaxStart:indRmaxStop));
        xc(kk) = mean(x(indOrbits{kk}));
        yc(kk) = mean(y(indOrbits{kk}));        
    end
    %%
    vx = diff(xc)./diff(tOrbits)';
    vy = diff(yc)./diff(tOrbits)';
    th = atan2d(xc(1:(end-1))+diff(xc)/2,yc(1:(end-1))+diff(yc)/2);
    vaztmp = +vy.*sind(th) - vx.*cosd(th); % azimuthal velocity     
    vc = sqrt(diff(xc).^2+diff(yc).^2)./diff(tOrbits)';
    %disp(['vaz = ' num2str(mean(vaztmp)) ' , vc = ' num2str(mean(vc))])
end
function plot_trajectories(x,y,lims)
strPlot = 'plot(x{1}*1e-3,y{1}*1e-3';
for ii = 2:numel(x)
    strPlot = [strPlot ',x{' num2str(ii) '}*1e-3,y{' num2str(ii) '}*1e-3'];
end
strPlot = [strPlot ')'];    
eval([strPlot ';'])
%lims = max(max([abs(get(gca,'xlim')) abs(get(gca,'ylim'))]));

axis equal
set(gca,'xlim',lims*[-1 1],'ylim',lims*[-1 1])
end
function  x_res = eom(t,x_vect,B0,lr,lz,phi0)
    e = 1.6022e-19;
    me = 9.10939999e-31;
    x0 = 0; y0 = 0;
    
    x = x_vect(1);
    y = x_vect(2);
    z = x_vect(3);
    vx = x_vect(4);
    vy = x_vect(5);
    vz = x_vect(6);

    x_res = zeros(6,1);

    x_res(1) = vx; % dx/dt = vx;
    x_res(2) = vy; % dy/dt = vy;
    x_res(3) = vz; % dz/dt = vz;
    x_res(4) = (-e/me)*((x/(lr*1e3)^2).*phi0.*exp(-0.5*((x-x0)/(lr*1e3)).^2-0.5*((y-y0)/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2) + vy*B0*1e-9);
    x_res(5) = (-e/me)*((y/(lr*1e3)^2).*phi0.*exp(-0.5*((x-x0)/(lr*1e3)).^2-0.5*((y-y0)/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2) - vx*B0*1e-9);
    x_res(6) = (-e/me)*((z/(lz*1e3)^2).*phi0.*exp(-0.5*((x-x0)/(lr*1e3)).^2-0.5*((y-y0)/(lr*1e3)).^2-0.5*(z/(lz*1e3)).^2));                                              
end                
    
