%% What to do
%nParticles = 1000;
rzbinning = 1;
xyzbinning = 1;
doPlot = 1;
toPlot = [2];
doSave = 1;
isInvisible = 0;
doClear = 1;

saveIn = '/Users/Cecilia/Research/EH/TestParticleSimulation/'; % on local
saveIn = '/home/cecilia/Research/ElectronHoleSimulation/SavedData/'; % on spis

global B0 lr lz phi0 L rlim zlim


% Choose and loop between models
for modelnr = modelNr
%try
ExB.model;

% Define time to run each particle
L = 8*lz;  % length of box
T = veh/L; % time for eh to pass box !!! This is waaay too long.
T = 0.01; % this gives 97975 timesteps for B0 = ?

% Set up simulation timestepping
%units=irf_units;
units.e = 1.6022e-19;
units.me = 9.10939999999e-31;
fce=units.e*B0*1e-9/units.me/2/pi; % gyro freq. Hz
tce=1/fce; % gyro period. s
nce = T/tce; % number of gyroperiods
ntpertce = 8000; % number of timesteps per gyroperiod
nt = round(ntpertce*nce); % number of timesteps
t=linspace(0,T,nt)'; % time vector
dt=diff(t(1:2)); % time step

%% Set up grid.
% number of bins
nx = 80;
ny = 80;
nz = 80; 
nr = 40; 
% edges of box [km]
zlim = L/2; 
rlim = 4*lr; 
% centers of bins
xGrid = linspace(-rlim,rlim,nx+1); 
yGrid = linspace(-rlim,rlim,ny+1); 
zGrid = linspace(-zlim,zlim,nz+1); 
rGrid = linspace(0,rlim,nr+1); 
% bin size
dx = diff(xGrid(1:2)); 
dy = diff(yGrid(1:2)); 
dz = diff(zGrid(1:2)); 
dr = diff(rGrid(1:2)); 
% edges between bins
xSurf = xGrid(1:end-1)+dx/2; 
ySurf = yGrid(1:end-1)+dy/2;  
zSurf = zGrid(1:end-1)+dz/2; 
rSurf = rGrid(1:end-1)+dr/2;

%% Initialize variables
% ask for number of particles if not already defined
if ~exist('nParticles','var'); nParticles = input('nParticles = '); end
nOver = zeros(nParticles,1); % number of overshoots the electron makes
nBounce = zeros(nParticles,1); % number of overshoots the electron makes

vaz1 = cell(nr,nz); % azimuthal velocity
vaz2 = cell(nx,ny,nz); % azimuthal velocity
numbrz = zeros(nr,nz); % number of particles per bin
numbxyz = zeros(nx,ny,nz); % number of particles per bin

% Initialize starting position
%x0 = rand(nParticles,1)*rlim.*sign(rand(nParticles,1)-0.5); 
%y0 = rand(nParticles,1)*rlim.*sign(rand(nParticles,1)-0.5); 
%z0 = rand(nParticles,1)*zlim.*sign(rand(nParticles,1)-0.5); 

x0 = -rlim + 2*rlim*rand(nParticles,1); 
y0 = -rlim + 2*rlim*rand(nParticles,1); 
z0 = -zlim + 2*zlim*rand(nParticles,1); 

% Extra starting position array for particles that overshoot. They will be
% sent to the other side of the box but with new random starting positions
% in x and y.
% x01 = nans(nParticles,1); % Save starting positions here after first overshoot
% x02 = nans(nParticles,1); % second overshoot
% x03 = nans(nParticles,1); % third overshoot
% x04 = nans(nParticles,1); % fourth overshoot, skip fifth etc.
% y01 = nans(nParticles,1); % Save starting positions here after first overshoot
% y02 = nans(nParticles,1); % second overshoot
% y03 = nans(nParticles,1); % third overshoot
% y04 = nans(nParticles,1); % fourth overshoot

% Or just make an array that is built up after hand
x0overshoot = [];
y0overshoot = [];
z0overshoot = [];
vx0overshoot = [];
vy0overshoot = [];
vz0overshoot = [];

%x0 = abs(x0);
%y0 = abs(y0);

% Initialize starting velocity
vtpar = cn_eV2v(Tpar,'ev'); % eV -> km/s (parallel thermal velocity)
vtper = cn_eV2v(Tper,'ev'); % eV -> km/s (perpendicular thermal velocity)
vx0 = vtper*randn(nParticles,1)/sqrt(2);
vy0 = vtper*randn(nParticles,1)/sqrt(2);
vz0 = vtpar*randn(nParticles,1)/sqrt(2)-veh; % shift parallel velocity by electron hole velocity

% Initial acceleration, for when using leapfrog
ax0 = (-units.e/units.me)*((x0/(lr*1e3)^2).*phi0.*exp(-0.5*(x0/(lr*1e3)).^2-0.5*(y0/(lr*1e3)).^2-0.5*(z0/(lz*1e3)).^2) + vy0*B0*1e-9);
ay0 = (-units.e/units.me)*((y0/(lr*1e3)^2).*phi0.*exp(-0.5*(x0/(lr*1e3)).^2-0.5*(y0/(lr*1e3)).^2-0.5*(z0/(lz*1e3)).^2) - vx0*B0*1e-9);
az0 = (-units.e/units.me)*((z0/(lz*1e3)^2).*phi0.*exp(-0.5*(x0/(lr*1e3)).^2-0.5*(y0/(lr*1e3)).^2-0.5*(z0/(lz*1e3)).^2));
                
if 0%leapfrog
    vx0 = vx0 + ax0*dt/2;
    vy0 = vy0 + ay0*dt/2;
    vz0 = vz0 + az0*dt/2;
end

% Initializing variable arrays
if 0
x=zeros(nt,1);
y=zeros(nt,1);
z=zeros(nt,1);
vx=zeros(nt,1);
vy=zeros(nt,1);
vz=zeros(nt,1);
dvx=zeros(nt,1);
dvy=zeros(nt,1);
dvz=zeros(nt,1);
ax=zeros(nt,1);
ay=zeros(nt,1);
az=zeros(nt,1);
end

mVaz=zeros(nParticles,1);
    
%% Run through particles
%wait = waitbar(0,'Calculating trajectories, please wait...');
%nwait = 10;

% Terminate integration when the particle passes outside the box in
% z-direction
options = odeset('Events',@ExB.events);%,'OutputSel',1,'Refine',refine);

tic
for jj=1:nParticles 
    tint = [0 T];
 %   try
        %if mod(jj,round(nParticles/nwait)) == 0; waitbar(jj/nParticles,wait); end % waittick each 100 particle            

        % Initial positions and velocities
        % (input is in km and km/s)
        x(1) = x0(jj)*1e3;   % m
        y(1) = y0(jj)*1e3;   % m 
        z(1) = z0(jj)*1e3;   % m
        vx(1) = vx0(jj)*1e3; % m/s
        vy(1) = vy0(jj)*1e3; % m/s
        vz(1) = vz0(jj)*1e3; % m/s                       
        
        nOverInd = 0; % number of overshoots
        nBounceInd = 0; % number of bounces
        
        % Integrate trajectory
        if 0 % Straightforward integration
            for ii = 1:(nt-1)              

                    vx(ii+1)=vx(ii)+dt*(-units.e/units.me)*((x(ii)/(lr*1e3)^2).*phi0.*exp(-0.5*(x(ii)/(lr*1e3)).^2-0.5*(y(ii)/(lr*1e3)).^2-0.5*(z(ii)/(lz*1e3)).^2) + vy(ii)*B0*1e-9);
                    vy(ii+1)=vy(ii)+dt*(-units.e/units.me)*((y(ii)/(lr*1e3)^2).*phi0.*exp(-0.5*(x(ii)/(lr*1e3)).^2-0.5*(y(ii)/(lr*1e3)).^2-0.5*(z(ii)/(lz*1e3)).^2) - vx(ii)*B0*1e-9);
                    vz(ii+1)=vz(ii)+dt*(-units.e/units.me)*((z(ii)/(lz*1e3)^2).*phi0.*exp(-0.5*(x(ii)/(lr*1e3)).^2-0.5*(y(ii)/(lr*1e3)).^2-0.5*(z(ii)/(lz*1e3)).^2));

                    x(ii+1)=x(ii)+vx(ii)*dt;    
                    y(ii+1)=y(ii)+vy(ii)*dt;   
                    z(ii+1)=z(ii)+vz(ii)*dt;                        

                % Move to other side of box and generate new x y if electron 
                % overshoots
                if z(ii+1)>L*1e3/2 
                    x(ii+1) = -rlim + 2*rlim *rand(1,1); % new x position
                    y(ii+1) = -rlim + 2*rlim *rand(1,1); % new y position
                    z(ii+1) = z(ii+1)-L*1e3;                
                    nOverInd = nOverInd  + 1;
                    x0overshoot = [x0overshoot x(ii+1)];
                    y0overshoot = [y0overshoot y(ii+1)];
                    z0overshoot = [z0overshoot z(ii+1)];
                    vx0overshoot = [vx0overshoot vx0(jj)];
                    vy0overshoot = [vy0overshoot vy0(jj)];
                    vz0overshoot = [vz0overshoot vz0(jj)];
                elseif z(ii+1)<-L*1e3/2 
                    x(ii+1) = -rlim + 2*rlim *rand(1,1); % new x position
                    y(ii+1) = -rlim + 2*rlim *rand(1,1); % new y position
                    z(ii+1) = z(ii+1)+L*1e3;
                    nOverInd = nOverInd  + 1;                
                    x0overshoot = [x0overshoot x(ii+1)];
                    y0overshoot = [y0overshoot y(ii+1)];
                    z0overshoot = [z0overshoot z(ii+1)];
                    vx0overshoot = [vx0overshoot vx0(jj)];
                    vy0overshoot = [vy0overshoot vy0(jj)];
                    vz0overshoot = [vz0overshoot vz0(jj)];
                end
            end
        end
        if 1 % Matlab ode solver
            x = [];
            y = [];
            z = [];
            vx = [];
            vy = [];
            vz = [];
            
             % Initial positions and velocities
            % (input is in km and km/s)
            x(1) = x0(jj)*1e3;   % m
            y(1) = y0(jj)*1e3;   % m 
            z(1) = z0(jj)*1e3;   % m
            vx(1) = vx0(jj)*1e3; % m/s
            vy(1) = vy0(jj)*1e3; % m/s
            vz(1) = vz0(jj)*1e3; % m/s   
            
            tend = 0;
            
            % If the integration terminats beforehand, due to the passing
            % of the particle outside the box, run the integration again 
            % with the same initialv velocities, but new position at other 
            % edge of box. This is done until the whole time-interval is 
            % used.
            
            while tend < T
            % Initial positions and velocities
            % (input is in km and km/s)
            [t,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x(end);y(end);z(end);vx(1);vy(1);vz(1)],options);
            x = [x; x_sol(2:end,1)]; 
            y = [y; x_sol(2:end,2)];
            z = [z; x_sol(2:end,3)];
            vx = [vx; x_sol(2:end,4)];
            vy = [vy; x_sol(2:end,5)];
            vz = [vz; x_sol(2:end,6)];
            tend = t(end);
            if tend >= T; break; end            
            newx = (-rlim + 2*rlim*rand(1,1))*1e3;
            newy = (-rlim + 2*rlim*rand(1,1))*1e3;
            if z(end)>0; newz = z(end)-0.5*L*1e3; % move to bottom of box                
            else newz = z(end)+0.5*L*1e3; % move to top of box
            end
            x = [x; newx];
            y = [y; newy];
            z = [z; newz];
            
            vx = [vx; vx(1)];
            vy = [vy; vy(1)];
            vz = [vz; vz(1)];
            
            nOverInd = nOverInd  + 1;
            x0overshoot = [x0overshoot newx];
            y0overshoot = [y0overshoot newy];
            z0overshoot = [z0overshoot newz];            
            
            tint = [tend T];
            end
             
        end  
       % plot3(x*1e-3,y*1e-3,z*1e-3)
       % set(gca,'xlim',[-rlim rlim],'ylim',[-rlim rlim])
        1;
%             % Count the number of bounces (= vz reversals) the electron does
%             if sign(vz(ii+1)) ~= sign(vz(ii)) 
%                 nBounceInd = nBounceInd + 1;
%             end        
        
        % Save bounce and overshoot statistics
        %[~,~,zbounce] = find(diff(sign(vz))~=0);        
        %nBounceInd = numel(zbounce);
        
        nOver(jj) = nOverInd;
        %nBounce(jj) = nBounceInd;        
        nBounce(jj) = numel(find(diff(sign(vz))~=0));        

        % Get azimuthal velocity and radial position
        xg0 = find(x>=0); % right half xy-plane
        xs0 = find(x<0); % left half xy-plane
        th = zeros(size(x));
        th(xg0) = atand(y(xg0)./x(xg0));
        th(xs0) = atand(y(xs0)./x(xs0)) + 180;        
        vaztmp = -vx.*sind(th) + vy.*cosd(th); % azimuthal velocity 
        r = sqrt(x.^2 + y.^2);
        
        % Save average azimuthal velocity over an even number of
        % gyroperiods
%         restGyroperiod = mod(nce,1);
%         evenNumberOfGyroperiods = nce - restGyroperiod;
%         timestepsGyroperiods = evenNumberOfGyroperiods*ntpertce;
%         restTimesteps = nt-timestepsGyroperiods;
%         gyroIndex = (1:timestepsGyroperiods) + fix(restTimesteps/2);
%         mVaz(jj) = mean(vaztmp(gyroIndex));
        
        if rzbinning % rz binning
            % Put values into discrete bins
            [nzz,binz] = histc(z,zSurf*1e3);
            [nrr,binr] = histc(r,rSurf*1e3);
            [arz,brz,crz] = unique([binr binz],'rows');

            % Remove indices when the particle is outside the box
            remr = find(arz(:,1)==0);
            remz = find(arz(:,2)==0);
            remrz = unique([remr;remz]);
            arz(remrz,:) = [];
            brz(remrz) = [];                        

            nun = max(crz);            
            for ui = 1:nun                  
                ic = find(crz==ui);
                ir = max(binr(ic)); if ir == 0; continue; end
                iz = max(binz(ic)); if iz == 0; continue; end
                vaz1{ir,iz} = [vaz1{ir,iz}  mean(vaztmp(ic))];            
            end   
        end 
        
        if xyzbinning % xyz binning            
            [nxx,binx] = histc(x,xSurf*1e3);
            [nyy,biny] = histc(y,ySurf*1e3);
            [nzz,binz] = histc(z,zSurf*1e3);
            [axyz,bxyz,cxyz] = unique([binx biny binz],'rows');
            
            % Remove indices when the particle is outside the box
            remxyzx = find(axyz(:,1)==0); % find indices
            remxyzy = find(axyz(:,2)==0);
            remxyzz = find(axyz(:,3)==0);
            remxyz = unique([remxyzx;remxyzy;remxyzz]); % find unique pairs
            axyz(remxyz,:) = []; % remove them
            bxyz(remxyz) = [];
         
            nun = max(cxyz);         
            for ui = 1:nun      
                ic = find(cxyz==ui);
                ix = max(binx(ic)); if ix == 0; continue; end
                iy = max(biny(ic)); if iy == 0; continue; end
                iz = max(binz(ic)); if iz == 0; continue; end
                vaz2{ix,iy,iz} = [vaz2{ix,iy,iz}  mean(vaztmp(ic))];            
            end                    
        end  
        
        if mod(jj,100) == 0; clc; disp([num2str(jj) '/' num2str(nParticles)]); end 
    %catch me
    %   disp(['Error. Skipping particle #' num2str(jj) '.'])        
    %   continue;
    %end    
end

% Change saved overshoot velocities to km/s
vx0overshoot = vx0overshoot*1e-3;
vy0overshoot = vy0overshoot*1e-3;
vz0overshoot = vz0overshoot*1e-3;

%close(wait)
tTot = toc;
disp(['elapsed time = ' num2str(tTot) ' s'])
disp(['time per particle = ' num2str(tTot/nParticles) ' s'])

%% Make average azimuthal velocity
if rzbinning % rz binning system
    vazmat = cellfun(@mean,vaz1);
    numbrz = cellfun(@numel,vaz1);
end

if xyzbinning % xyz binning system
    vazmat2 = cellfun(@mean,vaz2);
    numbxyz = cellfun(@numel,vaz2);
    
    % Average over whole column
    vazmat2yz = squeeze(nanmean(vazmat2,1));
    vazmat2xz = squeeze(nanmean(vazmat2,2));
    vazmat2xy = squeeze(nanmean(vazmat2,3));
    numbyz = squeeze(nanmean(numbxyz,1));
    numbxz = squeeze(nanmean(numbxyz,2));
    numbxy = squeeze(nanmean(numbxyz,3));
    
    % Average over a couple of bins centered about x,y,z=0    
    if mod(nx+1,2)==0 % uneven number of xbins
        nxc = ceil(nx/2);
        xind = (nxc-3):(nxc+4);
    else % even number of bins
        nxc = ceil(nx/2);
        xind = (nxc-2):(nxc+4);
    end
    if mod(ny+1,2)==0 % uneven number of xbins
        nyc = ceil(ny/2);
        yind = (nyc-3):(nyc+4);
    else % even number of bins
        nyc = ceil(ny/2);
        yind = (nyc-2):(nyc+4);
    end
    if mod(nz+1,2)==0 % uneven number of xbins
        nzc = ceil(nz/2);
        zind = (nzc-3):(nzc+4);
    else % even number of bins
        nzc = ceil(nz/2);
        zind = (nzc-2):(nzc+4);        
    end
    zCenter = zGrid(zind([1 end]))*1e-3;    
    vazmat2yzCenter = squeeze(nanmean(vazmat2(xind,:,:),1));
    vazmat2xzCenter = squeeze(nanmean(vazmat2(:,yind,:),2));
    vazmat2xyCenter = squeeze(nanmean(vazmat2(:,:,zind),3));       
end

%% Plot
if doPlot;
    hh=1;
    for plotNr = 1:numel(toPlot);                
        irf_colormap('poynting');
        close
        toplot = toPlot(plotNr);
        ExB.plot; 
        %cn.print([name '_' num2str(plotNr)]);        
        cn.print([name '_' datestr(now,'yyyymmddTHHMMSS') '_' num2str(plotNr)]);        
        hh=hh+1;        
    end
end

%% Save 
if doSave; ExB.save; end

%clear nParticles
%catch me
%    continue;
%end
end
