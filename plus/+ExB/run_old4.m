%% What to do
%nParticles = 400000;
rzbinning = 1;
xyzbinning = 1;
doPlot = 0;
toPlot = [2 3 4];
doSave = 0;
isInvisible = 0;
doClear = 1;
isMother = 1;
initialization = 'r';


global B0 lr lz phi0 L rlim zlim


%% Choose and loop between models
for modelnr = 2
%try
ExB.model;

% Define time to run each particle
L = 8*lz;  % length of box
T = 0.01; % this gives 97975 timesteps for B0 = ?

%% Set up grid.
% number of bins
nx = 120;
ny = 120;
nz = 120; 
nr = 60; 
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

% Initialize starting velocity
vtpar = cn_eV2v(Tpar,'ev'); % eV -> km/s (parallel thermal velocity)
vtper = cn_eV2v(Tper,'ev'); % eV -> km/s (perpendicular thermal velocity)
vx0 = vtper*randn(nParticles,1)/sqrt(2);
vy0 = vtper*randn(nParticles,1)/sqrt(2);
vz0 = vtpar*randn(nParticles,1)/sqrt(2)-veh; % shift parallel velocity by electron hole velocity

if 0
vperp = sqrt(vx0.^2+vy0.^2); % km/s
me = 9.10939999999e-31;
e = 1.6022e-19;
rg0 = me*vperp/e/(B0*1e-9); % km, the starting position should be shifted 
                           % with this distance
vxg0 = find(vx0>=0); % right half xy-plane
vxs0 = find(vx0<0); % left half xy-plane
vth = zeros(size(vx0));
vth(vxg0) = atand(vy0(vxg0)./vx0(vxg0));
vth(vxs0) = atand(vy0(vxs0)./vx0(vxs0)) + 180; 
vth = vth-90;



% Initialize starting position
%x0 = -rlim + 2*rlim*rand(nParticles,1); 
%y0 = -rlim + 2*rlim*rand(nParticles,1); 
%z0 = -zlim + 2*zlim*rand(nParticles,1); 

%x0 = 0.5*rlim*randn(nParticles,1); 
%y0 = 0.5*rlim*randn(nParticles,1);
%r0 = sqrt(x0.^2+y0.^2);


switch initialization
    case 'normal'
        r0 = rlim*randn(nParticles,1);
        angle = 180*rand(nParticles,1);
        x0 = r0.*cosd(angle);
        y0 = r0.*sind(angle);
        z0 = -zlim + 2*zlim*rand(nParticles,1); 
        m = 2*pi*rlim*exp(r0.^2/rlim.^2/2);
    case 'r'
        r0 = rlim*rand(nParticles,1);
        angle = 360*rand(nParticles,1);
        x0 = r0.*cosd(angle);
        y0 = r0.*sind(angle);
        z0 = -zlim + 2*zlim*rand(nParticles,1); 
        m = r0/rlim; 
end

x1 = x0 + rg0.*cosd(vth);
y1 = y0 + rg0.*sind(vth);
end

[xc,yc,x0,y0,m] = ExB.newposition(vx0,vy0);
z0 = -zlim + 2*zlim*rand(nParticles,1); 
z0 = -zlim + 2*zlim*rand(nParticles,1); 

if 0 % plot starting gyrocenters particle position and velocity in xy plane
    plot(xc,yc,'b.',x0,y0,'g.'); hold on;
    quiver(x1,y1,vx0,vy0)
    set(gca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])
    axis equal
end
if 1
    circle = 0:pi/10:2*pi;
    for ip=1:nParticles
        hs=patch(sin(circle)+xc(ip),cos(circle)+yc(ip),'b','edgecolor','none');
        alpha(hs,m(ip)/max(m)*0.1)
    end
end
% Initialize particle mass, normalized to me
%m = ones(nParticles,1); % all have same weight, this is for uniform particle initialization
%m = 0.5./rlim./randn(nParticles,1); 
%m = 1-2*rlim*randn(nParticles,1); 


% Initialize mean velocity sum and mass sum
sumMVrz = zeros(nr,nz);
sumMrz = zeros(nr,nz);
sumMVxyz = zeros(nx,ny,nz);
sumMxyz = zeros(nx,ny,nz);

% Position array for overhooting particles that is built up after hand
x0overshoot = [];
y0overshoot = [];
z0overshoot = [];
vx0overshoot = [];
vy0overshoot = [];
vz0overshoot = [];


mVaz=zeros(nParticles,1);
    
%% Run through particles
%wait = waitbar(0,'Calculating trajectories, please wait...');
nwait = 100;

% Terminate integration when the particle passes outside the box in
% z-direction
options = odeset('Events',@ExB.events,'InitialStep',2.5e-5);%,'OutputSel',1,'Refine',refine);
fprintf('%s\n',['Integration started at ' datestr(now,'yyyymmddTHHMMSS')])
fprintf('%s\n',['Estimated stop time:   ' datestr(now+nParticles*T*0.08/0.01/60/60/24,'yyyymmddTHHMMSS')])
n = 0;
msgProgress = ['Progress: 0/' num2str(nParticles)];
fprintf('%s',msgProgress)
tic
mm = [];
for jj=1:nParticles     
    tint = [0 T];
 %   try
        %if mod(jj,round(nParticles/nwait)) == 0; waitbar(jj/nParticles,wait); end % waittick each 100 particle                                     
        
        nOverInd = 0; % number of overshoots
        nBounceInd = 0; % number of bounces
                    
        % Initial positions and velocities
        % (input is in km and km/s)
        %clear x y x vx vy vz
        x = x0(jj)*1e3;   % m
        y = y0(jj)*1e3;   % m 
        z = z0(jj)*1e3;   % m
        vx = vx0(jj)*1e3; % m/s
        vy = vy0(jj)*1e3; % m/s
        vz = vz0(jj)*1e3; % m/s   
        newm = m(jj);
                
        % If the integration terminats beforehand, due to the passing
        % of the particle outside the box, run the integration again 
        % with the same initialv velocities, but new position at other 
        % edge of box. This is done until the whole time-interval is 
        % used.
        
        % Integrate trajectory
        % Matlab ode solver        
        tend = 0;
        while tend < T            
            [t,x_sol] = ode45(@ExB.EquationOfMotion,tint,[x(end);y(end);z(end);vx(1);vy(1);vz(1)],options);
            x = [x; x_sol(2:end,1)]; 
            y = [y; x_sol(2:end,2)];
            z = [z; x_sol(2:end,3)];
            vx = [vx; x_sol(2:end,4)];
            vy = [vy; x_sol(2:end,5)];
            vz = [vz; x_sol(2:end,6)];
            mm = [mm; repmat(newm,numel(t),1)];
            tend = t(end);
            if tend >= T; break; end  
            % The electron has overshoot, so new random x,y, z moves to
            % other end of box and vx,vy,vz is the same as the initial
            % values
            %newx = (-rlim + 2*rlim*rand(1,1))*1e3;
            %newy = (-rlim + 2*rlim*rand(1,1))*1e3;
            [newxc,newyc,newx0,newy0,newm] = ExB.newposition(vx0(jj),vy0(jj));
            %newr = rlim*randn(1,1);
            %newangle = 360*rand(1,1);
            %newx = newr*cosd(newangle)*1e3;
            %newy = newr*sind(newangle)*1e3;
            if z(end)>0; newz = z(end)-L*1e3; % move to bottom of box                
            else newz = z(end)+L*1e3; % move to top of box
            end
            x = [x; newx0*1e3];
            y = [y; newy0*1e3];
            z = [z; newz];
            vx = [vx; vx(1)];
            vy = [vy; vy(1)];
            vz = [vz; vz(1)];
            %m(jj) = [mm;newm];

            nOverInd = nOverInd  + 1;
            
            % Save new positions in array to see afterwards
            x0overshoot = [x0overshoot newxc*1e3];
            y0overshoot = [y0overshoot newyc*1e3];
            z0overshoot = [z0overshoot newz];          
            
            vx0overshoot = [vx0overshoot vx(1)];
            vy0overshoot = [vy0overshoot vy(1)];
            vz0overshoot = [vz0overshoot vz(1)];
            % The rest of the time interval to integrate over again
            tint = [tend T];
        end
                
        nOver(jj) = nOverInd;               
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
                mass = 2*pi*rlim*1e3*exp(rSurf(ir).^2/rlim.^2/2);
                %mass = m(jj);
                %mass = mean(mm(ic));
                %mass = 1;
                mass = rSurf(ir)/rlim;
                sumMVrz(ir,iz) = sumMVrz(ir,iz) + mean(vaztmp(ic)).*mass;
                sumMrz(ir,iz) = sumMrz(ir,iz) + mass;
                %vaz1{ir,iz} = [vaz1{ir,iz}  mean(vaztmp(ic))];            
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
                mass = 2*pi*rlim*1e3*exp((xSurf(ix).^2+ySurf(iy).^2)/rlim.^2/2);
                mass = m(jj);
                mass = mean(mm(ic));
                %mass = 1;
                mass = sqrt(xSurf(ix).^2+ySurf(iy).^2)/rlim;
                sumMVxyz(ix,iy,iz) = sumMVxyz(ix,iy,iz) + mean(vaztmp(ic)).*mass;
                sumMxyz(ix,iy,iz) = sumMxyz(ix,iy,iz) + mass;                
                %vaz2{ix,iy,iz} = [vaz2{ix,iy,iz}  mean(vaztmp(ic))];            
            end                    
        end  
        
        if mod(jj,round(nParticles/nwait)) == 0; 
            n = numel(msgProgress);
            msgProgress = ['Progress: ' num2str(jj) '/' num2str(nParticles)]; 
            msgReverse = repmat(sprintf('\b'), 1, n);
            fprintf('%s',msgReverse)
            fprintf('%s',[msgProgress]); 
            
        end 
    %catch me
    %   disp(['Error. Skipping particle #' num2str(jj) '.'])        
    %   continue;
    %end    
end
%fprintf([msgReverse '\b'])
fprintf('\n%s\n',['Integration stopped at ' datestr(now,'yyyymmddTHHMMSS')])
%close(wait)
tTot = toc;
fprintf('%s\n',['Elapsed time: ' num2str(tTot) ' s'])
fprintf('%s\n',['Time per particle: ' num2str(tTot/nParticles) ' s'])
tic
fprintf('%s\n',['Finish binning...'])

%% Make average azimuthal velocity
if rzbinning % rz binning system
    %vazmat = cellfun(@mean,vaz1);
    %numbrz = cellfun(@numel,vaz1);
    meanVrz = sumMVrz./sumMrz;
    vazmat = meanVrz;
    numbrz = sumMrz;
    
end

if xyzbinning % xyz binning system
    %vazmat2 = cellfun(@mean,vaz2);
    %numbxyz = cellfun(@numel,vaz2);
    meanVxyz = sumMVxyz./sumMxyz;
    vazmat2 = meanVxyz;
    numbxyz = sumMxyz;
    
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
tBin = toc;
fprintf('%s\n',['Elapsed time: ' num2str(tBin) ' s'])

%% Save 
if doSave; ExB.save; end

%% Plot
if doPlot;
    hh=1;
    try
    for plotNr = 1:numel(toPlot);                
        irf_colormap('poynting');
        close
        toplot = toPlot(plotNr);
        ExB.plot; 
        %cn.print([name '_' num2str(plotNr)]);        
        cn.print([name '_' datestr(now,'yyyymmddTHHMMSS') '_' num2str(plotNr)]);        
        hh=hh+1;        
    end
    catch me
        disp('Could not plot.')
    end
end


%clear nParticles
%catch me
%    continue;
%end
end
