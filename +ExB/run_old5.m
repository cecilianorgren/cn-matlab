%% What to do
%nParticles = 400000;
rzbinning = 1;
xyzbinning = 1;
doPlot = 1;
toPlot = [2 3 4 5];
doSave = 1;
isInvisible = 0;
doClear = 1;
isMother = 1;
initialization = 'r';
cellMassgen = {'bin','hm','pre','uni'};
massgen = cellMassgen{3};


global B0 lr lz phi0 L rlim zlim


%% Choose and loop between models
for modelnr = modelNr
ExB.model;
% Define time to run the particles
T = 0.01;

%% Set up grid.
L = 8*lz;  % length of box
W = 8*lr;  % width of box
% number of bins
nx = 120;
ny = 120;
nz = 120; 
nr = 60; 
% edges of box [km]
zlim = L/2; 
rlim = W/2; 
% centers of bins
xGrid = linspace(-rlim,rlim,nx+1); 
yGrid = linspace(-rlim,rlim,ny+1); 
zGrid = linspace(-zlim,zlim,nz+1); 
rGrid = linspace(0,rlim,nr+1); 
% edges between bins
xSurf = xGrid(1:end-1)+diff(xGrid(1:2))/2; 
ySurf = yGrid(1:end-1)+diff(yGrid(1:2))/2;  
zSurf = zGrid(1:end-1)+diff(zGrid(1:2))/2; 
rSurf = rGrid(1:end-1)+diff(rGrid(1:2))/2;

%% Initialize variables
% ask for number of particles if not already defined
if ~exist('nParticlesBox','var'); nParticles = input('nParticles = '); end

% Initialize starting velocity
vtpar = cn_eV2v(Tpar,'ev'); % eV -> km/s (parallel thermal velocity)
vtper = cn_eV2v(Tper,'ev'); % eV -> km/s (perpendicular thermal velocity)

% shift parallel velocity by electron hole velocity
vz0box = vtpar*randn(nParticlesBox,1)/sqrt(2)-veh; 

% See how many particles of a certain velocity should pass by the box,
% given that they are one at a time inside the box.
turns = vz0box*1e3*T/(L*1e3);
direction = sign(turns);
partialTurns = direction.*mod(abs(turns),1);
completeTurns = turns - partialTurns;
nParticlesFlow = sum(abs(completeTurns))+numel(partialTurns);

% Add that number of particles with that velocity to a new vector 
if 1 % Preallocated array
    vz0flow = zeros(nParticlesFlow,1);
    tintFlow = zeros(nParticlesFlow,1);    
    tic
    index = 1;
    for kk = 1:nParticlesBox        
        vz0flow(index:index+abs(completeTurns(kk))) = ones(abs(completeTurns(kk))+1,1)*vz0box(kk);         
        zPartial = partialTurns(kk)*L;
        tPartial = zPartial/vz0box(kk);
        tintFlow(index:index+abs(completeTurns(kk))-1) = ones(abs(completeTurns(kk)),1)*T; 
        tintFlow(index+abs(completeTurns(kk))) = tPartial;
        index = index + abs(completeTurns(kk)) + 1;
    end
    toc   
end

tint = [ones(nParticlesBox,1)*T; tintFlow];

% The total number of particles
nParticles = nParticlesBox + nParticlesFlow;

% Velocities
vx0 = vtper*randn(nParticles,1)/sqrt(2);
vy0 = vtper*randn(nParticles,1)/sqrt(2);
vz0 = [vz0box; vz0flow];

nBounce = zeros(nParticles,1); % number of bounces the electrons make

% Generate z-positions for all the particles
z0flow = ones(nParticlesFlow,1)*zlim;
% particles with positive velocity starts at bottom of box
z0flow(find(vz0flow>0)) = -z0flow(find(vz0flow>0));
z0box = -zlim + 2*zlim*rand(nParticlesBox,1); 
z0flow = ones(nParticlesFlow,1)*zlim;
z0 = [z0box; z0flow];

% Generate xy-positions and mass for all the particles
[xc,yc,x0,y0,m] = ExB.newposition(vx0,vy0,initialization);

if 0 % plot 'density'
    circle = 0:pi/10:2*pi;
    for ip=1:nParticles
        hs=patch(sin(circle)+xc(ip),cos(circle)+yc(ip),'b','edgecolor','none');
        alpha(hs,m(ip)/max(m)*0.1)
    end
end
if 0 % plot starting gyrocenters particle position and velocity in xy plane
    plot(xc,yc,'b.',x0,y0,'g.'); hold on;
    quiver(x1,y1,vx0,vy0)
    set(gca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])
    axis equal
end

% Initialize mean velocity sum and mass sum
sumMVrz = zeros(nr,nz);
sumMrz = zeros(nr,nz);
sumMVxyz = zeros(nx,ny,nz);
sumMxyz = zeros(nx,ny,nz);
sumMVXxyz = zeros(nx,ny,nz);
sumMVYxyz = zeros(nx,ny,nz);
sumMVZxyz = zeros(nx,ny,nz);

%% Run through particles
nwait = 100;
fprintf('%s\n',['Integration started at ' datestr(now,'yyyymmddTHHMMSS')])
fprintf('%s\n',['Estimated stop time:   ' datestr(now+nParticles*T*0.03/0.01/60/60/24,'yyyymmddTHHMMSS')])
msgProgress = ['Progress: 0/' num2str(nParticles)];
fprintf('%s',msgProgress)
tic

tend = zeros(nParticles,1);

% Terminate integration when the particle passes outside the box in
% z-direction
options = odeset('Events',@ExB.events,'InitialStep',2.5e-5);%,'OutputSel',1,'Refine',refine);
% Integration time

for jj=1:nParticles           
    nOverInd = 0; % number of overshoots
    nBounceInd = 0; % number of bounces

    % Initial positions and velocities                                   
    x_init = [x0(jj);y0(jj);z0(jj);vx0(jj);vy0(jj);vz0(jj)]*1e3; % m, m/s

    % Integrate trajectory
    % Matlab ode solver :             
    % If the integration terminats beforehand, due to the passing
    % of the particle outside the box, the integration stops, this is
    % defined in options.
    [t,x_sol] = ode45(@ExB.EquationOfMotion,[0 tint(jj)],x_init,options);
    x = x_sol(:,1);
    y = x_sol(:,2);
    z = x_sol(:,3);
    vx = x_sol(:,4);
    vy = x_sol(:,5);
    vz = x_sol(:,6);
    tend(jj) = t(end);
    
    % See if parallel velocity changes sign = bounce
    nBounce(jj) = numel(find(diff(sign(vz))~=0)); 

    % Get azimuthal velocity and radial position
    xg0 = find(x>=0); % right half xy-plane
    xs0 = find(x<0); % left half xy-plane
    th = zeros(size(x));
    th(xg0) = atand(y(xg0)./x(xg0));
    th(xs0) = atand(y(xs0)./x(xs0)) + 180;        
    vaztmp = -vx.*sind(th) + vy.*cosd(th); % azimuthal velocity 
    r = sqrt(x.^2 + y.^2);              

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

            switch massgen
                case 'pre'; mass = m(jj);
                case 'hm';  mass = 2*pi*rlim*1e3*exp(rSurf(ir).^2/rlim.^2/2);
                case 'uni'; mass = 1;
                case 'bin'; mass = rSurf(ir)/rlim;
            end 
            sumMVrz(ir,iz) = sumMVrz(ir,iz) + mean(vaztmp(ic)).*mass;
            sumMrz(ir,iz) = sumMrz(ir,iz) + mass;            
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
            switch massgen
                case 'pre'; mass = m(jj);
                case 'hm';  mass = 2*pi*rlim*1e3*exp((xSurf(ix).^2+ySurf(iy).^2)/rlim.^2/2);
                case 'uni'; mass = 1;
                case 'bin'; mass = sqrt(xSurf(ix).^2+ySurf(iy).^2)/rlim;
            end                                
            sumMVxyz(ix,iy,iz) = sumMVxyz(ix,iy,iz) + mean(vaztmp(ic)).*mass;
            sumMxyz(ix,iy,iz) = sumMxyz(ix,iy,iz) + mass;                            
            sumMVXxyz(ix,iy,iz) = sumMVXxyz(ix,iy,iz) + mean(vx(ic)).*mass;
            sumMVYxyz(ix,iy,iz) = sumMVYxyz(ix,iy,iz) + mean(vy(ic)).*mass;
            sumMVZxyz(ix,iy,iz) = sumMVZxyz(ix,iy,iz) + mean(vz(ic)).*mass;
        end                    
    end  

    % status update
    if mod(jj,round(nParticles/nwait)) == 0; 
        n = numel(msgProgress);
        msgReverse = repmat(sprintf('\b'), 1, n);
        msgProgress = ['Progress: ' num2str(jj) '/' num2str(nParticles)];         
        fprintf('%s',msgReverse)
        fprintf('%s',msgProgress); 
    end 

end

fprintf('\n%s\n',['Integration stopped at ' datestr(now,'yyyymmddTHHMMSS')])
tTot = toc;
fprintf('%s\n',['Elapsed time: ' num2str(tTot) ' s'])
fprintf('%s\n',['Time per particle: ' num2str(tTot/nParticles) ' s'])

%% Make average azimuthal velocity
if rzbinning % rz binning system
    meanVrz = sumMVrz./sumMrz;
    vazmat = meanVrz;
    numbrz = sumMrz;    
end
if xyzbinning % xyz binning system
    meanVxyz = sumMVxyz./sumMxyz;
    meanVXxyz = sumMVXxyz./sumMxyz;
    meanVYxyz = sumMVYxyz./sumMxyz;
    meanVZxyz = sumMVZxyz./sumMxyz;
    vazmat2 = meanVxyz;
    numbxyz = sumMxyz;
    
    % Average over whole column
    vazmat2yz = squeeze(nanmean(vazmat2,1));
    vazmat2xz = squeeze(nanmean(vazmat2,2));
    vazmat2xy = squeeze(nanmean(vazmat2,3));
    vxmat2xy  = squeeze(nanmean(meanVXxyz,3));
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
    numbyzCenter = squeeze(nanmean(numbxyz(xind,:,:),1));
    numbxzCenter = squeeze(nanmean(numbxyz(:,yind,:),2));
    numbxyCenter = squeeze(nanmean(numbxyz(:,:,zind),3));
end

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

end