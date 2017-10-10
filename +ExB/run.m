%% What to do
%nParticles = 400000;
%parpool local
rzbinning = 1;
xyzbinning = 1;
%doPlot = 1;
toPlot = [2 3 4 5];
%doSave = 0;
isInvisible = 0;
doClear = 1;
doMax = 1;
doMean = 1;
%isMother = 1;
initialization = 'uniform'; % one of 'r', 'uniform', etc...
cellMassgen = {'bin','hm','pre','uni'};
massgen = cellMassgen{3};
doODE45 = 1;
noForBin=0;
if ~exist('doPool','var'); doPool=1; end
if ~exist('noGyrocenter','var'); noGyrocenter=0; end
if ~exist('gyrocenterMass','var'); gyrocenterMass=1; end

%% Choose and loop between models
modeliter = 0;
for modelnr = modelNr
    modeliter = modeliter + 1;
    try
        nParticlesBox = vecParticles(modeliter);
    catch
        disp('Error.')
        nParticlesBox = 100;
    end
ExB.model;
% Define time to run the particles
T = 0.03;

%% Set up grid.
% Width of box should be either decided by the perpendicular thermal
% velocity or gyroradius, or the width of the electron hole
re = sqrt(10*Tper)/B0;
W = 12*lr + 0*re; % width of box, km
L = 14*lz;        % length of box, km
dr = 1; % km
%dr = 0.5;

% displace center of electron hole
%eh_x0 = 2*re*1e3; % m
%eh_y0 = 2*re*1e3; % m

% number of bins
nx = fix(W/dr); if mod(nx,2) ~= 1; nx = nx + 1; end
ny = fix(W/dr); if mod(ny,2) ~= 1; ny = ny + 1; end
nz = 161; 
nr = fix(W/dr/2); 90; 
disp(['nx=' num2str(nx) ', ny=' num2str(ny) ', nz=' num2str(nz)])
% edges of box [km]
zlim = L/2; 
rlim = W/2; 
disp(['x,y=[-' num2str(rlim) ', ' num2str(rlim) ', z=[-' num2str(zlim) ', ' num2str(zlim) ']'])
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
if ~exist('nParticlesBox','var'); nParticlesBox = input('nParticles = '); end

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
% Preallocated array
vz0flow = zeros(nParticlesFlow,1);
tintFlow = zeros(nParticlesFlow,1);        
index = 1;
for kk = 1:nParticlesBox        
    vz0flow(index:index+abs(completeTurns(kk))) = ones(abs(completeTurns(kk))+1,1)*vz0box(kk);         
    zPartial = partialTurns(kk)*L;
    tPartial = zPartial/vz0box(kk);
    tintFlow(index:index+abs(completeTurns(kk))-1) = ones(abs(completeTurns(kk)),1)*T; 
    tintFlow(index+abs(completeTurns(kk))) = tPartial;
    index = index + abs(completeTurns(kk)) + 1;
end   

tint = [ones(nParticlesBox,1)*T; tintFlow];

% The total number of particles
nParticles = nParticlesBox + nParticlesFlow;

% Velocities
vx0 = vtper*randn(nParticles,1)/sqrt(2);
vy0 = vtper*randn(nParticles,1)/sqrt(2);
vz0 = [vz0box; vz0flow];

% Generate z-positions for all the particles
z0flow = ones(nParticlesFlow,1)*zlim;
% particles with positive velocity starts at bottom of box
z0flow(find(vz0flow>0)) = -z0flow(find(vz0flow>0));
z0box = -zlim + 2*zlim*rand(nParticlesBox,1); 
z0 = [z0box; z0flow];

% Generate xy-positions and mass for all the particles
[x1,y1,x2,y2,m1,m2] = ExB.newposition(vx0,vy0,initialization,rlim+3*re,B0);
if noGyrocenter(modeliter)    
    x0 = x2; y0 = y2;
    disp('Moving starting position out to gyroradius.')
else    
    x0 = x1; y0 = y1;
    disp('Using directly generated starting positions.')    
end
if gyrocenterMass;
    m = m2; 
    disp('Using the gyrocenter for deciding mass.')
else
    m = m1; 
    disp('Using the directly generated starting positions for deciding mass.')
end


% Plot starting gyrocenters particle position and velocity in xy plane
if 0 
    plot(xc,yc,'b.',x0,y0,'g.'); hold on;
    quiver(x0,y0,vx0,vy0)
    set(gca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])
    legend('c','0')
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

if doMax
sumMVrzMax = zeros(nr,nz);
sumMrzMax = zeros(nr,nz);
sumMVxyzMax = zeros(nx,ny,nz);
sumMxyzMax = zeros(nx,ny,nz);
sumMVXxyzMax = zeros(nx,ny,nz);
sumMVYxyzMax = zeros(nx,ny,nz);
sumMVZxyzMax = zeros(nx,ny,nz);
end
    

%% Set up parpool
% See if pool already exist.

if ~exist('npl','var'); npl = 2; end % default is two workers
if doPool
checkpool = gcp('nocreate');
isNoPool = isempty(checkpool);
if isNoPool
    pool = parpool(npl);        
    disp(['Opened a new pool with ' num2str(pool.NumWorkers) ' workers.'])
else
    disp(['Pool already exist with ' num2str(checkpool.NumWorkers) ' workers.'])
    if npl == checkpool.NumWorkers; 
        disp('Going ahead with this pool.'); 
    else
        delete(checkpool); disp('Deleted existing pool.')
        pool = parpool(npl);        
        disp(['Opened a new pool with ' num2str(pool.NumWorkers) ' workers.'])
    end
end
end
% Set up the defined number of parallel pools
%[~,hostname]=system('hostname');
%if ~exist('pool') && ~strcmp(hostname(1:4),'moth')
%    pool = parpool(npl);
%end

% nwait = 100;
fprintf('%s\n',['Integration started at ' datestr(now,'yyyymmddTHHMMSS')])
% fprintf('%s\n',['Estimated stop time:   ' datestr(now+nParticles*T*0.03/0.01/60/60/24,'yyyymmddTHHMMSS')])
% msgProgress = ['Progress: 0/' num2str(nParticles)];
% fprintf('%s',msgProgress)
tic
% Terminate integration when the particle passes outside the box in
% z-direction
stopfunction = @(t,y) ExB.events(t,y,L);
options = odeset('Events',stopfunction);%,'InitialStep',2.5e-5,'OutputSel',1,'Refine',refine);
% Integration time

% Divide the total number of particles into 12 groups

indGroup = 1:fix(nParticles/npl):(nParticles+1);
indGroup(end) = nParticles+1;
disp('Number of particles in respective group')
disp(num2str(diff(indGroup)))

cellnBounce = cell(npl,1); % number of bounces the electrons make
celltend = cell(npl,1); % how long each particle runs

% Initialize end position
cellx1= cell(npl,1);
celly1= cell(npl,1);
cellz1= cell(npl,1);
cellx = cell(npl,1);

%% Run through particles
% Parfor-loop
warning off
%if doPool, eval('parfor np = 1:npl'),
%else eval('for np = 1:npl'), end
disp('Progress: ')
parfor np = 1:npl
    parsumMVrz   = zeros(nr,nz);
    parsumMrz    = zeros(nr,nz);
    parsumMVxyz  = zeros(nx,ny,nz);
    parsumMxyz   = zeros(nx,ny,nz);
    parsumMVXxyz = zeros(nx,ny,nz);
    parsumMVYxyz = zeros(nx,ny,nz);
    parsumMVZxyz = zeros(nx,ny,nz);
    
    if doMax
    parsumMVrzMax   = zeros(nr,nz);
    parsumMrzMax    = zeros(nr,nz);
    parsumMVxyzMax  = zeros(nx,ny,nz);
    parsumMxyzMax   = zeros(nx,ny,nz);
    parsumMVXxyzMax = zeros(nx,ny,nz);
    parsumMVYxyzMax = zeros(nx,ny,nz);
    parsumMVZxyzMax = zeros(nx,ny,nz);
    end
    
    nfor = indGroup(np):(indGroup(np+1)-1);
    nParParticles = nfor(end) - nfor(1) + 1;
    
    parnBounce = zeros(nParParticles,1);
    partend = zeros(nParParticles,1);
    partx1= zeros(nParParticles,1);
    party1= zeros(nParParticles,1);
    partz1= zeros(nParParticles,1);
    nProg = 20;
    nPart = 0;    
for jj = nfor 
    nPart = nPart + 1;    
    progg=0;
    if jj/nParParticles>progg;
        fprintf(['w' num2str(np) ':' num2str(progg,'%.f') '% ']);
        progg=progg+0.1;
    end
    %if mod(jj,(nParticles-mod(nParticles,nProg))/nProg) == 0; disp(['Progress: particle ' num2str(nPart) '/' num2str(nParticles)]); end
    %if mod(jj,(nParParticles-mod(nParParticles,nProg))/nProg) == 0; 
    %    progress = 100*(jj-nfor(1))/nParParticles;
    %    fprintf(['w' num2str(np) ':' num2str(progress,'%.f') '% ']); 
    %end
    nBounceInd = 0; % number of bounces 
    
    
    if doODE45
        % Initial positions and velocities                                   
        x_init = [x0(jj);y0(jj);z0(jj);vx0(jj);vy0(jj);vz0(jj)]*1e3; % m, m/s

        % Integrate trajectory
        % Matlab ode solver :             
        % If the integration terminates beforehand, due to the passing
        % of the particle outside of the box, the integration stops, this is
        % defined in options and ExB.events.
        eh_x0 = 0;
        eh_y0 = 0;
        EoM = @(ttt,xxx) ExB.EquationOfMotion(ttt,xxx,B0,lr,lz,phi0,L,rlim,zlim,eh_x0,eh_y0);
        [t,x_sol] = ode45(EoM,[0 tint(jj)],x_init,options);
        x = x_sol(:,1);
        y = x_sol(:,2);
        z = x_sol(:,3);
        vx = x_sol(:,4);
        vy = x_sol(:,5);
        vz = x_sol(:,6);
    else % Runge-Kutta integration
        %h = 2.4e-5;
        %t = 0;
        % solve dx/dt = v;
        %x_init = [x0(jj);y0(jj);z0(jj)]*1e3; % m        
        %x  = x_init(1); y  = x_init(2); z  = x_init(3);         
        %nstep = 0;
        %while ~any([z(end)<-L*1e3/2 z(end)>L*1e3/2])
        %    ns = ns + 1;
        %    xk1 = vx(ns); yk1 = vy(ns); zk1 = vz(ns); 
        %    k2 = 
        
        % solve dv/dt = -(e/me)*(E+vxB)
        %v_init = [vx0(jj);vy0(jj);vz0(jj)]*1e3; %  m/s
        %vx = v_init(1); vy = v_init(2); vz = v_init(3);
        
        %end
    end
    % Save end positions
    partx1(nPart) = x(end);
    party1(nPart) = y(end);
    partz1(nPart) = z(end);
    % Save end time for diagnostics
    partend(nPart) = t(end);
    
    % See if parallel velocity changes sign = bounce
    parnBounce(nPart) = numel(find(diff(sign(vz))~=0)); 

    % Get azimuthal velocity and radial position
    th = atan2d((y-eh_y0),(x-eh_x0));
    %xg0 = find(x>=0); % right half xy-plane
    %xs0 = find(x<0); % left half xy-plane
    %th = zeros(size(x));
    %th(xg0) = atand(y(xg0)./x(xg0));
    %th(xs0) = atand(y(xs0)./x(xs0)) + 180;        
    vaztmp = -vx.*sind(th) + vy.*cosd(th); % azimuthal velocity 
    %vaztmp = -vy.*sind(th) + vx.*cosd(th); % azimuthal velocity 
    r = sqrt(x.^2 + y.^2);              

    if rzbinning % rz binning
        % Put values into discrete bins
        [nzz,binz] = histc(z,zGrid*1e3);
        [nrr,binr] = histc(r,rGrid*1e3);
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
            ir = max(binr(ic)); if ir == 0; continue; end; if ir > nr; continue; end
            iz = max(binz(ic)); if iz == 0; continue; end; if iz > nz; continue; end
            
            
            switch massgen
                case 'pre'; mass = m(jj);
                case 'hm';  mass = 2*pi*rlim*1e3*exp(rSurf(ir).^2/rlim.^2/2);
                case 'uni'; mass = 1;
                case 'bin'; mass = rSurf(ir)/rlim;
            end 
            if doMax
                parsumMVrzMax(ir,iz) = parsumMVrzMax(ir,iz) + cn.max(vaztmp(ic)).*m(jj);
                parsumMrz(ir,iz) = parsumMrz(ir,iz) + m(jj);            
            end
            if doMean            
                parsumMVrz(ir,iz) = parsumMVrz(ir,iz) + mean(vaztmp(ic)).*m(jj);
                parsumMrz(ir,iz) = parsumMrz(ir,iz) + m(jj);            
            end
        end   
    end 

    if noForBin
        vbin = ExB.bindata3([vaztmp vx vy vz ones(size(vaztmp))],x,y,z,xGrid*1e3,yGrid*1e3,zGrid*1e3);
        vazb = vbin{1};
        vxb  = vbin{2};
        vyb  = vbin{3};
        vzb  = vbin{4};
        mb   = vbin{5};
        parsumMVxyz  = parsumMVxyz  + vazb.*m(jj);
        parsumMxyz   = parsumMxyz   + mb*m(jj);                            
        parsumMVXxyz = parsumMVXxyz + vxb.*m(jj);
        parsumMVYxyz = parsumMVYxyz + vyb.*m(jj);
        parsumMVZxyz = parsumMVZxyz + vzb.*m(jj);
    else
    if xyzbinning % xyz binning            
        [nxx,binx] = histc(x,xGrid*1e3);
        [nyy,biny] = histc(y,yGrid*1e3);
        [nzz,binz] = histc(z,zGrid*1e3);
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
            ix = max(binx(ic)); if ix == 0; continue; end; if ix > nx; continue; end
            iy = max(biny(ic)); if iy == 0; continue; end; if iy > ny; continue; end
            iz = max(binz(ic)); if iz == 0; continue; end; if iz > nz; continue; end
                                   
            switch massgen
                case 'pre'; mass = m(jj);
                case 'hm';  mass = 2*pi*rlim*1e3*exp((xSurf(ix).^2+ySurf(iy).^2)/rlim.^2/2);
                case 'uni'; mass = 1;
                case 'bin'; mass = sqrt(xSurf(ix).^2+ySurf(iy).^2)/rlim;
            end
            
            parsumMVxyz(ix,iy,iz)  = parsumMVxyz(ix,iy,iz)  + mean(vaztmp(ic)).*m(jj);
            parsumMxyz(ix,iy,iz)   = parsumMxyz(ix,iy,iz)   + m(jj);                            
            parsumMVXxyz(ix,iy,iz) = parsumMVXxyz(ix,iy,iz) + mean(vx(ic)).*m(jj);
            parsumMVYxyz(ix,iy,iz) = parsumMVYxyz(ix,iy,iz) + mean(vy(ic)).*m(jj);
            parsumMVZxyz(ix,iy,iz) = parsumMVZxyz(ix,iy,iz) + mean(vz(ic)).*m(jj);                                     
            
            if doMax
                parsumMVxyzMax(ix,iy,iz)  = parsumMVxyzMax(ix,iy,iz)  + cn.max(vaztmp(ic)).*m(jj);
                parsumMxyzMax(ix,iy,iz)   = parsumMxyzMax(ix,iy,iz)   + m(jj);                            
                parsumMVXxyzMax(ix,iy,iz) = parsumMVXxyzMax(ix,iy,iz) + cn.max(vx(ic)).*m(jj);
                parsumMVYxyzMax(ix,iy,iz) = parsumMVYxyzMax(ix,iy,iz) + cn.max(vy(ic)).*m(jj);
                parsumMVZxyzMax(ix,iy,iz) = parsumMVZxyzMax(ix,iy,iz) + cn.max(vz(ic)).*m(jj);            
            end
        end                    
    end  
    end

    % status update
%     if mod(jj,round(nParticles/nwait)) == 0; 
%         n = numel(msgProgress);
%         msgReverse = repmat(sprintf('\b'), 1, n);
%         msgProgress = ['Progress: ' num2str(jj) '/' num2str(nParticles)];         
%         fprintf('%s',msgReverse)
%         fprintf('%s',msgProgress); 
%     end 

end

% Add up the results for each parfor-loop
sumMVrz   = sumMVrz   + parsumMVrz;
sumMrz    = sumMrz    + parsumMrz;
sumMVxyz  = sumMVxyz  + parsumMVxyz;
sumMxyz   = sumMxyz   + parsumMxyz;
sumMVXxyz = sumMVXxyz + parsumMVXxyz;
sumMVYxyz = sumMVYxyz + parsumMVYxyz;
sumMVZxyz = sumMVZxyz + parsumMVZxyz;

if doMax
sumMVrzMax   = sumMVrzMax   + parsumMVrzMax;
sumMrzMax    = sumMrzMax    + parsumMrzMax;
sumMVxyzMax  = sumMVxyzMax  + parsumMVxyzMax;
sumMxyzMax   = sumMxyzMax   + parsumMxyzMax;
sumMVXxyzMax = sumMVXxyzMax + parsumMVXxyzMax;
sumMVYxyzMax = sumMVYxyzMax + parsumMVYxyzMax;
sumMVZxyzMax = sumMVZxyzMax + parsumMVZxyzMax;
end

cellx1{np} = partx1; 
celly1{np} = party1;
cellz1{np} = partz1;
celltend{np} = partend;
cellnBounce{np} = parnBounce;

end % end of parfor loop
warning on
% Collect cells into arrays
tend = []; nBounce = []; x1 = []; y1 = []; z1 = [];
for uu = 1:npl    
    tend = [tend; celltend{uu}];
    nBounce = [nBounce; cellnBounce{uu}];
    x1 = [x1; cellx1{uu}];
    y1 = [y1; celly1{uu}];
    z1 = [z1; cellz1{uu}];
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
    
    if doMax
        meanVxyzMax = sumMVxyzMax./sumMxyzMax;
        meanVXxyzMax = sumMVXxyzMax./sumMxyzMax;
        meanVYxyzMax = sumMVYxyzMax./sumMxyzMax;
        meanVZxyzMax = sumMVZxyzMax./sumMxyzMax;
        vazmat2Max = meanVxyzMax;
        numbxyzMax = sumMxyzMax;
    end
    
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
    if mod(ny+1,2)==0 % uneven number of ybins
        nyc = ceil(ny/2);
        yind = (nyc-3):(nyc+4);
    else % even number of bins
        nyc = ceil(ny/2);
        yind = (nyc-2):(nyc+4);
    end
    if mod(nz+1,2)==0 % uneven number of zbins
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

%% Calculate B for both model and simulation
try if doB; 
        % Subtract noise level from simulated ExB-drift in xyz-coordinates
        % Take noise to be the mean between 10 and 20 and 80 and 90%.
        limits = [0.02 0.25 0.75 0.98];
        zlimits = -zlim+zlim*limits*2;
        int1 = find(zSurf>(zlimits(1)) & zSurf<(zlimits(2)));
        int2 = find(zSurf>(zlimits(3)) & zSurf<(zlimits(4)));
        VXxyz = sumMVXxyz./sumMxyz;
        VYxyz = sumMVYxyz./sumMxyz;
        VZxyz = sumMVZxyz./sumMxyz;
        VXnoise = repmat(mean(VXxyz(:,:,[int1, int2]),3),1,1,nz);
        VYnoise = repmat(mean(VYxyz(:,:,[int1, int2]),3),1,1,nz);
        VZnoise = repmat(mean(VZxyz(:,:,[int1, int2]),3),1,1,nz);
        VXfix = VXxyz-VXnoise;
        VYfix = VYxyz-VYnoise;
        VZfix = VZxyz-VZnoise;
        rbin=35;
        thbin=24;
        zbin=35;
        %rbin=17; thbin=15; zbin=15;
        for sim = [0 1 2]; ExB.calculateB_rz; end
    end
catch meB; disp('Could not calculate B. Bummer.')
end

%% Save 
if doSave
    saveIn = {'/home/cecilia/Research/ElectronHoleSimulation/SavedData/',... % on spis
              '/Users/Cecilia/Research/EH/TestParticleSimulation/'}; % on local
    for ss = 1:2
        try
            saveFilepath = [saveIn{ss} datestr(now,'yyyymmddTHHMMSS') '-' name '-' num2str(nParticles)];
            save(saveFilepath);
            disp(['Saved ' saveFilepath])
        catch meSave
            disp(['Error. Could not save ' saveFilepath])
            continue; 
        end 
    end
end
% if doSave; ExB.save; end

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
if 0
    try % only shut down pool if working on spis
        [~,hostname]=system('hostname');
        if ~strcmp(hostname(1:4),'moth')        
            delete(gcp)
        end
    catch
    end
end
end
clear doPool noGyrocenter