
% Define and choose model
model = 2;
switch model
    case 1 % 2007-08-31
        lr = 15;
        lz = 10;
        phi0 = 200;
        B0 = 25;
        Tpar = 1600;
        Tper = 1600;
        r0 = 15;
        a1 = 45;
        a2 = 0; 
        veh = 500;
    case 2 % Tao2011
        lr = 35;
        lz = 35;
        phi0 = 3*3000;
        B0 = 50;
        Tpar = 8000;
        Tper = 5000;
        r0 = lr;
        a1 = 45;
        a2 = 0; 
        veh = 100*1e6;
end  

% Set up grid
nz = 100; % number of z bins
nr = 100; % number of r bins
rlim = 4*lr; % km 
zlim = 4*lz; % km
rGrid = linspace(0,rlim,nr+1); % centers of bins
zGrid = linspace(-zlim,zlim,nz+1); % centers of bins
dr = diff(rGrid(1:2)); % bin size
dz = diff(zGrid(1:2)); % bin size
rSurf = rGrid(1:end-1)+dr/2; % edges between bins
zSurf = zGrid(1:end-1)+dz/2; % edges between bins

% Initialize variables
nParticles = 1;
rIndex = zeros(nParticles,1);
zIndex = zeros(nParticles,1);
distribution = zeros(nr,nz);
vaz = cell(nr,nz);
numb = zeros(nr,nz);
x0s = NaN(nParticles,1);
y0s = NaN(nParticles,1);
vx0s = NaN(nParticles,1);
vy0s = NaN(nParticles,1);

    vtz = cn_eV2v(Tpar,'ev');    % eV -> km/s (parallel thermal velocity)
    vtperp = cn_eV2v(Tper,'ev');   % eV -> km/s (perpendicular thermal velocity)
    veh = 0; % km/s, electron hole velocity    
    
    vtperp=vtperp.*randn(nParticles,1);%*ones(nParticles,1);
    
nTot = 0; % We throw away a number of particles to not overpopulate the outer rbin
wait = waitbar(0,'Calculation trajectories, please wait...');

for ii=1:nParticles     
    %tic;
    try
        if mod(ii,100) == 0
            waitbar(ii/nParticles,wait)
        end
    % Initialize starting position
    x0 = rand*(rlim+dr); 
    y0 = rand*(rlim+dr);
    z0 = rand*zlim*sign(rand-0.5);        
    
    %if sqrt(x0^2+y0^2) > rlim; continue; end           
    nTot = nTot + 1;
    
    
    % Calculate trajectory
    %tic
    [x y z vx vy vz] = art2.ExB(lr,lz,phi0,B0,x0,y0,vtperp(ii),vtz-veh);    
    %toc
    % Get azimuthal velocity and radial position
    th = atan(y./x);
    vaztmp = -vx.*sin(th) + vy.*cos(th);
    r = sqrt(x.^2+y.^2);
    str1 = [num2str(max(r)) '    ' num2str(max(rlim*1e3))];    
    %if max(r) > rlim*1e3; continue; end % second security
        
    % Put values into discrete bins    
    [nzz,binz]=histc(z,zSurf*1e3);
    [nrr,binr]=histc(r,rSurf*1e3);
        %%
    [a,b,c]=unique([binz binr],'rows');
    if 1
    remz = find(a(:,1)==0);
    remr = find(a(:,2)==0);
    rem = unique([remz;remr]);
    %%
    a(rem,:)= [];
    b(rem)= [];
    end
    %%
    % a is the pair of binz which are unique
    % b is the stopping index for each unique pair.
    % c is an array with indexing indicating the indices for each pair. 
    % example: 
    % a(1:3,:) = [6 30;7 30;8 30]
    % b(1:3) = [42;207;303]  
    % c(1:42) are all ones
    % c(43:207) are all twos
    % c(208:303) are all threes
    
    % loop over unique pairs
    nun = numel(b);
    b = [0;b]; % add a zero for first loop
    for ui = 1:nun        
        vaz{a(ui,1),a(ui,2)} = [vaz{a(ui,1),a(ui,2)}  mean(vaztmp(b(ui)+1:b(ui+1)))];
        %numb(a(ui,1),a(ui,2)) = numb(a(ui,1),a(ui,2)) + 1; % one more particles passed by this bin
    end
    %disp(['Particle #' num2str(ii) '. ' str1 '. Time: ' num2str(toc) '.'])
    % Save starrting positions for diagnostics
    x0s(ii) = x(1);
    y0s(ii) = y(1);
    vx0s(ii) = vx(1);
    vy0s(ii) = vy(1);
    catch me
        disp(['Error. Skipping particle #' num2str(ii) '.'])        
        continue;
    end    
end
close(wait)

% average azimuthal velocities
vazmat = zeros(size(vaz));
for iz = 1:nz % loop over z bins    
        for ir = 1:nr % loop over r bins       
            vazmat(ir,iz) = nanmean(vaz{ir,iz});
            numb(ir,iz) = numel(vaz{ir,iz});
        end
end
%%

if 1 % plot
    h(1) = subplot(2,1,1);
    hca = h(1);
    surf(hca,rGrid',zGrid,zeros(nr+1,nz+1),vazmat*1e-3)
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    cha1 = colorbar;
    grid(hca,'off')
    view([0 0 1])
    %strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,  r_0 = ' num2str(r0*1e-3,'%.f') ' km'],...
    %        ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtz,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtperp,'%.f') ' km/s,   n_{particles} = ' num2str(nTot)  '/' num2str(nParticles)],...
    %        };
    %title(hca,strTitle)
    
    h(2) = subplot(2,1,2);
    hca = h(2);
    surf(hca,rGrid',zGrid,zeros(nr+1,nz+1),numb)
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')    
    cha2 = colorbar;
    grid(hca,'off')
    view([0 0 1])
end
%%
%surf(rGrid',zGrid,zeros(numel(zGrid),numel(rGrid)),distribution)
% surf(rSurf',zSurf,distribution*0,distribution)
%xlabel('r')
%ylabel('z')








%    [h,x,y,z,vx,vy,vz,ax,ay,az]=art2.ExB(lr,lz,phi0,B0,Tpar,Tper,r0,a1,a2);
    