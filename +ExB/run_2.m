% Define and choose model
modelnr = 2;
switch modelnr
    case 1 % 2007-08-31        
        name = '2007-08-31';
        lr   = 15;
        lz   = 10;
        phi0 = 200;
        B0   = 25;
        Tpar = 1600;
        Tper = 1600;        
        veh  = 500;
    case 2 % Tao2011
        name = 'Tao2011';
        lr   = 35;
        lz   = 35;
        phi0 = 3000;
        B0   = 50;
        Tpar = 8000;
        Tper = 5000;        
        veh  = 100000;
    case 3 % No electric field
        name = 'Basic';
        lr   = 20;
        lz   = 20;
        phi0 = 1;
        B0   = 20;
        Tpar = 1000;
        Tper = 1000;        
        veh  = 0;
        % Tper = 1000; B0 = 10 gives rho_e ~ 10.        
    case 4 % Exaggerated Tao
        name = 'ExagTao2011';
        lr   = 35;
        lz   = 35;
        phi0 = 10000;
        B0   = 50;
        Tpar = 8000;
        Tper = 5000;        
        veh  = 0;
end  

% Define time to run each particle
L = 6*lz;  % length of box
T = veh/L; % time for eh to pass box !!! This is waaay too long.
T = 0.01; % this gives 97975 timesteps

% Set up simulation timestepping
units=irf_units;
fce=units.e*B0*1e-9/units.me/2/pi; % gyro freq. Hz
tce=1/fce; % gyro period. s
nce = T/tce; % number of gyroperiods
ntpertce = 7000; % number of timesteps per gyroperiod
nt = round(ntpertce*nce); % number of timesteps
t=linspace(0,T,nt)'; % time vector
dt=diff(t(1:2)); % time step

%% Set up grid.
% number of bins
nx = 200;
ny = 200;
nz = 100; 
nr = 200; 
% edges of box [km]
zlim = L/2; 
rlim = 3*lr; 
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
nParticles = input('nParticles = '); % number of particles
%nParticles = str2num(nParticlesStr);
nOver = zeros(nParticles,1); % number of overshoots the electron makes
nBounce = zeros(nParticles,1); % number of overshoots the electron makes

vaz1 = cell(nr,nz); % azimuthal velocity
vaz2 = cell(nx,ny,nz); % azimuthal velocity
numbrz = zeros(nr,nz); % number of particles per bin
numbxyz = zeros(nx,ny,nz); % number of particles per bin

% Initialize starting position
x0 = rand(nParticles,1)*rlim.*sign(rand(nParticles,1)-0.5); 
y0 = rand(nParticles,1)*rlim.*sign(rand(nParticles,1)-0.5); 
z0 = rand(nParticles,1)*zlim.*sign(rand(nParticles,1)-0.5); 

% Initialize starting velocity
vtpar = cn_eV2v(Tpar,'ev'); % eV -> km/s (parallel thermal velocity)
vtper = cn_eV2v(Tper,'ev'); % eV -> km/s (perpendicular thermal velocity)
%vx0 = ExB.randv(nParticles,vtper).*sign(rand(nParticles,1)-0.5); 
%vy0 = ExB.randv(nParticles,vtper).*sign(rand(nParticles,1)-0.5); 
%vz0 = ExB.randv(nParticles,vtpar).*sign(rand(nParticles,1)-0.5); 
vx0 = vtper*randn(nParticles,1);
vy0 = vtper*randn(nParticles,1);
vz0 = vtpar*randn(nParticles,1)-veh;

% Initializing variable arrays
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

% hist(vx0,100) % plot velocity distributions
    
%% Run through particles
nTot = 0; % We throw away a number of particles to not overpopulate the outer rbin
wait = waitbar(0,'Calculating trajectories, please wait...');
nwait = 10;
for jj=1:nParticles     
%    try
        if mod(jj,round(nParticles/nwait)) == 0; waitbar(jj/nParticles,wait); end % waittick each 100 particle            

        % Initial positions and velocities
        % (input is in km and km/s)
        x(1) = x0(jj)*1e3;   % m
        y(1) = y0(jj)*1e3;   % m 
        z(1) = z0(jj)*1e3;   % m
        vx(1) = vx0(jj)*1e3; % m/s
        vy(1) = vy0(jj)*1e3; % m/s
        vz(1) = vz0(jj)*1e3; % m/s                
        
        nOverInd = 0;
        nBounceInd = 0;
        % Integrate trajectory
        for ii = 1:(nt-1)  
            vx(ii+1)=vx(ii)+dt*(-units.e/units.me)*((x(ii)/(lr*1e3)^2).*phi0.*exp(-0.5*(x(ii)/(lr*1e3)).^2-0.5*(y(ii)/(lr*1e3)).^2-0.5*(z(ii)/(lz*1e3)).^2) + vy(ii)*B0*1e-9);
            vy(ii+1)=vy(ii)+dt*(-units.e/units.me)*((y(ii)/(lr*1e3)^2).*phi0.*exp(-0.5*(x(ii)/(lr*1e3)).^2-0.5*(y(ii)/(lr*1e3)).^2-0.5*(z(ii)/(lz*1e3)).^2) - vx(ii)*B0*1e-9);
            vz(ii+1)=vz(ii)+dt*(-units.e/units.me)*((z(ii)/(lz*1e3)^2).*phi0.*exp(-0.5*(x(ii)/(lr*1e3)).^2-0.5*(y(ii)/(lr*1e3)).^2-0.5*(z(ii)/(lz*1e3)).^2));

            x(ii+1)=x(ii)+vx(ii)*dt;    
            y(ii+1)=y(ii)+vy(ii)*dt;   
            z(ii+1)=z(ii)+vz(ii)*dt;

            % Move to other side of box if electron overshoots
            if z(ii+1)>L*1e3/2 
                z(ii+1) = z(ii+1)-L*1e3;
                nOverInd = nOverInd  + 1;
            elseif z(ii+1)<-L*1e3/2 
                z(ii+1) = z(ii+1)+L*1e3;
                nOverInd = nOverInd  + 1;
            end
            % Count the number of bounces (= vz reversals) the electron does
            if sign(vz(ii+1)) ~= sign(vz(ii)) 
                nBounceInd = nBounceInd + 1;
            end
        end
        
        % Save bounce and overshoot statistics
        nOver(jj) = nOverInd;
        nBounce(jj) = nBounceInd;

        % Get azimuthal velocity and radial position
        xg0 = find(x>=0); % top half xy-plane
        xs0 = find(x<0); % bottom half xy-plane
        th = zeros(nt,1);
        th(xg0) = atand(y(xg0)./x(xg0));
        th(xs0) = - atand(y(xs0)./x(xs0));
        %th = atand(y./x); % angle of particle position
        vaztmp = -vx.*sind(th) + vy.*cosd(th); % azimuthal velocity
        r = sqrt(x.^2+y.^2); % radius of particle position        

        if 1 % rz binning
            % Put values into discrete bins
            [nzz,binz] = histc(z,zSurf*1e3);
            [nrr,binr] = histc(r,rSurf*1e3);

            [arz,brz,crz] = unique([binr binz],'rows');

            % Remove indices when the particle is outside the box
            % a is the pair of binz which are unique
            % b is the stopping index for each unique pair.
            % c is an array with indexing indicating the indices for each pair. 
            % example: 
            % a(1:3,:) = [6 30;7 30;8 30]
            % b(1:3) = [42;207;303]  
            % c(1:42) are all ones
            % c(43:207) are all twos
            % c(208:303) are all threes             

            remr = find(arz(:,1)==0);
            remz = find(arz(:,2)==0);
            remrz = unique([remr;remz]);
            arz(remrz,:) = [];
            brz(remrz) = [];                        

            % Loop over unique pairs and save into bins         
    %         nun = numel(brz);
    %         brz = [0;brz]; % add a zero for first loop
    %         for ui = 1:nun        
    %             vaz1{arz(ui,1),arz(ui,2)} = [vaz1{arz(ui,1),arz(ui,2)}  mean(vaztmp(brz(ui)+1:brz(ui+1)))];            
    %         end

            nun = max(crz);
            %brz = [0;brz]; % add a zero for first loop
            for ui = 1:nun                  
                ic = find(crz==ui);
                ir = max(binr(ic)); if ir == 0; continue; end
                iz = max(binz(ic)); if iz == 0; continue; end
                vaz1{ir,iz} = [vaz1{ir,iz}  mean(vaztmp(ic))];            
            end   
        end 
        
        if 1 % xyz binning            
            [nxx,binx] = histc(x,xSurf*1e3);
            [nyy,biny] = histc(y,ySurf*1e3);
            [nzz,binz] = histc(z,zSurf*1e3);
            [axyz,bxyz,cxyz] = unique([binx biny binz],'rows');
            remxyzx = find(axyz(:,1)==0); % find indices
            remxyzy = find(axyz(:,2)==0);
            remxyzz = find(axyz(:,3)==0);
            remxyz = unique([remxyzx;remxyzy;remxyzz]); % find unique pairs
            axyz(remxyz,:) = []; % remove them
            bxyz(remxyz) = [];

    %         if 0
    %         nun = numel(bxyz);
    %         bxyz = [0;bxyz]; % add a zero for first loop
    %         for ui = 1:nun        
    %             vaz2{axyz(ui,1),axyz(ui,2),axyz(ui,3)} = [vaz2{axyz(ui,1),axyz(ui,2),axyz(ui,3)}  mean(vaztmp(bxyz(ui)+1:bxyz(ui+1)))];            
    %         end
    %         end           
            nun = max(cxyz);
            %brz = [0;brz]; % add a zero for first loop
            for ui = 1:nun      
                %ui
                ic = find(cxyz==ui);
                ix = max(binx(ic)); if ix == 0; continue; end
                iy = max(biny(ic)); if iy == 0; continue; end
                iz = max(binz(ic)); if iz == 0; continue; end
                vaz2{ix,iy,iz} = [vaz2{ix,iy,iz}  mean(vaztmp(ic))];            
            end                    
        end  
              
%    catch me
%       disp(['Error. Skipping particle #' num2str(jj) '.'])        
%       continue;
%    end    
end
close(wait)
tTot = toc;
disp(['elapsed time = ' num2str(tTot) ' s'])
disp(['time per particle = ' num2str(tTot/nParticles) ' s'])
%% Make average azimuthal velocity in rz binning system
vazmat = cellfun(@nanmean,vaz1);
numbrz = cellfun(@numel,vaz1);
%% Make average azimuthal velocity in xyz binning system
if 0 % average azimuthal velocities (squeezed over vy)
vazmat2 = cellfun(@nanmean,vaz2);
numbxyz = cellfun(@numel,vaz2);
vazmat2yz = squeeze(nanmean(vazmat2,1));
vazmat2xz = squeeze(nanmean(vazmat2,2));
vazmat2xy = squeeze(nanmean(vazmat2,3));
numbyz = squeeze(nanmean(numbxyzb,1));
numbxz = squeeze(nanmean(numbxyzb,2));
numbxy = squeeze(nanmean(numbxyzb,3));
end

toplot = 2;
switch toplot
    case 1
%% Test plot what to expect, debug
if 1
for k = 1:6
    h(k) = subplot(3,2,k);
end
isub = 1;
irf_colormap('poynting');
if 1
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,vazmat'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    cha1 = colorbar('peer',hca);
    caxis(hca,20*[-1 1]*1e3)    
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])   
end

[RS,ZS] = meshgrid(rSurf,zSurf);
Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Etot = sqrt(Er.^2+Ez.^2);

if 1 % Illustrate the amplitude of the radial electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Er*1e3) % mV/m    
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha1 = colorbar('peer',hca); ylabel(cha1,'E_r [mV/m]')
    %caxis(hca,20*[-1 1]*1e3)    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Illustrate the amplitude of the parallel electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Ez*1e3)    
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
    %caxis(hca,20*[-1 1]*1e3)    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Illustrate the amplitude of the parallel electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Etot*1e3)    
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
    %caxis(hca,20*[-1 1]*1e3)    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Illustrate the expected ExB drift       
    ExBdrift = - Er / (B0*1e-9); % m/s 
    ExBdrift = ExBdrift*1e-3; % km/s
    
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);
    surf(hca,rGrid,zGrid,surfa,ExBdrift)
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha3 = colorbar('peer',hca); ylabel(cha3,'km/s')
    %caxis(hca,20*[-1 1]*1e3)
    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 0 % Illustrate the grid           
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,RS,ZS,RS)
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha4 = colorbar('peer',hca); ylabel(cha4,'r [km]')    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 0 % Illustrate the grid           
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);
    surf(hca,RS,ZS,ZS)
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha5 = colorbar('peer',hca); ylabel(cha5,'z [km]')       
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
 

strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
            ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
            };
title(h(1),strTitle') 
end

    case 2
%% Plot results
if 1
for k = 1:9
    h(k) = subplot(3,3,k);
end
isub = 1;
irf_colormap('poynting');

[RS,ZS] = meshgrid(rSurf,zSurf);
Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
%Etot = sqrt(Er.^2+Ez.^2);
ExBdrift = - Er / (B0*1e-9); % m/s 
ExBdrift = ExBdrift*1e-3; % km/s
    
if 1 % Illustrate the amplitude of the radial electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Er*1e3) % mV/m  
    Elim = max(max(Er*1e3));
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha1 = colorbar('peer',hca); ylabel(cha1,'E_r [mV/m]')
    caxis(hca,Elim*[-1 1])    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Illustrate the amplitude of the parallel electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Ez*1e3)    
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
    %caxis(hca,20*[-1 1]*1e3)    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Illustrate the expected ExB drift           
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);
    surf(hca,rGrid,zGrid,surfa,ExBdrift)
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')
    vlim = max(max(abs(ExBdrift)));
    caxis(hca,vlim*[-1 1])
    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
    title(hca,'ExB-drift')
end
if 1 % Un normalized number of passes per bin
    hca = h(isub); isub = isub + 1;
    surf(hca,rGrid',zGrid,surfa,numbrz')
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')    
    cha2 = colorbar('peer',hca);
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
end
if 1 % Normalized number of passes per bin
    hca = h(isub); isub = isub + 1;
    normnumbrz = repmat(rSurf,numel(zSurf),1);
    surf(hca,rGrid',zGrid,surfa,numbrz'./normnumbrz)
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')    
    cha2 = colorbar('peer',hca);
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,vazmat'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    cha1 = colorbar('peer',hca);
    vlim = max(max(abs(ExBdrift)));
    %caxis(hca,20*[-1 1]*1e3)
    caxis(hca,vlim*[-1 1])    
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Starting positions
    hca = h(isub); isub = isub + 1;    
    plot3(hca,x0,y0,z0,'.'); hold(hca,'on');
    quiver3(hca,x0,y0,z0,vx0,vy0,vz0); hold(hca,'off');
    title(hca,'Starting positions')
    view(hca,[0 0 1])
    xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
    set(hca,'xlim',[-rlim rlim],'ylim',[-rlim rlim])
end  
if 1 % 1D starting velocities, vx vy vz
    hca = h(isub); isub = isub + 1;
    hist(hca,[vx0*1e-3,vy0*1e-3,vz0*1e-3],30);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v [10^3 km/s]'); ylabel(hca,'# of particles');
    legend(hca,'v_x','v_y','v_z')
end  
if 0 % v_perp
    hca = h(isub); isub = isub + 1;
    hist(hca,sqrt(vx0.^2+vy0.^2)*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Perpendicular starting velocities')
    xlabel(hca,'v_{\perp} [10^3 km/s]'); ylabel(hca,'# of particles');
end 
if 0 % v_tot
    hca = h(isub); isub = isub + 1;
    hist(hca,sqrt(vx0.^2+vy0.^2+vz0.^2)*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Total starting velocities')
    xlabel(hca,'v_{tot} [10^3 km/s]'); ylabel(hca,'# of particles');
end 
if 0 % v_x
    hca = h(isub); isub = isub + 1;
    hist(hca,vx0*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_x [10^3 km/s]'); ylabel(hca,'# of particles');
end    
if 0 % v_y
    hca = h(isub); isub = isub + 1;
    hist(hca,vy0*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_y [10^3 km/s]'); ylabel(hca,'# of particles');
end    
if 0 % v_z
    hca = h(isub); isub = isub + 1;
    hist(hca,vz0*1e-3,100); hold(hca,'on');
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_z [10^3 km/s]'); ylabel(hca,'# of particles');
    f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
    f = @(v,vt) (1/vt).*exp(-(v/vt).^2);
    %vt = 50e3; % km/s
    [bincounts,binpositions] = hist(hca,vz0,100); hold(hca,'on');
    binwidth = binpositions(2) - binpositions(1);
    histarea = binwidth*sum(bincounts);
    v = linspace(0,4*vtpar,100)*1e-3;
    plot(hca,v,f(v,vtpar*1e-3)*histarea*2e-5,'r'); hold(hca,'off');
    
    %plot(hca,v,nParticles/10*f(v,vt*1e-3));
     hold off
end 
if 0 % bounce statistics
    hca = h(isub); isub = isub + 1;
    maxBounce = max(nBounce);
    [nnn,~] = histc(nBounce,0:maxBounce);    
    bar(hca,0:maxBounce,nnn);
    fractionBounce = sum(nnn(2:end))/sum(nnn);
    bounceStr = ['n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')];    
    title(hca,['Bouncing particles: ' bounceStr ])
    xlabel(hca,'# of bounces'); ylabel(hca,'# of particles');    
end   
if 0 % overshoot statistics
    hca = h(isub); isub = isub + 1;    
    maxOver = max(nOver);
    [nnn,~] = histc(nOver,0:maxOver);    
    bar(hca,0:maxOver,nnn);    
    title(hca,'Overshoot')
    xlabel(hca,'# of overshoots'); ylabel(hca,'# of particles');
end   

strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
            ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
            };
title(h(2),strTitle')  
%cn.print('run100000_4');
end
    case 3
%% % Plot xyz binning
for k = 1:4
    h(k) = subplot(2,2,k);
end
isub = 1;
irf_colormap('poynting');

[RS,ZS] = meshgrid(rSurf,zSurf);
Er = RS./(lr.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(RS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
%Etot = sqrt(Er.^2+Ez.^2);
ExBdrift = - Er / (B0*1e-9); % m/s 
ExBdrift = ExBdrift*1e-3; % km/s
    
if 0 % Illustrate the amplitude of the radial electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Er*1e3) % mV/m  
    Elim = max(max(Er*1e3));
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha1 = colorbar('peer',hca); ylabel(cha1,'E_r [mV/m]')
    caxis(hca,Elim*[-1 1])    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 0 % Illustrate the amplitude of the parallel electric field               
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,Ez*1e3)    
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha2 = colorbar('peer',hca); ylabel(cha2,'E_z [mV/m]')
    %caxis(hca,20*[-1 1]*1e3)    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 0 % Illustrate the expected ExB drift           
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);
    surf(hca,rGrid,zGrid,surfa,ExBdrift)
    shading(hca,'flat')
    xlabel(hca,'r'); ylabel(hca,'z')
    cha = colorbar('peer',hca); ylabel(cha,'km/s')
    vlim = max(max(abs(ExBdrift)));
    caxis(hca,vlim*[-1 1])
    
    grid(hca,'off'); view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim]) 
    title(hca,'ExB-drift')
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nr+1);    
    surf(hca,rGrid,zGrid,surfa,vazmat'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')
    cha1 = colorbar('peer',hca);
    %caxis(hca,20*[-1 1]*1e3)
    %caxis(hca,vlim*[-1 1])   
    axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])    
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(ny+1,nx+1);    
    surf(hca,xGrid,yGrid,surfa,vazmat2xy'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'x')
    ylabel(hca,'y')
    cha1 = colorbar('peer',hca);
    %caxis(hca,20*[-1 1]*1e3)
    %caxis(hca,vlim*[-1 1])  
    axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',rlim*[-1 1])    
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,nx+1);    
    surf(hca,xGrid,zGrid,surfa,vazmat2xz'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'x')
    ylabel(hca,'z')
    cha1 = colorbar('peer',hca);
    %caxis(hca,20*[-1 1]*1e3)
    %caxis(hca,vlim*[-1 1]) 
    axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
end
if 1 % Plot average azimuthal velocity
    hca = h(isub); isub = isub + 1;        
    surfa = zeros(nz+1,ny+1);    
    surf(hca,yGrid,zGrid,surfa,vazmat2yz'*1e-3)    
    shading(hca,'flat')
    xlabel(hca,'y')
    ylabel(hca,'z')
    cha1 = colorbar('peer',hca);
    %caxis(hca,20*[-1 1]*1e3)
    %caxis(hca,vlim*[-1 1])   
    axis(hca,'equal')
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',rlim*[-1 1],'ylim',zlim*[-1 1])    
end

if 0 % Un normalized number of passes per bin
    hca = h(isub); isub = isub + 1;
    surf(hca,rGrid',zGrid,surfa,numbrz')
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')    
    cha2 = colorbar('peer',hca);
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
end
if 0 % Normalized number of passes per bin
    hca = h(isub); isub = isub + 1;
    normnumbrz = repmat(rSurf,numel(zSurf),1);
    surf(hca,rGrid',zGrid,surfa,numbrz'./normnumbrz)
    shading(hca,'flat')
    xlabel(hca,'r')
    ylabel(hca,'z')    
    cha2 = colorbar('peer',hca);
    grid(hca,'off')
    view(hca,[0 0 1])
    set(hca,'xlim',[0 rlim],'ylim',[-zlim zlim])
end
if 0 % Starting positions
    hca = h(isub); isub = isub + 1;    
    plot3(hca,x0,y0,z0,'.'); hold(hca,'on');
    quiver3(hca,x0,y0,z0,vx0,vy0,vz0); hold(hca,'off');
    title(hca,'Starting positions')
    view(hca,[0 0 1])
    xlabel(hca,'x [km]'); ylabel(hca,'y [km]'); zlabel(hca,'z [km]')
    set(hca,'xlim',[-rlim rlim],'ylim',[-rlim rlim])
end  
if 0 % 1D starting velocities, vx vy vz
    hca = h(isub); isub = isub + 1;
    hist(hca,[vx0*1e-3,vy0*1e-3,vz0*1e-3],30);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v [10^3 km/s]'); ylabel(hca,'# of particles');
    legend(hca,'v_x','v_y','v_z')
end  
if 0 % v_perp
    hca = h(isub); isub = isub + 1;
    hist(hca,sqrt(vx0.^2+vy0.^2)*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Perpendicular starting velocities')
    xlabel(hca,'v_{\perp} [10^3 km/s]'); ylabel(hca,'# of particles');
end 
if 0 % v_tot
    hca = h(isub); isub = isub + 1;
    hist(hca,sqrt(vx0.^2+vy0.^2+vz0.^2)*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Total starting velocities')
    xlabel(hca,'v_{tot} [10^3 km/s]'); ylabel(hca,'# of particles');
end 
if 0 % v_x
    hca = h(isub); isub = isub + 1;
    hist(hca,vx0*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_x [10^3 km/s]'); ylabel(hca,'# of particles');
end    
if 0 % v_y
    hca = h(isub); isub = isub + 1;
    hist(hca,vy0*1e-3,100);
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_y [10^3 km/s]'); ylabel(hca,'# of particles');
end    
if 0 % v_z
    hca = h(isub); isub = isub + 1;
    hist(hca,vz0*1e-3,100); hold(hca,'on');
    set(hca,'xtick',[-1.5 -1 -0.5 0 0.5 1 1.5]*1e2)
    title(hca,'Starting velocities')
    xlabel(hca,'v_z [10^3 km/s]'); ylabel(hca,'# of particles');
    f = @(v,vt) (v/vt).*exp(-(v/vt).^2);
    f = @(v,vt) (1/vt).*exp(-(v/vt).^2);
    %vt = 50e3; % km/s
    [bincounts,binpositions] = hist(hca,vz0,100); hold(hca,'on');
    binwidth = binpositions(2) - binpositions(1);
    histarea = binwidth*sum(bincounts);
    v = linspace(0,4*vtpar,100)*1e-3;
    plot(hca,v,f(v,vtpar*1e-3)*histarea*2e-5,'r'); hold(hca,'off');
    
    %plot(hca,v,nParticles/10*f(v,vt*1e-3));
     hold off
end 
if 0 % bounce statistics
    hca = h(isub); isub = isub + 1;
    maxBounce = max(nBounce);
    [nnn,~] = histc(nBounce,0:maxBounce);    
    bar(hca,0:maxBounce,nnn);
    fractionBounce = sum(nnn(2:end))/sum(nnn);
    bounceStr = ['n_b/n_{tot} = ' num2str(fractionBounce,'%.2f')];    
    title(hca,['Bouncing particles: ' bounceStr ])
    xlabel(hca,'# of bounces'); ylabel(hca,'# of particles');    
end   
if 0 % overshoot statistics
    hca = h(isub); isub = isub + 1;    
    maxOver = max(nOver);
    [nnn,~] = histc(nOver,0:maxOver);    
    bar(hca,0:maxOver,nnn);    
    title(hca,'Overshoot')
    xlabel(hca,'# of overshoots'); ylabel(hca,'# of particles');
end   

strTitle = {['B_0 = ' num2str(B0) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr) ' km,  l_z = ' num2str(lz) ' km,   n_{particles} = ' num2str(nParticles)],...
            ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtpar,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtper,'%.f') ' km/s,   v_{eh} = ' num2str(veh,'%.f') ' km/s'],...
            };
title(h(2),strTitle')  
end
%%
% save run
save_run.name = name;
save_run.nOver = nOver;
save_run.nBounce = nBounce;
save_run.x0 = x0;
save_run.y0 = y0;
save_run.z0 = z0;
save_run.vx0 = vx0;
save_run.vy0 = vy0;
save_run.vz0 = vz0;
save_run.vaz1 = vaz1;
save_run.vazmat = vazmat;
save_run.numbrz = numbrz;

if 0
try
save('run_2014-03-01_a','name','nOver','nBounce','x0','y0','z0','vx0','vy0','vz0','vaz1','vazmat','numbrz')
catch
end
try
save('run_2014-03-01_b','save_run')
catch
end
end
