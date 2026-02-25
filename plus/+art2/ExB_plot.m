function varargout = ExB_plot(lr,lz,phi0,B0,Tpar,Tper,r0,a1,a2) 
% Plots trajectory of electron as it goes through an electron hole.
% - vxB and E acceleration (only in perpendicular direction, vz is constant)
%
%   out = ExB_plot(lr,lz,phi0,B0,Tpar,Tper,r0,a1,a2);
%
%       lr - half perpendicular width [km]
%       lz - half parallel width [km]
%       phi0 - center potential [V]
%       B0 - background magnetic field [nT]
%       Tpar - parallel electron tempreature [eV]
%       Tper - perpendicular electron tempreature [eV]
%       r0 - initial radial position of particle from center of hole [km]
%       a1 - defines starting postition [degrees]
%       a2 - defines starting direction [degrees]
%
%   Examples:
%       [h,x,y,z,vx,vy,vz,ax,ay,az]=art2.ExB(20,35,3000,50,8000,5000,15,0,0); % Tao2011
%       [h,x,y,z,vx,vy,vz,ax,ay,az]=art2.ExB(15,10,300 ,25,2000,2000,15,0,0); % 2007-08-31

units=irf_units; % Physical units

% Constants and initial values
lr = lr*1e3; % km -> m
lz = lz*1e3; % km -> m
if exist('r0','var'); r0 = r0*1e3; % km -> m
else r0 = 2*lr*rand; end
phi0 = phi0; % V
B0 = B0*1e-9; % nT -> T
vtz = cn_eV2v(Tpar,'ev')*1e3;    % eV -> m/s
vtxy = cn_eV2v(Tper,'ev')*1e3;   % eV -> m/s
L = lz*7; % Length of box
T = L/vtz; % Time to pass the box
fce=units.e*B0/units.me/2/pi; % gyro freq. Hz
tce=1/fce; % gyro period. s
nce = T/tce; % number of gyroperiods
ntpertce = 7000; % number of timesteps per gyroperiod
nt = round(ntpertce*nce); % number of timesteps
t=linspace(0,T,nt)'; % time vector
dt=diff(t(1:2)); % time step

% Initializing variable vector
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

% Initial positions and velocities
if ~exist('a1','var'); a1 = rand*360; end % degrees
x(1) = cosd(a1)*r0;  % m
y(1) = sind(a1)*r0;  % m 
z(1) = -L/2;            % m
if ~exist('a2','var'); a2 = rand*360; end % degrees
vx(1) = cosd(a2)*vtxy;    % m/s
vy(1) = sind(a2)*vtxy;    % m/s
vz(1) = vtz;                 % m/s

% Model electric field
Ex0 = @(x,y,z,lr,lz,phi0)(x/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ey0 = @(x,y,z,lr,lz,phi0)(y/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ez0 = @(x,y,z,lr,lz,phi0)(z/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ex = @(x,y,z) Ex0(x,y,z,lr,lz,phi0); % V/m
Ey = @(x,y,z) Ey0(x,y,z,lr,lz,phi0); % V/m
Ez = @(x,y,z) Ez0(x,y,z,lr,lz,phi0); % V/m

% Calculating trajectory
wait = waitbar(0,'Calculation trajectory, please wait...');
for ii=1:(nt-1)
    %waitbar(ii/nt,wait)
    ax(ii) = (-units.e/units.me)*(Ex(x(ii),y(ii),z(ii)) + vy(ii)*B0);
    ay(ii) = (-units.e/units.me)*(Ey(x(ii),y(ii),z(ii)) - vx(ii)*B0);   
    az(ii) = (-units.e/units.me)*(Ez(x(ii),y(ii),z(ii)));   
    
    dvx(ii)=ax(ii)*dt;
    dvy(ii)=ay(ii)*dt;
    dvz(ii)=az(ii)*dt;
        
    vx(ii+1)=vx(ii)+dvx(ii);
    vy(ii+1)=vy(ii)+dvy(ii);
    vz(ii+1)=vz(ii)+dvz(ii);
    
    x(ii+1)=x(ii)+vx(ii)*dt;    
    y(ii+1)=y(ii)+vy(ii)*dt;   
    z(ii+1)=z(ii)+vz(ii)*dt;       
end
close(wait)

% Calculate angular deplacement of gyrocenter
gind = ntpertce; % number of indices for one gyroperiod
x1 = mean(x(1:gind));
x2 = mean(x((end-gind):end));
y1 = mean(y(1:gind));
y2 = mean(y((end-gind):end));

alfa1 = atand(y1/x1);
alfa2 = atand(y2/x2);
dalfa = alfa2-alfa1;

% Make figure
xylim1 = 1.1*max(abs([x;y]))*1e-3;
xylim2 = 1.1*(r0+vtxy/fce/2/pi)*1e-3;
xylim = max([xylim1 xylim2]);
%xylim = 2*lr*1e-3;
if xylim < lr*1e-3; xylim = lr*1e-3; end

plot3(x*1e-3,y*1e-3,z*1e-3); axl=gca;
xlabel(axl,'x [km]'); ylabel(axl,'y [km]'); zlabel(axl,'z [km]');
axis equal
axis square
set(axl,'xlim',xylim*[-1 1],'ylim',xylim*[-1 1])
view([0 0 1])
hold(axl,'on')

strTitle = {['B_0 = ' num2str(B0*1e9) ' nT,  \phi_0 = ' num2str(phi0) ' V, l_r = ' num2str(lr*1e-3) ' km,  l_z = ' num2str(lz*1e-3) ' km,  r_0 = ' num2str(r0*1e-3,'%.f') ' km'],...
            ['T_{||} = ' num2str(Tpar) ' eV,  T_{\perp} = ' num2str(Tper) ' eV,  v_{t||} = ' num2str(vtz*1e-3,'%.f') ' km/s,  v_{t\perp} = ' num2str(vtxy*1e-3,'%.f') ' km/s'],...
            ['T = ' num2str(T*1e3) ' ms,  t_{ce} = ' num2str(tce*1e3) ' ms,  n_{gyroperiod} = ' num2str(nce) '  ms,  d_{\alpha} = ' num2str(dalfa) ' deg']};
title(axl,strTitle)

% Draw ellipses to indicate eh
if 1
    ellipse(lr*1e-3,lr*1e-3,0,0,0,[1 0 0],100,'xy');
    ellipse(lr*1e-3,lz*1e-3,0,0,0,[1 0 0],100,'xz');
    ellipse(lz*1e-3,lr*1e-3,0,0,0,[1 0 0],100,'yz');
end

% Plot electric field vectors
if 1
    nEx = 5;
    nEy = 5;
    nEz = 3;
    xE = linspace(-xylim,xylim,nEx);
    yE = linspace(-xylim,xylim,nEy);
    zE = linspace(-2*lz*1e-3,2*lz*1e-3,nEz);
    [XE,YE,ZE] = meshgrid(xE,yE,zE);
    qXE = zeros(nEx,nEy,nEz);
    qYE = zeros(nEx,nEy,nEz);
    qZE = zeros(nEx,nEy,nEz);

    for ii=1:nEx
        for jj=1:nEy
            for kk=1:nEz
                qXE(ii,jj,kk) = Ex(XE(ii,jj,kk)*1e3,YE(ii,jj,kk)*1e3,ZE(ii,jj,kk)*1e3)*1e3; % mV/m
                qYE(ii,jj,kk) = Ey(XE(ii,jj,kk)*1e3,YE(ii,jj,kk)*1e3,ZE(ii,jj,kk)*1e3)*1e3; % mV/m
                qZE(ii,jj,kk) = Ez(XE(ii,jj,kk)*1e3,YE(ii,jj,kk)*1e3,ZE(ii,jj,kk)*1e3)*1e3; % mV/m
            end
        end
    end
    quiver3(axl,XE,YE,ZE,qXE,qYE,qZE,0.5)
end
hold(axl,'off')


%%
% Prepare output
varargout{1} = axl;
varargout{2} = x;
varargout{3} = y;
varargout{4} = z;
varargout{5} = vx;
varargout{6} = vy;
varargout{7} = vz;
varargout{8} = ax;
varargout{9} = ay;
varargout{10} = az;