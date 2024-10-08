function varargout = integrate(lr,lz,phi0,B0,t,L,x0,y0,z0,vx0,vy0,vz0) 
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
%tic
% Constants and initial values
lr = lr*1e3;  % km -> m
lz = lz*1e3;  % km -> m
B0 = B0*1e-9; % nT -> T
L = L*1e3;    % km -> m

if 0 % this is done outside now
    L = lz*7; % Length of box
    T = L/vtz; % Time to pass the box
    fce=units.e*B0/units.me/2/pi; % gyro freq. Hz
    tce=1/fce; % gyro period. s
    nce = T/tce; % number of gyroperiods
    ntpertce = 7000; % number of timesteps per gyroperiod
    nt = round(ntpertce*nce); % number of timesteps
    t=linspace(0,T,nt)'; % time vector    
end
nt = numel(t);
dt = diff(t(1:2)); % time step

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

% Initial positions and velocities
% (input is in km and km/s)
x(1) = x0*1e3;   % m
y(1) = y0*1e3;   % m 
z(1) = z0*1e3;   % m
vx(1) = vx0*1e3; % m/s
vy(1) = vy0*1e3; % m/s
vz(1) = vz0*1e3; % m/s

%toc
%tic
% Calculating trajectory
nTurn = 0;
nBounce = 0;
for ii=1:(nt-1)  
    vx(ii+1)=vx(ii)+dt*(-units.e/units.me)*((x(ii)/lr^2).*phi0.*exp(-0.5*(x(ii)/lr).^2-0.5*(y(ii)/lr).^2-0.5*(z(ii)/lz).^2) + vy(ii)*B0);
    vy(ii+1)=vy(ii)+dt*(-units.e/units.me)*((y(ii)/lr^2).*phi0.*exp(-0.5*(x(ii)/lr).^2-0.5*(y(ii)/lr).^2-0.5*(z(ii)/lz).^2) - vx(ii)*B0);
    vz(ii+1)=vz(ii)+dt*(-units.e/units.me)*((z(ii)/lz^2).*phi0.*exp(-0.5*(x(ii)/lr).^2-0.5*(y(ii)/lr).^2-0.5*(z(ii)/lz).^2));

    x(ii+1)=x(ii)+vx(ii)*dt;    
    y(ii+1)=y(ii)+vy(ii)*dt;   
    z(ii+1)=z(ii)+vz(ii)*dt;
    
    % Move to other side of box if electron overshoots
    if z(ii+1)>L/2 
        z(ii+1) = z(ii+1)-L;
        nTurn = nTurn + 1;
    elseif z(ii+1)<-L/2 
        z(ii+1) = z(ii+1)+L;
        nTurn = nTurn + 1;
    end
    % Count the number of bounces the electron does
    if sign(vz(ii+1)) ~= sign(vz(ii)) 
        nBounce = nBounce + 1;
    end
end

% Prepare output
varargout{1} = x;
varargout{2} = y;
varargout{3} = z;
varargout{4} = vx;
varargout{5} = vy;
varargout{6} = vz;
varargout{7} = nTurn;
varargout{8} = nBounce;