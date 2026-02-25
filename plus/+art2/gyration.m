%% 1D:Linear acceleration
units=irf_units;
nt=500;
t=linspace(0,1,nt)';
v=zeros(nt,1);
x=zeros(nt,1);
E=10; % V/m
a=units.e*E/units.me;
dt=diff(t(1:2));
for ii=1:(nt-1)
    dv=a*dt;
    v(ii+1)=v(ii)+dv;
    x(ii+1)=x(ii)+v(ii)*dt;    
end

plot(t,v,t,x)
xlabel('t [s]')
legend('v','x')

%% 2D: vxB and E acceleration in infinite column electric field
vx=zeros(nt,1);
vy=zeros(nt,1);

dvx=zeros(nt,1);
dvy=zeros(nt,1);

x=zeros(nt,1);
y=zeros(nt,1);

ax=zeros(nt,1);
ay=zeros(nt,1);

units=irf_units;
lr = 4e3; % m
r = 10e3; % m
phi0 = 300; % V
B0 = 50e-9; % T
        
Ex0 = @(x,y,lr,phi0)(x/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2);
Ey0 = @(x,y,lr,phi0)(y/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2);
Ex = @(x,y) Ex0(x,y,lr,phi0); % V/m
Ey = @(x,y) Ey0(x,y,lr,phi0); % V/m

x(1) = lr/2;  % m
y(1) = 0;     % m  

vx(1) = 1e3*1e3;    % m/s
vy(1) = -5e3*1e3;    % m/s

fce=units.e*B0/units.me/2/pi; % Hz
tce=1/fce; % s
nt=10000;
ntg=10; % number of gyroperiods
t=linspace(0,ntg*tce,nt)';
dt=diff(t(1:2));

wait = waitbar(0,'Calculation trejectory, please wait...');
for ii=1:(nt-1)
    %waitbar(ii/nt,wait)
    ax(ii) = (-units.e/units.me)*(Ex(x(ii),y(ii)) + vy(ii)*B0);
    ay(ii) = (-units.e/units.me)*(Ey(x(ii),y(ii)) - vx(ii)*B0);   
    
    dvx(ii)=ax(ii)*dt;
    dvy(ii)=ay(ii)*dt;
    
    vx(ii+1)=vx(ii)+dvx(ii);
    vy(ii+1)=vy(ii)+dvy(ii);
    
    x(ii+1)=x(ii)+vx(ii)*dt;    
    y(ii+1)=y(ii)+vy(ii)*dt;   
    
end
close(wait)
%
%figure(108)
%ax=axis;
%
plot(0,0); %hold(ax,'on');
axl=gca;
axis(axl,'equal')

xlabel(gca,'x [km]');ylabel(gca,'y [km]');
%
if 0
for ii=1:10:nt
    pause(0.01)
    plot(axl,x*1e-3,y*1e-3,'-',x(ii)*1e-3,y(ii)*1e-3,'*'); 
    set(axl,'xlim',lr*[-1 1]*1e-3,'ylim',lr*[-1 1]*1e-3)
end
else
    plot(axl,x*1e-3,y*1e-3); 
    set(axl,'xlim',1*lr*[-1 1]*1e-3,'ylim',1*lr*[-1 1]*1e-3)
    %set(axl,'xlim',[0 -3]*1e-3,'ylim',[0 3]*1e-3)
    xlabel(axl,'x'); ylabel(axl,'y'); 
    axis equal
    %axis square
end
hold(axl,'off')
    %%
    if 0
figure(107)
for ii=1:6
    h(ii)=subplot(2,3,ii);
end
isub=1; ip=nt; % index to plot to

if 1
    hca=h(isub); isub=isub+1; 
    plot(hca,z(1:pi)*1e3,x(1:pi)*1e-3,z(1:pi)*1e-3,y(1:pi)*1e-3); 
    xlabel(hca,'t [s]')
end
if 1
    hca=h(isub); isub=isub+1; 
    plot(hca,t(1:pi),z(1:pi)); xlabel(hca,'t [s]')
end
if 1
    hca=h(isub); isub=isub+1; 
    plot(hca,z(1:pi),az(1:pi)); ylabel(hca,'a_z'); xlabel(hca,'t [s]')
end
if 1
    hca=h(isub); isub=isub+1; 
    plot(hca,z,Ex(2,0,z),z,Ey(2,0,z),z,Ez(2,0,z)); 
    xlabel(hca,'z'); ylabel(hca,'E [V/m]'); legend(hca,'x','y','z')
end
if 1
    hca=h(isub); isub=isub+1; plot(hca,z(1:pi),vz(1:pi))
    xlabel(hca,'z'); ylabel(hca,'v_z [km/s]');
end
if 1
    hca=h(isub); isub=isub+1;
end
    end

%% 2D: vxB and E acceleration in infinite column electric field
nt = 10000;
vx=zeros(nt,1);
vy=zeros(nt,1);
vz=zeros(nt,1);


dvx=zeros(nt,1);
dvy=zeros(nt,1);
dvz=zeros(nt,1);

x=zeros(nt,1);
y=zeros(nt,1);
z=zeros(nt,1);
%linspace(-30*1e3,30*1e3,nt)'*1;

ax=zeros(nt,1);
ay=zeros(nt,1);
az=zeros(nt,1);

units=irf_units;
lr = 9e3; % m
lz = 35e3;
r = 10e3; % m
phi0 = 3000; % V
B0 = 50e-9; % T
Tpar = 8000; % eV
        
Ex0 = @(x,y,z,lr,lz,phi0)(x/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ey0 = @(x,y,z,lr,lz,phi0)(y/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ex = @(x,y,z) Ex0(x,y,z,lr,lz,phi0); % V/m
Ey = @(x,y,z) Ey0(x,y,z,lr,lz,phi0); % V/m

fce=units.e*B0/units.me/2/pi; % Hz
tce=1/fce; % s
nt=20000;
ntg=8; % number of gyroperiods
T=ntg*tce;
t=linspace(0,T,nt)';
dt=diff(t(1:2));

vx(1) = -1e3*1e3;    % m/s
vy(1) = 5e3*1e3;    % m/s
vz(1) = cn_eV2v(Tpar,'ev')*1e3;    % m/s

x(1) = lr/4;  % m
y(1) = 0;     % m 
z(1) = -T*vz(1)/2;

wait = waitbar(0,'Calculation trajectory, please wait...');
for ii=1:(nt-1)
    %waitbar(ii/nt,wait)
    ax(ii) = (-units.e/units.me)*(Ex(x(ii),y(ii),1*z(ii)) + vy(ii)*B0);
    ay(ii) = (-units.e/units.me)*(Ey(x(ii),y(ii),1*z(ii)) - vx(ii)*B0);   
    az(ii) = 0;   
    
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
%
%figure(108)
%ax=axis;
%
plot(0,0); %hold(ax,'on');
axl=gca;
axis(axl,'equal')

xlabel(gca,'x [km]');ylabel(gca,'y [km]');
%
if 0
for ii=1:10:nt
    pause(0.01)
    plot3(axl,x*1e-3,y*1e-3,z*1e-3,'-',x(ii)*1e-3,y(ii)*1e-3,z*1e-3,'*'); 
    set(axl,'xlim',1.5*lr*[-1 1]*1e-3,'ylim',1.5*lr*[-1 1]*1e-3)
end
else
    plot3(axl,x*1e-3,y*1e-3,z*1e-3); 
    %set(axl,'xlim',1*lr*[-1 1]*1e-3,'ylim',1*lr*[-1 1]*1e-3)
    %set(axl,'xlim',[0 -3]*1e-3,'ylim',[0 3]*1e-3)
    xlabel(axl,'x'); ylabel(axl,'y'); zlabel(axl,'z');
    axis equal
    axis square
    view([0 0 1])
end
hold(axl,'off')





