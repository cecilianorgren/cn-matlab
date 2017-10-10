lr = 10e3;
lz = 10e3;
phi0 = 300;

Ex0 = @(x,y,z,lr,lz,phi0)(x/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ey0 = @(x,y,z,lr,lz,phi0)(y/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ez0 = @(x,y,z,lr,lz,phi0)(z/lr^2).*phi0.*exp(-0.5*(x/lr).^2-0.5*(y/lr).^2-0.5*(z/lz).^2);
Ex = @(x,y,z) Ex0(x,y,z,lr,lz,phi0); % V/m
Ey = @(x,y,z) Ey0(x,y,z,lr,lz,phi0); % V/m
Ez = @(x,y,z) Ez0(x,y,z,lr,lz,phi0); % V/m

xylim=20;

% Plot electric field vectors
nEx = 5;
nEy = 5;
nEz = 2;
xE = linspace(-xylim,xylim,nEx);
yE = linspace(-xylim,xylim,nEy);
zE = linspace(-lz*1e-3,lz*1e-3,nEz);
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

quiver3(XE,YE,ZE,qXE,qYE,qZE,0.5)
xlabel('x')
ylabel('y')
zlabel('z')