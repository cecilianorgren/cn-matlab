%% Make current matrix
% Integrate the magnetic field everywhere, using cartesian velocities
% Biot-Savart's law
% B = mu0/4pi*intdV(jxr/r^3)
mu0 = 1.2566e-6;
e = 1.6022e-19;

% Source grid/coordinates
nx = 100; ny = 100; %nz = 100;
xx = linspace(-lr*4,lr*4,nx); % km
yy = linspace(-lr*4,lr*4,ny); % km
zz = linspace(-lz*nlz,lz*nlz,nz); % km
volume = (xx(2)-xx(1))*(yy(2)-yy(1))*(zz(2)-zz(1))*1e9; % m^-3
[sX,sY,sZ] = cn.meshgrid(xx,yy,zz);
sR = sqrt(sX.^2 + sY.^2); % km

% The source velocities, model ExB-drift
EX = sX./(lr.^2).*phi0.*exp(-0.5*(sX/lr).^2-0.5*(sY/lr).^2-0.5*(sZ/lz).^2)*1e-3; % V
EY = sY./(lr.^2).*phi0.*exp(-0.5*(sX/lr).^2-0.5*(sY/lr).^2-0.5*(sZ/lz).^2)*1e-3; % V        
EZ = sZ./(lz.^2).*phi0.*exp(-0.5*(sX/lr).^2-0.5*(sY/lr).^2-0.5*(sZ/lz).^2)*1e-3; % V        
sVX = EY/(B0*1e-9); % m/s 
sVY = -EX/(B0*1e-9); % m/s 
sVZ = sVY * 0;

% The field coordinates
% make only along the z axis
[fX,fY,fZ] = cn.meshgrid(0,0,zz); 

% The magnetic field matrices
Bx = zeros(size(fX)); 
By = zeros(size(fY));
Bz = zeros(size(fZ));

tic
progress = 0.1;
for xxx = 1:size(fX,1)    
    for yyy = 1:size(fX,2)
        for zzz = 1:size(fX,3)            
            if (xxx*yyy*zzz)/numel(fX) > progress                
                fprintf([num2str(100*progress,'%.f') ' %']);
                progress = progress + 0.1;
            end
            X = (fX(xxx,yyy,zzz) - sX)*1e3; % m
            Y = (fY(xxx,yyy,zzz) - sY)*1e3; % m
            Z = (fZ(xxx,yyy,zzz) - sZ)*1e3; % m         
            R = sqrt(X.^2 + Y.^2 + Z.^2); % m
            if R < 1e3; continue; end % smaller then 1 km
            Bx(xxx,yyy,zzz) = nansum(nansum(nansum(-e*n*1e6*(sVY.*Z-sVZ.*Y)./R./R./R*volume*mu0/4/pi)));            
            By(xxx,yyy,zzz) = nansum(nansum(nansum(-e*n*1e6*(sVZ.*X-sVX.*Z)./R./R./R*volume*mu0/4/pi)));
            Bz(xxx,yyy,zzz) = nansum(nansum(nansum(-e*n*1e6*(sVX.*Y-sVY.*X)./R./R./R*volume*mu0/4/pi)));            
        end        
    end
end
toc

% Get radial magnetic field
%xg0 = find(fX>=0); % right half xy-plane
%xs0 = find(fX<0);  % left half xy-plane
% th = zeros(size(fX));
% th(xg0) = atand(fY(xg0)./fX(xg0));
% th(xs0) = atand(fY(xs0)./fX(xs0)) + 180;        
% Br = Bx.*cosd(th) + By.*cosd(th); % azimuthal velocity     
% 
% % Save magnetic field as either simulation field or model field
% if sim; 
%     simBx = Bx; simBy = By; simBz = Bz; simBr = Br;
%     simRr = fR; simRx = fX; simRy = fY; simRz = fZ;
% else
%     modBx = Bx; modBy = By; modBz = Bz; modBr = Br;
%     modRr = fR; modRx = fX; modRy = fY; modRz = fZ;
% end
% end