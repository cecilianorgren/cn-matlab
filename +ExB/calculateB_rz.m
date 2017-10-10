%% Make current matrix
% Integrate the magnetic field everywhere, using cartesian velocities
% Biot-Savart's law
% B = mu0/4pi*intdV(jxr/r^2)
re = sqrt(10*Tper)/B0;

mu0 = 1.2566e-6;
e = 1.6022e-19;
rlimproc = (rlim-5*re)/rlim;

% Grid volume
volume = (xGrid(2)-xGrid(1))*(yGrid(2)-yGrid(1))*(zGrid(2)-zGrid(1))*1e9;

% The source velocities
switch sim
    case 0  % model ExB drift
        [XS,YS,ZS] = cn.meshgrid(xSurf,ySurf,zSurf);
        Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
        Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
        Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
        VXmod = Ey/(B0*1e-9); % m/s 
        VYmod = -Ex/(B0*1e-9); % m/s 
        VZmod = VYmod * 0;
        sVX = VXmod; % m/s 
        sVY = VYmod; % m/s 
        sVZ = VZmod;
    case 1 % simulation
        %sVAZ = nansumMVxyz./nansumMxyz; % m/s
        VXsim = sumMVXxyz./sumMxyz; % m/s
        VYsim = sumMVYxyz./sumMxyz; % m/s
        VZsim = sumMVZxyz./sumMxyz; % m/s
        sVX = VXsim; % m/s
        sVY = VYsim; % m/s
        sVZ = VZsim*0; % m/s
    case 2 % subtracted noise       
        sVX = VXfix;
        sVY = VYfix;
        sVZ = VZfix;
end

% The source coordinates
% This has to correspond to the simulation xyz binning system
[sX,sY,sZ] = cn.meshgrid(xSurf,ySurf,zSurf); % km
sR = sqrt(sX.^2 + sY.^2);                    % km

% Put to NaN all values outside a certain radius
rmax = rlimproc*rlim;   % km
rind = find(sR>rmax);
%rind = find(sR>45);
doTrimR = 1;
if doTrimR
    trimR = 70;
    sVX(sR>trimR) = NaN;
    sVY(sR>trimR) = NaN;
    sVZ(sR>trimR) = NaN;
    sVX(sR<1) = NaN;
    sVY(sR<1) = NaN;
    sVZ(sR<1) = NaN;
    rind = find(sR>45);
end
doTrimZ = 1;
if doTrimZ
    zi = 5*lz;
    sVX(sZ>zi) = NaN;
    sVX(sZ<-zi) = NaN;
    sVY(sZ>zi) = NaN;
    sVY(sZ<-zi) = NaN;
    sVZ(sZ>zi) = NaN;
    sVZ(sZ<-zi) = NaN;
end

%rindmin = find(sR<3.1);
%sVX(rind) = 0;
%sVY(rind) = 0;
%sVZ(rind) = 0;

% The field coordinates, r-theta-z generated
if ~exist('rbin','var'); rbin = 7; end
if ~exist('thbin','var'); thbin = 7; end
if ~exist('zbin','var'); zbin = 7; end
%rbin=1; thbin=1; zbin=100;
rg = linspace(0,rlim,rbin+1); % km
thg = linspace(0,360,thbin+1); % deg
zg = linspace(-zlim,zlim,zbin+1); % km
if exist('doOneR','var')
    if doOneR
        rg = [16 18];
        zg = linspace(-1.5*zlim,1.5*zlim,zbin+1); % km
        rbin=1;
    end
end
[fR,fTH,fZ] = cn.meshgrid(rg(1:end-1) + 0.5*diff(rg(1:2)),...
                          thg(1:end-1)+ 0.5*diff(thg(1:2)),...
                          zg(1:end-1) + 0.5*diff(zg(1:2))); 
if rbin==1; fR = fR*0; end
fX = fR.*cosd(fTH); 
fY = fR.*sind(fTH);   

% The magnetic field matrices
Bx = zeros(size(fX)); 
By = zeros(size(fY));
Bz = zeros(size(fZ));

% Exclude some source points, in order to see which region contributes to
% what
if ~exist('z_maxes')
    z_maxes = [min(zSurf) max(zSurf)];
    r_maxes = [min(rSurf) max(rSurf)];
end

tic
if sim == 0; disp('modB: '); elseif sim == 1; disp('simB:'); elseif sim == 2; disp('fixB:'); end
for rr = 1:rbin     
    try fprintf([num2str(100*rr/rbin,'%.1f'),' ']), end    
    
    for th = 1:thbin                
        for zz = 1:zbin                           
            % there was no minus in front before, but it seemed to work, 
            % now it doesn?t anymore so i put the minus there, correctly
            X = -(fX(rr,th,zz) - sX)*1e3; % m 
            Y = -(fY(rr,th,zz) - sY)*1e3; % m
            Z = -(fZ(rr,th,zz) - sZ)*1e3; % m         
            R = sqrt(X.^2 + Y.^2 + Z.^2); % m
            if R < 1e3; continue; end % smaller then 1 km
            Bx(rr,th,zz) = nansum(nansum(nansum(-e*n*1e6*(sVY.*Z-sVZ.*Y)./R./R./R*volume*mu0/4/pi)));            
            By(rr,th,zz) = nansum(nansum(nansum(-e*n*1e6*(sVZ.*X-sVX.*Z)./R./R./R*volume*mu0/4/pi)));
            Bz(rr,th,zz) = nansum(nansum(nansum(-e*n*1e6*(sVX.*Y-sVY.*X)./R./R./R*volume*mu0/4/pi)));            
        end        
    end
end
toc

% Get radial magnetic field
xg0 = find(fX>=0); % right half xy-plane
xs0 = find(fX<0); % left half xy-plane
th = zeros(size(fX));
th(xg0) = atand(fY(xg0)./fX(xg0));
th(xs0) = atand(fY(xs0)./fX(xs0)) + 180;        
Br = Bx.*cosd(th) + By.*cosd(th); % azimuthal velocity     

% Save magnetic field as either simulation field or model field
switch sim
    case 1 % simulation
        simBx = Bx; simBy = By; simBz = Bz; simBr = Br;
        simRr = fR; simRx = fX; simRy = fY; simRz = fZ;
        simVX = sVX;
    case 0 % model
        modBx = Bx; modBy = By; modBz = Bz; modBr = Br; 
        modRr = fR; modRx = fX; modRy = fY; modRz = fZ; 
        modVX = sVX;
    case 2 % noise fixed
        fixBx = Bx; fixBy = By; fixBz = Bz; fixBr = Br; 
        fixRr = fR; fixRx = fX; fixRy = fY; fixRz = fZ; 
        fixVX = sVX;
end