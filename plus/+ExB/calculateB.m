%% Make current matrix
% Integrate the magnetic field everywhere, using cartesian velocities
% Biot-Savart's law
% B = mu0/4pi*intdV(jxr/r^2)
for sim = [0 1]
mu0 = 1.2566e-6;
e = 1.6022e-19;
rlimproc = 0.7;

% Grid volume
volume = (xGrid(2)-xGrid(1))*(yGrid(2)-yGrid(1))*(zGrid(2)-zGrid(1))*1e9;

% The source velocities
if sim % simulation
    sVAZ = sumMVxyz./sumMxyz; % m/s
    sVX = sumMVXxyz./sumMxyz; % m/s
    sVY = sumMVYxyz./sumMxyz; % m/s
    sVZ = (sumMVZxyz./sumMxyz)*0; % m/s

    sVX = permute(sVX,[2 1 3]);
    sVY = permute(sVY,[2 1 3]);
    sVZ = permute(sVZ,[2 1 3]);
else  % model ExB drift
    [XS,YS,ZS] = meshgrid(xSurf,ySurf,zSurf);
    Ex = XS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V
    Ey = YS./(lr.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    Ez = ZS./(lz.^2).*phi0.*exp(-0.5*(XS/lr).^2-0.5*(YS/lr).^2-0.5*(ZS/lz).^2)*1e-3; % V        
    sVX = Ey/(B0*1e-9); % m/s 
    sVY = -Ex/(B0*1e-9); % m/s 
    sVZ = sVY * 0;
end

% The source coordinates
[sX,sY,sZ] = meshgrid(xSurf,ySurf,zSurf); % km
sX = sX*1e3; sY = sY*1e3; sZ = sZ*1e3;    % m
sR = sqrt(sX.^2 + sY.^2);                 % m

% The field coordinates
fX = sX; 
fY = sY;
fZ = sZ;

% Put to NaN all values outside a certain radius
rmax = rlimproc*rlim;   % km
rind = find(sR>rmax*1e3);
sVX(rind) = NaN;
sVY(rind) = NaN;
sVZ(rind) = NaN;

% Try to take away some noise
if 0
sV = sqrt(sVX.^2 + sVY.^2 + sVZ.^2);
vhigh = find(sV>2000*1e3);
sVX(vhigh) = NaN;
sVY(vhigh) = NaN;
sVZ(vhigh) = NaN;
end

% The magnetic field matrices
Bx = zeros(size(sX)); 
By = zeros(size(sY));
Bz = zeros(size(sZ));

% The source currents
sJX = -e*n*1e6*sVX; % A/m^2 ?
sJY = -e*n*1e6*sVY;
sJZ = -e*n*1e6*sVZ;

if 1 % Without for-loops. Use VX, VY, VZ
xind = 10:2:numel(xSurf)-10;  
yind = 10:2:numel(ySurf)-10;  
zind = 60;%20:20:numel(zSurf); % fix(numel(zSurf)/2);
xn = 0;
tic
for xx = xind  
    xn = xn + 1; disp(['xx = ' num2str(xn) '/' num2str(numel(xind))])
    for yy = yind        
        if sqrt(sX(xx,yy,1).^2+sY(xx,yy,1).^2)*1e-3<rlim*rlimproc
        for zz = zind                        
            fR = [sX(xx,yy,zz) sY(xx,yy,zz) sZ(xx,yy,zz)];
            X = fX(xx,yy,zz) - sX;
            Y = fY(xx,yy,zz) - sY;
            Z = fZ(xx,yy,zz) - sZ;            
            R = sqrt(X.^2 + Y.^2 + Z.^2);
            if R == 0; continue; end
            % JxR_x = sJY.*Z-sJZ.*Y;
            % JxR_y = sJZ.*X-sJX.*Z;
            % JxR_z = sJX.*Y-sJY.*X;
            newBx = nansum(nansum(nansum((sJY.*Z-sJZ.*Y)./R./R./R*volume*mu0/4/pi)));
            newBy = nansum(nansum(nansum((sJZ.*X-sJX.*Z)./R./R./R*volume*mu0/4/pi)));
            newBz = nansum(nansum(nansum((sJX.*Y-sJY.*X)./R./R./R*volume*mu0/4/pi)));
            Bx(xx,yy,zz) = newBx;            
            By(xx,yy,zz) = newBy;
            Bz(xx,yy,zz) = newBz;            
        end
        end
    end
end
toc
end

if sim; simBx = Bx; simBy = By; simBz = Bz; 
else modBx = Bx; modBy = By; modBz = Bz; end
end