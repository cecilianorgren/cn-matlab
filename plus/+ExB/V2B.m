function [Bx By Bz] = V2B(sJX,sJY,sJZ,sX,sY,sZ,fX,fY,fZ)
% ExB.V2B Calculate magnetic field from Biot-Savart's law.   
%   [Bx By Bz] = V2B(sVX,sVY,sVZ,fVX,fVY,fVZ);
%   
%   Biot-Savart's law: B = mu0/4pi*intdV(jxr/r^2)
%   sJX, sJY, sJZ - source velocitites
%   sX, sY, sZ - source coordinates, has to correspond to the 
%                simulation xyz binning system
%   fX, fY, fZ - field coordinates

% Physical constants
mu0 = 1.2566e-6;
e = 1.6022e-19;

% Grid volume
dsX = sX(1,1,1)-sX(2,1,1); % m
dsY = sY(1,1,1)-sY(1,2,1); % m
dsZ = sZ(1,1,1)-sZ(1,1,2); % m
volume = dsX*dsY*dsZ*1e9; % km^3

% Radial source coordinates
sR = sqrt(sX.^2 + sY.^2); % km

% The magnetic field matrices
Bx = zeros(size(fX)); 
By = zeros(size(fY));
Bz = zeros(size(fZ));

[n1 n2 n3] = size(fX);
nn = n1*n2*n3;
nprog = 0;
nprint = 0;
disp('Progress: ');
for rr = 1:n1         
    for th = 1:n2                
        for zz = 1:n3    
            nprog=nprog+1;
            if 100*nprog/nn>nprint
                fprintf([num2str(100*nprog/nn,'%.0f'),' ']);
                nprint=nprint+10;
            end
            X = (sX - fX(rr,th,zz))*1e3; % m
            Y = (sY - fY(rr,th,zz))*1e3; % m
            Z = (sZ - fZ(rr,th,zz))*1e3; % m         
            R = sqrt(X.^2 + Y.^2 + Z.^2); % m
            if R < 1e3; continue; end % smaller then 1 km, then skip
            Bx(rr,th,zz) = nansum(nansum(nansum((sJY.*Z-sJZ.*Y)./R./R./R*volume*mu0/4/pi)));            
            By(rr,th,zz) = nansum(nansum(nansum((sJZ.*X-sJX.*Z)./R./R./R*volume*mu0/4/pi)));
            Bz(rr,th,zz) = nansum(nansum(nansum((sJX.*Y-sJY.*X)./R./R./R*volume*mu0/4/pi)));            
        end        
    end
end