function [Bx By] = j2B(sJZ,sX,sY,sZ,fX,fY)
% tool.j2B Calculate magnetic field from Biot-Savart's law.   
%           (See also ExB.V2B for 3D fields.)
%   [Bx By] = V2B(jz,sX,sY,fX,fY);
%   
%   Biot-Savart's law: B = mu0/4pi*intdV(jxr/r^2)
%   js     - source current
%   sX, sY - source coordinates
%   fX, fY - field coordinates

% Physical constants
mu0 = 1.2566e-6;
e = 1.6022e-19;

% Grid volume
dsX = sX(1,1,1)-sX(1,2,1); % m
dsY = sY(1,1,1)-sY(2,1,1); % m
dsZ = sZ(1,1,1)-sZ(1,1,2); % m; min([dsX dsY]); % m
volume = dsX*dsY*dsZ*1e6; % km^3

%sZ = max(max(max(sX)));

% Radial source coordinates
sR = sqrt(sX.^2 + sY.^2); % km

% The magnetic field matrices
Bx = zeros(size(fX)); 
By = zeros(size(fY));

[n1 n2 n3] = size(fX);
nn = n1*n2*n3;
nprog = 0;
nprint = 0;
disp('Progress: ');
for xx = 1:n1         
    for yy = 1:n2     
        for zz = 1:n3       
            nprog=nprog+1;
            if 100*nprog/nn>nprint
                fprintf([num2str(100*nprog/nn,'%.0f'),' ']);
                nprint=nprint+10;
            end
            X = (sX - fX(xx,yy))*1e3; % m
            Y = (sY - fY(xx,yy))*1e3; % m   
            Z = (sZ - 0)*1e3; % m
            R = sqrt(X.^2 + Y.^2 + Z.^2); % m
            if R < 1e3; continue; end % smaller then 1 km, then skip
            Bx(xx,yy) = nansum(nansum(nansum((-sJZ.*Y)./R./R./R*volume*mu0/4/pi)));
            By(xx,yy) = nansum(nansum(nansum((sJZ.*X)./R./R./R*volume*mu0/4/pi)));
        end
    end
end