%% Make current matrix
% Integrate the magnetic field everywhere, using cartesian velocities
% Biot-Savart's law
% B = mu0/4pi*intdV(jxr/r^2)

e = 1.6022e-19;

toLoad = '20140802T231112-2007-08-31-500V_1600eV_9_5-5074743'; 
loadPath = '/Users/Cecilia/Research/EH/TestParticleSimulation/Spis/'; 
load([loadPath toLoad],'sumMVXxyz','sumMVYxyz','sumMVZxyz',...
    'sumMxyz','zGrid','rGrid','xSurf','ySurf','zSurf','phi0',...
    'lr','lz','B0','Tper','zlim','rlim','n','nz','nr');

re = sqrt(10*Tper)/B0;
rlimproc = (rlim-5*re)/rlim;

% The source velocities
%sim = 0;
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
        limits = [0.02 0.20 0.80 0.98];
        zlimits = -zlim+zlim*limits*2;
        zlimits = [-zlim -11 11 zlim];
        int1 = find(zSurf>(zlimits(1)) & zSurf<(zlimits(2)));
        int2 = find(zSurf>(zlimits(3)) & zSurf<(zlimits(4)));
        VXxyz = sumMVXxyz./sumMxyz;
        VYxyz = sumMVYxyz./sumMxyz;
        VZxyz = sumMVZxyz./sumMxyz;
        VXnoise = repmat(mean(VXxyz(:,:,[int1, int2]),3),1,1,nz);
        VYnoise = repmat(mean(VYxyz(:,:,[int1, int2]),3),1,1,nz);
        VZnoise = repmat(mean(VZxyz(:,:,[int1, int2]),3),1,1,nz);
        VXfix = VXxyz-VXnoise;
        VYfix = VYxyz-VYnoise;
        VZfix = VZxyz-VZnoise;
        sVX = VXfix;
        sVY = VYfix;
        sVZ = VZfix;
end

sJX = -e*n*1e6*sVX;
sJY = -e*n*1e6*sVY;
sJZ = -e*n*1e6*sVZ;

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

% The field coordinates, r-theta-z generated
if ~exist('rbin','var'); rbin = 7; end
if ~exist('thbin','var'); thbin = 7; end
if ~exist('zbin','var'); zbin = 7; end

zbin = 80;
rg = linspace(0,rlim,rbin+1); % km
thg = linspace(0,360,thbin+1); % deg
zg = linspace(-zlim,zlim,zbin+1)*3; % km
doOneR = 1;
if exist('doOneR','var')
    if doOneR
        rg = [16 18];        
    end
end
[fR,fTH,fZ] = cn.meshgrid(rg(1:end-1) + 0.5*diff(rg(1:2)),...
                          thg(1:end-1)+ 0.5*diff(thg(1:2)),...
                          zg(1:end-1) + 0.5*diff(zg(1:2))); 
%if rbin==1; fR = fR*0; end
fX = fR.*cosd(fTH); 
fY = fR.*sind(fTH);   

% The magnetic field matrices
%Bx = zeros(size(fX)); 
%By = zeros(size(fY));
%Bz = zeros(size(fZ));

tic
[Bx By Bz] = ExB.V2B(sJX,sJY,sJZ,sX,sY,sZ,fX,fY,fZ);
toc

% Get radial magnetic field
%xg0 = find(fX>=0); % right half xy-plane
%xs0 = find(fX<0); % left half xy-plane
%th = zeros(size(fX));
%th(xg0) = atand(fY(xg0)./fX(xg0));
%th(xs0) = atand(fY(xs0)./fX(xs0)) + 180;        
%Br = Bx.*cosd(th) + By.*cosd(th); % azimuthal velocity     

% Save magnetic field as either simulation field or model field
switch sim
    case 1 % simulation
        simBx = Bx; simBy = By; simBz = Bz; %simBr = Br;
        simRr = fR; simRx = fX; simRy = fY; simRz = fZ;
        simVX = sVX;
    case 0 % model
        modBx = Bx; modBy = By; modBz = Bz; %modBr = Br; 
        modRr = fR; modRx = fX; modRy = fY; modRz = fZ; 
        modVX = sVX;
    case 2 % noise fixed
        fixBx = Bx; fixBy = By; fixBz = Bz; %fixBr = Br; 
        fixRr = fR; fixRx = fX; fixRy = fY; fixRz = fZ; 
        fixVX = sVX;
end