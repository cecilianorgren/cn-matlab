%tint_psbl = irf.tint('2017-07-06T13:54:05.520Z/2017-07-06T13:54:05.640Z');
t = irf_time('2017-07-06T00:54:19.020108000Z','utc>EpochTT');
t = irf_time('2017-07-06T00:54:10.020108000Z','utc>EpochTT');
%ePDist = ePDist1(500);%.convertto('s^3/m^6');

u = irf_units;
ePDist = ePDist1;

% get index of measurement nearest t
it = interp1(ePDist.time.epochUnix,1:ePDist.length,t.epochUnix,'nearest');

M = u.me;
%M = u.mp;

% set velocity grid (same in all directions)
vmax = 5e7; % m/s
nvg = 100;
vg = linspace(-vmax,vmax,nvg);

%[VXG,VYG,VZG] = PD.v;
[VXG,VYG,VZG,Fg] = get3Ddist(ePDist1.convertto('s^3/m^6'),it,vg,M);

%%
Fg(Fg==0) = NaN;
fsurf = prctile(Fg(:),70) %
%fsurf = 2e-16;
% Fg = squeeze(PD.data);
% VXG = squeeze(VXG);
% VYG = squeeze(VYG);
% VZG = squeeze(VZG);

% isosurface does not allow for axis input, this might be an issue if you
% have many panels. It seems better to let isosurface initiate the axis
%close
isosurface(VXG*1e-4,VYG*1e-4,VZG*1e-4,Fg,fsurf); 
axis equal;
%%
hf.FaceAlpha = 0.1;
%hf.EdgeColor = 'none';
%hf.FaceColor = 'interp';
axis equal;
hold on
hf2 = patch(isosurface(VXG*1e-4,VYG*1e-4,VZG*1e-4,Fg,fsurf*1e-1)); 
hf2.FaceAlpha = 0.1;
%hf2.EdgeColor = 'none';
%hf2.FaceColor = 'interp';
hold off
%5
%colormap(prism(28))
%%
eDist = ePDist1(it);
vgr = -70e3:1e3:70e3; % km/s

%ef3D = eDist.elim([40 inf]).rebin('cart',{vgr,vgr,vgr},orient); % Rebins skymap into 3D cartesian grid

[VX,VY,VZ] = meshgrid(ef3D.depend{1},ef3D.depend{2},ef3D.depend{3});
F = permute(squeeze(ef3D.data),[2 1 3]);

%% F(F==0) = NaN;
fsurf = prctile(F(:),70); %

% isosurface does not allow for axis input, this might be an issue if you
% have many panels. It seems better to let isosurface initiate the axis

isosurface(VX*1e-4,VY*1e-4,VZ*1e-4,F,fsurf); 
axis equal;
%%
function [VXG,VYG,VZG,Fg] = get3Ddist(dist,it,vg,M)
% GET3DDIST get 3D distribution in cartesian grid (DMPA coordinates)
%  [VXG,VYG,VZG,F] = get3Ddist(dist,vg,M) get velocity mesh grids VXG, VYG,
%  VZG, given PDist dist, velocity array vg, and mass M.

% elementary charge
qe = 1.6022e-19;

emat = double(dist.depend{1}); % in eV
energy = emat(it,:);
v = sqrt(2*energy*qe/M); % m/s

% azimuthal angle
phi = double(dist.depend{2}(it,:)); % in degrees
%phi = phi+180;
%phi(phi>360) = phi(phi>360)-360;
phi = phi-180;
phi = phi*pi/180; % in radians

% elevation angle
th = double(dist.depend{3}); % polar angle in degrees
th = th-90; % elevation angle in degrees
th = th*pi/180; % in radians


nAz = length(phi);
nEle = length(th);
nV = length(v);

%
% 3D matrices for instrumental bin centers
TH = repmat(th,nV,1,nAz);       % [phi,th,v]
TH = permute(TH,[1,3,2]);       % [v,phi,th]
PHI = repmat(phi,nV,1,nEle);    % [v,phi,th]
VEL = repmat(v,nAz,1,nEle);     % [phi,v,th]
VEL = permute(VEL,[2,1,3]);     % [v,phi,th]

F = double(squeeze(dist.data(it,:,:,:)));

% instrument grid in DMPA
[VX,VY,VZ] = sph2cart(PHI,TH,VEL);

% make grid

% allow for different grids in future ?
vxg = vg; vyg = vg; vzg = vg;

[VXG,VYG,VZG] = meshgrid(vxg,vyg,vzg);

% get the cartesian distribution, this is the slow part
Fg = griddata(VX,VY,VZ,F,VXG,VYG,VZG,'linear');

end
