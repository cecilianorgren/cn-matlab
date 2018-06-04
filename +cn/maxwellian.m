function out = maxwellian(v,T,n,vd,species,varargin)
% Generates a maxwellian distribution.
%   f = cn.maxwellian(v,T,n,vd,species,optional);
%       T - eV
%       n - cc
%       v - km/s
%       vd - drift speed, km/s
%       species - 'e', 'i' or 'p'
%       optional - 1 for 3D distribution, leave empty or 0 for 1D
%                   distribution
%
%       f is a matrix with dimensions (#v)x(#n) 

e = 1.6022e-19;
kB = 1.38e-23;
me = 9.1094e-31;
m = 9.1094e-31;
mp = 1.6726e-27;

do3D = 0;
if ~isempty(varargin)
    do3D = varargin{1};
end
%do3D = 0;
disp(['do3D = ' num2str(do3D)])

if numel(v)>1
    v = torow(v)*1e3; % m/s
else
    v = v*1e3; % m/s
end
vd = vd*1e3; % m/s
if numel(n)>1  
  n = tocolumn(n)*1e6; % m^-3
else
  n = n*1e6; % m^-3
end
T = T*e/kB; % K

if any([strcmp(species,'i') strcmp(species,'p')]); % ions
    w = @(T) sqrt(2*kB*T/(m*1836)); % m/s
else % electron
    w = @(T) sqrt(2*kB*T/(m*1)); % m/s
end
if do3D; f = @(v,T,n,vd) n/((pi)^(3/2)*w(T).^3)*exp(-(v-vd).^2./w(T)./w(T));
else;    f = @(v,T,n,vd) n/((pi)^(1/2)*w(T))*exp(-(v-vd).^2./w(T)./w(T)); end


out = tocolumn(f(v,T,n,vd)); 