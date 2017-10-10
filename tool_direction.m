function [x y z correlation PhiE Bz Ek En ufEn ufEk] = tool_direction(B,E,angles,f)
% IRF_MATCH_PHIB_DIR Get propagation direction by matching dBpar and phi.
%   Tries different propagation directions and finds the direction 
%   perpendicular to the magnetic field that gives the best correlation
%   between the electrostatic potential and the parallel wave magnetic field
%   according to phi_E = waveB*B0/n*e*mu0.
%
%   [x y z correlation phiE Bz Ek En ufEn ufEk] = IRF_MATCH_PHIB_DIR(B,E,angles,f)
%
%   Input 
%       B - magnetic field (to be filtered if f is given)
%       E - electric field (to be filtered if f is given)
%       angles - the angles in degrees to try (1-180 default)
%       f - filter frequency
% 
%   Returns 
%       x - normal direction (size: n_tries x 3)
%       y - propagation direction
%       z - magnetic field direction.
%       correlation - correlation vector
%       phiE - potential
%       Bz - wave magnetic field in parallel direction
%       Ek - wave electric field in propagation direction
%       En - wave electric field in propagation normal direction
%       ufEk - electric field in propagation direction
%       ufEn - electric field in propagation normal direction
% 
%   See also IRF_MATCH_PHIB_V, IRF_MATCH_PHIB_VIS

% Resample B to E if they have different size
if size(B,1) ~= size(E,1)
    B=irf_resamp(B,E);
end

% Filter if f is given, otherwise assume it is filtered
if exist('f','var') 
    BAC=irf_filt(B,f,0,450,5);
    EAC=irf_filt(E,f,0,450,5);
else
    BAC=B;
    EAC=E;
end

% If no angles are specified, set 1:180 as default
if ~exist('angles','var') 
    angles=1:180;
end
na=length(angles); % number of angles

% Set up coordinate systems
if 1
    z=repmat(irf_norm(mean(B(:,2:4),1)),na,1); % B/z direction, tries*3
    y=irf_norm(irf_cross(irf_cross(z,[1 0 0]),z)); % perp1
    x=irf_norm(irf_cross(y,z)); % perp2
    
    theta=(0:2*pi/tries:2*pi-pi/na)'; % angles
    xn=irf_norm(x.*repmat(cos(theta),1,3)+y.*repmat(sin(theta),1,3));
    y=irf_cross(z,xn);
    x=xn;
end
    
% Field aligned B
Bz=irf_dot(BAC,z(1,:));

% Allocate correlations
correlation=zeros(na,1);

% Allocate vectors, 4 first used for illustration
Ek=[EAC(:,1) zeros(size(E,1),na)]; % field to integrate
En=[EAC(:,1) zeros(size(E,1),na)]; % normal field
ufEk=[E(:,1) zeros(size(E,1),na)]; % unfiltered
ufEn=[E(:,1) zeros(size(E,1),na)]; % unfiltered
PhiE=[EAC(:,1) zeros(size(E,1),na)]; % potential

% Integrate E in all x-directions
for k=1:tries   
    Ek(:,k+1)=irf_dot(EAC(:,(2:4)),x(k,:));
    En(:,k+1)=irf_dot(EAC(:,(2:4)),y(k,:));    
    ufEk(:,k+1)=irf_dot(E(:,(2:4)),x(k,:));
    ufEn(:,k+1)=irf_dot(E(:,(2:4)),y(k,:));
    
    % Get Phi_E = int(Ek), there's no minus since the field is integrated
    % in the opposite direction of the wave propagation direction.
    intE=irf_integrate([EAC(:,1) Ek(:,k+1)]);    
    PhiE(:,k+1)=intE(:,2)-mean(intE(:,2));
       
    % Get correlation
    correlation(k,1)=xcorr(PhiE(:,k+1),Bz(:,2),0,'coeff');    
end