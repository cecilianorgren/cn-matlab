function [correlation,phi_E,phi_B] = irf_match_phib_v(B0,Bz,intEdt,n,v)
% IRF_MATCH_PHIB_V Get propagation velocity by matching dBpar and phi.
%   Used together with irf_match_phibe_dir.m. Finds best match in amplitude 
%   given, B0, dB_par, phi, propagation direction implied, for specified n 
%   and v given as vectors. Returns a matrix of correlations and the two
%   potentials that were correlated.
%
%   [correlation,phi_e,phi_B]=IRF_MATCH_PHIBE_V(B0,Bpar,intEdt,n,v)
%   
%   Input
%       B0 - average background magnetic field
%       B_par - parallel wave magnetic field
%       intEdt - int(E)dt (from highest correlation with irf_match_phib_dir.m)
%       n - vector of densities
%       v - vector of velocities
%
%   Output
%       correlation - correlation matrix (nn x nv)
%       phi_B = B0*dB_par/n_e*e*mu0
%       phi_E = int(E)dt*v (dl=-vdt => -dl=vdt)
%
%   Examples:
%       angles=1:180;
%       f_highpass=7;
%       [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
%       i_dir=find(corr(:,1)==max(corr(:,1)));
%       direction=x(i_dir,:);
%
%       n=linspace(0.01,0.1,10);
%       v=linspace(200,2000,10);
%       [corr_v,phi_E,phi_B] = irf_match_phib_v(B0,Bz,intEdt(:,1+i_dir),n,v)
%       i_v=find(corr_v(:,1)==min(corr_v(:,1)));
%       velocity=v(i_v);
%
%   See also IRF_MATCH_PHIBE_DIR, IRF_MATCH_PHIBE_VIS

% Define constants
mu0=4*pi*1e-7;
n=n*1e6; % density in #/m^3
e=1.6e-19;

% Allocate correlations matrix 
nn=length(n);
nv=length(v);
correlation=zeros(nn,nv); % rows: n, cols: v

% Setup potentials 
phi_E=[intEdt(:,1) intEdt(:,2:end)*torow(v)]; % depends on v
phi_B=[Bz(:,1) Bz(:,2)*B0*1e-18/mu0/e*(1./torow(n))]; % depends on n

% Get correlation
for k=1:nn;
    for p=1:nv;
        correlation(k,p)=sum((phi_E(:,1+p)-phi_B(:,1+k)).^2);
    end
end