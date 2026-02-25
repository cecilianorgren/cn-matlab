function [correlation,phi_E,phi_B] = irf_match_phibe_v(B0,Bz,intEdt,n,v)
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
%       % Direction
%       angles=1:3:360;
%       f_highpass=7;
%       [x y z corr_dir intEdt Bz B0 dEk dEn Ek En]=irf_match_phibe_dir(B,E,angles,f_highpass);
%       i_dir=find(corr_dir(:,1)==max(corr_dir(:,1)));
%       direction=x(i_dir,:);
%
%       % Velocity and density
%       n=linspace(0.01,0.1,100);
%       v=linspace(200,2000,100);
%       [corr_v,phi_E,phi_B]=IRF_MATCH_PHIBE_V(B0,Bz,intEdt(:,[1 1+i_dir]),n,v);
%       i_v=find(corr_v(:,1)==min(corr_v(:,1)));
%       velocity=v(i_v);
%   
%       % Figures
%       gif_stuff_dir = irf_match_phibe_vis('direction',x,y,z,corr_dir,intEdt,Bz,En,Ek);
%       imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);%       
%       imwrite(gif_stuff_dir.im,gif_stuff_dir.map,'mygif_dir.gif','DelayTime',0.01,'LoopCount',inf);
%
%       i_n=50; % if more than one densitiy, choose one by specifying index
%       gif_stuff_v = irf_match_phibe_vis('velocity',phi_E,phi_B(:,[1 i_n]),v,n(i_n));
%       imwrite(gif_stuff_v.im,gif_stuff_v.map,'mygif_v.gif','DelayTime',0.01,'LoopCount',inf);
%   
%       figure; h=axes;
%       axis_handle = irf_match_phibe_vis('velocity/density',h,n,v,corr_v);
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
phi_E=[intEdt(:,1) intEdt(:,2:end)*torow(double(v))]; % depends on v
phi_B=[Bz(:,1) Bz(:,2)*B0*1e-18/mu0/e*(1./torow(double(n)))]; % depends on n

% Get correlation
if 0
for k=1:nn;
    for p=1:nv;        
        correlation(k,p)=sum((phi_E(:,1+p)-phi_B(:,1+k)).^2);
    end
end
else % test wighted correlation
    nbins = 10;
    nint = ceil(numel(phi_B(:,1))/nbins);
for k=1:nn;
    for p=1:nv; 
        locE = phi_E(:,1+p);
        locB = phi_B(:,1+k);
        padsize = ceil((numel(locE)+1)/nbins)*nbins;
        padnum = padsize-numel(locE);
        padlocE = padarray(locE,[padnum 0],NaN,'post');
        matlocE = reshape(padlocE,padsize/nbins,nbins);
        padlocB = padarray(locB,[padnum 0],NaN,'post');
        matlocB = reshape(padlocB,padsize/nbins,nbins);
        
        for bin = 1:nbins                 
            loc_xcorr(bin) = xcorr(matlocE(:,bin),matlocB(:,bin),0,'coeff');
            loc_vcorr(bin) = sum((matlocE(:,bin)-matlocB(:,bin)).^2)*loc_xcorr(bin);
        end       
        correlation(k,p)=nansum(loc_vcorr);
    end
end
end
end