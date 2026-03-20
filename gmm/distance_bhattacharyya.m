function out = distance_bhattacharyya(mu1,sigma1,mu2,sigma2)
% distance_bhattacharyya: Bhattacharyya distance between 1D Gaussians.
%
% Inputs
%   sigma1, sigma2 : variances (i.e. covariance scalars, NOT standard deviations)
%
% For two Gaussians with covariances S1,S2 and mean difference dmu:
%   S = (S1+S2)/2
%   D_B = 1/8 * dmu^T S^{-1} dmu + 1/2 * ( log(det(S)) - 0.5*(log(det(S1))+log(det(S2))) )
%
% In terms of helper distances:
%   DM2 = dmu^T S^{-1} dmu
%   DS  = log(det(S)) - 0.5*(log(det(S1))+log(det(S2)))  (Stein divergence)
% so:
%   DB = 1/8 * DM2 + 1/2 * DS

DS = distance_stein(sigma1,sigma2);
DM2 = distance_mahalanobis(mu1,sigma1,mu2,sigma2);
out = (1/8)*DM2 + 0.5*DS;



