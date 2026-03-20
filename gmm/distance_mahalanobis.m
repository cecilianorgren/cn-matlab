function out = distance_mahalanobis(mu1,sigma1,mu2,sigma2)
% distance_mahalanobis: squared Mahalanobis distance between 1D Gaussians.
%
% Inputs
%   sigma1, sigma2 : variances (i.e. covariance scalars, NOT standard deviations)
%
% Output
%   out : DM2 = (mu1-mu2)^2 / ((sigma1 + sigma2)/2)
%

var_shared = (sigma1 + sigma2)/2;
DM2 = (mu1-mu2)*var_shared^(-1)*(mu1-mu2)';

out = DM2;



