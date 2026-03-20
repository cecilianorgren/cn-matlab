function out = distance_stein(sigma1,sigma2)
% distance_stein: Stein divergence between 1D Gaussians' covariance scalars.
%
% Inputs
%   sigma1, sigma2 : variances (i.e. covariance scalars, NOT standard deviations)
%
% For 1D, covariances are variances:
%   S = (S1+S2)/2
%
% Stein divergence:
%   D = log(det(S)) - 0.5*(log(det(S1)) + log(det(S2)))
% With scalars this reduces to:
%   D = log(S) - 0.5*(log(S1) + log(S2))

var1 = sigma1;
var2 = sigma2;
var_shared = (var1 + var2)/2;

out = log(var_shared) - 0.5*(log(var1) + log(var2));



