function [wS,muS,SigS] = merge_gaussians(w0,mu0,Sig0,S)

wS = sum(w0(S));
muS = (w0(S) * mu0(S,:)) / wS; % 1x3

SigS = zeros(3,3);
for ii = 1:numel(S)
  k = S(ii);
  dmu = (mu0(k,:) - muS)';
  SigS = SigS + w0(k) * (Sig0(:,:,k) + (dmu*dmu')); % law of shared covariances
end
SigS = SigS / wS;