function out = fit_ntrunc(X,K,v_max)

maxIter = 100;

heaviside_up = @(x) 1.0*(x<=0);
pdf = @(X, mu, sigma) (normpdf(X , mu, sigma)./normcdf(v_min, mu, sigma) .* heaviside_l(sqrt(X(:,1).^2 + X(:,2).^2 + X(:,3).^2) - v_max));

converged = false;
for iIter = 1:maxIter
  [phat, phat_ci]  = mle(X , 'pdf', norm_trunc,'start', [mean(dat_normal), std(dat_normal)]);
  

end