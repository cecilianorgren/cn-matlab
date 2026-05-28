function out = gmm_monte_carlo_sampling(g,N,group)
% Generate N samples from the specified groups of the given Gaussian
% Mixture Model

mu = g.mu(group,:); % km/s
Sigma = g.Sigma(:,:,group); % (km/s)^2
w = g.ComponentProportion(group)/sum(g.ComponentProportion((group)));
%numComp = g.NumComponents;

p = zeros(N,3);
Ncum = 0;
for ik = 1:numel(group)
  Nk = round(N*w(ik)); % Generate particles based on component proportion
  p(Ncum+[1:Nk],:) = mvnrnd(mu(ik,:), Sigma(:,:,ik),Nk);
  %histogram(p(Ncum+[1:Nk],1),-2500:50:2500)
  %hold on
  Ncum = Ncum + Nk;
end

out = p;