w = [1 1 1 1];
mu = [-500 0 900 950];     % km/s (or arbitrary units)
std = [500 500 200 300];  
sig = std.^2;

w = w(:); mu = mu(:); sig = sig(:);
w = w./sum(w);

nComp = numel(w);

threshold = 0.5;
% Calculate distances
sigVar = sig.^2; % distance helpers expect covariance/variance (NOT std)
for i = 1:nComp
  for j = 1:nComp
    d.m(i,j) = distance_mahalanobis(mu(i),sig(i),mu(j),sig(j));
    d.s(i,j) = distance_stein(sig(i),sig(j));
    d.b(i,j) = distance_bhattacharyya(mu(i),sig(i),mu(j),sig(j));
  end
end

% Check overlap

% Compute component pdfs and mixture
f_comp = zeros(K,numel(x));
for k = 1:K
  f_comp(k,:) = normpdf1(x, mu(k), std(k));
end
f_mix = (w.' * f_comp); % 1 x Nx

% Law of total variance: merge subset S
wS = sum(w(S));
muS = sum(w(S).*mu(S)) / wS;
% Total variance for subset-mixture:
% Var = E[Var(X|Z)] + Var(E[X|Z]) where Z is component index
var_within = sum(w(S).*(sig(S).^2)) / wS;
var_between = sum(w(S).*(mu(S)-muS).^2) / wS;
varS = var_within + var_between;
stdS = sqrt(varS);

% The subset-mixture density (normalized within subset) and merged Gaussian
f_subset_mix = ((w(S).'/wS) * f_comp(S,:)); % 1 x Nx
f_subset_merged = normpdf1(x, muS, stdS);

% Plot
figure(801); clf

h = setup_subplots(2,3);
isub = 1;

% Panel 1: full mixture + components

hca = h(isub); isub = isub + 1;
hold(hca,'on')
cols = lines(K);
for k = 1:K
  plot(hca,x, w(k)*f_comp(k,:), 'Color', cols(k,:), 'LineWidth', 1.2);
end
plot(hca,x, f_mix, 'k', 'LineWidth', 2.0);
hold(hca,'off')
grid(hca,'on')
xlabel(hca,'x')
ylabel(hca,'pdf')
title(hca,'Components (weighted) and full mixture')
legend(hca,[arrayfun(@(k) sprintf('w%g N(\\mu=%g,\\sigma=%g)',k,mu(k),sig(k)),1:K,'UniformOutput',false), ...
        {'mixture'}], 'Location','best')

% Panel 2: subset-mixture vs merged Gaussian (both normalized to area 1)
hca = h(isub); isub = isub + 1;
hold(hca,'on')
plot(hca,x, f_subset_mix, 'k', 'LineWidth', 2.0);
plot(hca,x, f_subset_merged, 'r--', 'LineWidth', 2.0);
hold(hca,'off')
grid(hca,'on')
xlabel(hca,'x')
ylabel(hca,'pdf')
title(hca,sprintf('Subset S=%s (normalized) vs merged Gaussian (moment-matched)', mat2str(S)))
legend(hca,{sprintf('subset-mixture (components %s)', mat2str(S)), ...
        sprintf('merged N(\\mu=%.2f, \\sigma=%.2f)', muS, sigS)}, 'Location','best')


hca = h(isub); isub = isub + 1;
imagesc(hca,d.m2)
hcb = colorbar(hca);
hcb.YLabel.String = 'Mahalanobis distance squared';

hca = h(isub); isub = isub + 1;
imagesc(hca,d.s)
hcb = colorbar(hca);
hcb.YLabel.String = 'Stein distance';


hca = h(isub); isub = isub + 1;
imagesc(hca,d.b)
hcb = colorbar(hca);
hcb.YLabel.String = 'Bhattacharyya distance';

% Print merged parameters (analog of ComponentProportion/mu/Sigma)
fprintf('Merged subset S=%s:\n', mat2str(S));
fprintf('  ComponentProportion (weight) wS = %.6f\n', wS);
fprintf('  muS = %.6f\n', muS);
fprintf('  SigmaS = %.6f (variance), sigmaS = %.6f\n', varS, sigS);