% Demo: law of total variance for merging Gaussian components.
%
% Goal:
% - Define several 1D Gaussian components with weights.
% - Plot individual components and the full mixture.
% - Pick a subset of components to "merge" into a single Gaussian by matching
%   0th/1st/2nd central moments (law of total variance).
% - Overlay the merged Gaussian against the subset-mixture it approximates.
%
% Works without any toolboxes (uses only exp/sqrt/pi).

clear; clc

%% Define components (1D)
% weights should sum to 1
w = [0.30 0.20 0.35 0.15];
mu = [-500  -50   400  900];     % km/s (or arbitrary units)
sig = [ 120   60   250  180];    % standard deviation (same units as mu)

%w = [0.5 0.5];
%mu = [-500  500];     % km/s (or arbitrary units)
%sig = [ 120   120];    % standard deviation (same units as mu)

w = [1 1 1 1]/3;
mu = [-500 0 900 950];     % km/s (or arbitrary units)
sig = [500 500 200 300];    % standard deviation (same units as mu)

K = numel(w);
assert(numel(mu)==K && numel(sig)==K, 'w/mu/sig sizes must match.');
w = w(:); mu = mu(:); sig = sig(:);
w = w./sum(w);

% Choose which components to merge (subset S)
S = [1 2 3 4];   % merge these into one equivalent Gaussian

% Bhattacharyya overlap diagnostics (1D)
% Bhattacharyya distance between N(mu1,s1^2) and N(mu2,s2^2):
%   DB = 1/4 * ( (mu1-mu2)^2 / (s1^2+s2^2) ) + 1/2 * ln( (s1^2+s2^2) / (2*s1*s2) )
%
% Smaller DB => more overlap. Use DB thresholding to suggest which components
% are "similar enough" to merge (including transitive overlap).
doDB = true;
DBthresh = 0.6;                  % try 0.2..1.5
DBthreshSweep = [0.2 0.4 0.6 0.8 1.2];

% Helpers
normpdf1 = @(x,m,s) exp(-0.5*((x-m)./s).^2)./(sqrt(2*pi).*s);

% Grid for plotting
xmin = min(mu - 5*sig);
xmax = max(mu + 5*sig);
x = linspace(xmin,xmax,2000);

% Compute component pdfs and mixture
f_comp = zeros(K,numel(x));
for k = 1:K
  f_comp(k,:) = normpdf1(x, mu(k), sig(k));
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
sigS = sqrt(varS);

% The subset-mixture density (normalized within subset) and merged Gaussian
f_subset_mix = ((w(S).'/wS) * f_comp(S,:)); % 1 x Nx
f_subset_merged = normpdf1(x, muS, sigS);

% Plot
figure(801); clf
tiledlayout(2,1,'TileSpacing','compact')

% Panel 1: full mixture + components
nexttile
hold on
cols = lines(K);
for k = 1:K
  plot(x, w(k)*f_comp(k,:), 'Color', cols(k,:), 'LineWidth', 1.2);
end
plot(x, f_mix, 'k', 'LineWidth', 2.0);
hold off
grid on
xlabel('x')
ylabel('pdf')
title('Components (weighted) and full mixture')
legend([arrayfun(@(k) sprintf('w%g N(\\mu=%g,\\sigma=%g)',k,mu(k),sig(k)),1:K,'UniformOutput',false), ...
        {'mixture'}], 'Location','best')

% Panel 2: subset-mixture vs merged Gaussian (both normalized to area 1)
nexttile
hold on
plot(x, f_subset_mix, 'k', 'LineWidth', 2.0);
plot(x, f_subset_merged, 'r--', 'LineWidth', 2.0);
hold off
grid on
xlabel('x')
ylabel('pdf')
title(sprintf('Subset S=%s (normalized) vs merged Gaussian (moment-matched)', mat2str(S)))
legend({sprintf('subset-mixture (components %s)', mat2str(S)), ...
        sprintf('merged N(\\mu=%.2f, \\sigma=%.2f)', muS, sigS)}, 'Location','best')

% Print merged parameters (analog of ComponentProportion/mu/Sigma)
fprintf('Merged subset S=%s:\n', mat2str(S));
fprintf('  ComponentProportion (weight) wS = %.6f\n', wS);
fprintf('  muS = %.6f\n', muS);
fprintf('  SigmaS = %.6f (variance), sigmaS = %.6f\n', varS, sigS);

%% Bhattacharyya distance diagnostics + suggested merges
if doDB
  DB = nan(K,K);
  for i = 1:K
    for j = i+1:K
      DB(i,j) = bhattacharyya_1d(mu(i), sig(i), mu(j), sig(j));
      DB(j,i) = DB(i,j);
    end
  end
  DB(1:K+1:end) = 0;

  figure(802); clf
  imagesc(DB);
  axis xy equal tight
  colorbar
  title('Bhattacharyya distance DB between components (1D)')
  xlabel('Component index')
  ylabel('Component index')

  groups = overlap_groups_from_distance(DB, DBthresh);
  fprintf('\nSuggested overlap groups (DB <= %.3f):\n', DBthresh);
  disp(groups(:));

  S_best = most_overlapping_group(DB, DBthresh);
  fprintf('Most-overlapping group (connected via DB <= %.3f): S_best = %s\n', DBthresh, mat2str(S_best));

  % Sweep thresholds to build intuition
  for ith = 1:numel(DBthreshSweep)
    th = DBthreshSweep(ith);
    gs = overlap_groups_from_distance(DB, th);
    fprintf('  DB <= %.3f -> groups: ', th);
    disp(gs);
  end
end

% Monte Carlo sanity check (optional)
doMC = true;
if doMC
  N = 200000;
  % Sample from subset-mixture (normalized within subset)
  wS_norm = w(S)./wS;
  comp = randsample(numel(S), N, true, wS_norm);
  xs = zeros(N,1);
  for ii = 1:numel(S)
    mask = comp==ii;
    k = S(ii);
    xs(mask) = mu(k) + sig(k)*randn(sum(mask),1);
  end
  mu_mc = mean(xs);
  var_mc = var(xs,1); % population variance
  fprintf('MC check (subset-mixture, N=%d): mu=%.6f, var=%.6f\n', N, mu_mc, var_mc);
end

%% Local helper functions (script-local)
function DB = bhattacharyya_1d(mu1, s1, mu2, s2)
  v1 = s1.^2;
  v2 = s2.^2;
  v = 0.5*(v1+v2);
  DB = 0.125*((mu1-mu2).^2)./v + 0.5*log(v./sqrt(v1.*v2));
end

function groups = overlap_groups_from_distance(D, thresh)
  % Connected components of graph where edge exists if D(i,j) <= thresh.
  K = size(D,1);
  Adj = (D <= thresh);
  Adj(1:K+1:end) = false;
  unvisited = true(1,K);
  groups = {};
  for i = 1:K
    if ~unvisited(i); continue; end
    stack = i;
    comp = [];
    unvisited(i) = false;
    while ~isempty(stack)
      v = stack(end); stack(end) = [];
      comp(end+1) = v; %#ok<AGROW>
      nbrs = find(Adj(v,:) & unvisited);
      unvisited(nbrs) = false;
      stack = [stack nbrs]; %#ok<AGROW>
    end
    groups{end+1} = sort(comp); %#ok<AGROW>
  end
end

function Sbest = most_overlapping_group(D, thresh)
  % Return the single "most overlapping" group under the threshold.
  % Definition:
  % - Build overlap groups as connected components of edges D(i,j) <= thresh.
  % - Pick the largest group (transitive overlap counts).
  % - Tie-breaker: smallest mean pairwise distance within the group.
  groups = overlap_groups_from_distance(D, thresh);
  if isempty(groups)
    Sbest = [];
    return
  end

  sizes = cellfun(@numel, groups);
  maxSize = max(sizes);
  cand = find(sizes == maxSize);
  if numel(cand) == 1
    Sbest = groups{cand};
    return
  end

  meanDB = nan(size(cand));
  for ic = 1:numel(cand)
    S = groups{cand(ic)};
    if numel(S) < 2
      meanDB(ic) = inf;
      continue
    end
    sub = D(S,S);
    mask = triu(true(numel(S)),1);
    vals = sub(mask);
    meanDB(ic) = mean(vals,'omitnan');
  end
  [~,iMin] = min(meanDB);
  Sbest = groups{cand(iMin)};
end

