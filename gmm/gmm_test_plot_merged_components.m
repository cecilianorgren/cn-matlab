% Plot original vs merged GMM components for one (it,iK).
%
% Prereqs in workspace:
%   gm : cell array of gmdistribution from fitgmdist, size (nt x nK)
%
% This script:
% - picks g = gm{it,iK}
% - computes Bhattacharyya distance between components
% - makes overlap groups with a chosen threshold (transitive)
% - merges groups by moment matching (law of total variance)
% - plots 1D marginals (e.g. along vz) for original + merged

%% Choose which distribution
it = 41;
iK = 3;
g = gm{it,iK};

% Settings
doSort = true;          % sort by temperature proxy before grouping/merging
DBthresh = 0.4;         % overlap threshold (tune)
dim = 'z';              % 'x' | 'y' | 'z' for 1D marginal plot
nSigmaPlot = 6;         % plot range in +/- nSigmaPlot around means
nGrid = 2000;

% Compute overlap groups and merge
bd = bhattacharyya_distance(g, 'sort', doSort);
sd = stein_distance(gm{it,iK}, 'sort', true);
%groups = overlap_groups_from_distance(bd.DB, DBthresh);
groups = overlap_groups_from_distance(sd.D, DBthresh);
merged = gmm_merge_components(g, groups, 'isort', bd.isort);
gM = merged.gmMerged;

% Extract 1D marginal parameters
switch lower(dim)
  case 'x'
    idim = 1;
  case 'y'
    idim = 2;
  case 'z'
    idim = 3;
  otherwise
    error('Unknown dim: %s', dim);
end

% Original (possibly sorted) for plotting
K = g.NumComponents;
if doSort
  isort = bd.isort;
else
  isort = 1:K;
end
mu0 = g.mu(isort,idim);
sig0 = squeeze(sqrt(g.Sigma(idim,idim,isort)));
w0 = g.ComponentProportion(isort);

% Merged
KM = gM.NumComponents;
mu1 = gM.mu(:,idim);
sig1 = squeeze(sqrt(gM.Sigma(idim,idim,:)));
w1 = gM.ComponentProportion(:).';

% Build x-grid covering both sets
xmin0 = min(mu0 - nSigmaPlot*sig0);
xmax0 = max(mu0 + nSigmaPlot*sig0);
xmin1 = min(mu1 - nSigmaPlot*sig1);
xmax1 = max(mu1 + nSigmaPlot*sig1);
xmin = min([xmin0 xmin1]);
xmax = max([xmax0 xmax1]);
xmin = -3000;
xmax = 3000;
x = linspace(xmin, xmax, nGrid);

normpdf1 = @(x,m,s) exp(-0.5*((x-m)./s).^2)./(sqrt(2*pi).*s);

f0 = zeros(K, numel(x));
for k = 1:K
  f0(k,:) = normpdf1(x, mu0(k), sig0(k));
end
f0mix = (w0(:).')*f0;

f1 = zeros(KM, numel(x));
for k = 1:KM
  f1(k,:) = normpdf1(x, mu1(k), sig1(k));
end
f1mix = (w1(:).')*f1;

% Plot
figure(811); clf
tiledlayout(2,1,'TileSpacing','compact')

nexttile
hold on
cols = lines(max(K,KM));
for k = 1:K
  plot(x, w0(k)*f0(k,:), 'Color', cols(k,:), 'LineWidth', 1.2);
end
plot(x, f0mix, 'k', 'LineWidth', 2.0);
plot(vdf_fz.depend{1}(it,:),vdf_fz.data(it,:)*1e-3*4,'k--')
hold off
grid on
xlabel(sprintf('v_%s', dim))
ylabel('pdf')
title(sprintf('Original components (weighted) + mixture, it=%d iK=%d, DBthresh=%.2f', it, iK, DBthresh))
legs = arrayfun(@(x) sprintf('%g',x),1:K,'UniformOutput',false);
irf_legend(gca,legs',[0.98 0.98])

nexttile
hold on
for k = 1:KM
  plot(x, w1(k)*f1(k,:), 'Color', cols(k,:), 'LineWidth', 1.2);
end
plot(x, f1mix, 'k', 'LineWidth', 2.0);
plot(vdf_fz.depend{1}(it,:),vdf_fz.data(it,:)*1e-3*4,'k--')
hold off
grid on
xlabel(sprintf('v_%s', dim))
ylabel('pdf')
title(sprintf('Merged groups (weighted) + mixture (KM=%d)', KM))
legs = cellfun(@(x) sprintf('%g',x),groups,'UniformOutput',false);
irf_legend(gca,legs',[0.98 0.98])

%% Print groups
fprintf('Overlap groups (indices in %s space) with DBthresh=%.3f:\n', ternary(doSort,'sorted','original'), DBthresh);
disp(groups(:));
fprintf('Mapped to original component indices used for merging:\n');
disp(merged.groups_unsorted(:));

%% small helper
function out = ternary(cond, a, b)
if cond; out = a; else; out = b; end
end

