% Illustrate overlap diagnostics between GMM components.
%
% Prereqs in workspace:
%   times : EpochTT or similar time array (nt x 1)
%   gm    : cell array (nt x nK) or (nt x 1) of gmdistribution objects
%   iK    : which column of gm to use if gm is (nt x nK)
%
% The overlap metric uses a symmetric Mahalanobis distance between means:
%   d_ij^2 = (mu_i-mu_j)' * inv((Sigma_i+Sigma_j)/2) * (mu_i-mu_j)
% Small d_ij => strong overlap.
%
% This script also computes the Bhattacharyya distance between Gaussians:
%   DB = 1/8 * dmu' * inv(S) * dmu + 1/2 * ln( det(S) / sqrt(det(Si)*det(Sj)) )
% where S = (Si+Sj)/2. This captures "embedded" narrow-in-broad overlap better
% than mean-separation alone.

%% Select GMM time series
if ~exist('iK','var'); iK = 1; end
if iscell(gm) && size(gm,2) > 1
  gm_ts = gm(:,iK);
else
  gm_ts = gm(:);
end

nt = numel(gm_ts);
K = gm_ts{1}.NumComponents;

%% Settings
doSort = true;           % sort components by temperature proxy before comparing
dOverlap = 2.0;          % threshold for "overlapping-ish" (for dij)
dBOverlap = 0.5;         % threshold for "overlapping-ish" (for DB) - tune
itPlot = max(1, round(nt*0.3)); % which time index to show heatmap for
itPlot = 80;
metricForGrouping = 'DB'; % 'dij' or 'DB'

%% Compute pairwise distances per time
D = nan(K,K,nt);
dmin = nan(nt,1);
DB = nan(K,K,nt);
DBmin = nan(nt,1);

for it = 1:nt
  g = gm_ts{it};
  if doSort
    isort = gmm_sort(g);
  else
    isort = 1:K;
  end
  mu = g.mu(isort,:);          % Kx3
  Sigma = g.Sigma(:,:,isort);  % 3x3xK

  for i = 1:K
    for j = i+1:K
      dmu = (mu(i,:)-mu(j,:))';
      S = 0.5*(Sigma(:,:,i) + Sigma(:,:,j));
      d2 = dmu'*(S\dmu);
      dij = sqrt(max(d2,0));
      D(i,j,it) = dij;
      D(j,i,it) = dij;

      % Bhattacharyya distance
      % logdet via chol for numerical stability (S should be SPD)
      Si = Sigma(:,:,i);
      Sj = Sigma(:,:,j);
      % Regularize lightly if needed
      reg = 0;
      for tries = 1:3
        [Ri,pS] = chol(S + reg*eye(3));
        [R1,p1] = chol(Si + reg*eye(3));
        [R2,p2] = chol(Sj + reg*eye(3));
        if pS==0 && p1==0 && p2==0
          break
        end
        reg = max(1e-10, 10^(tries-10));
      end
      logdetS = 2*sum(log(diag(Ri)));
      logdetSi = 2*sum(log(diag(R1)));
      logdetSj = 2*sum(log(diag(R2)));
      term1 = 0.125 * dmu'*(S\dmu);
      term2 = 0.5 * (logdetS - 0.5*(logdetSi + logdetSj));
      db = term1 + term2;
      DB(i,j,it) = db;
      DB(j,i,it) = db;
    end
  end
  tmp = D(:,:,it);
  tmp(1:K+1:end) = NaN;
  dmin(it) = min(tmp,[],'all','omitnan');

  tmp = DB(:,:,it);
  tmp(1:K+1:end) = NaN;
  DBmin(it) = min(tmp,[],'all','omitnan');
end

%% Plot: heatmap at one time
figure(701); clf
imagesc(D(:,:,itPlot))
axis xy equal tight
colorbar
caxis([0 5])
title(sprintf('Pairwise d_{ij} at it=%d (doSort=%d)', itPlot, doSort))
xlabel('Component index')
ylabel('Component index')

figure(703); clf
imagesc(DB(:,:,itPlot))
axis xy equal tight
colorbar
title(sprintf('Pairwise Bhattacharyya DB at it=%d (doSort=%d)', itPlot, doSort))
xlabel('Component index')
ylabel('Component index')

%% Plot: closest-pair distance vs time
figure(702); clf
plot(dmin,'k','LineWidth',1.5); hold on
yline(dOverlap,'k--','LineWidth',1.2);
hold off
grid on
xlabel('it')
ylabel('min_{i<j} d_{ij}')
title('Closest component-pair separation vs time')

figure(704); clf
plot(DBmin,'k','LineWidth',1.5); hold on
yline(dBOverlap,'k--','LineWidth',1.2);
hold off
grid on
xlabel('it')
ylabel('min_{i<j} DB_{ij}')
title('Closest component-pair Bhattacharyya distance vs time')

%% Find overlap groups (connected components) for one time index
itGroup = itPlot;
switch lower(metricForGrouping)
  case 'dij'
    Adj = D(:,:,itGroup) < dOverlap;
  case 'db'
    Adj = DB(:,:,itGroup) < dBOverlap;
  otherwise
    error('Unknown metricForGrouping: %s', metricForGrouping);
end
Adj(1:K+1:end) = false;

unvisited = true(1,K);
groups = {};
for i = 1:K
  if ~unvisited(i); continue; end
  % BFS/DFS
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

fprintf('Overlap groups at it=%d with metric=%s (doSort=%d):\n', itGroup, metricForGrouping, doSort);
if strcmpi(metricForGrouping,'dij')
  fprintf('  threshold dOverlap = %.2f\n', dOverlap);
else
  fprintf('  threshold dBOverlap = %.2f\n', dBOverlap);
end
disp(groups(:));

