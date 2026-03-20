function out = gmm_merge_components(gm, groups, varargin)
%GMM_MERGE_COMPONENTS Merge GMM components using moment matching (3D).
%
% out = gmm_merge_components(gm, groups, 'Name', value, ...)
%
% Inputs
%   gm     : either
%            - gmdistribution, or
%            - cell array of gmdistribution objects (nt x nK)
%   groups : either
%            - 1xNg cell array of vectors of component indices (for single gm)
%            - cell array matching size(gm), where each entry is 1xNg cell array
%              of vectors (output of overlap_groups_from_distance on a cell array)
%
% Each group is merged into one new Gaussian component by matching 0th/1st/2nd
% moments:
%   w   = sum_i w_i
%   mu  = (1/w) * sum_i w_i * mu_i
%   Sig = (1/w) * sum_i w_i * (Sig_i + (mu_i-mu)(mu_i-mu)')
%
% Options (name-value)
%   'isort' : permutation mapping sorted->unsorted indices. Use this if your
%             groups were computed in sorted index space (e.g. after calling
%             bhattacharyya_distance(...,'sort',true)). Can be:
%             - numeric 1xK (for single gm), or
%             - cell array matching size(gm) with numeric vectors per entry.
%
% Output (struct)
%   out.gmMerged         : merged gmdistribution (or cell array)
%   out.mu               : merged mu (M x 3) (or cell array)
%   out.Sigma            : merged Sigma (3 x 3 x M) (or cell array)
%   out.ComponentProportion : merged weights (1 x M) (or cell array)
%   out.groups_unsorted  : groups after mapping to original gm indexing
%   out.groupWeights     : total weight per merged component (or cell array)
%

isortOpt = [];
if ~isempty(varargin)
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'isort'
        isortOpt = varargin{iarg+1};
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end

% Normalize inputs to cell arrays
if isa(gm, 'gmdistribution')
  gm_cell = {gm};
  groups_cell = {groups};
  if isempty(isortOpt)
    isort_cell = {[]};
  else
    isort_cell = {isortOpt};
  end
elseif iscell(gm)
  gm_cell = gm;
  groups_cell = groups;
  if isempty(isortOpt)
    isort_cell = cell(size(gm_cell));
  else
    isort_cell = isortOpt;
  end
else
  error('gm must be a gmdistribution or cell array.');
end

if ~iscell(groups_cell)
  error('groups must be a cell array (or cell array of cell arrays).');
end

[nt, nK] = size(gm_cell);

gmMerged = cell(nt,nK);
mu_out = cell(nt,nK);
Sigma_out = cell(nt,nK);
w_out = cell(nt,nK);
groups_unsorted_out = cell(nt,nK);
groupWeights_out = cell(nt,nK);

for it = 1:nt
  for iK = 1:nK
    g = gm_cell{it,iK};
    if isempty(g)
      continue
    end
    if ~isa(g, 'gmdistribution')
      error('gm{%d,%d} is not a gmdistribution.', it, iK);
    end

    gs = groups_cell{it,iK};
    if isempty(gs)
      % default: no merging, keep as-is
      gs = arrayfun(@(k) k, 1:g.NumComponents, 'UniformOutput', false);
    end
    if ~iscell(gs)
      error('groups{%d,%d} must be a cell array of index vectors.', it, iK);
    end

    isort = [];
    if ~isempty(isort_cell)
      isort = isort_cell{it,iK};
    end
    if ~isempty(isort)
      % Map sorted indices -> original (unsorted) indices
      gs_unsorted = cellfun(@(S) isort(S), gs, 'UniformOutput', false);
    else
      gs_unsorted = gs;
    end
    groups_unsorted_out{it,iK} = gs_unsorted;

    % Extract original parameters
    mu0 = g.mu;                 % Kx3
    Sigma0 = g.Sigma;           % 3x3xK
    w0 = g.ComponentProportion; % 1xK
    K0 = g.NumComponents;
    if numel(w0) ~= K0
      error('Unexpected ComponentProportion size at it=%d iK=%d.', it, iK);
    end

    M = numel(gs_unsorted);
    muM = zeros(M,3);
    SigmaM = zeros(3,3,M);
    wM = zeros(1,M);

    for m = 1:M
      S = gs_unsorted{m};
      S = unique(S(:))';
      if any(S < 1) || any(S > K0)
        error('Group index out of range at it=%d iK=%d.', it, iK);
      end
      wS = sum(w0(S));
      if wS <= 0
        % Degenerate: keep placeholder with tiny weight
        wS = 0;
        muS = [0 0 0];
        SigS = eye(3);
      else
        % Weighted mean
        muS = (w0(S) * mu0(S,:)) / wS; % 1x3
        % Total covariance (within + between)
        SigS = zeros(3,3);
        for ii = 1:numel(S)
          k = S(ii);
          dmu = (mu0(k,:) - muS)';
          SigS = SigS + w0(k) * (Sigma0(:,:,k) + (dmu*dmu')); % law of shared covariances
        end
        SigS = SigS / wS;
      end

      wM(m) = wS;
      muM(m,:) = muS;
      SigmaM(:,:,m) = SigS;
    end

    % Error if ComponentProportion == 0 so put it to a very small value
    wM(wM==0) = 1e-12;

    % Renormalize weights to sum to 1 for gmdistribution
    if sum(wM) > 0
      wMnorm = wM./sum(wM);
    else
      wMnorm = ones(1,M)./M;
    end
    
    mu_out{it,iK} = muM;
    Sigma_out{it,iK} = SigmaM;
    w_out{it,iK} = wMnorm;
    groupWeights_out{it,iK} = wM;
    gmMerged{it,iK} = gmdistribution(muM, SigmaM, wMnorm);
  end
end

out = struct();
out.gmMerged = gmMerged;
out.mu = mu_out;
out.Sigma = Sigma_out;
out.ComponentProportion = w_out;
out.groups_unsorted = groups_unsorted_out;
out.groupWeights = groupWeights_out;

% Collapse convenience outputs for single gm input
if isa(gm, 'gmdistribution')
  out.gmMerged = out.gmMerged{1,1};
  out.mu = out.mu{1,1};
  out.Sigma = out.Sigma{1,1};
  out.ComponentProportion = out.ComponentProportion{1,1};
  out.groups_unsorted = out.groups_unsorted{1,1};
  out.groupWeights = out.groupWeights{1,1};
end

end

