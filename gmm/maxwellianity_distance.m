function out = maxwellianity_distance(gm)
%
% Merge all permutations (or component-wise) of populations, see what under
% what combination they remain more or less maxwellian.
%

% out = stein_distance(gm, 'Name', value, ...)
%
% Inputs
%   gm : either
%        - gmdistribution object, or
%        - cell array of gmdistribution objects (nt x nK) or (nt x 1)
%
% Options (name-value)
%   'sort'        : true/false (default false). If true, sort components by
%                   temperature proxy using gmm_sort() before computing distance.
%
% Output (struct)
%   out.D{it,iK}     : KxK Stein divergence matrix (0 on diagonal)
%   out.isort{it,iK} : component sort order used (sorted->unsorted)
%   out.sort         : sort flag
%   out.minD(it,iK)  : minimum off-diagonal distance
%   out.meanD(it,iK) : mean off-diagonal distance
%
% Definition for maxwellianity:
%   ---
%
% Notes
% - This penalizes covariance/temperature mismatch strongly and ignores mean differences.
% - Use with overlap_groups_from_distance(out.D, threshold).

% Defaults
doSort = false;
reg0 = 1e-10;
maxRegTries = 6;

if ~isempty(varargin)
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'sort'
        doSort = varargin{iarg+1};
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end

% Normalize input to cell array gm_cell (nt x nK)
if isa(gm, 'gmdistribution')
  gm_cell = {gm};
elseif iscell(gm)
  gm_cell = gm;
else
  error('gm must be a gmdistribution or a cell array of gmdistribution objects.');
end

[nt, nK] = size(gm_cell);
out = struct();
out.sort = doSort;
out.D = cell(nt, nK);
out.isort = cell(nt, nK);
out.minD = nan(nt, nK);
out.meanD = nan(nt, nK);

for it = 1:nt
  for iK = 1:nK
    g = gm_cell{it, iK};
    if isempty(g)
      continue
    end
    if ~isa(g, 'gmdistribution')
      error('gm{%d,%d} is not a gmdistribution.', it, iK);
    end
    K = g.NumComponents;

    if doSort
      isort = gmm_sort(g);
    else
      isort = 1:K;
    end
    %isort = reshape(isort,permute(size(isort),[2 1]));
    out.isort{it,iK} = isort;

    %Sigma = g.Sigma(:,:,isort); % 3x3xK

    D = nan(K,K);
    D(1:K+1:end) = 0;
    for i = 1:K
      Si = Sigma(:,:,i);
      for j = i+1:K
        group = {i,j};
        g_merged = gmm_merge_components()
        Sj = Sigma(:,:,j);
        






        D(i,j) = dij;
        D(j,i) = dij;
      end
    end

    out.D{it,iK} = D;
    tmp = D;
    tmp(1:K+1:end) = NaN;
    out.minD(it,iK) = min(tmp,[],'all','omitnan');
    out.meanD(it,iK) = mean(tmp(~isnan(tmp)),'omitnan');
  end
end

% If input was a single gmdistribution, collapse outputs for convenience
if isa(gm, 'gmdistribution')
  out.D = out.D{1,1};
  out.isort = out.isort{1,1};
  out.minD = out.minD(1,1);
  out.meanD = out.meanD(1,1);
end

end

