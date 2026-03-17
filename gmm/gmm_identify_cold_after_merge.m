function out = gmm_identify_cold_after_merge(gm, steinThreshold, varargin)
%GMM_IDENTIFY_COLD_AFTER_MERGE Merge similar components, then identify cold ones.
%
% out = gmm_identify_cold_after_merge(gm, steinThreshold, 'Name', value, ...)
%
% Workflow per time:
%   1) Compute Stein distance between components (optionally sorted by gmm_sort)
%   2) Group components by connected components under steinThreshold
%   3) Merge each group by moment matching (law of total variance)
%   4) Identify cold components on the merged GMM (gmm_identify_cold_components)
%   5) Keep only those cold components whose merge-group size is 1 (unmerged)
%
% Inputs
%   gm            : cell array (nt x 1) of gmdistribution objects
%   steinThreshold: scalar threshold passed to overlap_groups_from_distance
%
% Options (name-value)
%   'steinSort' : true/false (default true). If true, stein_distance uses gmm_sort.
%   'coldArgs'  : cell array of name-value pairs forwarded to gmm_identify_cold_components
%                 (default {}), e.g. {'minGapLog10',0.3,'maxColdK',2,'minColdK',0}
%
% Output (struct)
%   out.merged            : output of gmm_merge_components (contains gmMerged, groups_unsorted, ...)
%   out.coldMerged        : output of gmm_identify_cold_components on merged.gmMerged
%   out.coldMerged_kept{it}   : cold merged-component indices (in merged gm indexing) that were kept
%   out.coldOriginal{it}      : corresponding original component indices (unmerged singletons)
%   out.isColdOriginal_mat    : (nt x Korig) logical mask in original component indexing
%

steinSort = true;
coldArgs = {};

if ~isempty(varargin)
  for iarg = 1:2:numel(varargin)
    switch lower(varargin{iarg})
      case 'steinsort'
        steinSort = varargin{iarg+1};
      case 'coldargs'
        coldArgs = varargin{iarg+1};
      otherwise
        error('Unknown option: %s', varargin{iarg});
    end
  end
end

if ~iscell(gm)
  error('gm must be a cell array (nt x 1) of gmdistribution objects.');
end
if size(gm,2) ~= 1
  error('gm must be size (nt x 1). Pass gm(:,iK) if needed.');
end

nt = size(gm,1);
Korig = gm{1}.NumComponents;

% 1-3) merge
sd = stein_distance(gm, 'sort', steinSort);
groups_sorted = overlap_groups_from_distance(sd.D, steinThreshold);
merged = gmm_merge_components(gm, groups_sorted, 'isort', sd.isort);

% 4) cold identification on merged GMM
coldMerged = gmm_identify_cold_components(merged.gmMerged, coldArgs{:});

% 5) keep only cold merged components that were not merged (singleton groups)
coldMerged_kept = cell(nt,1);
coldOriginal = cell(nt,1);
isColdOriginal_mat = false(nt, Korig);

for it = 1:nt
  if isempty(merged.gmMerged{it}); continue; end

  % coldMerged.cold_unsorted are indices in merged component indexing
  coldIdxMerged = coldMerged.cold_unsorted{it};
  if isempty(coldIdxMerged)
    coldMerged_kept{it} = [];
    coldOriginal{it} = [];
    continue
  end

  kept = [];
  orig = [];
  for ii = 1:numel(coldIdxMerged)
    im = coldIdxMerged(ii);
    grp = merged.groups_unsorted{it}{im};
    if numel(grp) == 1
      kept(end+1) = im; %#ok<AGROW>
      orig(end+1) = grp; %#ok<AGROW>
    end
  end

  coldMerged_kept{it} = kept;
  coldOriginal{it} = orig;
  if ~isempty(orig)
    isColdOriginal_mat(it, orig) = true;
  end
end

out = struct();
out.merged = merged;
out.coldMerged = coldMerged;
out.coldMerged_kept = coldMerged_kept;
out.coldOriginal = coldOriginal;
out.isColdOriginal_mat = isColdOriginal_mat;

end

