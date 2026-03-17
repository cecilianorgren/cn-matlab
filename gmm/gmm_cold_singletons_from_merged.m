function out = gmm_cold_singletons_from_merged(merged, coldMerged)
%GMM_COLD_SINGLETONS_FROM_MERGED Map cold components back to original singletons.
%
% out = gmm_cold_singletons_from_merged(merged, coldMerged)
%
% Inputs
%   merged     : struct output from gmm_merge_components (cell arrays)
%                must contain merged.groups_unsorted{it,1} where each entry is
%                a 1xKM cell array of original component index vectors.
%   coldMerged : struct output from gmm_identify_cold_components run on merged.gmMerged
%                must contain coldMerged.cold_unsorted{it} in merged-component indexing.
%
% Output (struct)
%   out.coldMerged_kept{it} : merged-component indices that are cold AND singleton groups
%   out.coldOriginal{it}    : corresponding original component indices (singletons)
%   out.isColdOriginal_mat  : (nt x Korig) logical mask in original component indexing
%

if ~isstruct(merged) || ~isfield(merged,'groups_unsorted')
  error('merged must be output from gmm_merge_components (needs groups_unsorted).');
end
if ~isstruct(coldMerged) || ~isfield(coldMerged,'cold_unsorted')
  error('coldMerged must be output from gmm_identify_cold_components (needs cold_unsorted).');
end

nt = size(merged.groups_unsorted,1);

% infer Korig from first non-empty group list
Korig = [];
for it = 1:nt
  if isempty(merged.groups_unsorted{it,1}); continue; end
  allIdx = [merged.groups_unsorted{it,1}{:}];
  if ~isempty(allIdx)
    Korig = max(allIdx);
    break
  end
end
if isempty(Korig)
  Korig = 0;
end

coldMerged_kept = cell(nt,1);
coldOriginal = cell(nt,1);
isColdOriginal_mat = false(nt, Korig);

for it = 1:nt
  if isempty(merged.groups_unsorted{it,1})
    coldMerged_kept{it} = [];
    coldOriginal{it} = [];
    continue
  end

  coldIdxMerged = coldMerged.cold_unsorted{it};
  kept = [];
  orig = [];
  for ii = 1:numel(coldIdxMerged)
    im = coldIdxMerged(ii);
    grp = merged.groups_unsorted{it,1}{im};
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
out.coldMerged_kept = coldMerged_kept;
out.coldOriginal = coldOriginal;
out.isColdOriginal_mat = isColdOriginal_mat;

end

