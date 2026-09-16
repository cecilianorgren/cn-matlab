function out = gmm_choose_grouping_from_threshold(groups,quality,threshold)
% gmm_choose_grouping_from_threshold Get the most merged gm that satisfy
% the threshold for the given quality measure.
%
%
% new_gm = gmm_choose_grouping_from_threshold(gm,quality,threshold)

1;
%gm_new = 

[nt, nK,~] = size(quality);
iG = nan(nt,nK);
groups_selected = cell(nt,nK);
quality_selected = nan(nt,nK);


for it = 1:nt
  for iK = 1:nK
    tmp = find(quality(it,iK,:)<threshold,1,'last');
    if ~isempty(tmp)
      iG(it,iK) = tmp;
      groups_selected(it,iK) = {groups{it,iK}{tmp}};
      quality_selected(it,iK) = quality(it,iK,tmp);
    end
  end
end

out.groups = groups_selected;
out.nGroups = iG;
out.quality = quality_selected;