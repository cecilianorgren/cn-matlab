function groups = ngroups_from_threshold(D,threshold_vec)

if isnumeric(D)
  groups = local_groups_one(D, threshold);
elseif iscell(D)
  groups = cell(size(D));
  for it = 1:size(D,1)
    for iK = 1:size(D,2)
      if isempty(D{it,iK})
        groups{it,iK} = {};
      else
        groups{it,iK} = local_groups_one(D{it,iK}, threshold_vec);
      end
    end
  end
else
  error('D must be numeric or cell array.');
end

end

function n_groups = local_groups_one(D, threshold)

nThresh = numel(threshold);
n_groups = zeros(1,nThresh);
for iThresh = 1:nThresh 
  groups = overlap_groups_from_distance(D, threshold(iThresh));
  n_groups(iThresh) = numel(groups);
end
end
