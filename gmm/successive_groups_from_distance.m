function groups = successive_groups_from_distance(D)


if isnumeric(D)
  groups = local_groups_one(D);
elseif iscell(D)
  groups = cell(size(D));
  for it = 1:size(D,1)
    for iK = 1:size(D,2)
      if isempty(D{it,iK})
        groups{it,iK} = {};
      else
        groups{it,iK} = local_groups_one(D{it,iK});
      end
    end
  end
else
  error('D must be numeric or cell array.');
end
end


function groups = local_groups_one(D)

D_unique = unique(D);
threshold_vec = D_unique(1:end-1) + 0.5*diff(D_unique);

nThresh = numel(threshold_vec);
%n_groups = zeros(1,nThresh);
groups{1} = overlap_groups_from_distance(D, 0);
for iThresh = 2:nThresh 
  new_grouping = overlap_groups_from_distance(D, threshold_vec(iThresh));
  if ~isequal(groups{end},new_grouping)
    groups{end+1} = new_grouping;
  end  
end
end
