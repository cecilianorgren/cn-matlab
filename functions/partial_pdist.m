function varargout = partial_pdist(pdist,MP,idx)
% PARTIAL_PDIST Separates PDist based on group index.

nGroups = numel(unique(idx));
pdist_group = cell(1,nGroups);
Dep1_edges = 0.5:1:(size(pdist.depend{1},2)+0.5);
Dep2_edges = 0.5:1:(size(pdist.depend{2},2)+0.5);
Dep3_edges = 0.5:1:(size(pdist.depend{3},2)+0.5);

for iGroup = 1:nGroups  
  % Bin the macroparticles into their original PDist bins
  [average_groupid,~,mid,~] = histcn([MP.iDep1, MP.iDep2, MP.iDep3], Dep1_edges, Dep2_edges, Dep3_edges, 'AccumData', idx, 'Fun', @mean);
  
  % Find mean group id for each bin
  average_groupid = round(average_groupid); 
  
  % Copy PDist and then put all bins that are not in this group to zero.
  pdist_tmp = pdist;
  id_exclude = find(average_groupid ~= iGroup);  
  pdist_tmp.data(:,id_exclude) = 0;  

  % Collect for output
  pdist_group{iGroup} = pdist_tmp;
end

varargout{1} = pdist_group;
