file = strsplit(printpath,'/'); file = strrep([strjoin(file(1:5),filesep) filesep 'Research/GMM/GMM_BFFs/db_gmm.mat'],'''','"');

DB.metadata.reduction = [];
DB.metadata.list_events = [];
for K = 1:5
  [partitions,unique_groups,map_part2group] = allPartitions(K);
  DB.metadata.reduction(K).numPopulations = K;
  for iGroup = 1:numel(unique_groups)
    DB.metadata.reduction(K).groups(iGroup).id = iGroup;
    DB.metadata.reduction(K).groups(iGroup).members = unique_groups(iGroup);
  end
  for iPartition = 1:numel(partitions)
    DB.metadata.reduction(K).partitions(iPartition).id = iPartition;
    DB.metadata.reduction(K).partitions(iPartition).members = map_part2group(iPartition);
  end
  %DB.metadata_reduction(K).partitions = partitions;
  %DB.metadata_reduction(K).map_partitions_to_groups = map_part2group;
end
save(file,'DB')
%%
DB = load(file); DB = DB.DB;
id_event = iDF;
id_run = 'run_00'; 

nMacroParticles = size(R,1);
DB = gmm_save_to_db(DB, id_event, id_run, 'gmm',gm, nMacroParticles);


save(file,'DB')
% DB
% ├── metadata.reduction
% │   └── K
% │       ├── groups
% │       │   ├── id
% │       │   └── members
% │       ├── partitions
% │           ├── id
% │           ├── members % sufficient with group_ids here
% │           └── % map_partitions_to_groups
% └── events
%     ├── event_001
%     │   ├── metadata
%     │   │   ├── event_id
%     │   │   ├── time_start
%     │   │   ├── time_stop
%     │   │   ├── 
%     │   │   ├── 
%     │   │   ├── 
%     │   │   ├── 
%     │   │   └── 
%     │   ├── data
%     │   └── runs
%     │       ├── run_001
%     │       │   ├── settings
%     │       │   ├── GMM
%     │       │   ├── reduction
%     │       │   └── cold
%     │       └── run_002
%     └── event_002