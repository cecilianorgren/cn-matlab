function id = sql_parse_flag(flag)

flag_ = cellfun(@(x) [x ','],flag,'UniformOutput', false);
uniqueFlags = unique(strsplit([flag_{:}],','));


allFlags = cellfun(@(x) strsplit(x,','),flag,'UniformOutput', false);
nFlags = numel(allFlags{1});
nEvents = numel(allFlags);
%aa=[flag_{:}]

indAllFlags = zeros(nEvents,nFlags);
labAllFlags = cell(nFlags,1);
for iEvent = 1:nEvents
  
  %tmpFlags = allFlags{iEvent};
  labAllFlags;
end
%%



aa = cellfun(@(x) strsplit(x,','),uniqueFlags,'UniformOutput', false);


'split'

aa = cellfun(@(x) unique(x),flag,'UniformOutput', false);
id = find(~cellfun(@isempty,aa));