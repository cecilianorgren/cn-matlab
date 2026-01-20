function out = find_flag(data,flag)

nEvents = numel(data);

idx = zeros(nEvents,1);
for iEvent = 1:nEvents
  str = data{iEvent};
  data_tmp = strsplit(data{iEvent},',');
  encounters = cellfun(@(x)strcmp(x,flag),data_tmp);
  if any(encounters==1)
    idx(iEvent) = 1;
  end
end

out = idx;