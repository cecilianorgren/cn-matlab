function id = sql_find_event(eventId,findId)

aa = cellfun(@(x) strfind(x,'20151112_071854'),eventId,'UniformOutput', false);
id = find(~cellfun(@isempty,aa));