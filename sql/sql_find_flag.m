function id = sql_find_flag(flag)
1;
uniqueFlags = unique(flag(:));

str = uniqueFlags{1};
pat = '\s,';
s = regexp(str, pat, 'split')
             

'split'

aa = cellfun(@(x) unique(x),flag,'UniformOutput', false);
id = find(~cellfun(@isempty,aa));