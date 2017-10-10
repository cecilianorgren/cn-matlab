% bm_make_plots

[t1,t2,tint_comments]=textread('/Users/Cecilia/Exjobb/BM/BM.txt','%s%s%[^\n]');
for j=1:size(t1,1)
    j
    tint{j}=[iso2epoch(t1{j}) iso2epoch(t1{j})];
    datestring{j}=datestr(fromepoch(tint{j}(1)),'yyyymmdd');
end
%clear t1 t2 tint_comments j