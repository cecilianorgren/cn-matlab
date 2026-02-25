function tint = bm_tint
[t1,t2,tint_comments]=textread('/Users/Cecilia/Research/BM.txt','%s%s%[^\n]');
for j=1:size(t1,1),
    tint(j,1)=iso2epoch(t1{j});
    tint(j,2)=iso2epoch(t2{j});
end