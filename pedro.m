nInd = 20;
a = round(rand(nInd,1));
b = round(rand(nInd,1));

ab = [a b];

a1 = find(a==1);
b0 = find(b==0);

a1b0 = intersect(a1,b0);

c = ab(a1b0,:);
%%
[a b]

