%rng(1)
N = 1e4;
tmp_r = rand(N,1);

vmin = 0; vmax = 1;
edges = linspace(vmin,vmax,11);
centers = edges(1:end-1) + 0.5*diff(edges);


v1 = (vmin^1 + tmp_r*(vmax^1-vmin^1)).^(1/1);
v2 = (vmin^2 + tmp_r*(vmax^2-vmin^2)).^(1/2);
v3 = (vmin^3 + tmp_r*(vmax^3-vmin^3)).^(1/3);


[N1,edges] = histcounts(v1,edges);
[N2,edges] = histcounts(v2,edges);
[N3,edges] = histcounts(v3,edges);

bar(centers,[N1' N2' N3'])
