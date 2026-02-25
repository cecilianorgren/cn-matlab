x = 0:100;
edges = 0:10:100;
[N,bin] = histc(x,edges);
bar(edges,N,'histc')

