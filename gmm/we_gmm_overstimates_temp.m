ffun = @(x,t) exp(-x.^2/t^2);

x_edges = linspace(0,1,11);
%x_edges = logspace(-2,0,11);
x_centers = x_edges(1:end-1) + 0.5*diff(x_edges);
t = 0.5;
f = ffun(x_centers-0.0,t);

Ntot = 1e4;
x1 = generate_mp_from_vf(x_edges,f,1,Ntot);
x2 = generate_mp_from_vf(x_edges,f,2,Ntot);
x3 = generate_mp_from_vf(x_edges,f,3,Ntot);

x_edges_up = linspace(min(x_edges),max(x_edges),100);
%x_edges_up = x_edges;
x_centers_up = x_edges_up(1:end-1) + 0.5*diff(x_edges_up);
N1 = histcounts(cat(1,x1{:}),x_edges_up);
N2 = histcounts(cat(1,x2{:}),x_edges_up);
N3 = histcounts(cat(1,x3{:}),x_edges_up);

m1 = mean(cat(1,x1{:}));
m2 = mean(cat(1,x2{:}));
m3 = mean(cat(1,x3{:}));

m1_2 = mean(cat(1,x1{:}).^2);
m2_2 = mean(cat(1,x2{:}).^2);
m3_2 = mean(cat(1,x3{:}).^2);

h = setup_subplots(4,1,'vertical');
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,x_centers,f,'*-')
hca.XTick = x_edges;
hca.XGrid = 'on';

hca = h(isub); isub = isub + 1;
bar(hca,x_centers_up,N1)
hca.XTick = x_edges;
irf_legend(hca,sprintf('mean(x) = %g, mean(x^2) = %g',m1,m1_2),[0.98 0.98])
hca.Title.String = 'Initialising uniform in v^1';

hca = h(isub); isub = isub + 1;
bar(hca,x_centers_up,N2)
hca.XTick = x_edges;
irf_legend(hca,sprintf('mean(x) = %g, mean(x^2) = %g',m2,m2_2),[0.98 0.98])
hca.Title.String = 'Initialising uniform in v^2';

hca = h(isub); isub = isub + 1;
bar(hca,x_centers_up,N3)
hca.XTick = x_edges;
irf_legend(hca,sprintf('mean(x) = %g, mean(x^2) = %g',m2,m3_2),[0.98 0.98])
hca.Title.String = 'Initialising uniform in v^3';

c_eval('h(?).XGrid = ''on'';',1:numel(h))
c_eval('h(?).XTick = x_edges;',1:numel(h))
c_eval('h(?).FontSize = 12;',1:numel(h))

%%
hca = h(isub); isub = isub + 1;
plot(hca,x_centers,f)

hca = h(isub); isub = isub + 1;
%bar(hca,x_centers_up,N1)
histogram(hca,cat(1,x1{:}),x_edges,'Normalization','pdf')
%hca.XTick = x_edges;

hca = h(isub); isub = isub + 1;
%bar(hca,x_centers_up,N2)
histogram(hca,cat(1,x2{:}),x_edges,'Normalization','pdf')
hca.XTick = x_edges;

hca = h(isub); isub = isub + 1;
%bar(hca,x_centers_up,N3)
histogram(hca,cat(1,x3{:}),x_edges,'Normalization','pdf')
hca.XTick = x_edges;

c_eval('h(?).XGrid = ''on'';',1:numel(h))
c_eval('h(?).XTick = x_edges;',1:numel(h))
%%
h = setup_subplots(3,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,x_centers.^1,f,'*-')
hca.XLabel.String = 'x^1';


hca = h(isub); isub = isub + 1;
plot(hca,x_centers.^2,f,'*-')
hca.XLabel.String = 'x^2';


hca = h(isub); isub = isub + 1;
plot(hca,x_centers.^3,f,'*-')
hca.XLabel.String = 'x^3';