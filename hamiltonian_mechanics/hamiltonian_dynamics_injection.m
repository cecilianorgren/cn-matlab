h5filepath = '/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5';
no02m = PIC(h5filepath);


%%
pic = no02m.twpelim(20000);

Ay = pic.A;
Ez = pic.Ez;
By = pic.By;

Ez_edges = linspace(-0.2,0.2,100);
By_edges = linspace(-0.02,0.02,101);

[N,edges,mid,locs] = histcn([Ez(:) By(:)],Ez_edges,By_edges);



nRows = 1;
nCols = 1;
h = setup_subplots(nRows,nCols);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,mid{1},mid{2},log10(N)');
shading(hca,'flat')
hca.YLabel.String = 'B_y';
hca.XLabel.String = 'E_z';
