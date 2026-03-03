%% With TSeries implemented
en = iPDist_counts.find_low_counts('counts',8,'nMovMean',[3 3],'output','energy');
mat = iPDist_counts.find_low_counts('counts',8,'nMovMean',[3 3],'output','mat');

pd_mat = iPDist.mask('energy','mat',mat);
pd_en = iPDist.mask('energy','max',en);

h = irf_plot(3);

hca = irf_panel('pd');
irf_spectrogram(hca,iPDist.omni.deflux.specrec)
%hca.YScale = 'log';

hca = irf_panel('pd max E');
irf_spectrogram(hca,pd_en.omni.deflux.specrec)
%hca.YScale = 'log';

hca = irf_panel('pd mask E');
irf_spectrogram(hca,pd_mat.omni.deflux.specrec)
%hca.YScale = 'log';

irf_plot_axis_align(h)
hlinks = linkprop(h,{'CLim','XLim','YLim'});
colormap(flipdim(irf_colormap(hca,"spectral"),1))
h(end).XTickLabelRotation = 0;

%%
%tsElow = iPDist_counts.find_noise_energy_limit_counts(5,nMovMean);
tsElow = iPDist_counts.find_low_counts('counts',5,'nMovMean',3);
en = iPDist_counts.find_low_counts('counts',5,'nMovMean',3,'output','energy');
mat = iPDist_counts.find_low_counts('counts',5,'nMovMean',3,'output','mat');

%data = nansum(iPDist_counts.data(:,:,:),3);
%data = movmean(data,[5 5],1);
%mask = zeros(size(data));
%mask(find(data<5)) = 1;
%limit = find(data<5,1,'first');
pd_en = iPDist.mask('energy','mat',en);
pd_mat = iPDist.mask('energy','max',mat);

h = setup_subplots(4,1);
isub = 1;

hca = h(isub); isub = isub + 1;
pcolor(hca,data')
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,mask')
shading(hca,'flat')

hca = h(isub); isub = isub + 1;
pcolor(hca,log(pd.omni.deflux.data)')
shading(hca,'flat')
hcb = colorbar(hca);

hca = h(isub); isub = isub + 1;
pcolor(hca,log(iPDist.omni.deflux.data)')
shading(hca,'flat')
hcb = colorbar(hca);

irf_plot_axis_align(h)
hlinks = linkprop(h([3 4]),{'CLim','XLim','YLim'});
irf_plot_axis_align(h)