%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/dists.h5');

twpe = 24000;
xlim = mean(no02m.xi) + [-50 10];
zlim = [-8 8];

cmap = pic_colors('thermal');
clabels = {'Ion temperature','Electron temperature'};
clabels = {'',''};
clims = {[0 0.7],[0 0.3]};
no02m.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map({'ti','te'}','cbarlabels',clabels,'clim',clims,'cmap',cmap,'A',0.5,'sep');

ha = findobj(gcf,'type','axes'); ha = ha(end:-1:1);

hc = findobj(gcf,'type','contour');
c_eval('hc(?).Color = 0.3*[1 1 1];',1:numel(hc))

c_eval('ha(?).FontSize = 12;',1:numel(ha))
ha(1).Title.String = 'Plasma temperature';

irf_legend(ha(1),{'Ions'},[0.95 0.92],'fontsize',12,'fontweight','bold','color',0.8*[1 1 1])
irf_legend(ha(2),{'Electrons'},[0.95 0.92],'fontsize',12,'fontweight','bold','color',0.8*[1 1 1])

%% Select boxes
inds_ = [48 378];
inds_maybe = [148 248 348];
inds_maybe = [148 228 348];

twpe = 24000;
xlim = mean(no02m.xi) + [-50 10];
zlim = [-8 8];

cmap = pic_colors('thermal');
clabels = {'Ion temperature (xx)','Electron temperature (xx)'};
clabels = {'',''}; 
clims = {[0 0.7],[0 0.3]};
%no02m.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map({'txx(3)','txx(4)'}','cbarlabels',clabels,'clim',clims,'cmap',cmap,'A',0.5,'sep');

nrows = 2;
ncols = 2;
h = setup_subplots(nrows,ncols,'horizontal');
isub = 1;
ele = [2 4 6];
ion = [1 3 5];
sumdim = 3;
inds = inds_;
for ii = inds  
  hca = h(isub); isub = isub + 1;
  f = ds.update_inds({ii}).fxyz(1,1,ele,sumdim);
  fxyz = ds.update_inds({ii}).f(1,1,ele);
  imagesc(hca,f.v,f.v,log10(f.f)')
  
  hca = h(isub); isub = isub + 1;
  f = ds.update_inds({ii}).fxyz(1,1,ion,sumdim);
  imagesc(hca,f.v,f.v,log10(f.f)')
end
colormap(pic_colors('candy4'))