%no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

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

