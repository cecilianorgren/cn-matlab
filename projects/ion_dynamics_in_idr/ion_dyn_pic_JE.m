no02m = PIC('/Users/cecilianorgren/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

%%

pic = no02m.twpelim(21000).xlim([90 110]).zlim([-5 5]);

varstrs = {'-jex.*Ex','-jey.*Ey','-jez.*Ez','jix.*Ex','jiy.*Ey','jiz.*Ez'}';
clims = {0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1]};
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br,cmap_br,cmap_br,cmap_br,cmap_br,cmap_br}; 
pic.plot_map(varstrs,'clim',clims,'cmap',cmaps,'A',1,'sep')
%pic.plot_map(varstrs,'A',1)

%%
pic = no02m.twpelim(21000).xlim([101.5 102]).zlim([-2 2]);

varstrs = {{'ni','ne'},...
  {'vex','vey','vez','vExBy','vExBz'},...
  {'vix','viy','viz','vExBy','vExBz'},...
  {'-jex.*Ex','-jey.*Ey','-jez.*Ez','jix.*Ex','jiy.*Ey','jiz.*Ez'},...
  {'-jex.*Ex','-jey.*Ey','-jez.*Ez'},...
  {'jix.*Ex','jiy.*Ey','jiz.*Ez'},...
  {'Jx.*Ex','Jy.*Ey','Jz.*Ez'},...
  {'Jx.*Ex+Jy.*Ey+Jz.*Ez','jix.*Ex+jiy.*Ey+jiz.*Ez','-jex.*Ex-jey.*Ey-jez.*Ez'}}';
clims = {0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1],0.1*[-1 1]};
cmap_br = pic_colors('blue_red');
cmaps = {cmap_br,cmap_br,cmap_br,cmap_br,cmap_br,cmap_br}; 
pic.plot_line('z',varstrs)
%pic.plot_map(varstrs,'A',1)

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = 16;',1:numel(h))

linkprop(h(4:5),'YLim')