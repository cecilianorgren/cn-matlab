pic1 = PIC('/Users/cecilia/Data/PIC/rec_onset_4/data_h5/fields_F10_E01.h5');
%pic1 = PIC('/Users/cecilia/Data/PIC/rec_onset_4/data_h5/fields_F10_E01.h5');
tite05 = PIC('/Users/cecilia/Data/PIC/varying_tite/tite_05/fields.h5');
tite10 = PIC('/Users/cecilia/Data/PIC/varying_tite/tite_10/fields.h5');

%%
pic = tite05;
pic = pic(pic.length).xlim(mean(pic.xi) + 1 + [-10  10]).zlim([-3 3]);
pic.plot_map({'Ez','By'}','A',1,'smooth',1)
colormap(pic_colors('blue_re'))
colormap(pic_colors('blue_red'))
%%

pic = tite10;
pic = pic(pic.length).xlim(mean(pic.xi) + 1 + [-15  15]).zlim([-5 5]);
pic.plot_map({'Ez','By'}','A',1,'smooth',1)
colormap(pic_colors('blue_re'))
colormap(pic_colors('blue_red'))

%%
pic = pic1;
pic = pic(pic.length).xlim(mean(pic.xi) + 1 + [-10  10]).zlim([-3 3]);
pic.plot_map({'Ez','By'}','A',1,'smooth',1)
colormap(pic_colors('blue_re'))
colormap(pic_colors('blue_red'))