%tr100 = PICTraj('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/trajectories.h5');
no02m = PIC('/Volumes/Fountain/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

tr = tr100.find([tr100.z0]==0,[tr100.x0]>93.5);

tr.plot_all_xz;