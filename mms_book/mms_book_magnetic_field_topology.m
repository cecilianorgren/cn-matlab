no02m = PIC('/Users/cno062/Data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');

twpe = 24000;
pic = no02m.twpelim(twpe);

Ay = pic.A;


%%
hca = subplot(1,1,1);

Alev = min(Ay(:)):1:max(Ay(:));
Alev = -8:1:-5;
hc = contour(hca,pic.xi,pic.zi,Ay',Alev,'k','LineWidth',1);

hca.XLim = [105 140];

hca.Visible = "off";