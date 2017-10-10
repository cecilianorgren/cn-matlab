% irf_plot({j,je1+0*ji1.resample(je1.time),je2+0*ji2.resample(je2.time),je3+0*ji3.resample(je3.time),je4+0*ji4.resample(je4.time)},'comp')

h = irf_plot(10);

ic = 1;
hca = irf_panel('B?');
c_eval('irf_plot(hca,dmpaB?);',ic);
hca.YLabel.String = {irf_ssub('B_?',ic),'[nT]'};
irf_legend(hca,{'B_x','B_y','B_z'},[0.95 0.95]);

hca = irf_panel('jx e');
hlines = irf_plot(hca,{je1.x,je2.x,je3.x,je4.x,j.x},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{x}','[nA/m^2]'};
irf_legend(hca,{'j_e mms1','j_e mms2','j_e mms3','j_e mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jy e');
hlines = irf_plot(hca,{je1.y,je2.y,je3.y,je4.y,j.y},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{y}','[nA/m^2]'};
irf_legend(hca,{'j_e mms1','j_e mms2','j_e mms3','j_e mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jz e');
hlines = irf_plot(hca,{je1.z,je2.z,je3.z,je4.z,j.z},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{z}','[nA/m^2]'};
irf_legend(hca,{'j_e mms1','j_e mms2','j_e mms3','j_e mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jx i');
hlines = irf_plot(hca,{ji1.x,ji2.x,ji3.x,ji4.x,j.x},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{x}','[nA/m^2]'};
irf_legend(hca,{'j_i mms1','j_i mms2','j_i mms3','j_i mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jy i');
hlines = irf_plot(hca,{ji1.y,ji2.y,ji3.y,ji4.y,j.y},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{y}','[nA/m^2]'};
irf_legend(hca,{'j_i mms1','j_i mms2','j_i mms3','j_i mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jz i');
hlines = irf_plot(hca,{ji1.z,ji2.z,ji3.z,ji4.z,j.z},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{z}','[nA/m^2]'};
irf_legend(hca,{'j_i mms1','j_i mms2','j_i mms3','j_i mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jx tot');
hlines = irf_plot(hca,{je1.x+ji1.x.resample(je1.time),je2.x+ji2.x.resample(je2.time),je3.x+ji3.x.resample(je3.time),je4.x+ji4.x.resample(je4.time),j.x},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{x}','[nA/m^2]'};
irf_legend(hca,{'j_{tot} mms1','j_{tot} mms2','j_{tot} mms3','j_{tot} mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jy tot');
hlines = irf_plot(hca,{je1.y+ji1.y.resample(je1.time),je2.y+ji2.y.resample(je2.time),je3.y+ji3.y.resample(je3.time),je4.y+ji4.y.resample(je4.time),j.y},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{y}','[nA/m^2]'};
irf_legend(hca,{'j_{tot} mms1','j_{tot} mms2','j_{tot} mms3','j_{tot} mms4','4sc'},[0.95 0.95]);

hca = irf_panel('jz tot');
hlines = irf_plot(hca,{je1.z+ji1.z.resample(je1.time),je2.z+ji2.z.resample(je2.time),je3.z+ji3.z.resample(je3.time),je4.z+ji4.z.resample(je4.time),j.z},'comp');
hlines.Children(1).LineWidth = 2;
hca.YLabel.String = {'j_{z}','[nA/m^2]'};
irf_legend(hca,{'j_{tot} mms1','j_{tot} mms2','j_{tot} mms3','j_{tot} mms4','4sc'},[0.95 0.95]);

irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align

for ii = 2:10
  h(ii).YLim = 2000*[-1 1];
end


