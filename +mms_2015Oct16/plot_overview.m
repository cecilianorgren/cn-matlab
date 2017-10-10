tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
h = irf_plot(9);

%ic = 1;

hca = irf_panel('B');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {'B','(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[1.01 0.5]);

hca = irf_panel('brst n');
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
c_eval('irf_plot(hca,ne?_lowres.tlim(tint));',ic);
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
irf_legend(hca,{'n_e'},[1.01 0.5]);

hca=irf_panel('edist');
c_eval('irf_spectrogram(hca,eDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{e,OMNI}','(eV)'},'Interpreter','tex');

hca=irf_panel('idist');
c_eval('irf_spectrogram(hca,iDEFomni?,''log'',''donotfitcolorbarlabel'');',ic)
set(hca,'yscale','log');
set(hca,'ytick',[1e1 1e2 1e3 1e4]);
ylabel(hca,{'E_{i,OMNI}','(eV)'},'Interpreter','tex');

hca = irf_panel('brst Te');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{Te?_lowres.tlim(tint).x,Te?_lowres.tlim(tint).y,Te?_lowres.tlim(tint).z},''comp'');',ic);
hca.YLabel.String = {'T','(eV)'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'T_{xx}','T_{yy}','T_{zz}'},[1.01 0.5]);

hca = irf_panel('brst vi');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{vi?_lowres.tlim(tint).x,vi?_lowres.tlim(tint).y,vi?_lowres.tlim(tint).z},''comp'');',ic);
hca.YLabel.String = {'v_i','(km/s)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'v_x','v_y','v_z'},[1.01 0.5]);

hca = irf_panel('brst ve');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{ve?brst.tlim(tint).x,ve?brst.tlim(tint).y,ve?brst.tlim(tint).z},''comp'');',ic);
hca.YLabel.String = {'v_e','(km/s)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'v_x','v_y','v_z'},[1.01 0.5]);

hca = irf_panel('j');
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{j.x,j.y,j.z},'comp');
%hca.YLabel.String = {'v_i','(km/s)'};
ylabel(hca,{'J','(nA/m^2)'},'interpreter','tex')
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'J_x','J_y','J_z'},[1.01 0.5]);

hca = irf_panel('brst E');
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y'},[1.01 0.5]);



irf_zoom(h,'x',tint)
irf_plot_axis_align
irf_zoom(h([1:2 5:9]),'y')
h(1).Title.String = irf_ssub('MMS ?',ic);
h(1).Title.FontSize = 14;
for ii = 1:9, h(ii).YLabel.FontSize = 14; end
