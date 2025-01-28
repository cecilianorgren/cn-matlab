tint = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z'); % magnetosphere-magnetosheath-magnetosphere
h = irf_plot(10,'newfigure');

ic  = 3;
hca = irf_panel(irf_ssub('B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?brst.tlim(tint).x,dmpaB?brst.tlim(tint).y,dmpaB?brst.tlim(tint).z,dmpaB?brst.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {irf_ssub('B',ic),'(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z','|B|'},[0.95 0.95]);

hca = irf_panel(irf_ssub('brst E?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dslE?brst.tlim(tint).x,dslE?brst.tlim(tint).y,dslE?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'E','(mV/m)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'E_x','E_y','E_z'},[0.95 0.95]);

hca = irf_panel(irf_ssub('brst B?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{dmpaB?scm.tlim(tint).x,dmpaB?scm.tlim(tint).y,dmpaB?scm.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'B','(nT)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'B_x','B_y','B_z'},[0.95 0.95]);

hca = irf_panel('grad P');
set(hca,'ColorOrder',mms_colors('xyza'))
irf_plot(hca,{gradPe.tlim(tint).x,gradPe.tlim(tint).y,gradPe.tlim(tint).z},'comp');
hca.YLabel.String = {'\nabla P_e ','(nP/m^2)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'\nabla P_e_x','\nabla P_e_y','\nabla P_e_z'},[0.95 0.95]);

hca = irf_panel(irf_ssub('ve?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{ve?brst.tlim(tint).x,ve?brst.tlim(tint).y,ve?brst.tlim(tint).z},''comp'');',ic)
hca.YLabel.String = {'v_e','(km/s)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'v_x','v_y','v_z'},[0.95 0.95]);

hca = irf_panel(irf_ssub('ve? z',ic));
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{ve1brst.tlim(tint).z,ve2brst.tlim(tint).z,ve3brst.tlim(tint).z,ve4brst.tlim(tint).z},'comp');
hca.YLabel.String = {'v_e_z','(km/s)'};
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('brst n');
%c_eval('irf_plot(hca,{ne?_lowres.tlim(tint),ni?_lowres.tlim(tint)},''comp'');',ic);
set(hca,'ColorOrder',mms_colors('1234'))
irf_plot(hca,{ne1_lowres.tlim(tint),ne2_lowres.tlim(tint),ne3_lowres.tlim(tint),ne4_lowres.tlim(tint)},'comp');
hca.YLabel.String = {'n','(cm^{-3})'};
hca.YScale = 'lin';
set(hca,'ColorOrder',mms_colors('1234'))
irf_legend(hca,{'mms1','mms2','mms3','mms4'},[0.95 0.95]);

hca = irf_panel('E spectrogram');
c_eval('specrec = pfftE?ts; specrec.p = specrec.p{2}; specrec.f_label = ''f (Hz)''; specrec.p_label = ''(mV/m)^2/Hz'';',ic)
irf_spectrogram(hca,specrec)
hca.YScale = 'log';
hca.YTick = 10.^[1 2 3 4];
if 0
hca = irf_panel('E spectrogram hmfe');
c_eval('specrec = pfftE?hmfe; specrec.p = specrec.p{1}; specrec.f_label = ''f (Hz)''; specrec.p_label = ''(mV/m)^2/Hz'';',ic)
irf_spectrogram(hca,specrec)
hca.YScale = 'log';
hca.YTick = 10.^[1 2 3 4 5];
ylim = hca.YLim; 
hold(hca,'on')
c_eval('irf_plot(hca,flh?brst,''color'',[0 0 0])',ic)
c_eval('irf_plot(hca,fce?brst,''color'',[0 1 1])',ic)
c_eval('irf_plot(hca,fpe?brst,''color'',[1 0 1])',ic)
hold(hca,'off')
hca.YLim = ylim;
grid(hca,'off')
end

hca = irf_panel('B spectrogram');
c_eval('specrec = pfftB?ts; specrec.p = specrec.p{2}; specrec.f_label = ''f (Hz)''; specrec.p_label = ''(nT)^2/Hz'';',ic)
irf_spectrogram(hca,specrec)
hca.YScale = 'log';
hca.CLim = [-13 1];
hca.YTick = 10.^[1 2 3 4];

irf_zoom(h,'x',tint)
%irf_zoom(h,'y')
irf_plot_axis_align
h(1).Title.String = irf_ssub('MMS ?',ic);


%%
hca = irf_panel(irf_ssub('grad P/P?',ic));
set(hca,'ColorOrder',mms_colors('xyza'))
c_eval('irf_plot(hca,{gradPe.tlim(tint).x/Pe?_lowres.abs,gradPe.tlim(tint).y/Pe?_lowres.tlim(tint).abs,gradPe.tlim(tint).abs/Pe?_lowres.tlim(tint).abs},''comp'');',ic)
hca.YLabel.String = {'\nabla P_e/P ','(nP/m^2)'};
set(hca,'ColorOrder',mms_colors('xyza'))
irf_legend(hca,{'\nabla P_e_x','\nabla P_e_y','\nabla P_e_z'},[0.95 0.95]);

