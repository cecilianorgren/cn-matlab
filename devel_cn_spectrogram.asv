h = irf_plot(3);

hca = irf_panel('irf');
irf_spectrogram(hca,iPDist1.deflux.omni.specrec,'log')
hca.YScale = 'log';


hca = irf_panel('cn');
cn_spectrogram(hca,iPDist1.deflux.omni.specrec,'log')
hca.YScale = 'log';

linkprop(h,{'XLim','YLim','CLim'})
h(end).XTickLabelRotation = 0;


