h = irf_plot(4);

hca = irf_panel('orig');
irf_spectrogram(hca,ef1D1_nobg.specrec)

hca = irf_panel('smooth 3');
irf_spectrogram(hca,ef1D1_nobg.smooth(3).specrec)

hca = irf_panel('smooth 6');
irf_spectrogram(hca,ef1D1_nobg.smooth(6).specrec)

hca = irf_panel('smooth 3 3');
irf_spectrogram(hca,ef1D1_nobg.smooth(3).smooth(3).specrec)

%hca = irf_panel('smooth 6');
%irf_spectrogram(hca,ef1D1_orig.smooth(6).specrec)
