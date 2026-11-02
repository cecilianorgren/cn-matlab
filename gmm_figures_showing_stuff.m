  
%% plot difference between original and cleaned/smoothed VDF
h = irf_plot(2);
c_eval('h(?).Position(2) = h(?).Position(2) + 0.03;',1:numel(h))

hca = irf_panel('original deflux');
irf_spectrogram(hca,iPDist.omni.specrec,'log');

hca = irf_panel('cleand and smoothed deflux');
irf_spectrogram(hca,PD.omni.specrec,'log');

%hlinks = linkprop(h(end:-1:1),{'CLim'});
hlinks = linkprop(h,{'CLim'});
irf_timeaxis(h)
h(end).XTickLabelRotation = 0;
c_eval('h(?).FontSize = 12;',1:numel(h))
