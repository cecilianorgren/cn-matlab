pdist = obsPDist.tlim(tintObs);


npanels = 8;
h1 = irf_plot(npanels);
elim = [000 10000];
clim = [-36 -26];
colors = mms_colors('matlab');

nangles = 15;
isub = 1;

if 0 % .omni.specrec('energy')
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.tlim(tintObs).omni.specrec('energy'),'log');
  hca.YScale = 'log';
  hca.YLim = elim;
  %hca.CLim = clim;
end
if 0 % .e64.omni.specrec('energy')
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.tlim(tintObs).e64.omni.specrec('energy'),'log');
  hca.YScale = 'log';
  hca.YLim = elim;
  %hca.CLim = clim;
end
if 0 % .einterp('linear').omni.specrec('energy')
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.tlim(tintObs).einterp('linear').omni.specrec('energy'),'log');
  hca.YScale = 'log';
  hca.YLim = elim;
  %hca.CLim = clim;
end
if 0 % .einterp('pchip').omni.specrec('energy')
  hca = h1(isub); isub = isub + 1;
  irf_spectrogram(hca,obsPDist.tlim(tintObs).einterp('pchip').omni.specrec('energy'),'log');
  hca.YScale = 'log';
  hca.YLim = elim;
  %hca.CLim = clim;
end

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).pitchangles(gseB1,nangles,0).specrec('pa'),'log');
irf_legend(hca,'pdist.pitchangles(gseB1,nangles,0).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).pitchangles(gseB1,nangles,1).specrec('pa'),'log');
irf_legend(hca,'pdist.pitchangles(gseB1,nangles,1).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).e64.pitchangles(gseB1,nangles,0).specrec('pa'),'log');
irf_legend(hca,'pdist.e64.pitchangles(gseB1,nangles,0).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).e64.pitchangles(gseB1,nangles,1).specrec('pa'),'log');
irf_legend(hca,'pdist.e64.pitchangles(gseB1,nangles,1).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).einterp('pchip').pitchangles(gseB1,nangles,0).specrec('pa'),'log');
irf_legend(hca,'pdist.einterp(''pchip'').pitchangles(gseB1,nangles,0).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).einterp('pchip').pitchangles(gseB1,nangles,1).specrec('pa'),'log');
irf_legend(hca,'pdist.einterp(''pchip'').pitchangles(gseB1,nangles,1).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).einterp('linear').pitchangles(gseB1,nangles,0).specrec('pa'),'log');
irf_legend(hca,'pdist.einterp(''linear'').pitchangles(gseB1,nangles,0).specrec(''pa'')',[0.99 0.99])

hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.tlim(tintObs).einterp('linear').pitchangles(gseB1,nangles,1).specrec('pa'),'log');
irf_legend(hca,'pdist.einterp(''linear'').pitchangles(gseB1,nangles,1).specrec(''pa'')',[0.99 0.99])

for ip = 1:npanels
  h1(ip).CLim = [-26.7 -26.3];
end
%%
tic; obsPDist.pitchangles(gseB1,nangles,1); toc
tic; obsPDist.pitchangles(gseB1,nangles,0); toc
%%


hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.omni.tlim(tintObs).elim(elim).specrec('energy'),'log'); hca.YScale = 'log';
hca.YLim = [10 1000];
hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.pitchangles(gseB1,nangles).tlim(tintObs).elim(elim).specrec('pa'),'log');
hca = h1(isub); isub = isub + 1;
irf_spectrogram(hca,obsPDist.einterp.pitchangles(gseB1,nangles).tlim(tintObs).elim(elim).specrec('pa'),'log');

% h1(1).CLim = [7.3 8.1];
% h1(2).CLim = [7.0 8.1];
% h1(3).CLim = [-26.5 -26.0];
% h1(4).CLim = [-26.5 -24];
irf_plot_axis_align
