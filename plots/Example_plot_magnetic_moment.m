%% Load data

%% Construct mu distributions
mu2i = iPDist2.elim([100 40000]).mu(dmpaB2);
mu2e_elim000 = ePDist2.elim([000 40000]).mu(dmpaB2);
mu2e_elim050 = ePDist2.elim([050 40000]).mu(dmpaB2);
mu2e_elim100 = ePDist2.elim([100 40000]).mu(dmpaB2);

%% Plot results
nPanels = 6;
h = irf_plot(nPanels);

hca = irf_panel('B abs');
set(hca,'colororder',mms_colors('1'));
irf_plot(hca,gseB2.abs);
hca.YLabel.String = '|B| (nT)';
hca = irf_panel('scpot');
set(hca,'colororder',mms_colors('1'));
irf_plot(hca,scPot2);
hca.YLabel.String = 'V_{sc} (V)';

hca = irf_panel('mu i');
irf_spectrogram(hca,mu2i.specrec,'log');
hca.YLabel.String = '\mu (eV/T)';
hca.YScale = 'log';

hca = irf_panel('mu e elim 000');
irf_spectrogram(hca,mu2e_elim000.specrec,'log');
hca.YLabel.String = '\mu (eV/T)';
hca.YScale = 'log';

hca = irf_panel('mu e elim 050');
irf_spectrogram(hca,mu2e_elim050.specrec,'log');
hca.YLabel.String = '\mu (eV/T)';
hca.YScale = 'log';

hca = irf_panel('mu e elim 100');
irf_spectrogram(hca,mu2e_elim100.specrec,'log');
hca.YLabel.String = '\mu (eV/T)';
hca.YScale = 'log';

hlinks_e = linkprop(h(4:end),{'YLim','CLim'});
h(4).CLim = [-10 -4];

hca = irf_panel('mu i'); hca.CLim = [-5 -1];

irf_plot_axis_align(h)