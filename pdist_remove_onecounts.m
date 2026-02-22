

% Compare density with and without removing one counts
%c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);

time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT') +- 0;
nMovMean = 7;
elim = [3000 Inf];

ipdist = iPDist3.movmean(nMovMean);
ipdist_clean = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts);
ipdist_clean_elim = iPDist3.movmean(nMovMean,'RemoveOneCounts',iPDist3_counts).elim(elim);
ipdist_elim = iPDist3.elim(elim);

%epdist = ePDist3.movmean(nMovMean).tlim(time_xline + 10*[-1 1]);
%%
% nan2zero
moms_e = mms.psd_moments(epdist,scPot3);
moms_i = mms.psd_moments(ipdist,scPot3);
moms_i_clean = mms.psd_moments(ipdist_clean,scPot3);
moms_i_clean_elim = mms.psd_moments(ipdist_clean_elim,scPot3);
moms_i_elim = mms.psd_moments(ipdist_elim,scPot3);
%%
h = irf_plot(7);

ylim_e = [1.2850 3.1765e+04];

hca = irf_panel('DEF i');
irf_spectrogram(hca,iPDist3.deflux.omni.specrec);
hca.YScale = 'log';
hca.YLim = ylim_e;

hca = irf_panel('DEF i movmean');
irf_spectrogram(hca,ipdist.deflux.omni.specrec);
hca.YScale = 'log';
hca.YLim = ylim_e;

hca = irf_panel('DEF i one counts removed');
irf_spectrogram(hca,ipdist_clean.deflux.omni.specrec);
hca.YScale = 'log';
hca.YLim = ylim_e;

hca = irf_panel('DEF i one counts removed elim');
irf_spectrogram(hca,ipdist_clean_elim.deflux.omni.specrec);
hca.YScale = 'log';
hca.YLim = ylim_e;

hca = irf_panel('DEF i elim');
irf_spectrogram(hca,ipdist_elim.deflux.omni.specrec);
hca.YScale = 'log';
hca.YLim = ylim_e;

hca = irf_panel('density');
irf_plot(hca,{ne3,ni3,moms_e.n_psd,moms_i.n_psd,moms_i_clean.n_psd,moms_i_clean_elim.n_psd,moms_i_elim.n_psd},'comp')
irf_legend(gca,{'ne3','ni3','moms_e','moms_i','moms_i onecounts removed','moms_i onecounts removed elim','moms_i elim'},[0.02 0.98])
hca.YLabel.String = 'n (cc)';

hca = irf_panel('density 2');
irf_plot(hca,{ne3,moms_i_clean.n_psd,moms_i_clean_elim.n_psd,moms_i_elim.n_psd},'comp')
irf_legend(gca,{'ne3','moms_i onecounts removed','moms_i onecounts removed elim','moms_i elim'},[0.02 0.98])
hca.YLabel.String = 'n (cc)';

irf_pl_mark(h,time_xline,'k')
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
irf_plot_axis_align


