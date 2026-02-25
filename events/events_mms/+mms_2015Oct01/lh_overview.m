
tlimit = irf.tint(B3fg.time.start.utc,B3fg.time.stop.utc);
tlimit = irf.tint(dslE3.time.start.utc,dslE3.time.stop.utc);


h = irf_plot(7);

hca = irf_panel('B');
irf_plot(hca,B3fg);
hca.YLabel.String = ['B_{FG} [' B3fg.units ']'];

hca = irf_panel('Bscm');
irf_plot(hca,B3sc.tlim(tlimit));
hca.YLabel.String = 'B_{SC} [nT]';

hca = irf_panel('E');
irf_plot(hca,dslE3);
hca.YLabel.String = ['E [' dslE3.units ']'];

hca = irf_panel('P_{SC}');
irf_plot(hca,P3*(-1));
hca.YLabel.String = ['-V_{SC} [' P3.units ']'];

hca = irf_panel('n');
irf_plot(hca,{ne3_lowres,ni3_lowres},'comp')
hca.YLabel.String = ['n [' ne3.units ']'];
irf_legend(hca,{'n_e','n_i'},[0.98 0.95])

hca = irf_panel('T');
irf_plot(hca,Te3_lowres)
hca.YLabel.String = ['T_e [' Te3.units ']'];



hca = irf_panel('vi');
irf_plot(hca,vi3)
hca.YLabel.String = ['v_i [' vi3.units ']'];

hca = irf_panel('ve');
irf_plot(hca,ve3)
hca.YLabel.String = ['v_e [' vi3.units ']'];

irf_zoom(h,'x',tlimit)
irf_zoom(h,'y')