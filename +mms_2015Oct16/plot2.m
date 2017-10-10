tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T10:50:00.000Z');
%tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
h = irf_plot(6);
hca = irf_panel('B1');
irf_plot(hca,dfg1);
hca.YLabel.String = {'B_{1,DMPA}','[nT]'};
irf_legend(hca,{'B_x','B_y','B_z'},[0.99 0.95]);

hca = irf_panel('vi1');
irf_plot(hca,vi1fast);
hca.YLabel.String = {'v_{i,DSC}','[km/s]'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.99 0.95]);

hca = irf_panel('ve1');
irf_plot(hca,ve1fast);
hca.YLabel.String = {'v_{e,DSC}','[km/s]'};
irf_legend(hca,{'v_x','v_y','v_z'},[0.99 0.95]);

hca = irf_panel('n1');
irf_plot(hca,{ne1fast,ni1fast},'comp');
hca.YLabel.String = {'n_{1}','[cm^{-1}]'};
irf_legend(hca,{'n_e','n_i'},[0.99 0.95]);
hca.YScale = 'log';
hca.YTick = 10.^[-1:1:3];
hca.YLim = [5e-1 1e2];

hca = irf_panel('E1');
irf_plot(hca,edp1);
hca.YLabel.String = {'E_{1,DSL}','[mV/m]'};
irf_legend(hca,{'E_x','E_y','E_z'},[0.99 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J_{DMPA}','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'J_x','J_y','J_z'},[0.99 0.95]);


irf_zoom(h,'x',tint)
irf_zoom(h,'y')

hca = irf_panel('n1');
hca.YLim = [5e-1 1e2];