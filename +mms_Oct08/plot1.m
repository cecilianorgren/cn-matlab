tint = irf.tint('2015-10-16T10:30:00.000Z/2015-10-16T10:50:00.000Z');
tint = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');

%tint_zoom = [edp3.time.start.epochUnix+60*60 edp3.time.stop.epochUnix];
h = irf_plot(6);
hca = irf_panel('Bx');
irf_plot(hca,{dfg1.x.tlim(tint),dfg2.x.tlim(tint),dfg3.x.tlim(tint),dfg4.x.tlim(tint)},'comp');
hca.YLabel.String = {'B_{x,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('By');
irf_plot(hca,{dfg1.y.tlim(tint),dfg2.y.tlim(tint),dfg3.y.tlim(tint),dfg4.y.tlim(tint)},'comp');
hca.YLabel.String = {'B_{y,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Bz');
irf_plot(hca,{dfg1.z.tlim(tint),dfg2.z.tlim(tint),dfg3.z.tlim(tint),dfg4.z.tlim(tint)},'comp');
hca.YLabel.String = {'B_{z,DMPA}','[nT]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ex');
irf_plot(hca,{edp1.x.tlim(tint),edp2.x.tlim(tint),edp3.x.tlim(tint),edp4.x.tlim(tint)},'comp');
hca.YLabel.String = {'E_{x,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('Ey');
irf_plot(hca,{edp1.y.tlim(tint),edp2.y.tlim(tint),edp3.y.tlim(tint),edp4.y.tlim(tint)},'comp');
hca.YLabel.String = {'E_{y,DSL}','[mV/m]'};
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);

hca = irf_panel('J');
irf_plot(hca,j);
%hca.YLabel.String = 'J [nAm^{-2}]';
ylabel(hca,{'J_{DMPA}','[nA m^{-2}]'},'Interpreter','tex');
irf_legend(hca,{'M1','M2','M3','M4'},[0.95 0.95]);


irf_zoom(h,'x',tint)
irf_zoom(h,'y')

