% Load B
tint_J = irf.tint('2015-10-16T08:00:00.000Z/2015-10-16T14:00:00.000Z');
c_eval('[Bxyz?,dobj_Bxyz?] = cn_get_ts(''mms?_dfg_srvy_ql'',''mms?_dfg_srvy_dmpa'',tint_J(1));',sc)
c_eval('Bxyz? = Bxyz?.resample(Bxyz1);',1:4);

R  = mms.get_data('R_gse',tint_J);
c_eval('Rxyz? = TSeries(R.time,R.gseR?,''to'',1);',1:4);
c_eval('Rxyz? = Rxyz?.resample(Bxyz1);',1:4);

% Assuming GSE and DMPA are the same coordinate system.
[j,divB,B,jxB,divTshear,divPb] = c_4_j('Rxyz?','Bxyz?');

%% Make plot
h = irf_plot(3,'newfigure');

hca = irf_panel('B1');
irf_plot(hca,Bxyz1);
ylabel(hca,{'B_{DMPA}','[nT]'},'Interpreter','tex');
irf_legend(hca,{'B_{x}','B_{y}','B_{z}'},[0.88 0.10])
irf_legend(hca,'(a)',[0.99 0.98],'color','k')

hca = irf_panel('R1');
irf_plot(hca,Rxyz1);
ylabel(hca,{'B_{DMPA}','[nT]'},'Interpreter','tex');
irf_legend(hca,{'M1','M2','M3'},[0.88 0.10])
irf_legend(hca,'(b)',[0.99 0.98],'color','k')

hca = irf_panel('J');
j.data = j.data*1e9;
irf_plot(hca,j);
ylabel(hca,{'J_{DMPA}','(nA m^{-2})'},'Interpreter','tex');
irf_legend(hca,{'J_{x}','J_{y}','J_{z}'},[0.88 0.10])
irf_legend(hca,'(c)',[0.99 0.98],'color','k')

%%
irf_zoom(h,'x',tint_J)
irf_zoom(h,'y')

% Roy Torberts requested intervals
tint1 = irf.tint('2015-10-16T13:00:00',8*60);
tint2 = irf.tint('2015-10-16T11:25:00',60*7);
irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')
