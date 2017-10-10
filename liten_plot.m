figure('name','E and Ne');
h=irf_plot(4);
tint2=toepoch([t1,t2]);
irf_plot(h(1),diE3);
irf_zoom('x',tint2)
irf_plot(h(2),gsmB3);
irf_zoom('x',tint2)
irf_plot(h(3),Ne3);
irf_zoom('x',tint2)
irf_plot(h(4),Te);
irf_zoom(h,'x',tint2)