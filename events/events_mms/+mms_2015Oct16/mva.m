tint = irf.tint('2015-10-16T10:33:00.000Z/2015-10-16T10:36:00.000Z');
irf_minvar_gui(dmpaB3.tlim(tint))
%%
tint = irf.tint('2015-10-16T10:33:24.000Z/2015-10-16T10:36:32.000Z');
%irf_minvar_gui(dmpaB3scm.tlim(tint))
%irf_minvar_gui(irf_filt(dmpaB4scm,100,4000,8192,3))
irf_minvar_gui(irf_filt(dslE4brst,100,4000,8192,3))
