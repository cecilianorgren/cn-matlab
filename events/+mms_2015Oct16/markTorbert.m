% Roy Torberts requested intervals
tint1 = irf.tint('2015-10-16T13:00:00',8*60);
tint2 = irf.tint('2015-10-16T11:25:00',60*7);
irf_pl_mark(h,tint1.epochUnix')
irf_pl_mark(h,tint2.epochUnix')