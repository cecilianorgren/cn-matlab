% 
tintOverview = irf.tint('2015-10-16T10:33:20.00Z/2015-10-16T10:34:00.00Z');
tintDist = tintOverview(1)+[1 1.1];
for ii = 1:38  
  tintDist = tintDist+1;
  tintDist.disp
  mms_2015Oct16.plot_psds_4
  cn.print(['plot_psds_4_' irf_time(tintDist(1),'epochtt>utc_yyyymmddTHHMMSS.mmm')]);
end