function f = f(in)
% CN.F takes time series as input and returns average frequency. Assumes
% first column is time.

dt = diff(in(:,1));
fall = 1./dt;
f = irf.nanmean(fall);