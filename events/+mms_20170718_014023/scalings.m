%% Egedal2015
% eq. 6
units = irf_units;
tref = irf_time('2017-07-18T01:40:35.633701660Z','utc>EpochTT');
betaref = beta1e.filt(0,0.5,[],3).resample(tref).data;
vteref = vte1.filt(0,0.5,[],3).resample(tref).data;
vepar = vteref*4*sqrt(units.me/units.mp)/sqrt(betaref);
