% Wal?n test
units = irf_units;
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?); facTe? = irf.ts_tensor_xyz(facTe?.time,facTe?.data);',ic)

%% Magnetosheat reference point.
timeMSH = irf_time('2015-10-16T10:33:35.000Z','utc>epochTT');
c_eval('dP? = irf.ts_scalar(facPe?.time,facPe?.xx.data-0.5*(facPe?.yy.data+facPe?.zz.data));',ic)
c_eval('Alfa? = dP?.resample(mvaB?)/units.mu0/mvaB?.abs2;',ic)
c_eval('Ve? = mvaVe?;',ic)
c_eval('B? = mvaB?;',ic)

c_eval('mshAlfa? = Alfa?.resample(timeMSH);',ic)
c_eval('mshVe? = mvaVe?.resample(timeMSH);',ic)
c_eval('mshB? = mvaB?.resample(timeMSH);',ic)

c_eval('dV?obs = Ve?-mshVe?.data;',ic)
c_eval('dV?th = (mvaB?*(-Alfa?+1)-mshB?.data*(1-mshAlfa?.data))/ne?.sqrt/units.mu0*1e-6;',ic)

c_eval('walAngle? = irf.ts_scalar(dV?obs.time,acosd(dV?obs.dot(dV?th.resample(dV?obs)/dV?th.abs/dV?obs.abs).data));',ic)
c_eval('walR? = dV?obs.abs/dV?th.resample(dV?obs).abs;',ic)

h = irf_plot({dV1obs,dV1th,walAngle1,walR1});
mshRefMark = irf_pl_mark(h,timeMSH.epochUnix');

