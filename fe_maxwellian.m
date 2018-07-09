function fout = fe_maxwellian(v,n,vt,vd)

v = v-vd; % move into drifting reference frame
fout = n*(1/pi./vt.^2)^(1/2)*exp(-v.^2./vt.^2);
