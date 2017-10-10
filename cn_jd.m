function j = cn_jd(vph,n)
% Calculates the current needed to account for the measured electron holes
% phase velocity, if it is the Buneman instability that is responsible.

units=irf_units;
vd=vph*(16*units.mp/units.me)^(1/3);
j=units.e*n*1e6*vd*1e3; % A
j=j*1e9; % nA