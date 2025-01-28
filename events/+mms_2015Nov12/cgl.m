
tref = irf_time('2015-11-12T07:19:21.00Z','utc>epochtt');
tref = irf_time('2015-11-12T07:19:30.00Z','utc>epochtt');
c_eval('cglPe?par = pi/6*(ne?/ne?.resample(tref)).^3.*(gseB?.abs.resample(tref).data/gseB?.resample(ne?).abs).^2;',ic)
%c_eval('Pr?par'
%cglPe1par = 
%aa = pi/6*(ne1/ne1.resample(tref)).^3;
%bb = (gseB1.abs.resample(tref).data/gseB1.resample(ne1).abs).^2;
