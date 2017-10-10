% vD = -grad(p)xB/qsnsB^2 - opposite directions for e- and i+

vd = 150; % km/s
n_loc_cm = n_loc*1e6; % m^-3

units=irf_units;
Te_loc_MK = irf_resamp(parTe3MK,mean(tint),'nearest'); Te_loc_MK = Te_loc_MK(2);
%Ti_loc_MK = irf_resamp(Ti3MK,mean(tint),'nearest'); Ti_loc_MK = Ti_loc_MK(2);
Te_loc = Te_loc_MK*units.kB/units.e*1e6;
%Ti_loc = Ti_loc_MK*units.kB/units.e*1e6;
Ti_loc=3*Te_loc;
Ln = units.kB*Te_loc_MK*1e6/units.e/(B0*1e-9)/vd*1e-6; % km/s
ri_loc = re_loc*sqrt(1836)*sqrt(Ti_loc/Te_loc);
disp(['for vD = ' num2str(vd) 'km/s  B0 = ' num2str(B0,'%.0f') 'nT  n = ' num2str(n_loc,'%.0f') 'cm-3  Te = ' num2str(Te_loc,'%.f') 'eV  Ti = ' num2str(Ti_loc,'%.f') 'eV: Ln = ' num2str(Ln,'%.0f') 'km (rop = ' num2str(ri_loc,'%.0f') 'km => Ln/rop = ' num2str(Ln/ri_loc,'%.1f') ')'])