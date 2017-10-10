units = irf_units;

planetDistance = [0.39 0.723 1 1.524 5.203 9.539 19.18 30.06 39.53]*units.AU;
planetRadius = [4878 12104 12756 6787 142796 120660 51118 48600 2274]*1e3/2;
planetMP = [0.5 0 11 0 80 20 20 25 0].*planetRadius;

vsw = 400*1e3; % m/s
nswE = 5e6; % m-3
nSW = nswE*units.AU^2./(planetDistance.^2);
flux = nSW*vsw.*planetMP.^2;