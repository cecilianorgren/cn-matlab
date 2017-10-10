ic = 4;

%% Load data
interval = '20151030051444';
yyyy = interval(1:4);
mm = interval(5:6);
dd = interval(7:8);
starttime = interval(9:end);

tmp = ['tmpDataObj = dataobj(''/Volumes/DansHD2/data/mms?/fpi/brst/l2/des-dist/' yyyy '/' mm '/' dd '/mms?_fpi_brst_l2_des-dist_' interval '_v2.1.0.cdf'');'];
c_eval(tmp,ic);
c_eval('diste = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_dist_brst''));',ic);
c_eval('energy0 = get_variable(tmpDataObj,''mms?_des_energy0_brst'');',ic);
energy0 = energy0.data;
c_eval('energy1 = get_variable(tmpDataObj,''mms?_des_energy1_brst'');',ic);
energy1 = energy1.data;
c_eval('stepTable = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_steptable_parity_brst''));',ic);
c_eval('phi = mms.variable2ts(get_variable(tmpDataObj,''mms?_des_phi_brst''));',ic);
c_eval('theta = get_variable(tmpDataObj,''mms?_des_theta_brst'');',ic);

tmp = ['tmpDataObjB = dataobj(''/Volumes/DansHD2/data/mms?/dfg/brst/l2pre/' yyyy '/' mm '/' dd '/mms?_dfg_brst_l2pre_' interval '_v3.11.0.cdf'');'];
c_eval(tmp,ic);
c_eval('Bxyz = mms.variable2ts(get_variable(tmpDataObjB,''mms?_dfg_brst_l2pre_dmpa''));',ic);

tmp = ['tmpDataObj = dataobj(''/Volumes/DansHD2/data/mms?/edp/fast/l2/scpot/' yyyy '/' mm '/mms?_edp_fast_l2_scpot_' yyyy mm dd '000000_v1.0.0.cdf'');'];
c_eval(tmp,ic);
c_eval('SCpot = mms.variable2ts(get_variable(tmpDataObj,''mms?_edp_psp''));',ic);
offset1 = 1.3; offset2 = 1.5; offset3 = 1.2; offset4 = 0.0; %For v1 data
c_eval('SCpot.data = -SCpot.data*1.2+offset?;',ic);

tint = irf.tint('2015-10-30T05:15:40.00Z/2015-10-30T05:15:50.00Z');

diste.data = diste.data*1e12;

%% Compute PAD
[paddist,thetap,energymat,tint] = mms.get_pitchangledist(diste,phi,theta,stepTable,energy0,energy1,Bxyz,tint);
SCpot = SCpot.resample(paddist);

SCpotmat = SCpot.data*ones(1,32);
energymat = energymat-SCpotmat;

PSDpar = squeeze(paddist.data(:,:,1));
PSDapar = squeeze(paddist.data(:,:,12));
PSDperp = squeeze(paddist.data(:,:,[6 7]));
PSDperp = irf.nanmean(PSDperp,3);

%% Find fits and estimate
Units = irf_units; % Use IAU and CODATA values for fundamental constants.
qe = Units.e;
me = Units.me;

enrange = [3:15];
%enrange2 = [3:15];


phipar = zeros(length(paddist.time),1);
Teperp = zeros(length(paddist.time),1);
neback = zeros(length(paddist.time),1);

for ii=1:length(paddist.time)

PSDperps = double(log10(PSDperp(ii,enrange)));
PSDpars = double(log10(PSDpar(ii,enrange)));
PSDapars = double(log10(PSDapar(ii,enrange)));
vvec =  double(sqrt(2*qe*energymat(ii,enrange)/me));

guess = [5e6 sqrt(2*qe*50/me)];

fun = fittype('log10(a/(sqrt(pi)*b)^3*exp(-x^2/b^2))');
options = fitoptions(fun);
options.StartPoint = guess;
options.Lower = [1e6 1e6];
options.Upper = [5e7 1e7];
fit1 = fit(vvec',PSDperps',fun,options);


nefit = fit1.a;
Tefit = me*fit1.b^2/(2*qe);

phipar(ii) = fitparpotential(nefit,Tefit,PSDpars,PSDapars,double(energymat(ii,enrange)));
Teperp(ii) = Tefit;
phipar(ii)
neback(ii) = nefit/1e6;
neback(ii);

end

phipar = TSeries(paddist.time,phipar);
Teperp = TSeries(paddist.time,Teperp);
neback = TSeries(paddist.time,neback);