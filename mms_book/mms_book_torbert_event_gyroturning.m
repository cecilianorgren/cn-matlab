%% Load data
ic = [3];
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');

units = irf_units;

% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);

%c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

% Load spacecraft position
disp('Loading spacecraft position...')
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

% Particle moments
% Skymap distributions
if 1
  %%
  ic_ = ic;
  %ic = 3;
disp('Loading skymaps...')
c_eval('ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?);',ic)
%c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
%c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
%c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
%c_eval('iPDist?_onecount = iPDist?; iPDist?_onecount.data = (iPDist?_onecount.data./iPDistErr?.data).^2;',ic)
ic = ic_;
end
% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',1:4) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',1:4);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% Velocity
disp('Loading bulk velocities...'); tic
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',1:4)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc

% Calculate stress tensor from pressure, density and speed
% P_psd(:,1,1) = P_psd(:,1,1)-pmass*n_psd.*V_psd(:,1).*V_psd(:,1);
%
comps = ['x','y','z'];
c_eval('Se? = gsePe?.data*0;',1:4);
%c_eval('Si? = gsePi?.data*0;',1:4);
c_eval('nvve? = gsePe?.data*0;',1:4);
%c_eval('nvvi? = gsePi?.data*0;',1:4);
for ic1 = 1:3
  for ic2 = 1:3
    c_eval('nvve?(:,ic1,ic2) = units.me*1e6*1e3*1e3*ne?.data.*gseVe?.data(:,ic1).*gseVe?.data(:,ic2);');
    %c_eval('nvvi?(:,ic1,ic2) = units.mp*1e6*1e3*1e3*ni?.data.*gseVi?.data(:,ic1).*gseVi?.data(:,ic2);',ic);    
    c_eval('gseNVVe? = irf.ts_tensor_xyz(gsePe?.time,nvve?*1e9); gseNVVe?.name = ''Inertia''; gseNVVe?.units = ''nPa'';')
    %c_eval('gseNVVi? = irf.ts_tensor_xyz(gsePi?.time,nvvi?*1e9); gseNVVi?.name = ''Iertia''; gseNVVi?.units = ''nPa'';',ic)
%    
    c_eval('Se?(:,ic1,ic2) = gsePe?.data(:,ic1,ic2)*1e-9 + units.me*1e6*1e3*1e3*ne?.data.*gseVe?.data(:,ic1).*gseVe?.data(:,ic2);');
    %c_eval('Si?(:,ic1,ic2) = gsePi?.data(:,ic1,ic2)*1e-9 + units.mp*1e6*1e3*1e3*ni?.data.*gseVi?.data(:,ic1).*gseVi?.data(:,ic2);',ic);
    c_eval('gseSe? = irf.ts_tensor_xyz(gsePe?.time,Se?*1e9); gseSe?.name = ''Stress tensor''; gseSe?.units = ''nPa'';')
    %c_eval('gseSi? = irf.ts_tensor_xyz(gsePi?.time,Si?*1e9); gseSi?.name = ''Stress tensor''; gseSi?.units = ''nPa'';',ic)
  end
end
%

c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('gseVexB? = cross(gseVe?,gseB?.resample(gseVe?.time))*1e-3; gseVExB?.units = ''mV/m'';',ic) % km/s

% EDR signatures
c_eval('facPepp? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''pp'');',ic); % Peperp1 = Peperp2
c_eval('facPeqq? = mms.rotate_tensor(gsePe?,''fac'',gseB?,''qq'');',ic); % Peperp1 and Peperp2 are most unequal
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic);
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic);

c_eval('PB? = 1e9*1e-18*gseB?.abs2/2/units.mu0;',ic)

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs;',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

c_eval('Te?par = facTe?.xx;',ic);
c_eval('Te?perp = (facTe?.yy+facTe?.yy)/2;',ic);

c_eval('vte? = sqrt(2*units.eV*(facTe?.trace)/3/units.me)*1e-3;',ic)

if 0 % Gradient from 4 sc
  %%
  ne = (ne1.resample(ne1) + ne2.resample(ne1) + ne3.resample(ne1) + ne4.resample(ne1))/4; % gseE1 not there?
  
  nsm_pre = 1;
  c_eval('R? = gseR?.resample(gsePe!);',1:4,ic)
  c_eval('Pxy? = gsePe?.xy.resample(gsePe!).smooth(nsm_pre);',1:4,ic)
  gradPexy_4sc = c_4_grad('R?','Pxy?','grad');
  gradPexy_ne_4sc = gradPexy_4sc*(1e-9*1e-3)/(ne*1e6)/units.e*1e3; % mV/m
  gradPexy_ne_4sc.name = 'div Pexy/ne';

  c_eval('Sxy? = gseSe?.xy.resample(gseSe!).smooth(nsm_pre);',1:4,ic)
  gradSexy_4sc = c_4_grad('R?','Sxy?','grad');
  gradSexy_ne_4sc = gradSexy_4sc*(1e-9*1e-3)/(ne*1e6)/units.e*1e3; % mV/m
  gradSexy_ne_4sc.name = 'div Sexy/ne';


  % Gradient terms
  %c_eval('gseGrad_Se? = irf.ts_scalar(gseSe?.time,gradient(gseSe?.xy.data,gseSe?.time-gseSe?.time(1))); gseGrad_Se?.name = ''grad(n*vx*vy + Pxy)'';',ic)
  %c_eval('gseGrad_NVve? = irf.ts_scalar(gseNVVe?.time,gradient(gseNVVe?.xy.data,gseSe?.time-gseSe?.time(1))); gseGrad_NVve?.name = ''grad(n*vx*vy)'';',ic)
  %c_eval('gseNVgradVe? = ne?*gseVe?.x*irf.ts_scalar(gseVe?.time,gradient(gseVe?.y.data,gseSe?.time-gseSe?.time(1))); gseNVgradVe?.name = ''n*vx*grad(vy)'';',ic)
  %c_eval('gseVgradNVe? = gseVe?.y*irf.ts_scalar(gseVe?.time,gradient(ne?.data.*gseVe?.x.data,gseSe?.time-gseSe?.time(1))); gseVgradNVe?.name = ''n*vy*grad(n*vx)'';',ic)
end

% Rotate things into new coordinate system 
% Nakamura?
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn = [L;M;N];

% Torbert
%lmn = [0.971, 0.216, -0.106; -0.234, 0.948, -0.215; 0.054, 0.233, 0.971]; % gse


c_eval('tsLgse? = irf.ts_vec_xyz(ePDist?.time,repmat(L,ePDist?.length,1));',ic)
c_eval('tsMgse? = irf.ts_vec_xyz(ePDist?.time,repmat(M,ePDist?.length,1));',ic)
c_eval('tsNgse? = irf.ts_vec_xyz(ePDist?.time,repmat(N,ePDist?.length,1));',ic)

c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)
c_eval('tsLdsl? = mms_dsl2gse(tsLgse?,defatt?,-1);',ic)
c_eval('tsMdsl? = mms_dsl2gse(tsMgse?,defatt?,-1);',ic)
c_eval('tsNdsl? = mms_dsl2gse(tsNgse?,defatt?,-1);',ic)

doDefatt = 0;
if doDefatt
tt = irf_time('2017-07-11T22:34:02.00Z','utc>EpochTT'); %20151112071854
tint_defatt = tt + [-2 2];
zra = mms.db_get_variable('mms2_ancillary_defatt','zra',tint_defatt);
zdec = mms.db_get_variable('mms2_ancillary_defatt','zdec',tint_defatt);
defatt = zra; defatt.zdec = zdec.zdec;
tsL = irf.ts_vec_xyz(tt,L);
tsM = irf.ts_vec_xyz(tt,M);
tsN = irf.ts_vec_xyz(tt,N);
L_dsl = mms_dsl2gse(tsL, defatt, -1); L_dsl = L_dsl.data;
M_dsl = mms_dsl2gse(tsM, defatt, -1); M_dsl = M_dsl.data;
N_dsl = mms_dsl2gse(tsN, defatt, -1); N_dsl = N_dsl.data;
lmn_dsl = [L_dsl; M_dsl; N_dsl];
end

% Rotate data
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)

c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''VExB LMN'';',ic)
c_eval('mvaVexB? = gseVexB?*lmn''; mvaVexB?.name = ''VExB LMN'';',ic)
%c_eval('mvaJ? = gseJ?*lmn''; mvaJ?.name = ''J LMN'';',ic)
%c_eval('mvaJe? = gseJe?*lmn''; mvaJe?.name = ''Je LMN'';',ic)
%c_eval('mvaJi? = gseJi?*lmn''; mvaJi?.name = ''Ji LMN'';',ic)
%mvaJcurl = gseJcurl*lmn'; mvaJcurl.name = 'J LMN CURL';
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;')
c_eval('mvaTi? = lmn*gseTi?*lmn''; mvaTi?.units = gseTi?.units;',ic)
c_eval('mvaTe? = lmn*gseTe?*lmn''; mvaTe?.units = gseTe?.units;',ic)
%
try
c_eval('mvaNVVi? = lmn*gseNVVi?*lmn''; mvaNVVi?.units = gseNVVi?.units;',ic)
end
c_eval('mvaNVVe? = lmn*gseNVVe?*lmn''; mvaNVVe?.units = gseNVVe?.units;',ic)

%c_eval('mvaSi? = lmn*gseSi?*lmn''; mvaSi?.units = gseSi?.units;',ic)
c_eval('mvaSe? = lmn*gseSe?*lmn''; mvaSe?.units = gseSe?.units;',ic)
c_eval('mvaPe? = lmn*gsePe?*lmn''; mvaPe?.units = gsePe?.units;',ic)

try
mvaGradPexy = gradPexy*lmn''; mvaGradPexy.units = mvaGradPexy.units;
mvaGradSexy = gradSexy*lmn''; mvaGradSexy.units = mvaGradSexy.units;
end


if doDefatt*0
  c_eval('mvaVe?_dsl = mms_dsl2gse(mvaVe?,defatt,-1); mvaVe?_dsl.name = ''Ve LMN dsl'';',ic)
  c_eval('mvaE?_dsl = mms_dsl2gse(mvaE?,defatt,-1); mvaE?_dsl.name = ''E LMN dsl'';',ic)
  c_eval('mvaB?_dsl = mms_dsl2gse(mvaB?,defatt,-1); mvaB?_dsl.name = ''B LMN dsl'';',ic)
  c_eval('mvaE?_dsl2 = dslE?*lmn_dsl''; mvaE?_dsl2.name = ''E LMN DSL2'';',ic)
elseif doDefatt
  c_eval('mvaVe?_dsl = mms_dsl2gse(mvaVe?,defatt,-1); mvaVe?_dsl.name = ''Ve LMN dsl'';',ic)
  c_eval('mvaE?_dsl = mms_dsl2gse(mvaE?,defatt,-1); mvaE?_dsl.name = ''E LMN dsl'';',ic)
  c_eval('mvaB?_dsl = mms_dsl2gse(mvaB?,defatt,-1); mvaB?_dsl.name = ''B LMN dsl'';',ic)
  c_eval('mvaE?_dsl2 = dslE?*lmn_dsl''; mvaE?_dsl2.name = ''E LMN DSL2'';',ic)
end

%mvaGradSe_4sc = gradPexy_4sc*lmn'; mvaGradSe_4sc.units = gradSexy_4sc.units;
%mvaGradPe_4sc = gradPexy_4sc*lmn'; mvaGradPe_4sc.units = gradPexy_4sc.units;


%vsc = 10e3;
c_eval('tt = mvaVe?.time-mvaVe?.time(1);',ic)
%tt = tt*vsc;
c_eval('mvaGrad_Se? = irf.ts_scalar(mvaSe?.time,gradient(mvaSe?.xy.data,tt)); mvaGrad_Se?.name = ''grad(n*vx*vy + Pxy)'';',ic)
c_eval('mvaGrad_Pe? = irf.ts_scalar(mvaPe?.time,gradient(mvaPe?.xy.data,tt)); mvaGrad_Pe?.name = ''grad(Pxy)'';',ic)
c_eval('mvaGrad_NVve? = irf.ts_scalar(mvaNVVe?.time,gradient(mvaNVVe?.xy.data,tt)); mvaGrad_NVve?.name = ''grad(n*vx*vy)'';',ic)
c_eval('mvaNVgradVe? = units.me*ne?*1e6*gseVe?.x*1e3*irf.ts_scalar(mvaVe?.time,gradient(mvaVe?.y.data*1e3,tt))*1e9; mvaNVgradVe?.name = ''n*vx*grad(vy,x)'';',ic)
c_eval('mvaVgradNVe? = units.me*gseVe?.y*1e3*irf.ts_scalar(mvaVe?.time,gradient(ne?.data*1e6.*mvaVe?.x.data*1e3,tt))*1e9; mvaVgradNVe?.name = ''vy*grad(n*vx)'';',ic)
c_eval('mvaGradNVe? = irf.ts_scalar(mvaVe?.time,gradient(ne?.data*1e6.*mvaVe?.x.data*1e3,tt)); mvaGradNVe?.name = ''grad_x(n*vx)'';',ic)

np_smooth = 30;
c_eval('mvaGrad_Se?_smooth = irf.ts_scalar(mvaSe?.time,gradient(mvaSe?.smooth(np_smooth).xy.data)); mvaGrad_Se?.name = ''grad(n*vx*vy + Pxy)'';',ic)
c_eval('mvaGrad_NVve?_smooth = irf.ts_scalar(mvaNVVe?.time,gradient(mvaNVVe?.smooth(np_smooth).xy.data)); mvaGrad_NVve?.name = ''grad(n*vx*vy)'';',ic)
c_eval('mvaNVgradVe?_smooth = units.me*ne?*1e6*gseVe?.x*1e3*irf.ts_scalar(mvaVe?.time,gradient(mvaVe?.smooth(np_smooth).y.data*1e3))*1e9; mvaNVgradVe?.name = ''n*vx*grad(vy,x)'';',ic)


c_eval('mvaR? = gseR?*lmn''; mvaR?.name = ''R LMN'';',1:4)

%% Plot: Short overview and 2D distribution f(vx,vy) with locations shown
ic = 3;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
%times = times + 0.25;
times = times + 0.30*1;
%times = times + 3;

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z';...
  '2017-07-11T22:34:07.300000000Z']);
times = times + 0.30*1;
%times = times + 3;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;

%[h1,h2] = initialize_combined_plot('topbottom',2,2,times.length/2,0.5,'horizontal');
[h1,h2] = initialize_combined_plot('topbottom',2,1,times.length,0.5,'horizontal');


comps = 'LM';
comps_ts = 'LM';
comps_str = ["L","M"];
%fig = gcf; fig.Position = [343     1   938   496];

isub = 1;
hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{mvaB3.x,mvaB3.y,mvaB3.z},'comp')
hca.YLabel.String = sprintf('B (nT)');
hca.YLabel.String = {'B','(nT)'};
hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.98])


if 0 % E
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaE3.x,mvaE3.y,mvaE3.z},'comp')
  hca.YLabel.String = sprintf('E (mV/m)');
  hca.YLabel.Interpreter = 'tex';
end

hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{mvaVe3.x*1e-3,mvaVe3.y*1e-3,mvaVe3.z*1e-3},'comp')
hca.YLabel.String = {'v_{e}','(10^3 km/s)'};
hca.YLabel.Interpreter = 'tex';
set(hca,'ColorOrder',mms_colors('xyz'))
irf_legend(hca,{'L','M','N'},[0.98 0.5])

irf_zoom(h1,'x',irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:09.000Z'))
irf_zoom(h1,'y')
irf_plot_axis_align(h1)
c_eval('xtickangle(h1(?),0)',1:numel(h1))
%c_eval('h1(?).Position(1) = h1(?).Position(1) + 0.1; h1(?).Position(3) = h1(?).Position(3) - 0.2;',1:numel(h1))



fontsize_B_amp = 13;
markersize = 5;
contour_levels = 10.^([-14 -13 -12 -11 -10 -9]+0.5);

times_exact = {};
isub = 1;

% f
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('E = mvaE?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  

  c_eval('B = mvaB?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic) 
  c_eval('ve = mvaVe?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)   


  t_dist_center = dist.time.start + (dist.time.stop-dist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data;tsMdsl?.resample(t_dist_center).data;tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';
  
  vdf = dist.reduce('2D',Ldsl,Mdsl,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'contour',contour_levels);
  hc_.Colorbar.YLabel.String = sprintf('log_{10}f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  %hc_.Colorbar.YLabel.String = {sprintf('log_{10}f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));

  if 0 % plot E+vxB=0
    %%
    hold(hca,'on')
    B__ = mean(B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]).data,1);
    E__ = mean(E.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]).data,1);
    vlimit_x = (E__(2)*1e-3)/(B__(3)*1e-9)*1e-6;
    vlimit_y = (E__(2)*1e-3)/(B__(1)*1e-9)*1e-6;
    plot(hca,hca.XLim,vlimit_y*[1 1],'r--')
    plot(hca,vlimit_x*[1 1],hca.YLim,'b--')    
    hold(hca,'off')
  end
  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B__ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B__.data,1);
    B_std = std(B__.data,1);
    b = B_/norm(B_);
    B_std_inplane = std(B__.data(:,1:2),1);
    B_inplane = sqrt(sum(B_(1:2).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(2)/b(1);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);
  end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.02 0.07],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hbulk = plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+w','MarkerFaceColor','w','markersize',markersize,'linewidth',2);
    hold(hca,'off')    
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1:times.length,1:numel(h1))

c_eval('h1(?).Position(3) = -1*h1(?).Position(1) + h2(end).Position(1) + h2(end).Position(3); ',1:numel(h1))

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
%hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);

%ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(1:numel(h2)-1))
%ih = nt;
%hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
%hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%ih = nt;
%hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
%hb(ih).Position(4) = hb(ih).Position(4)*0.9; 

c_eval('h2(?).YLabel.String = [];',2:numel(h2))
c_eval('h2(?).YTickLabel = [];',2:numel(h2))
c_eval('h2(?).XTickLabelRotation = 0;',1:numel(h2))

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).XLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).YLim = 0.99*[-60 60];',1:numel(h2))
%c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)
c_eval('h2(?).Layer = ''top'';',1:numel(h2))

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';


c_eval('h2(?).Color = 1*[1 1 1];',1:numel(h2))
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:numel(h2))
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:numel(h2))


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
nInd = 1;
for ii = 1:numel(h1)
  irf_legend(h1(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h1(ii).FontSize = 14;
end
for ii = 1:numel(h2)
  irf_legend(h2(ii),legends{nInd},[0.04 0.97],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h2(ii).FontSize = 14;
end


c_eval('h1(?).Position(1) = h1(?).Position(1) + 0.1; h1(?).Position(3) = h1(?).Position(3) - 0.2;',1:numel(h1))
%c_eval('h1(?).Position(1) = h2(1).Position(1);',1:numel(h1))

c_eval('h2(?).YDir = ''reverse'';',1:numel(h2))
hb = findobj(gcf,'type','colorbar');
hb.Location = 'northoutside';
hb.Position = [ 0.6355    0.4495    0.1113    0.0227];
hb.Ticks = [-15 -13 -11];

h2(3).Title.String = 'X line';
h2(6).Title.String = 'Gyrotropic exhaust';
h2(1).Title.String = 'Towards inflow';
h2(2).Title.String = '<- Tailward';
h2(4).Title.String = 'Earthward ->';
c_eval('h2(?).Title.FontWeight = ''light'';',1:numel(h2))

%% Plot: Short overview and 3D distribution with locations shown
ic = 3;
%t1 =  irf_time('2017-07-11T22:34:01.30Z','utc>EpochTT');
%t2 =  irf_time('2017-07-11T22:34:03.30Z','utc>EpochTT');
times = irf_time(['2017-07-11T22:34:01.30Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.03;
times = irf_time(['2017-07-11T22:34:01.00Z';'2017-07-11T22:34:03.30Z'],'utc>EpochTT')+0.4;
times = times(1):0.5:times(2);

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
%times = times + 0.25;
times = times + 0.30*1;
%times = times + 3;

times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z';...
  '2017-07-11T22:34:07.300000000Z']);
times = times + 0.30*1;
%times = times + 3;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;

%[h1,h2] = initialize_combined_plot('topbottom',2,2,times.length/2,0.5,'horizontal');
[h1,h2] = initialize_combined_plot('topbottom',2,1,times.length,0.5,'horizontal');


comps = 'xy';
comps_ts = 'xy';
comps_str = ["x","y"];
%fig = gcf; fig.Position = [343     1   938   496];

isub = 1;
hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{mvaB3.x,mvaB3.y,mvaB3.z},'comp')
hca.YLabel.String = sprintf('B (nT)');
hca.YLabel.String = {'B','(nT)'};
hca.YLabel.Interpreter = 'tex';


if 0 % E
  hca = h1(isub); isub = isub + 1;
  set(hca,'ColorOrder',mms_colors('xyz'))
  irf_plot(hca,{mvaE3.x,mvaE3.y,mvaE3.z},'comp')
  hca.YLabel.String = sprintf('E (mV/m)');
  hca.YLabel.Interpreter = 'tex';
end

hca = h1(isub); isub = isub + 1;
set(hca,'ColorOrder',mms_colors('xyz'))
irf_plot(hca,{mvaVe3.x*1e-3,mvaVe3.y*1e-3,mvaVe3.z*1e-3},'comp')
hca.YLabel.String = {'v_{e}','(10^3 km/s)'};
hca.YLabel.Interpreter = 'tex';

irf_zoom(h1,'x',irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:09.000Z'))
irf_zoom(h1,'y')
irf_plot_axis_align(h1)
c_eval('xtickangle(h1(?),0)',1:numel(h1))
%c_eval('h1(?).Position(1) = h1(?).Position(1) + 0.1; h1(?).Position(3) = h1(?).Position(3) - 0.2;',1:numel(h1))



fontsize_B_amp = 13;
markersize = 5;
contour_levels = 10.^([-14 -13 -12 -11 -10 -9]+0.5);

times_exact = {};
isub = 1;

% f
for itime = 1:times.length
  hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('E = mvaE?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  

  c_eval('B = mvaB?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic) 
  c_eval('ve = mvaVe?.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);',ic)   


  t_dist_center = dist.time.start + (dist.time.stop-dist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data;tsMdsl?.resample(t_dist_center).data;tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';
  
  iso_values = 1e-31;
  nSmooth = 1;
  hs = dist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',Tdsl);  
  
  hs.Patch(1).FaceAlpha = 1;
  view(hca,[0.2 -1 0.5])
  camlight(hca)
  
  %camlight%(hca,0,-45)
  lighting(hca,'gouraud')

  set(hca,'colororder',[0.5 1 0.5; 1 0.5 0.5])
  irf_legend(hca,{sprintf('f = %g',Flev),sprintf('f = %g',Flev2)}',[0.98 0.98],'backgroundcolor',[1 1 1],'fontsize',12)




  times_exact{itime} = dist.time;
    
end
%%
c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1:times.length,1:numel(h1))

c_eval('h1(?).Position(3) = -1*h1(?).Position(1) + h2(end).Position(1) + h2(end).Position(3); ',1:numel(h1))

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
%hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

h2(1).XLim = 0.99*vg([1 end])*1e-3;
h2(1).YLim = 0.99*vg([1 end])*1e-3;

c_eval('h2(?).XTick = -60:20:60; h2(?).YTick = -60:20:60;',1:numel(h2))

c_eval('h2(?).LineWidth = 1;',1:numel(h2))
c_eval('h2(?).FontSize = 14;',1:numel(h2))
c_eval('axis(h2(?),''square'');',1:numel(h2))

c_eval('h1(?).LineWidth = 1;',1:numel(h1))
c_eval('h1(?).FontSize = 14;',1:numel(h1))

%compact_panels(h2,0.0,00)
compact_panels(h2,0.005,00.005)
hb = findobj(gcf,'type','colorbar'); 
c_eval('hb(?).FontSize = 14;',1:numel(h2))
hb = hb(end:-1:1);

%ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(1:numel(h2)-1))
%ih = nt;
%hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
%hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%ih = nt;
%hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
%hb(ih).Position(4) = hb(ih).Position(4)*0.9; 

c_eval('h2(?).YLabel.String = [];',2:numel(h2))
c_eval('h2(?).YTickLabel = [];',2:numel(h2))
c_eval('h2(?).XTickLabelRotation = 0;',1:numel(h2))

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).XLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).YLim = 0.99*[-60 60];',1:numel(h2))
%c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)
c_eval('h2(?).Layer = ''top'';',1:numel(h2))

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';


c_eval('h2(?).Color = 1*[1 1 1];',1:numel(h2))
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:numel(h2))
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:numel(h2))


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
nInd = 1;
for ii = 1:numel(h1)
  irf_legend(h1(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h1(ii).FontSize = 14;
end
for ii = 1:numel(h2)
  irf_legend(h2(ii),legends{nInd},[0.04 0.97],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h2(ii).FontSize = 14;
end


c_eval('h1(?).Position(1) = h1(?).Position(1) + 0.1; h1(?).Position(3) = h1(?).Position(3) - 0.2;',1:numel(h1))
