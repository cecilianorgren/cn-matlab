% US-Japan
%% Load data
ic = [3];
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
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
%
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

% Compute Q and Dng from facPepp
c_eval('Q? = (facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2)./(facPepp?.yy.data.^2+2*facPepp?.yy.data.*facPepp?.xx.data);',ic);
c_eval('Q? = irf.ts_scalar(ne?.time,sqrt(Q?));',ic);
c_eval('Dng? = sqrt(8*(facPepp?.xy.data.^2+facPepp?.xz.data.^2+facPepp?.yz.data.^2))./(facPepp?.xx.data+2*facPepp?.yy.data);',ic);
c_eval('Dng? = irf.ts_scalar(ne?.time,Dng?);',ic);

% Compute agyrotropy Aphi from facPeqq
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs',ic);

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
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn = [L;M;N];

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

%% Spacecraft constellation
tt = irf_time('2017-07-11T22:33:58.000000000Z','utc>EpochTT');
R0 = (mvaR1 + mvaR2.resample(mvaR1) + mvaR3.resample(mvaR1) + mvaR4.resample(mvaR1))/4;
c_eval('r? = mvaR?-R0;',1:4)
c_eval('r? = r?.resample(tt).data;',1:4)
%c_eval('dr(?,!) = dr? - dr!;',1:4,1:4)
r = [r1;r2;r3;r4];
rmax = max(r);
rmin = min(r);

h = subplot(1,1,1);
isub = 1;
hca = h(isub);
symbols = {'o','s','^','<'};
colors = mms_colors('1234');
colors(1,:) = [0.2 0.2 0.2];

colors = [0.2000    0.2000    0.2000;
          0.9000    0.2000         0;
               0    0.8000         0;
          0.1000    0.4000    0.9000];
linewidth = 2;
markersize = 20;

holdon = 0;
for isc = 1:4
  plot3(hca,r(isc,1),r(isc,2),r(isc,3),'s','markersize',markersize,'linewidth',linewidth,'color',colors(isc,:));
  if not(holdon)
    hold(hca,'on')
  end
end

hca.XLim = [-15 15];
hca.YLim = [-15 15];
hca.ZLim = [-15 15];

hca.XLim = [rmin(1) rmax(1)]+[-1 1]*2;
hca.YLim = [rmin(2) rmax(2)]+[-1 1]*2;
hca.ZLim = [rmin(3) rmax(3)]+[-1 1]*2;

xmin = hca.XLim(2);
ymin = hca.YLim(2);
zmin = hca.ZLim(1);
for isc = 1:4
  plot3(hca,xmin*[1 1],r(isc,2),r(isc,3),'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
  plot3(hca,r(isc,1),ymin*[1 1],r(isc,3),'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
  plot3(hca,r(isc,1),r(isc,2),zmin*[1 1],'s','markersize',10,'linewidth',1,'color',colors(isc,:).^0.125);  
end

for i1 = 1:4
  for i2 = i1:4
    plot3(hca,[r(i1,1) r(i2,1)],[r(i1,2) r(i2,2)],[r(i1,3) r(i2,3)],'linewidth',1,'color',[0 0 0]);
  end
end
for i1 = 1:4
  for i2 = i1:4
    plot3(hca,[xmin xmin],[r(i1,2) r(i2,2)],[r(i1,3) r(i2,3)],':','linewidth',1,'color',[1 1 1]*0.7);
    plot3(hca,[r(i1,1) r(i2,1)],[ymin ymin],[r(i1,3) r(i2,3)],':','linewidth',1,'color',[1 1 1]*0.7);
    plot3(hca,[r(i1,1) r(i2,1)],[r(i1,2) r(i2,2)],[zmin zmin],':','linewidth',1,'color',[1 1 1]*0.7);
  end
end

hold(hca,'off')

hca.XLabel.String = 'L (km)';
hca.YLabel.String = 'M (km)';
hca.ZLabel.String = 'N (km)';

c_eval('h(?).FontSize = 14;',1:numel(h))

hca.Box = 'on';
%axis(hca,'square');

%daspect(hca,[1 1 1])
hca.DataAspectRatio = [1 1 1];
%hca.XLim = [rmin(1) rmax(1)]+[-1 1]*2;
%hca.YLim = [rmin(2) rmax(2)]+[-1 1]*2;
%hca.ZLim = [rmin(3) rmax(3)]+[-1 1]*2;

%% Figure: Overview 1
ic = 3;

tint_edr = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:09.00Z'); %20151112071854

npanels = 7;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
fontsize = 12;
comps = 'xyz';
ylim_stress = [-4 3];

isub = 0;
zoomy = [];
leg_loc = [0.98 0.1];
%leg_loc = [1.02 0.9];
nsm = 1;

if 1 % B LMN
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},leg_loc,'fontsize',fontsize);
end 
if 1 % E LMN
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  fhigh = 10;
  c_eval('irf_plot(hca,{mvaE?.x.filt(0,fhigh,[],3),mvaE?.y.filt(0,fhigh,[],3),mvaE?.z.filt(0,fhigh,[],3)},''comp'');',ic)
  hca.YLabel.String = {'E (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},leg_loc,'fontsize',fontsize);
  irf_legend(hca,{sprintf('f<%g Hz',fhigh)},[0.98 0.98],'fontsize',fontsize,'color',[0 0 0]);
end 
if 1 % Tpar, Tperp
  hca = irf_panel('Tepar, Teperp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{Te?par,Te?perp},''comp'');',ic)  
  hca.YLabel.String = {'T_e (eV)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'T_{e,||}','T_{e,\perp}'},[0.98 0.15],'fontsize',fontsize);
end
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint)*1e-3,mvaVe?.y.tlim(tint)*1e-3,mvaVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  c_eval('irf_plot(hca,{-1*vte?.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  
  hca.YLabel.String = {'u_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{['u_' comps(1)],['u_' comps(2)],['u_' comps(3)],'-v_{te}'},[0.98,0.3],'fontsize',fontsize);
end
if 1 % NVVe-off LMN
  hca = irf_panel('mnvve off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaNVVe?.xy.tlim(tint).smooth(nsm)*1e3,mvaNVVe?.xz.tlim(tint).smooth(nsm)*1e3,mvaNVVe?.yz.tlim(tint).smooth(nsm)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'m_enu_eu_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1:2),comps([1 3]),comps([2 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % Pe-diag LMN
  hca = irf_panel('Pe diag LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xx.tlim(tint)*1e3,mvaPe?.yy.tlim(tint)*1e3,mvaPe?.zz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps([1 1]),comps([2 2]),comps([3 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LL','MM','NN'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Pe-off LMN
  hca = irf_panel('Pe off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xy.tlim(tint).smooth(nsm)*1e3,mvaPe?.xz.tlim(tint).smooth(nsm)*1e3,mvaPe?.yz.tlim(tint).smooth(nsm)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1:2),comps([1 3]),comps([2 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 1 % Se-off LMN
  hca = irf_panel('Se off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaSe?.xy.tlim(tint).smooth(nsm)*1e3,mvaSe?.xz.tlim(tint).smooth(nsm)*1e3,mvaSe?.yz.tlim(tint).smooth(nsm)*1e3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %c_eval('irf_plot(hca,{mvaSe?.xy.tlim(tint).smooth(nsm)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'X_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1:2),comps([1 3]),comps([2 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % Non-gyro
  hca = irf_panel('Non-gyrotropies');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{agyro?,Dng?,Q?.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'','Non-gyrotropy'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'A_{\phi}','Dng','Q^{1/2}'},leg_loc,'fontsize',fontsize);
end
if 0 % gradients, comparison  
  hca = irf_panel('gradients 1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaGrad_Se?.tlim(tint),mvaGrad_NVve?.tlim(tint),mvaNVgradVe?.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'Gradients (pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'grad(n*vx*vy+Pxy)','grad(n*vx*vy)','n*vx*grad(vy,x)'},leg_loc,'fontsize',fontsize);
end
if 0 % gradients, comparison    
  hca = irf_panel('gradients 1, smoothed after');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaGrad_Se?.tlim(tint).smooth(np_smooth),mvaGrad_NVve?.tlim(tint).smooth(np_smooth),mvaNVgradVe?.tlim(tint).smooth(np_smooth),mvaVgradNVe?.tlim(tint).smooth(np_smooth)},''comp'');',ic)  
  hca.YLabel.String = {'(pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\partial_x(mnv_xv_y+P_{xy})'},[0.01 0.15],'fontsize',fontsize);
  set(hca,'ColorOrder',mms_colors('yza'))
  irf_legend(hca,{'\partial_x(mnv_xv_y)','mnv_x\partial_x v_y'},[0.99 0.15],'fontsize',fontsize);
end
if 0 % gradients, comparison , smoothed
  hca = irf_panel('gradients 1 smoothed before');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaGrad_Se?_smooth.tlim(tint),mvaGrad_NVve?_smooth.tlim(tint),mvaNVgradVe?_smooth.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'Gradients (pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\partial_x(nv_xv_y + P_{xy})','\partial_x(nv_xv_y)','nv_x\partial_x v_y'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % Ve x B
  hca = irf_panel('Ve x B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVexB?.x.tlim(tint),mvaVexB?.y.tlim(tint),mvaVexB?.z.tlim(tint)},''comp'');',ic)
  hca.YLabel.String = {'v_e\times{B}','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},leg_loc,'fontsize',fontsize);
end
if 0 % Ve x B
  hca = irf_panel('Ve x B LMN, ycomp');
  set(hca,'ColorOrder',mms_colors('yxz'))
  c_eval('irf_plot(hca,{mvaVexB?.y.tlim(tint),mvaE?.y.tlim(tint).resample(mvaVe?),mvaVexB?.y.tlim(tint)+mvaE?.y.tlim(tint).resample(mvaVexB?)},''comp'');',ic)
  hca.YLabel.String = {'v_e\times{B}_y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('yxz'))
  irf_legend(hca,{'v_e\times{B}','E_y','v_e\times{B}+E_y'},leg_loc,'fontsize',fontsize);
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

time_reversal = irf_time('2017-07-11T22:34:02.641773681Z','utc>EpochTT');
c_eval('hm(?) = irf_pl_mark(h(?),time_reversal);',1:numel(h))
c_eval('hm(?).Color = [0.3 0.3 0.3];',1:numel(h))
c_eval('hm(?).LineStyle = ''--'';',1:numel(h))

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(h,'x',tint_edr)
%irf_zoom(h,'x',tint)
irf_zoom(h,'y')
irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
h(end).XTickLabelRotation = 0;
c_eval('h(?).YLim = ylim_stress;',(4:6)+1)

c_eval('hold(h(?),"on"); plot(h(?),h(?).XLim,[0 0],"color",[0.7 0.7 0.7],"linewidth",1); h(?).Children = circshift(h(?).Children,-1)',1:numel(h))

%% Figure: Gradients
ic = 3;

tint_edr = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:09.00Z'); %20151112071854

npanels = 6;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
fontsize = 14;
fontsize_leg = 11;
comps = 'xyz';

isub = 0;
zoomy = [];
leg_loc = [0.95 0.1];
%leg_loc = [1.02 0.9];
nsm_base = 10;
nsm = 30;

vsc = 130e3;10e3; % m/s, Nakamura 2019, from Lxline and timescale
vsc = 133e3;10e3; % m/s, Hasegawa 2018, 320 km - 2.4 s -> 133 km/s
vsc = 200e3;10e3; % m/s, Egedal 2019, 500 km, 2.5 s
vsc = 170e3; % Torbert
  
  
if 1 % ne
  hca = irf_panel('ne');
  set(hca,'ColorOrder',mms_colors('1'))
  c_eval('irf_plot(hca,{ne?.tlim(tint)},''comp'');',ic)        
  hca.YLabel.String = {'n_e','(cm^{-3})'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'n_e'},leg_loc,'fontsize',fontsize);
end
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint)*1e-3,mvaVe?.y.tlim(tint)*1e-3,mvaVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('b'))
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
    
  %yyaxis(hca,'right');
  %hca2 = gca;
  %c_eval('irf_plot(hca2,{ne?.tlim(tint)},''comp'',''--'');',ic)  
  
  hca.YLabel.String = {'u_e','(10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)},[0.98 0.15],'fontsize',fontsize_leg);
end
if 1 % Se, Pe, nvv, xy
  hca = irf_panel('stress x comp');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaSe?.xy.tlim(tint)*1e3,mvaPe?.xy.tlim(tint)*1e3,mvaNVVe?.xy.tlim(tint)*1e3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaSe?.xy.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'Stress','(pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'X_{xy}','P_{xy}','m_enu_xu_y'},[0.98 0.15],'fontsize',fontsize_leg);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % NVVe-off LMN
  hca = irf_panel('mnvve off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaNVVe?.xy.tlim(tint)*1e3,mvaNVVe?.xz.tlim(tint)*1e3,mvaNVVe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'m_env_ev_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1:2),comps([1 3]),comps([2 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % Pe-off LMN
  hca = irf_panel('Pe off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPe?.xy.tlim(tint)*1e3,mvaPe?.xz.tlim(tint)*1e3,mvaPe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'P_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1:2),comps([1 3]),comps([2 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % Se-off LMN
  hca = irf_panel('Se off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaSe?.xy.tlim(tint)*1e3,mvaSe?.xz.tlim(tint)*1e3,mvaSe?.yz.tlim(tint)*1e3},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaSe?.xy.tlim(tint)*1e3},''comp'');',ic)  
  hca.YLabel.String = {'S_e (pPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1:2),comps([1 3]),comps([2 3])},leg_loc,'fontsize',fontsize);
  %irf_legend(hca,{'LM','LN','MN'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % gradients, comparison
  hca = irf_panel('gradients 1');
  set(hca,'ColorOrder',mms_colors('xyzab'))
  c_eval('data1 = mvaGrad_Se?.tlim(tint);',ic)
  c_eval('data2 = mvaGrad_Pe?.tlim(tint);',ic)
  c_eval('data3 = mvaGrad_NVve?.tlim(tint);',ic)
  c_eval('data4 = mvaNVgradVe?.tlim(tint);',ic)
  c_eval('data5 = mvaVgradNVe?.tlim(tint);',ic)
  irf_plot(hca,{data1,data2,data3,data4,data5},'comp');
  hca.YLabel.String = {'Gradients','(pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyzab'))
  irf_legend(hca,{'\partial_x{S_{xy}}','\partial_x{P_{xy}}','\partial_x(nv_x*v_y)','n*vx*\partial_x{v_y}','v_y\partial(nv_x)'},leg_loc,'fontsize',fontsize);
end
if 0 % gradients, comparison
  nsm = 10;
  hca = irf_panel('gradients 1, smoothed 1');
  set(hca,'ColorOrder',mms_colors('xyzab'))
  c_eval('data1 = mvaGrad_Se?.tlim(tint).smooth(nsm);',ic)
  c_eval('data2 = mvaGrad_Pe?.tlim(tint).smooth(nsm);',ic)
  c_eval('data3 = mvaGrad_NVve?.tlim(tint).smooth(nsm);',ic)
  c_eval('data4 = mvaNVgradVe?.tlim(tint).smooth(nsm);',ic)
  c_eval('data5 = mvaVgradNVe?.tlim(tint).smooth(nsm);',ic)
  irf_plot(hca,{data1,data2,data3,data4,data5},'comp');
  hca.YLabel.String = {'Gradients','(pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyzab'))
  irf_legend(hca,{'\partial_x{S_{xy}}','\partial_x{P_{xy}}','\partial_x(nv_x*v_y)','n*vx*\partial_x{v_y}','v_y\partial(nv_x)'},leg_loc,'fontsize',fontsize);
  irf_legend(hca,{sprintf('N_{smooth} = %g',nsm)},[0.02 0.99],'fontsize',fontsize,'color','k');
end
if 1 % gradients, comparison
  
  hca = irf_panel('gradients 1, smoothed 2');
  set(hca,'ColorOrder',mms_colors('xyzap'))
  c_eval('data1 = mvaGrad_Se?.tlim(tint).smooth(nsm)*1e3;',ic) % nPa -> pPa
  c_eval('data2 = mvaGrad_Pe?.tlim(tint).smooth(nsm)*1e3;',ic)
  c_eval('data3 = mvaGrad_NVve?.tlim(tint).smooth(nsm)*1e3;',ic)
  c_eval('data4 = mvaNVgradVe?.tlim(tint).smooth(nsm)*1e3;',ic)
  c_eval('data5 = mvaVgradNVe?.tlim(tint).smooth(nsm)*1e3;',ic)
  irf_plot(hca,{data1,data2,data3,data4,data5},'comp');
  hca.YLabel.String = {'Gradients','(pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyzap'))
  %irf_legend(hca,{'\partial_x{S_{xy}}','\partial_x{P_{xy}}','\partial_x(nv_x*v_y)','n*vx*\partial_x{v_y}','v_y\partial(nv_x)'},leg_loc,'fontsize',fontsize);
  
  if 0
    set(hca,'ColorOrder',mms_colors('xy'))
    irf_legend(hca,{'\partial_t{X_{xy}}','  \partial_t{P_{xy}}'},[0.02 0.15],'fontsize',fontsize_leg);
    set(hca,'ColorOrder',mms_colors('z'))
    irf_legend(hca,{'m_e\partial_t(nu_xu_y)'},[0.98 0.98],'fontsize',fontsize_leg);  
    set(hca,'ColorOrder',mms_colors('ab'))
    irf_legend(hca,{'m_enu_x\partial_t{u_y}'},[0.98 0.35],'fontsize',fontsize_leg);
    set(hca,'ColorOrder',mms_colors('b'))
    irf_legend(hca,{'m_eu_y\partial_t(nu_x)'},[0.98 0.15],'fontsize',fontsize_leg);
  else
    set(hca,'ColorOrder',mms_colors('xy'))
    irf_legend(hca,{'\Delta X_{xy}/\Delta{t}','  \Delta{P_{xy}}/\Delta{t}'},[0.02 0.15],'fontsize',fontsize_leg);
    set(hca,'ColorOrder',mms_colors('z'))
    irf_legend(hca,{'m_e\Delta(nu_xu_y)/\Delta{t}'},[0.98 0.98],'fontsize',fontsize_leg);  
    set(hca,'ColorOrder',mms_colors('ap'))
    irf_legend(hca,{'m_enu_x\Delta{u_y}/\Delta{t}'},[0.98 0.35],'fontsize',fontsize_leg);
    set(hca,'ColorOrder',mms_colors('p'))
    irf_legend(hca,{'m_eu_y\Delta(nu_x)/\Delta{t}'},[0.98 0.15],'fontsize',fontsize_leg);
  end

  %irf_legend(hca,{'m_env_x\partial_t{v_y}','m_ev_y\partial_t(nv_x)'},[0.98 0.1],'fontsize',fontsize_leg);
  %irf_legend(hca,{sprintf('N_{smooth} = %g',nsm)},[0.5 0.99],'fontsize',fontsize,'color','k');
  %irf_legend(hca,{sprintf('N_{sm}=%g',nsm)},[1.01 0.99],'fontsize',fontsize_leg,'color','k');
  irf_legend(h(1),{sprintf('Gradients: N_{sm}=%g',nsm)},[0.1 1.1],'fontsize',fontsize_leg,'color','k');
end
if 1 % gradients, comparison, mV/m  
  hca = irf_panel('gradients 1, smoothed 2, mV/m');
  set(hca,'ColorOrder',mms_colors('xyzap'))
  c_eval('data1 =   mvaGrad_Se?*1e-9/(ne?*1e6*units.e*vsc);',ic)
  c_eval('data2 =   mvaGrad_Pe?*1e-9/(ne?*1e6*units.e*vsc);',ic)
  c_eval('data3 = mvaGrad_NVve?*1e-9/(ne?*1e6*units.e*vsc);',ic)
  c_eval('data4 =  mvaNVgradVe?*1e-9/(ne?*1e6*units.e*vsc);',ic)
  c_eval('data5 =  mvaVgradNVe?*1e-9/(ne?*1e6*units.e*vsc);',ic)
  
  data1 = data1.tlim(tint).smooth(nsm)*1e3; % mV/m ?
  data2 = data2.tlim(tint).smooth(nsm)*1e3;
  data3 = data3.tlim(tint).smooth(nsm)*1e3;
  data4 = data4.tlim(tint).smooth(nsm)*1e3;
  data5 = data5.tlim(tint).smooth(nsm)*1e3;
  
  irf_plot(hca,{-data1,-data2,-data3,-data4,-data5},'comp');
  hca.YLabel.String = {'Gradients','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyzap'))
  %irf_legend(hca,{'\partial_x{S_{xy}}','\partial_x{P_{xy}}','\partial_x(nv_x*v_y)','n*vx*\partial_x{v_y}','v_y\partial(nv_x)'},leg_loc,'fontsize',fontsize);
  set(hca,'ColorOrder',mms_colors('xy'))
  irf_legend(hca,{'-\partial_x{X_{xy}}/ne','  -\partial_x{P_{xy}}/ne'},[0.02 0.15],'fontsize',fontsize_leg);
  set(hca,'ColorOrder',mms_colors('z'))
  irf_legend(hca,{'-m_e\partial_x(nu_xu_y)/ne'},[0.98 0.98],'fontsize',fontsize_leg);  
  set(hca,'ColorOrder',mms_colors('ap'))
  %irf_legend(hca,{'-m_env_x\partial_x{v_y}/ne','-m_ev_y\partial_x(nv_x)/ne'},leg_loc,'fontsize',fontsize_leg);
  irf_legend(hca,{'-m_enu_x\partial_x{u_y}/ne'},[0.98 0.75],'fontsize',fontsize_leg);
  set(hca,'ColorOrder',mms_colors('p'))
  irf_legend(hca,{'-m_eu_y\partial_x(nu_x)/ne'},[0.98 0.15],'fontsize',fontsize_leg);
  %irf_legend(hca,{sprintf('N_{smooth} = %g',nsm)},[0.5 0.99],'fontsize',fontsize,'color','k');
  set(hca,'ColorOrder',mms_colors('1111'))
  %irf_legend(hca,{sprintf('N_{sm}=%g',nsm),sprintf('v=%g{km/s}',vsc*1e-3)}',[1.01 0.99],'fontsize',14,'color','k');
  
  %irf_legend(hca,{'\partial_x\rightarrow{v_x^{-1}}\partial_t'}',[0.05 0.99],'fontsize',14,'color','k');
  irf_legend(hca,{'\partial_x\rightarrow{v_{x,sc}^{-1}}/\Delta{t}'}',[0.05 0.99],'fontsize',14,'color','k');
  irf_legend(hca,{sprintf('v_{x,sc}=%g {km/s}',vsc*1e-3)}',[0.02 0.75],'fontsize',14,'color','k');

end
if 0 % gradients, comparison    
  hca = irf_panel('gradients 1, smoothed after');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaGrad_Se?.tlim(tint).smooth(np_smooth),mvaGrad_NVve?.tlim(tint).smooth(np_smooth),mvaNVgradVe?.tlim(tint).smooth(np_smooth),mvaVgradNVe?.tlim(tint).smooth(np_smooth)},''comp'');',ic)  
  hca.YLabel.String = {'(pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\partial_x(mnv_xv_y+P_{xy})'},[0.01 0.15],'fontsize',fontsize);
  set(hca,'ColorOrder',mms_colors('yza'))
  irf_legend(hca,{'\partial_x(mnv_xv_y)','mnv_x\partial_x v_y'},[0.99 0.15],'fontsize',fontsize);
end
if 0 % gradients, comparison , smoothed
  hca = irf_panel('gradients 1 smoothed before');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaGrad_Se?_smooth.tlim(tint),mvaGrad_NVve?_smooth.tlim(tint),mvaNVgradVe?_smooth.tlim(tint)},''comp'');',ic)  
  hca.YLabel.String = {'Gradients','(pPa/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\partial_x(nv_xv_y + P_{xy})','\partial_x(nv_xv_y)','nv_x\partial_x v_y'}',[1.02 0.9],'fontsize',fontsize);
end
if 0 % Ve x B
  hca = irf_panel('Ve x B LMN, ycomp');
  set(hca,'ColorOrder',mms_colors('yxz'))
  c_eval('irf_plot(hca,{mvaVexB?.y.tlim(tint),mvaE?.y.tlim(tint).resample(mvaVe?),mvaVexB?.y.tlim(tint)+mvaE?.y.tlim(tint).resample(mvaVexB?)},''comp'');',ic)
  hca.YLabel.String = {'v_e\times{B}_y','(mV/m)'};
  set(hca,'ColorOrder',mms_colors('yxz'))
  irf_legend(hca,{'v_e\times{B}','E_y','v_e\times{B}+E_y'},leg_loc,'fontsize',fontsize);
end
if 0 % n*Ve
  hca = irf_panel('N*Vex');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('data_pl = mvaVe?*1e3*ne?*1e6;',ic)
  irf_plot(hca,{data_pl.x,data_pl.x.smooth(nsm)},'comp');
  set(hca,'ColorOrder',mms_colors('xyza'))    
  hca.YLabel.String = {'nv_{ex}','(m^{-12}s^{-1})'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'original',sprintf('smoothed, %g',nsm)},leg_loc,'fontsize',fontsize);
end
if 0 % gradx(n*Ve)
  hca = irf_panel('gradx N*Vex');
  set(hca,'ColorOrder',mms_colors('xyza'))
  nsm2 = 10;
  c_eval('data_pl = mvaGradNVe?/vsc*1e-6;',ic)
  irf_plot(hca,{data_pl.smooth(10),data_pl.smooth(nsm)},'comp');
  set(hca,'ColorOrder',mms_colors('xyza'))    
  hca.YLabel.String = {'\partial_x(nv_{ex})','(cm^{-3}s^{-1})'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{sprintf('smoothed, %g',nsm2),sprintf('smoothed, %g',nsm)},leg_loc,'fontsize',fontsize);
end
if 0 % n*Ve, for paper, only one line
  hca = irf_panel('N*Vex');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('data_pl = mvaVe?*1e3*ne?*1e6;',ic)
  scaling = 1e6;
  irf_plot(hca,{data_pl.x*1e-4/scaling},'comp'); % 1e-4 to do: m-2  ->  cm-2
  set(hca,'ColorOrder',mms_colors('xyza'))    
  hca.YLabel.String = {'nv_{ex}',sprintf('(10^{%g} cm^{-2}s^{-1})',log10(scaling))};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'original',sprintf('smoothed, %g',nsm)},leg_loc,'fontsize',fontsize);
end
if 0 % gradx(n*Ve), less data
  hca = irf_panel('gradx N*Vex');
  set(hca,'ColorOrder',mms_colors('12'))
  nsm2 = 10;
  c_eval('data_pl = mvaGradNVe?/vsc*1e-6;',ic)
  irf_plot(hca,{data_pl.smooth(nsm)},'comp');
  set(hca,'ColorOrder',mms_colors('xyza'))    
  hca.YLabel.String = {'\partial_x(nv_{ex})','(cm^{-3}s^{-1})'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{sprintf('smoothed, %g',nsm2),sprintf('smoothed, %g',nsm)},leg_loc,'fontsize',fontsize);
  irf_legend(hca,{sprintf('N_{sm}=%g',nsm),sprintf('v=%g{km/s}',vsc*1e-3)}',[1.01 0.99],'fontsize',14,'color','k');
end
if 1 % n*Ve, grad_x(nVe) with two yaxes
  set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
  
  hca = irf_panel('N*Vex, grad nvex');
  set(hca,'ColorOrder',mms_colors('12'))
  c_eval('data_pl = mvaVe?*1e3*ne?*1e6;',ic)
  scaling = 1e6;
  irf_plot(hca,{data_pl.x*1e-4/scaling},'comp'); % 1e-4 to do: m-2  ->  cm-2
  set(hca,'ColorOrder',mms_colors('xyza'))    
  hca.YLabel.String = {'nu_{ex}',sprintf('(10^{%g} cm^{-2}s^{-1})',log10(scaling))};
  
  hca.YLim = [-29 39];
  
  %set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
  ax1 = hca;
  yyaxis(hca,'right')
  ax2 = gca;
  %ax2.Color = 'none';
  %ax2 = axes('position',hca.Position); 
  %ax2.YAxisLocation = 'right';   
  color_ax = [0    0.7000         0];
  %ax2.YAxis.Color = color_ax;
  c_eval('data_pl2 = mvaGradNVe?/vsc*1e-6;',ic)
  set(ax2,'ColorOrder',color_ax)
  irf_plot(ax2,{data_pl2.smooth(nsm)},'comp');
  %ax2.Visible = 'off';
  ax2.YLabel.String = {'\partial_x(nu_{ex}) (cm^{-3}s^{-1})'};
  
  set(hca,'ColorOrder',[mms_colors('1');color_ax])
  irf_legend(hca,{'nu_{x}','\partial_x(nu_{x})'},[0.98 0.95],'fontsize',fontsize_leg);
  %set(gcf,'defaultAxesColorOrder',[0 0 0; 0 0 0]);
  ax2.YLim = [-29 39]*4/39;
end

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

time_reversal = irf_time('2017-07-11T22:34:02.641773681Z','utc>EpochTT');
c_eval('hm(?) = irf_pl_mark(h(?),time_reversal);',1:numel(h))
c_eval('hm(?).Color = [0.3 0.3 0.3];',1:numel(h))
c_eval('hm(?).LineStyle = ''--'';',1:numel(h))

hh = findobj(gcf,'type','axes'); hh = hh(end:-1:1);

%irf_zoom(h(1:iisub),'x',fastTint)
irf_zoom(hh,'x',tint_edr)
%irf_zoom(h,'x',tint)
irf_zoom(h,'y')
%irf_plot_axis_align
%h(1).Title.String = irf_ssub('MMS ?',ic);
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('hh(?).LineWidth = 1;',1:numel(h))
c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
h(end).XTickLabelRotation = 0;

ylim_gradnvex = [-1.5 3.99];
ax2.YLim = ylim_gradnvex*vsc/130e3;
ax2.YLim = [-29 39]*ylim_gradnvex(2)/39;


hca = irf_panel('ne');
hca.YLim = [0.021 0.069];

hca = irf_panel('stress x comp');
hca.YLim = [-4.49 2.99];


hca = irf_panel('gradients 1, smoothed 2');
hca.YLim = [-5.99 3.99];

hca = irf_panel('gradients 1, smoothed 2, mV/m');
hca.YLim = [-3.99 6.99];

%% Figure 3. Pyz, 2D distribution and pressure and stress contributions, MN, X times, with locations shown
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
%times = times + 0.06;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;
markersize = 5;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(4,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(2);

h1(1).Position(1) = position(1); 
h1(2).Position(1) = 0.53; 
h1(2).Position(3) = 0.32; 
h1(1).Position(2) = position(2);
h1(2).Position(2) = position(2);
h1(1).Position(4) = position(4);
h1(2).Position(4) = position(4);
h1(1).Position(3) = 0.32; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*3)
  h2(ip) = subplot(4,nt,nt+ip);
end
comps = 'yz';
comps_ts = 'yz';
comps_str = ["y","z"];
fig = gcf; fig.Position = [343     1   938   696];

isub = 1;
if 1 % Pe
  hca = h1(isub); isub = isub + 1;
  c_eval('mvaPe = mvaPe?;',ic)
  %c_eval('mvaB = mvaB?;',ic)
  %irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
  hp = irf_patch(hca,{mvaPe.(comps_ts)*1e3,0},'color','k','linewidth',1);
  %hp = irf_patch(hca,{mvaB.y,0},'color','k','linewidth',1);
  %irf_plot(hca,{mvaB.x,mvaB.y},'comp')
  hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  hca.YLabel.Interpreter = 'tex';
  irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
  irf_zoom(hca,'y')
  %hmark = irf_pl_mark(hca,times.epochUnix,'k');
  %c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
  %c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
  xtickangle(h1,0)
  %h1.Position(4) = 0.16;
  ppsum = [];
  hp.EdgeColor = [0 0 0];
  hp.FaceColor = [0.5 0.5 0.5];
  hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  %h1(1).YLabel.String = 'P_{e,LM} (pPa)';
  h1(1).YLabel.Interpreter = 'tex';
end
if 1 % S
  hca = h1(isub); isub = isub + 1;
  c_eval('mvaSe = mvaSe?;',ic)  
  %irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
  hp = irf_patch(hca,{mvaSe.(comps_ts)*1e3,0},'color','k','linewidth',1);
  %irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
  hca.YLabel.String = ['X',sprintf('_{e%s} (pPa)',comps)];
  %hca.YLabel.String = 'S_{eLM} (pPa)';
  hca.YLabel.Interpreter = 'tex';
  irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
  irf_zoom(hca,'y')
  %hmark = irf_pl_mark(hca,times.epochUnix,'k');
  %c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
  %c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
  xtickangle(h1,0)
  %h1.Position(4) = 0.16;
  ppsum = [];
  hp.EdgeColor = [0 0 0];
  hp.FaceColor = [0.5 0.5 0.5];
  %h1(2).YLabel.String = 'S_{e,LM} (pPa)';
  h1(2).YLabel.Interpreter = 'tex';
  h1(2).YLim = h1(1).YLim;
end

ql = 35;
q_lw = 1.5;
dp_clim = 0.01999*[-1 1];
ds_clim = 0.01999*[-1 1];
contour_levels = 10.^([-14 -13 -12 -11 -10 -9]+0.75);

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
  Ldsl = Tdsl(2,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';
  
  vdf = dist.reduce('2D',Mdsl,Ndsl,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'contour',contour_levels);
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
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
  if 1 % plot B direction quiver
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(2)*ql,-bpl(3)*ql,bpl(2)*2*ql,bpl(3)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);    end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 0 % plot B direction
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(3)/b(2);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(2:3);
  end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    hbulk = plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
% S
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)   

  t_dist_center = dist.time.start + (dist.time.stop-dist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data;tsMdsl?.resample(t_dist_center).data;tsNdsl?.resample(t_dist_center).data];',ic)  
  %c_eval('Tdsl = [tsLdsl?.resample(dist).data;tsMdsl?.resample(dist).data;tsNdsl?.resample(dist).data];',ic)
  Ldsl = Tdsl(2,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';

  vdf = dist.reduce('2D',Mdsl,Ndsl,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress','contour',contour_levels);  
  hc_.Colorbar.YLabel.String = sprintf('f_e(v_%s,v_%s)v_%s v_%s (1/m^3)',comps(1),comps(2),comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)v_%s v_%s',comps(1),comps(2),comps(1),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Intgrate data to compare with moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    sstmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    sssum(itime) = sstmp;
    
    if sstmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    %irf_legend(hca,sprintf('p = %.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
    irf_legend(hca,sprintf('%.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
  end
  if 1 % change units
    newdata = ha_.CData;
    data_scale_num = 1;
    data_scale_let = 1e12;
    dv = (vg(2) - vg(1))*1e3; % m
    newdata = newdata*units.me*dv*dv*data_scale_num*data_scale_let;
    ha_.CData = newdata;
    %hc_.Colorbar.YLabel.String = {sprintf('m_ef_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)dv_%sdv_%s',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2),comps(1),comps(2)),'(1/m^3)'};
    %hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s} (10^{%.0f} pPa)',comps(1),comps(2),log10(data_scale_num))};
    hc_.Colorbar.YLabel.String = {['dX',sprintf('_{e%s%s}(v_%s,v_%s)',comps(1),comps(2),comps(1),comps(2))],'(pPa)'};
    hca.CLim = dp_clim;
    1;    
  end
  if 1 % plot B direction quiver
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(2)*ql,-bpl(3)*ql,bpl(2)*2*ql,bpl(3)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);    end
   % irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 0 % plot B direction
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
      B_std_inplane = std(B__.data(:,2:3),1);
      B_inplane = sqrt(sum(B_(2:3).^2));
      if B_inplane > 2*norm(B_std_inplane)
        k = b(3)/b(2);
        if k > 1
          plot(hca,xlim,xlim*k,'k')
        else
          plot(hca,xlim/k,ylim,'k')
        end
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(2:3);
    end
      %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    hbulk = plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
% P
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  t_dist_center = dist.time.start + (dist.time.stop-dist.time.start)/2;
  t_center{itime} = t_dist_center;

  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data;tsMdsl?.resample(t_dist_center).data;tsNdsl?.resample(t_dist_center).data];',ic)
  
  %c_eval('Tdsl = [tsLdsl?.resample(dist).data;tsMdsl?.resample(dist).data;tsNdsl?.resample(dist).data];',ic)
  Ldsl = Tdsl(2,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);

  c_eval('B = dmpaB?.resample(dist);',ic)  
  c_eval('ve = dbcsVe?.resample(dist);',ic)    
  ve = ve*Tdsl';
  B = B*Tdsl';
  

  %tint = irf.tint('2015-05-09T14:00:00Z/2015-05-09T17:59:59Z');
   %defatt = mms.db_get_variable('mms2_ancillary_defatt','zra',tint);
   %defatt.zdec = mms.db_get_variable('mms2_ancillary_defatt','zdec',tint).zdec;
   %gseB = mms_dsl2gse(B_dmpa,defatt);
   
  vdf = dist.reduce('2D',Mdsl,Ndsl,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.y.resample(dist),ve.z.resample(dist),'contour',contour_levels); 
 
  %hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  

  
  if 1 % Integrate data to compare to moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    %irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
    irf_legend(hca,sprintf('%.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end
  if 1 % change units
    newdata = ha_.CData;
    data_scale_num = 1;
    data_scale_let = 1e12;
    dv = (vg(2) - vg(1))*1e3; % m
    newdata = newdata*units.me*dv*dv*data_scale_num*data_scale_let;
    ha_.CData = newdata;
    %hc_.Colorbar.YLabel.String = {sprintf('m_ef_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)dv_%sdv_%s',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2),comps(1),comps(2)),'(1/m^3)'};
    %hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s} (10^{%.0f} pPa)',comps(1),comps(2),log10(data_scale_num))};
    %hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s} (pPa)',comps(1),comps(2))};
    %hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s}(v_%s,v_%s)',comps(1),comps(2),comps(1),comps(2)),'(pPa)'};
    hc_.Colorbar.YLabel.String = {['dP',sprintf('_{e%s%s}(v_%s,v_%s)',comps(1),comps(2),comps(1),comps(2))],'(pPa)'};
    hca.CLim = dp_clim;
    1;    
  end
  if 1 % plot B direction quiver
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(2)*ql,-bpl(3)*ql,bpl(2)*2*ql,bpl(3)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);    end
 %   irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 0 % plot B direction
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(3)/b(2);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(2:3);
    end

    %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    hbulk = plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);

    %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end

% Formatting
% c_eval('tmark = [times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1]; hmark = irf_pl_mark(h1(!),tmark,[0.5 0.5 0.5]); hmark.FaceAlpha = 0.4;',1:times.length,1:2)
for ip = 1:2
  for it = 1:numel(times_exact)
    axes(h1(ip))
    tmark = [times_exact{it}(1).epochUnix times_exact{it}(end).epochUnix] + 0.5*0.03*[-1 1]; 
    hmark = irf_pl_mark(h1(ip),tmark,[0.5 0.5 0.5]); 
    hmark.FaceAlpha = 0.5;
  end
end


hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
c_eval('hb(?).FontSize = 14;',1:numel(hb))
hb = hb(end:-1:1);
%
ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;

hb(5).Position(2) = h2(5).Position(2);
hb(5).Position(4) = h2(5).Position(4);

h1(1).Position(2) = 0.78;
h1(2).Position(2) = 0.78;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
%c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt+nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).XLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).YLim = 0.99*[-60 60];',1:numel(h2))
%c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)
c_eval('h2(?).Layer = ''top'';',1:numel(h2))

if 0 % Change units of data from m^-3 to kg*m^-3
for ip = 6:15 
  hsurf = findobj(h2(ip),'type','surface');
  hsurf.CData = hsurf.CData*units.me;
end
end

if 1 % make TSeries of integrated dP, dS
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  isub = 2;
  tsPP = irf.ts_scalar([t_center{:}],ppsum);
  hold(h1(isub),'on')
  hpp = irf_plot(h1(isub),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(isub),'off')
  
  isub = 1;
  tsSS = irf.ts_scalar([t_center{:}],sssum);
  hold(h1(isub),'on')
  hpp = irf_plot(h1(isub),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(isub),'off')
  
  h1(2).YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  h1(1).YLabel.String = ['X',sprintf('_{e%s} (pPa)',comps)];
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end



hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';


c_eval('h2(?).Color = 1*[1 1 1];',1:15)
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:15)
c_eval('h1(?).XGrid = ''off''; h1(?).YGrid = ''off'';',1:numel(h1))
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:15)

%
legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
nInd = 1;
for ii = 1:2
  irf_legend(h1(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h1(ii).FontSize = 14;
end
for ii = 1:numel(h2)
  irf_legend(h2(ii),legends{nInd},[0.04 0.97],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h2(ii).FontSize = 14;
end

%% Figure 4. Pxy, 2D distribution and pressure and stress contributions, LM, X times, with locations shown
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
%times = times + 0.06;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;
markersize = 5;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(4,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(2);

h1(1).Position(1) = position(1); 
h1(2).Position(1) = 0.53; 
h1(2).Position(3) = 0.32; 
h1(1).Position(2) = position(2);
h1(2).Position(2) = position(2);
h1(1).Position(4) = position(4);
h1(2).Position(4) = position(4);
h1(1).Position(3) = 0.32; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*3)
  h2(ip) = subplot(4,nt,nt+ip);
end
comps = 'xy';
comps_ts = 'xy';
comps_str = ["x","y"];
fig = gcf; fig.Position = [343     1   938   696];

isub = 2;
if 1 % P  
  hca = h1(isub);
  c_eval('mvaPe = mvaPe?;',ic)
  %c_eval('mvaB = mvaB?;',ic)
  %irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
  hp = irf_patch(hca,{mvaPe.(comps_ts)*1e3,0},'color','k','linewidth',1);
  %hp = irf_patch(hca,{mvaB.y,0},'color','k','linewidth',1);
  %irf_plot(hca,{mvaB.x,mvaB.y},'comp')
  hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  hca.YLabel.Interpreter = 'tex';
  irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
  irf_zoom(hca,'y')
  %hmark = irf_pl_mark(hca,times.epochUnix,'k');
  %c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
  %c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
  xtickangle(h1,0)
  %h1.Position(4) = 0.16;
  ppsum = [];
  hp.EdgeColor = [0 0 0];
  hp.FaceColor = [0.5 0.5 0.5];
  hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  %h1(1).YLabel.String = 'P_{e,LM} (pPa)';
  h1(1).YLabel.Interpreter = 'tex';
end
isub = 1;
if 1
  hca = h1(isub);
  c_eval('mvaSe = mvaSe?;',ic)  
  %irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
  hp = irf_patch(hca,{mvaSe.(comps_ts)*1e3,0},'color','k','linewidth',1);
  %irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
  %hca.YLabel.String = sprintf('S_{e%s} (pPa)',comps);
  hca.YLabel.String = ['X',sprintf('_{e%s} (pPa)',comps)];  
  %hca.YLabel.String = 'S_{eLM} (pPa)';
  hca.YLabel.Interpreter = 'tex';
  %hca.YLabel.String = [sprintf('\\mathcal{X}_{e%s} (pPa)',comps)]; 
  %hca.YLabel.String = ['$\mathcal{X}_{e' comps '}$ (pPa)']; 
  %hca.YLabel.Interpreter = 'latex';

  irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
  irf_zoom(hca,'y')
  %hmark = irf_pl_mark(hca,times.epochUnix,'k');
  %c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
  %c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
  xtickangle(h1,0)
  %h1.Position(4) = 0.16;
  ppsum = [];
  hp.EdgeColor = [0 0 0];
  hp.FaceColor = [0.5 0.5 0.5];
  %h1(2).YLabel.String = 'S_{e,LM} (pPa)';
  h1(2).YLabel.Interpreter = 'tex';
  h1(2).YLim = h1(1).YLim;
  dp_clim = 0.025*[-1 1];
  ds_clim = 0.025*[-1 1];
end


ql = 30;
q_lw = 1.5;
dp_clim = 0.01999*[-1 1];
ds_clim = 0.01999*[-1 1];
contour_levels = 10.^([-14 -13 -12 -11 -10 -9]+0.5);
%contour_levels = 10.^([-14 -13 -12 -11 -10 -9]+0.75);


times_exact = {};
isub = 1;
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
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
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
  if 1 % plot B direction quiver
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
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(1)*ql,-bpl(2)*ql,bpl(1)*2*ql,bpl(2)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      %plot(hca,[-bpl(1)*ql bpl(1)*ql],[-bpl(2)*ql,bpl(2)*ql],'color',[1 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);
    end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  
  if 0 % plot B direction
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
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hbulk = plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    hold(hca,'off')    
  end
end
%
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    


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
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress','contour',contour_levels);  
  hc_.Colorbar.YLabel.String = sprintf('f_e(v_%s,v_%s)v_%s v_%s (1/m^3)',comps(1),comps(2),comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)v_%s v_%s',comps(1),comps(2),comps(1),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Intgrate data to compare with moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    sstmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    sssum(itime) = sstmp;
    
    if sstmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    %irf_legend(hca,sprintf('p = %.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
    irf_legend(hca,sprintf('%.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
  end
  if 1 % change units
    newdata = ha_.CData;
    data_scale_num = 1;
    data_scale_let = 1e12;
    dv = (vg(2) - vg(1))*1e3; % m
    newdata = newdata*units.me*dv*dv*data_scale_num*data_scale_let;
    ha_.CData = newdata;
    %hc_.Colorbar.YLabel.String = {sprintf('m_ef_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)dv_%sdv_%s',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2),comps(1),comps(2)),'(1/m^3)'};
    %hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s} (10^{%.0f} pPa)',comps(1),comps(2),log10(data_scale_num))};
    %hc_.Colorbar.YLabel.String = {sprintf('dS_{e%s%s}(v_x,v_y)',comps(1),comps(2)),'(pPa)'};
    hc_.Colorbar.YLabel.String = {['dX',sprintf('_{e%s%s}(v_x,v_y)',comps(1),comps(2))],'(pPa)'};
    hca.CLim = dp_clim;
    1;    
  end
  if 1 % plot B direction quiver
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
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)      
      quiver(hca,-bpl(1)*ql,-bpl(2)*ql,bpl(1)*2*ql,bpl(2)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);
    end
  end
  
  if 0 % plot B direction
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
      %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hbulk = plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    hold(hca,'off')    
  end
end

for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  
  %tint = irf.tint('2015-05-09T14:00:00Z/2015-05-09T17:59:59Z');
   %defatt = mms.db_get_variable('mms2_ancillary_defatt','zra',tint);
   %defatt.zdec = mms.db_get_variable('mms2_ancillary_defatt','zdec',tint).zdec;
   %gseB = mms_dsl2gse(B_dmpa,defatt);



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
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.x.resample(dist),ve.y.resample(dist),'contour',contour_levels); 
 
  %hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  
  if 1 % Integrate data to compare to moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    %irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
    irf_legend(hca,sprintf('%.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end
  if 1 % change units
    newdata = ha_.CData;
    data_scale_num = 1;
    data_scale_let = 1e12;
    dv = (vg(2) - vg(1))*1e3; % m
    newdata = newdata*units.me*dv*dv*data_scale_num*data_scale_let;
    ha_.CData = newdata;
    %hc_.Colorbar.YLabel.String = {sprintf('m_ef_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)dv_%sdv_%s',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2),comps(1),comps(2)),'(1/m^3)'};
    %hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s} (10^{%.0f} pPa)',comps(1),comps(2),log10(data_scale_num))};
    hc_.Colorbar.YLabel.String = {sprintf('dP_{e%s%s} (pPa)',comps(1),comps(2))};
    
    hca.CLim = dp_clim;
    1;    
  end
  if 1 % plot B direction quiver
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
    bpl = B_;

    if B_inplane > 2*norm(B_std_inplane)
      quiver(hca,-bpl(1)*ql,-bpl(2)*ql,bpl(1)*2*ql,bpl(2)*2*ql,0,'color',[0 0 0],'linewidth',q_lw)      
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(1:2);
    end
    %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 0 % plot B direction
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

    %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    %plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hbulk = plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ok','MarkerFaceColor','w','markersize',markersize);
    hold(hca,'off')    
  end
end


% Formatting

% for ip = 1:2
%   for it = 1:numel(times_exact)
%     axes(h1(ip))
%     tmark = [times_exact{it}(1).epochUnix times_exact{it}(end).epochUnix] + 0.5*0.03*[-1 1]; 
%     hmark = irf_pl_mark(h1(ip),tmark,[0.5 0.5 0.5]); 
%     %hmark.FaceAlpha = 0.5;
%   end
% end
c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]); hmark(?).FaceAlpha = 0.5;',1:times.length,1:2)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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

ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;

hb(5).Position(2) = h2(5).Position(2);
hb(5).Position(4) = h2(5).Position(4);

h1(1).Position(2) = 0.78;
h1(2).Position(2) = 0.78;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt+nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).XLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).YLim = 0.99*[-60 60];',1:numel(h2))
%c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)
c_eval('h2(?).Layer = ''top'';',1:numel(h2))

if 0 % Change units of data from m^-3 to kg*m^-3
for ip = 6:15 
  hsurf = findobj(h2(ip),'type','surface');
  hsurf.CData = hsurf.CData*units.me;
end
end

if 1 % make TSeries of integrated dP, dS
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  isub = 2;
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(isub),'on')
  hpp = irf_plot(h1(isub),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(isub),'off')
  
  isub = 1;
  tsSS = irf.ts_scalar(times,sssum);
  hold(h1(isub),'on')
  hpp = irf_plot(h1(isub),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(isub),'off')
  
  h1(2).YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  h1(1).YLabel.String = ['X',sprintf('_{e%s} (pPa)',comps)];
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';

c_eval('h1(?).XGrid = ''off''; h1(?).YGrid = ''off'';',1:2)

c_eval('h2(?).Color = 1*[1 1 1];',1:15)
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:15)
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:15)


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)','n)','o)','p)','q)','r)'};
nInd = 1;
for ii = 1:2
  irf_legend(h1(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h1(ii).FontSize = 14;
end

for ii = 1:numel(h2)
  irf_legend(h2(ii),legends{nInd},[0.04 0.97],'color',[0 0 0],'fontsize',14)
  nInd = nInd + 1;
  h2(ii).FontSize = 14;
end

%% 2D distribution and pressure and stress contributions, MN, adapted from LM, X times, with locations shown
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
%times = times + 0.06;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(4,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(2);

h1(1).Position(1) = position(1); 
h1(2).Position(1) = 0.53; 
h1(2).Position(3) = 0.32; 
h1(1).Position(2) = position(2);
h1(2).Position(2) = position(2);
h1(1).Position(4) = position(4);
h1(2).Position(4) = position(4);
h1(1).Position(3) = 0.32; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*3)
  h2(ip) = subplot(4,nt,nt+ip);
end
comps = 'yz';
comps_ts = 'yz';
comps_str = ["y","z"];
fig = gcf; fig.Position = [343     1   938   696];

hca = h1(1);
c_eval('mvaPe = mvaPe?;',ic)
%c_eval('mvaB = mvaB?;',ic)
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaPe.(comps_ts)*1e3,0},'color','k','linewidth',1);
%hp = irf_patch(hca,{mvaB.y,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaB.x,mvaB.y},'comp')
hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
%hp.EdgeColor = [0 0 0];
%hp.FaceColor = [0.5 0.5 0.5];
hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
%h1(1).YLabel.String = 'P_{e,LM} (pPa)';
h1(1).YLabel.Interpreter = 'tex';

hca = h1(2);
c_eval('mvaSe = mvaSe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaSe.(comps_ts)*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = sprintf('S_{e%s} (pPa)',comps);
%hca.YLabel.String = 'S_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];
%h1(2).YLabel.String = 'S_{e,LM} (pPa)';
h1(2).YLabel.Interpreter = 'tex';
h1(2).YLim = h1(1).YLim;


times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));

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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(3)/b(2);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(2:3);
  end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'+k')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.(comps(1)).resample(dist),ve.(comps(2)).resample(dist));  
  %hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Integrate data to compare to moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(3)/b(2);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(2:3);
    end

    %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
    end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'+k')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress');  
  hc_.Colorbar.YLabel.String = sprintf('f_e(v_%s,v_%s)v_%s v_%s (1/m^3)',comps(1),comps(2),comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)v_%s v_%s',comps(1),comps(2),comps(1),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Intgrate data to compare with moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    sstmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    sssum(itime) = sstmp;
    
    if sstmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
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
      B_std_inplane = std(B__.data(:,2:3),1);
      B_inplane = sqrt(sum(B_(2:3).^2));
      if B_inplane > 2*norm(B_std_inplane)
        k = b(3)/b(2);
        if k > 1
          plot(hca,xlim,xlim*k,'k')
        else
          plot(hca,xlim/k,ylim,'k')
        end
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(2:3);
    end
      %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'+k')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length,1:2)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)

if 0 % Change units of data from m^-3 to kg*m^-3
for ip = 6:15 
  hsurf = findobj(h2(ip),'type','surface');
  hsurf.CData = hsurf.CData*units.me;
end
end

if 1 % make TSeries of integrated dP, dS
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
  
  
  tsSS = irf.ts_scalar(times,sssum);
  hold(h1(2),'on')
  hpp = irf_plot(h1(2),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(2),'off')
  
  h1(1).YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  h1(2).YLabel.String = sprintf('S_{e%s} (pPa)',comps);
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';

c_eval('h2(?).Color = 1*[1 1 1];',1:15)
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:15)
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:15)


%% 1D distributions for the 5 times
times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
%times = times + 0.25;
times = times + 0.30;
%times = times + 0.06;
dt_dist = 4*0.062; % for X distributions

%elim = [100 Inf];

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -60000:2000:60000;

nrows = 3;
ncols = 2;
ip = 0;
for irow = 1:nrows
  for icol = 1:ncols
    ip = ip + 1;
    h(irow,icol) = subplot(nrows,ncols,ip);
  end
end

times_exact = {};

elim = [100 Inf];
for itime = 1:times.length  
  isub = 1;
  icol = 1;
  if itime == 2
    c_eval('hold(h(?,icol),''on'')',1:3)
  end
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  fL = dist.reduce('1D',L,'vint',vint,'scpot',scpot,'vg',vg);
  fM = dist.reduce('1D',M,'vint',vint,'scpot',scpot,'vg',vg);
  fN = dist.reduce('1D',N,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = fL.time;
    
  % Plot
  hca = h(isub,icol); isub = isub + 1;
  plot(hca,fL.depend{1}(1,:)*1e-3,mean(fL.data,1))
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)','x');
  hca.YLabel.String = sprintf('f_e (s/m^4)');
  
  hca = h(isub,icol); isub = isub + 1;
  plot(hca,fM.depend{1}(1,:)*1e-3,mean(fN.data,1))
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)','y');
  hca.YLabel.String = sprintf('f_e (s/m^4)');
  
  hca = h(isub,icol); isub = isub + 1;
  plot(hca,fM.depend{1}(1,:)*1e-3,mean(fN.data,1))
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)','z');
  hca.YLabel.String = sprintf('f_e (s/m^4)');

  if 0 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end

elim = [0 Inf];
for itime = 1:times.length
  icol = 2;
  isub = 1;
  if itime == 2
    c_eval('hold(h(?,icol),''on'')',1:3)
  end
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  fL = dist.reduce('1D',L,'vint',vint,'scpot',scpot,'vg',vg);
  fM = dist.reduce('1D',M,'vint',vint,'scpot',scpot,'vg',vg);
  fN = dist.reduce('1D',N,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = fL.time;
    
  % Plot
  hca = h(isub,icol); isub = isub + 1;
  plot(hca,fL.depend{1}(1,:)*1e-3,mean(fL.data,1))
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)','x');
  hca.YLabel.String = sprintf('f_e (s/m^4)');
  
  hca = h(isub,icol); isub = isub + 1;
  plot(hca,fM.depend{1}(1,:)*1e-3,mean(fN.data,1))
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)','y');
  hca.YLabel.String = sprintf('f_e (s/m^4)');
  
  hca = h(isub,icol); isub = isub + 1;
  plot(hca,fM.depend{1}(1,:)*1e-3,mean(fN.data,1))
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)','z');
  hca.YLabel.String = sprintf('f_e (s/m^4)');

  if 0 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end


c_eval('hold(h(?,1),''off'')',1:3)
c_eval('hold(h(?,2),''off'')',1:3)

%% 1D distributions for the 5 times, comparing elim with no elim
times = EpochTT(['2017-07-11T22:34:01.300000000Z';
  '2017-07-11T22:34:01.800000000Z';...
  '2017-07-11T22:34:02.300000000Z';...
  '2017-07-11T22:34:02.800000000Z';...
  '2017-07-11T22:34:03.300000000Z']);
%times = times([2 3 4])+0.25;
%times = times + 0.25;
times = times + 0.30;
%times = times + 0.06;
dt_dist = 4*0.062; % for X distributions
nMC = 500;
%elim = [100 Inf];

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -60000:1000:60000;

nrows = 3;
ncols = times.length;
ip = 0;
h = gobjects(0);
for irow = 1:nrows
  for icol = 1:ncols
    ip = ip + 1;
    h(irow,icol) = subplot(nrows,ncols,ip);
  end
end

times_exact = {};

elows = [0 50 100 150 200];
ehigh = Inf;

for itime = 1:times.length  
  icol = itime;
  time = times(itime);
  emin = {};
  % Reduce distributions
  for elow = elows
    isub = 1;
    elim = [elow ehigh];
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)  
    emin{end+1} = dist.depend{1}(1,1) - dist.ancillary.delta_energy_minus(1,1);
    c_eval('scpot = scPot?.resample(dist);',ic)  
    c_eval('B = mvaB?.resample(dist);',ic)  
    c_eval('ve = mvaVe?.resample(dist);',ic)  
    fL = dist.reduce('1D',L,'vint',vint,'scpot',scpot,'vg',vg,'nMC',nMC);
    fM = dist.reduce('1D',M,'vint',vint,'scpot',scpot,'vg',vg,'nMC',nMC);
    fN = dist.reduce('1D',N,'vint',vint,'scpot',scpot,'vg',vg,'nMC',nMC);
    times_exact{itime} = fL.time;
    
    % Plot
    hca = h(isub,icol); isub = isub + 1;
    plot(hca,fL.depend{1}(1,:)*1e-3,mean(fL.data,1))
    hca.XLabel.String = sprintf('v_%s (10^3 km/s)','x');
    hca.YLabel.String = sprintf('f_e(v_%s) (s/m^4)','x');

    hca = h(isub,icol); isub = isub + 1;
    plot(hca,fM.depend{1}(1,:)*1e-3,mean(fM.data,1))
    hca.XLabel.String = sprintf('v_%s (10^3 km/s)','y');
    hca.YLabel.String = sprintf('f_e(v_%s) (s/m^4)','y');

    hca = h(isub,icol); isub = isub + 1;
    plot(hca,fM.depend{1}(1,:)*1e-3,mean(fN.data,1))
    hca.XLabel.String = sprintf('v_%s (10^3 km/s)','z');
    hca.YLabel.String = sprintf('f_e(v_%s) (s/m^4)','z');
    drawnow
    c_eval('hold(h(?,icol),''on'')',1:3)     
  end 
  legends = cellfun(@(x)sprintf('E>%.0f eV',x),emin,'UniformOutput',false);
  irf_legend(h(1,1),legends',[0.98 0.98])
  
    if 0 % plot bulk speed
      %%
      hold(hca,'on')
      plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
      plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
      hold(hca,'off')    
    end
end
c_eval('h(?).XGrid = ''on'';',1:numel(h)) 
c_eval('h(?).YGrid = ''on'';',1:numel(h)) 
c_eval('h(?).XLim = [-30 30];',1:numel(h)) 
c_eval('hold(h(?),''off'')',1:numel(h))  
%compact_panels(0.01,0.01)

%% 2D distribution and pressure and stress contributions, LM, X times, with locations shown
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
%times = times + 0.06;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(4,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(2);

h1(1).Position(1) = position(1); 
h1(2).Position(1) = 0.53; 
h1(2).Position(3) = 0.32; 
h1(1).Position(2) = position(2);
h1(2).Position(2) = position(2);
h1(1).Position(4) = position(4);
h1(2).Position(4) = position(4);
h1(1).Position(3) = 0.32; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*3)
  h2(ip) = subplot(4,nt,nt+ip);
end
comps = 'xy';
comps_ts = 'xy';
comps_str = ["x","y"];
fig = gcf; fig.Position = [343     1   938   696];

hca = h1(1);
c_eval('mvaPe = mvaPe?;',ic)
%c_eval('mvaB = mvaB?;',ic)
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaPe.(comps_ts)*1e3,0},'color','k','linewidth',1);
%hp = irf_patch(hca,{mvaB.y,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaB.x,mvaB.y},'comp')
hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
%hp.EdgeColor = [0 0 0];
%hp.FaceColor = [0.5 0.5 0.5];
hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
%h1(1).YLabel.String = 'P_{e,LM} (pPa)';
h1(1).YLabel.Interpreter = 'tex';

hca = h1(2);
c_eval('mvaSe = mvaSe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaSe.(comps_ts)*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = sprintf('S_{e%s} (pPa)',comps);
%hca.YLabel.String = 'S_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];
%h1(2).YLabel.String = 'S_{e,LM} (pPa)';
h1(2).YLabel.Interpreter = 'tex';
h1(2).YLim = h1(1).YLim;


times_exact = {};
isub = 1;
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
  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
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
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
%
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress');  
  hc_.Colorbar.YLabel.String = sprintf('f_e(v_%s,v_%s)v_%s v_%s (1/m^3)',comps(1),comps(2),comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)v_%s v_%s',comps(1),comps(2),comps(1),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Intgrate data to compare with moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    sstmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    sssum(itime) = sstmp;
    
    if sstmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
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
      %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end

for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.x.resample(dist),ve.y.resample(dist));  
  %hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Integrate data to compare to moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
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

    %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
    end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'+k')
    plot(hca,mean(ve.x.data,1)*1e-3,mean(ve.y.data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length,1:2)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).XLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).YLim = 0.99*[-60 60];',1:numel(h2))
c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)
c_eval('h2(?).Layer = ''top'';',1:numel(h2))

if 0 % Change units of data from m^-3 to kg*m^-3
for ip = 6:15 
  hsurf = findobj(h2(ip),'type','surface');
  hsurf.CData = hsurf.CData*units.me;
end
end

if 1 % make TSeries of integrated dP, dS
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
  
  
  tsSS = irf.ts_scalar(times,sssum);
  hold(h1(2),'on')
  hpp = irf_plot(h1(2),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(2),'off')
  
  h1(1).YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  h1(2).YLabel.String = sprintf('S_{e%s} (pPa)',comps);
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';


c_eval('h2(?).Color = 1*[1 1 1];',1:15)
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:15)
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:15)

%% 2D distribution and pressure and stress contributions, MN, adapted from LM, X times, with locations shown
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
%times = times + 0.06;
dt_dist = 4*0.061; % for two distributions

fontsize_B_amp = 13;

vint = [-Inf Inf];
%vint = [-20000 20000];
vg = -100000:2000:100000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(4,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(2);

h1(1).Position(1) = position(1); 
h1(2).Position(1) = 0.53; 
h1(2).Position(3) = 0.32; 
h1(1).Position(2) = position(2);
h1(2).Position(2) = position(2);
h1(1).Position(4) = position(4);
h1(2).Position(4) = position(4);
h1(1).Position(3) = 0.32; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*3)
  h2(ip) = subplot(4,nt,nt+ip);
end
comps = 'yz';
comps_ts = 'yz';
comps_str = ["y","z"];
fig = gcf; fig.Position = [343     1   938   696];

hca = h1(1);
c_eval('mvaPe = mvaPe?;',ic)
%c_eval('mvaB = mvaB?;',ic)
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaPe.(comps_ts)*1e3,0},'color','k','linewidth',1);
%hp = irf_patch(hca,{mvaB.y,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaB.x,mvaB.y},'comp')
hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
%hp.EdgeColor = [0 0 0];
%hp.FaceColor = [0.5 0.5 0.5];
hca.YLabel.String = sprintf('P_{e%s} (pPa)',comps);
%h1(1).YLabel.String = 'P_{e,LM} (pPa)';
h1(1).YLabel.Interpreter = 'tex';

hca = h1(2);
c_eval('mvaSe = mvaSe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaSe.(comps_ts)*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = sprintf('S_{e%s} (pPa)',comps);
%hca.YLabel.String = 'S_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];
%h1(2).YLabel.String = 'S_{e,LM} (pPa)';
h1(2).YLabel.Interpreter = 'tex';
h1(2).YLim = h1(1).YLim;


times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = sprintf('log_{10} f_e(v_%s,v_%s) (s^2/m^5)',comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('log_{10} f_e(v_%s,v_%s)',comps(1),comps(2)),'(s^2/m^5)'};
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));

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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(3)/b(2);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(2:3);
  end
    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
    
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'+k')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.(comps(1)).resample(dist),ve.(comps(2)).resample(dist));  
  %hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)(v_%s-v_%s)(v_%s-v_%s)',comps(1),comps(2),comps(1),comps(1),comps(2),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Integrate data to compare to moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
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
    B_std_inplane = std(B__.data(:,2:3),1);
    B_inplane = sqrt(sum(B_(2:3).^2));
    if B_inplane > 2*norm(B_std_inplane)
      k = b(3)/b(2);
      if k > 1
        plot(hca,xlim,xlim*k,'k')
      else
        plot(hca,xlim/k,ylim,'k')
      end
      hold(hca,'off')
      hca.XLim = xlim;
      hca.YLim = ylim;
      B_ = B_(2:3);
    end

    %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
    end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'+k')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress');  
  hc_.Colorbar.YLabel.String = sprintf('f_e(v_%s,v_%s)v_%s v_%s (1/m^3)',comps(1),comps(2),comps(1),comps(2));
  hc_.Colorbar.YLabel.String = {sprintf('f_e(v_%s,v_%s)v_%s v_%s',comps(1),comps(2),comps(1),comps(2)),'(1/m^3)'};
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = sprintf('v_%s (10^3 km/s)',comps(1));
  hca.YLabel.String = sprintf('v_%s (10^3 km/s)',comps(2));
  
  if 1 % Intgrate data to compare with moments
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    sstmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    sssum(itime) = sstmp;
    
    if sstmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
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
      B_std_inplane = std(B__.data(:,2:3),1);
      B_inplane = sqrt(sum(B_(2:3).^2));
      if B_inplane > 2*norm(B_std_inplane)
        k = b(3)/b(2);
        if k > 1
          plot(hca,xlim,xlim*k,'k')
        else
          plot(hca,xlim/k,ylim,'k')
        end
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(2:3);
    end
      %irf_legend(hca,sprintf('B = %.1f nT', B_inplane),[0.02 0.98],'color','k','fontsize',fontsize_B_amp)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'+k')
    plot(hca,mean(ve.(comps(1)).data,1)*1e-3,mean(ve.(comps(2)).data,1)*1e-3,'ow')
    hold(hca,'off')    
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length,1:2)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;

c_eval('h2(?).CLim = 0.99*[-5000 5000];',6:15)

if 0 % Change units of data from m^-3 to kg*m^-3
for ip = 6:15 
  hsurf = findobj(h2(ip),'type','surface');
  hsurf.CData = hsurf.CData*units.me;
end
end

if 1 % make TSeries of integrated dP, dS
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
  
  
  tsSS = irf.ts_scalar(times,sssum);
  hold(h1(2),'on')
  hpp = irf_plot(h1(2),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(2),'off')
  
  h1(1).YLabel.String = sprintf('P_{e%s} (pPa)',comps);
  h1(2).YLabel.String = sprintf('S_{e%s} (pPa)',comps);
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

%hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
%hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';

c_eval('h2(?).Color = 1*[1 1 1];',1:15)
c_eval('h2(?).XGrid = ''off''; h2(?).YGrid = ''off'';',1:15)
color_grid = 0.5*[1 1 1];
c_eval('hold(h2(?),"on"); plot(h2(?),vg*0,vg,''color'',color_grid); plot(h2(?),vg,vg*0,''color'',color_grid);',1:15)

%% 2D distribution and pressure and stress contributions, MN, X times, with locations shown
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
times = times + 0.30;
%times = times + 0.06;
dt_dist = 4*0.062; % for two distributions

vint = [-Inf Inf];
%vint = [-20000 20000];
vint = [-10000 10000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(4,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(2);

h1(1).Position(1) = position(1); 
h1(2).Position(1) = 0.53; 
h1(2).Position(3) = 0.32; 
h1(1).Position(2) = position(2);
h1(2).Position(2) = position(2);
h1(1).Position(4) = position(4);
h1(2).Position(4) = position(4);
h1(1).Position(3) = 0.32; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*3)
  h2(ip) = subplot(4,nt,nt+ip);
end

hca = h1(1);
c_eval('mvaPe = mvaPe?;',ic)
c_eval('mvaB = mvaB?;',ic)
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaPe.yz*1e3,0},'color','k','linewidth',1);
%hp = irf_patch(hca,{mvaB.y,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaB.x,mvaB.y},'comp')
hca.YLabel.String = 'P_{eMN} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
%hp.EdgeColor = [0 0 0];
%hp.FaceColor = [0.5 0.5 0.5];
h1(1).YLabel.String = 'P_{e,LM} (pPa)';
h1(1).YLabel.Interpreter = 'tex';

hca = h1(2);
c_eval('mvaSe = mvaSe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaSe.yz*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = 'S_{eMN} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
%h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];
h1(2).YLabel.String = 'S_{e,LM} (pPa)';
h1(2).YLabel.Interpreter = 'tex';
h1(2).YLim = h1(1).YLim;


times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  %vint = [-1 1]*20e3;
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = 'v_M (10^3 km/s)';
  hca.YLabel.String = 'v_N (10^3 km/s)';

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(3)/b(2);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    B_ = B_(2:3);
    irf_legend(hca,sprintf('B = %.1f nT', sqrt(sum(B_.^2))),[0.02 0.98],'color','k','fontsize',14)
  end
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
    hold(hca,'off')    
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.y.resample(dist),ve.z.resample(dist));  
  hc_.Colorbar.YLabel.String = 'f_e(v_M,v_N)(v_M-v_M^{bulk})(v_N-v_N^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_M (10^3 km/s)';
  hca.YLabel.String = 'v_N (10^3 km/s)';
  
  if 1
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_dist = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(3)/b(2);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
  
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
    hold(hca,'off')    
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress');  
  hc_.Colorbar.YLabel.String = 'f_e(v_M,v_N)(v_M-v_M^{bulk})(v_N-v_N^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_M (10^3 km/s)';
  hca.YLabel.String = 'v_N (10^3 km/s)';
  
  if 1
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    sstmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    sssum(itime) = sstmp;
    
    if sstmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',sstmp),[0.98 0.98],'color',color,'fontsize',13)
  end

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_dist = 0.030;
    %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = B.tlim(dist.time([1 end]) + 0.5*0.03*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    disp(b)
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
  
  if 1 % plot bulk speed
    %%
    hold(hca,'on')
    plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
    hold(hca,'off')    
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length,1:2)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1) 2*nt+1:(3*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

%h1(1).Title.String = sprintf('MMS %g',ic);
%h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
%h1.Position(3) = h1.Position(3)*0.8;




if 1
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
  
  
  tsSS = irf.ts_scalar(times,sssum);
  hold(h1(2),'on')
  hpp = irf_plot(h1(2),tsSS,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(2),'off')
  
  h1(1).YLabel.String = 'P_{e,LM} (pPa)';
  h1(2).YLabel.String = 'S_{e,LM} (pPa)';
  h1(1).YLabel.Interpreter = 'tex';
  h1(2).YLabel.Interpreter = 'tex';
end

hbb = findobj(gcf,'type','colorbar'); hbb = hbb(end:-1:1);

hbb(2).YLabel.String = 'dP^e_{LM}/dv_Ldv_M (kg/m^3)';
hbb(3).YLabel.String = 'dS^e_{LM}/dv_Ldv_M (kg/m^3)';


%% 2D distribution and pressure contributions, MN, X times, with locations shown
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
times = times + 0.30;
%times = times + 0.06;
dt_dist = 2*0.062; % for two distributions

vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(3,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*2)
  h2(ip) = subplot(3,nt,nt+ip);
end

hca = h1(1);
c_eval('mvaPe = mvaPe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaPe.xy*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = 'P_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];

times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';

  if 0 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    %B_
    irf_legend(hca,sprintf('B = %.1f nT', sqrt(sum(B_.^2))),[0.02 0.98],'color','k','fontsize',14)
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.x.resample(dist),ve.y.resample(dist));  
  hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';
  
  if 1
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_dist = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(1),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

h1(1).Title.String = sprintf('MMS %g',ic);
h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
h1.Position(3) = h1.Position(3)*0.8;


if 1
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
end

%% 2D distribution and pressure contributions, LM, X times, comparing different satellites
ics = [3 4];
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
  '2017-07-11T22:34:03.300000000Z';...
  '2017-07-11T22:34:03.670000000Z';...
  '2017-07-11T22:34:04.300000000Z']);
%times = times([2 3 4])+0.25;
times = times + 0.05 + 0*0.25;
times = times + 0.05 + 2*0.03;


vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(5,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 

%delete(h1_pos);
for ip = 1:(nt*2*2)
  h2(ip) = subplot(5,nt,nt+ip);
end

hca = h1(1);
set(hca,'colororder',mms_colors('12'))
c_eval('mvape_1 = mvaPe?;',ics(1))  
c_eval('mvape_2 = mvaPe?;',ics(2))  
irf_plot(hca,{mvape_1.yz,mvape_2.yz},'comp')
hold(hca,'on')
fhigh = 2;
irf_plot(hca,{mvape_1.yz.filt(0,fhigh,[],5),mvape_2.yz.filt(0,fhigh,[],5)},'comp')
hold(hca,'off')

hca.YLabel.String = 'P_{eMN} (nPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
hmark = irf_pl_mark(hca,times.epochUnix,'k');
c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
%
isub = 1;
for ic = ics
  for itime = 1:times.length
    %hca = h2(isub); isub = isub + 1;
    time = times(itime);
    % Reduce distributions
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
    c_eval('scpot = scPot?.resample(dist);',ic)  
    vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
      
    % Plot
    hca = h2(isub); isub = isub + 1;
    [ha_,hb_,hc_] = vdf.plot_plane(hca);
    hc_.Colorbar.YLabel.String = 'log_{10} f_e(v_M,v_N) (s^2/m^5)';
    colormap(hca,pic_colors('candy4'))   
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
  end
  for itime = 1:times.length
    %hca = h2(isub); isub = isub + 1;
    time = times(itime);
    % Reduce distributions
    c_eval('dist = ePDist?.elim(elim).tlim(time+0.015*[-1 1]*1); dist = dist(:);',ic)
    c_eval('scpot = scPot?.resample(dist);',ic)  
    c_eval('ve = mvaVe?;',ic)  
    vdf = dist.reduce('2D',M,N,'vint',vint,'scpot',scpot,'vg',vg);
      
    % Plot
    hca = h2(isub); isub = isub + 1;
    [ha_,hb_,hc_] = vdf.plot_plane(hca,'off-diag-pres-cont',ve.y.resample(dist),ve.z.resample(dist));  
    hc_.Colorbar.YLabel.String = 'f_e(v_M,v_N)(v_M-v_M^{bulk})(v_N-v_N^{bulk}) (1/m^3)';
    colormap(hca,pic_colors('blue_red'))
    hca.CLim = max(abs(hca.CLim))*[-1 1];
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
  end
end

hlinks_all = linkprop(h2,{'XLim','YLim'});
%hlinks_f = linkprop(h2([1:nt]*[0 2]'),{'CLim'});
%hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});
hlinks_f = linkprop(h2([1:nt 2*nt+(1:nt)]),{'CLim'});
hlinks_p = linkprop(h2([(nt+1):2*nt 2*nt+((nt+1):2*nt)]),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1) (2*nt-1)+1:(3*nt-1) (3*nt-1)+1:(4*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
%ihsub = [2:nt nt+2:(2*nt)];
ihsub = [2:nt nt+2:(2*nt) 2*nt+2:(3*nt) 3*nt+2:(4*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)


%% 2D distribution and stress contributions, MN, X times, with locations shown
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
times = times + 0.30;
%times = times + 0.06;
dt_dist = 2*0.062; % for two distributions

vint = [-Inf Inf];
vint = [-20000 20000];
vg = -60000:2000:60000;
elim = [100 Inf];
%h = setup_subplots(2,times.length,'vertical'); 
nt = times.length;
h1_pos = subplot(3,nt,1:nt); position = h1_pos.Position;
h1 = irf_plot(1); 
h1.Position = position; 
%delete(h1_pos);
clear h2
for ip = 1:(nt*2)
  h2(ip) = subplot(3,nt,nt+ip);
end

hca = h1(1);
c_eval('mvaSe = mvaSe?;',ic)  
%irf_plot(hca,mvaPe.xy*1e3,'color','k','linewidth',1)
hp = irf_patch(hca,{mvaSe.xy*1e3,0},'color','k','linewidth',1);
%irf_plot(hca,{mvaPe1.xy*1e3,mvaB1.abs},'comp')
hca.YLabel.String = 'S_{eLM} (pPa)';
hca.YLabel.Interpreter = 'tex';
irf_zoom(hca,'x',irf.tint('2017-07-11T22:33:58.30Z/2017-07-11T22:34:05.99Z'))
irf_zoom(hca,'y')
%hmark = irf_pl_mark(hca,times.epochUnix,'k');
%c_eval('hmark(?) = irf_pl_mark(hca,times(?).epochUnix + 0.5*0.03*[-1 1],''k'');',1:times.length)
%c_eval('hmark(?).LineWidth = 1;',1:numel(hmark))
xtickangle(h1,0)
h1.Position(4) = 0.16;
ppsum = [];
hp.EdgeColor = [0 0 0];
hp.FaceColor = [0.5 0.5 0.5];

times_exact = {};
isub = 1;
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)  
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
  times_exact{itime} = vdf.time;
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca);
  hc_.Colorbar.YLabel.String = 'log_{10} f_e (s^2/m^5)';
  colormap(hca,pic_colors('candy4'))   
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_distx = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    %B_
    irf_legend(hca,sprintf('B = %.1f nT', sqrt(sum(B_.^2))),[0.02 0.98],'color','k','fontsize',14)
  end
end
for itime = 1:times.length
  %hca = h2(isub); isub = isub + 1;
  time = times(itime);
  % Reduce distributions
  c_eval('dist = ePDist?.elim(elim).tlim(time+0.5*dt_dist*[-1 1]*1); dist = dist(:);',ic)
  c_eval('scpot = scPot?.resample(dist);',ic)  
  c_eval('B = mvaB?.resample(dist);',ic)  
  c_eval('ve = mvaVe?.resample(dist);',ic)    
  vdf = dist.reduce('2D',L,M,'vint',vint,'scpot',scpot,'vg',vg);
    
  % Plot
  hca = h2(isub); isub = isub + 1;
  [ha_,hb_,hc_] = vdf.plot_plane(hca,'stress');  
  hc_.Colorbar.YLabel.String = 'f_e(v_L,v_M)(v_L-v_L^{bulk})(v_M-v_M^{bulk}) (1/m^3)';
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'v_L (10^3 km/s)';
  hca.YLabel.String = 'v_M (10^3 km/s)';
  
  if 1
    %%
    cdata = ha_.CData;
    xdata = ha_.XData; 
    dv = xdata(2)-xdata(1);
    dv = dv*1e3;
    pptmp = nansum(cdata(:))*dv*dv*units.me*1e9*1e9;
    ppsum(itime) = pptmp;
    
    if pptmp < 0
      color = [0 0 1];
    else
      color = [1 0 0];
    end
    irf_legend(hca,sprintf('p = %.3f pPa',pptmp),[0.98 0.98],'color',color,'fontsize',13)
  end

  if 1 % plot B direction
    %%
    xlim = hca.XLim;
    ylim = hca.YLim;
    hold(hca,'on')
    %dt_dist = 0.030;
    B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
    B_ = mean(B_.data,1);
    b = B_/norm(B_);
    k = b(2)/b(1);
    if k > 1
      plot(hca,xlim,xlim*k,'k')
    else
      plot(hca,xlim/k,ylim,'k')
    end
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
  end
end

c_eval('hmark(?) = irf_pl_mark(h1(1),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]);',1:times.length)

hlinks_all = linkprop(h2,{'XLim','YLim'});
hlinks_f = linkprop(h2(1:nt),{'CLim'});
hlinks_p = linkprop(h2((nt+1):2*nt),{'CLim'});

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
%
ihsub = [1:nt-1 nt+1:(2*nt-1)];
delete(hb(ihsub))
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
ih = nt;
hb(ih).Position(2) = hb(ih).Position(2)+hb(ih).Position(4)*0.05; 
hb(ih).Position(4) = hb(ih).Position(4)*0.9; 
%hb(3).Position(4) = 0.22; hb(6).Position(4) = 0.22;
%ihsub = [2 3 5 6];
ihsub = [2:nt nt+2:(2*nt)];
c_eval('h2(?).YLabel.String = [];',ihsub)
c_eval('h2(?).YTickLabel = [];',ihsub)
c_eval('h2(?).XTickLabelRotation = 0;',(nt+1):2*nt)

h1(1).Title.String = sprintf('MMS %g',ic);
h1.Position(1) = h1.Position(1) + h1.Position(3)*0.2;
h1.Position(3) = h1.Position(3)*0.8;




if 1
  %%
  %for ip = 1:times.length
  %  hca = h2(ip);
  %  %cdata = hca.
  %end
  %diff(hc_.Surface.XData(1:2))
  %pp = [-1.483,-0.732,-0.042,0.128,0.513];
  tsPP = irf.ts_scalar(times,ppsum);
  hold(h1(1),'on')
  hpp = irf_plot(h1(1),tsPP,'or');
  hpp.MarkerFaceColor = 'k';
  hpp.MarkerEdgeColor = 'k';
  hpp.MarkerSize = 7;  
  hold(h1(1),'off')
end


%% PIC: example of plots, for explaining
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(21);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];

%xpick = 102.4;
%xpick = 106.0;

zpick = 0;
for xpick = 102.4:0.2:106.0

  %ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds);
  ds = ds100.twpelim(twpe).xfind(xpick).zfind(zpick);
  clear h
  h(1) = subplot(2,1,1);
  h(2) = subplot(2,1,2);
  isub = 1;
  
  
  %ispecies = [2 4 6];
  ispecies = [4];
  sumdim = 3;
  fontsize = 14;
  
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'bline',no02m);
  hca.CLim = [0 0.002];
  
  
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'off-diag','bline',no02m);
  colormap(hca,pic_colors('blue_red'));
  hca.CLim = 0.04*[-1 1];
  
  ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
  delete(ht)
  
  if 1
    hi = findobj(hca.Children,'type','image');
    sumf = sum(hi.CData(:));
    if sumf > 0
      color = [1 0 0];
    else
      color = [0 0 1];
    end
    irf_legend(hca,sprintf('%.3f',sumf),[0.02 0.98],'color',color,'fontsize',fontsize)
    hca.Color = color;
  end
  
  %c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))
  c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
  
  
  for ip = 1:numel(h)
    axis(h(ip),'square');
    h(ip).XTick = -10:5:10;
    h(ip).YTick = -10:5:10;
    h(ip).XTickLabelRotation = 0;
  end
  c_eval('h(?).FontSize = fontsize;',1:numel(h))
  
  compact_panels(h,0.02)
  
  drawnow

  if 1 % print
    h(1).Title.String = sprintf('x = %.1f, z = %.1f',xpick,zpick);
    cn.print(sprintf('ex_fxy_x=%.1f_e4',xpick))
  
    
    h(1).XLabel.String = 'v_L';
    h(2).XLabel.String = 'v_L';
    h(1).YLabel.String = 'v_M';
    h(2).YLabel.String = 'v_M';
    cn.print(sprintf('ex_fLM_x=%.1f_e4',xpick))
  end
end

%% PIC: example of plots, for explaining, plot location of picked boxes
figure(19);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];

xpick = 102.4:0.2:106.0;
zpick = 0;

pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);

h = pic.plot_map({'pexy'},'A',0.1,'sep');
h.CLim = 0.012*[-1 1];
colormap(pic_colors('blue_red'));
hold(h(1),'on')
[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(h(1));
hold(h(1),'off')
h(1).Position(2) = 0.18;
h(1).Position(4) = 0.7;
h(1).XLabel.String = 'L (d_i)';
h(1).YLabel.String = 'N (d_i)';
h(1).FontSize = 14;

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
hb.YLabel.String = 'P^e_{LM}';

for ii = 1:numel(hdl)
  c_eval('hdl(?).LineWidth = 1;',1:numel(hdl))
  c_eval('hdl(?).LineWidth = 2;',ii)
  cn.print(sprintf('ex_locationboxes_x=%.1f',xpick(ii)))
end

%% PIC: example of plots, for explaining, combined with location boxes
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(18);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];

%xpick = 102.4;
%xpick = 106.0;

doVertical = 1;
ispecies = [4 6];
zpick = 0;
sumdim = 3;
fontsize = 16;
pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-3 3]).zlim([-0.99 0.99]);

xpick_all = 102.4:0.2:106.0;
for xpick = xpick_all

  %ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds);
  ds = ds100.twpelim(twpe).xfind(xpick).zfind(zpick);
  clear h
  if not(doVertical)
    h(1) = subplot(1,4,[1 2]);
    h(2) = subplot(1,4,3);
    h(3) = subplot(1,4,4);
  
    h(1).Position(1) = 0.17;
    h(1).Position(3) = 0.3;
    
    %h(1).Position(1) = 0.08;
  else
    h(1) = subplot(3,4,[1 2 3 4]);
    h(2) = subplot(3,4,[5 6 9 10]);
    h(3) = subplot(3,4,[7 8 11 12]);
    h(3).Position(1) = h(3).Position(1) - 0.057;
    
    %h(1) = axes('Position',[0.1300    0.6593    0.7179    0.2157]);
    %h(2) = axes('Position',[0.1300    0.1100    0.3025    0.5154]);
    %h(3) = axes('Position',[0.5225    0.1100    0.3025    0.5154]);
    %h(1).Position(1) = 0.17;
    %h(1).Position(3) = 0.3;

  end
  isub = 1;

  % Location of boxes
  hca = h(isub); isub = isub + 1;
  hh = pic.plot_map(hca,{'pexy'},'A',0.1,'sep','smooth',3);
  hca.CLim = 0.012*[-1 1];
  colormap(hca,pic_colors('blue_red'));
  hold(hca,'on')
  %[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
  [hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
  hold(hca,'off')
  %hca.Position(2) = 0.18;
  %hca.Position(4) = 0.7;
  hca.XLabel.String = 'L (d_i)';
  hca.YLabel.String = 'N (d_i)';
  hca.FontSize = 14;
  
  hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
  hb.YLabel.String = 'P^e_{LM}';
  
  
  
  hdl.LineWidth = 2;
  
  
  % f
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'bline',no02m);
  hca.CLim = [0 0.002];
  
  % Pressure integrand
  hca = h(isub); isub = isub + 1;
  htmp = ds.plot_map(hca,ispecies,sumdim,'off-diag',1/100,'bline',no02m);
  colormap(hca,pic_colors('blue_red'));
  hca.CLim = 0.04*[-1 1];
  
  ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
  delete(ht)
  
  if 1
    hi = findobj(hca.Children,'type','image');
    sumf = sum(hi.CData(:));
    if sumf > 0
      color = [1 0 0];
    else
      color = [0 0 1];
    end
    irf_legend(hca,sprintf('%.3f',sumf*1/100),[0.98 0.98],'color',color,'fontsize',fontsize)
    hca.Color = color;
  end
  
  %c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))
  c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))
  
  
  for ip = 2:numel(h)
    axis(h(ip),'square');
    h(ip).XTick = -10:5:10;
    h(ip).YTick = -10:5:10;
    h(ip).XTickLabelRotation = 0;
  end
  c_eval('h(?).FontSize = fontsize;',1:numel(h))
  
  %compact_panels(h,0.02)
  
  drawnow


    h(2).XLabel.String = 'v_L';
    h(3).XLabel.String = 'v_L';
    h(2).YLabel.String = 'v_M';
    h(3).YLabel.String = 'v_M';

    if doVertical
      h(1).Position(2) = h(1).Position(2)-0.05;
      compact_panels(h(2:3),0,0.01)
      h(3).YLabel.String = [];
      h(3).YTickLabels = [];
    end
  if 1 % print
    %h(1).Title.String = sprintf('x = %.1f, z = %.1f',xpick,zpick);
  %  cn.print(sprintf('ex_fxy_x=%.1f_comb',xpick))
  
    
    cn.print(sprintf('ex_fLM_x=%.1f_comb_vertical',xpick))
  end
end

%% Map of off-diag PeLM
% for LM
xpicks = 102.4:0.4:105.0;
zpicks = -0.4:0.2:0;

% for MN
xpicks = 102.4:0.4:105.0;
zpicks = -0.2:0.1:0;

ds = ds100.twpelim(16000).xfind(xpicks).zfind(zpicks);


figure(103)
clear h
h(1) = subplot(1,1,1);
isub = 1;
hca = h(isub); isub = isub + 1;
hh = pic.plot_map(hca,{'peyz'},'A',0.1,'sep','smooth',3);
hca.CLim = 0.012*[-1 1];
colormap(hca,pic_colors('blue_red'));
hold(hca,'on')
%[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
[hd,hdl] = ds.plot_boxes(hca);
hold(hca,'off')
%hca.Position(2) = 0.18;
%hca.Position(4) = 0.7;
hca.XLabel.String = 'L (d_i)';
hca.YLabel.String = 'N (d_i)';
hca.FontSize = 14;

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
hb.YLabel.String = 'P^e_{MN}';


%%
figure(104)
sumdim = 1;

h = findobj(gcf,'type','axes'); h = h(end:-1:1); delete(h);
h = ds.plot_map([2 4 6],sumdim,'bline',pic);
compact_panels(h.ax,0.01,0.01)
%c_eval('axis(h.ax(?),''equal'');',1:numel(h.ax))
c_eval('axis(h.ax(?),''square'');',1:numel(h.ax))
c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))
hlinks = linkprop(h.ax,{'CLim'});

%hlinks.Targets(1).CLim = 0.04*[-1 1];
hlinks.Targets(1).CLim = 10e-3*[0 1];
c_eval('h.ax(?).XTick = -10:5:10;',1:numel(h.ax))
c_eval('h.ax(?).YTick = -10:5:10;',1:numel(h.ax))
c_eval('h.ax(?).XTickLabelRotation = 90;',1:numel(h.ax))

if sumdim == 3
  c_eval('h.ax(?).YLabel.String = ''v_M'';',1:numel(zpicks))
  c_eval('h.ax(?).XLabel.String = ''v_L'';',1:numel(zpicks):numel(h.ax))
  hb = findobj(gcf,'type','colorbar');
  hb.YLabel.String = 'f(v_M,v_L)';
elseif sumdim == 1
  c_eval('h.ax(?).YLabel.String = ''v_N'';',1:numel(zpicks))
  c_eval('h.ax(?).XLabel.String = ''v_M'';',1:numel(zpicks):numel(h.ax))
  hb = findobj(gcf,'type','colorbar');
  hb.YLabel.String = 'f(v_M,v_N)';


end
ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
delete(ht)

compact_panels(h.ax,0.003,0.003)

c_eval('h.ax(?).FontSize = 13;',1:numel(h.ax))


for ip = 1:numel(h.ax)
  hca = h.ax(ip);
  hi = findobj(hca.Children,'type','image');
  sumf = sum(hi.CData(:));
  if sumf > 0
    color = [1 0 0];
  else
    color = [0 0 1];
  end
  irf_legend(hca,sprintf('%.3f',sumf*1/100),[0.98 0.98],'color',color,'fontsize',fontsize)
  hca.Color = color;
  %colormap(hca,pic_colors('blue_red'))
  colormap(hca,pic_colors('candy4'))
  
end

%% Coordinate system
% illustrate_magnetic_reconnection

doVideo = 1;
doGif = 1;
fileName = 'illustration_magnetic_reconnection';

a = 5;
b = 1;
x = a*linspace(-10,10,200);
y = b*linspace(-10,10,100);
z = linspace(-10,10,5);
[X,Y] = meshgrid(x,y);
dx = x(2) - x(1);
dy = y(2) - y(1);
%dz = z(2) - z(1);
x_xline = x;
y_xline = x*b/a;

Ay = @(x,y) (x/a).^2 - (y/b).^2;
AY0 = Ay(X,Y);

%[FX,FY] = gradient(AY,dx,dy);
%Bx = -FX;
%By = FY;

colors = pic_colors('matlab');
colors = [colors; colors(end:-1:1,:)];

hca = subplot(1,1,1);
t = 0:30;
Astep = 20;
dA = Astep/numel(t);
AYlev0 = -100:Astep:(100 + Astep);

% Initiate

it = 1;
% Draw separatrix
if 0 % 2D
  plot(hca,x_xline,y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  hold(hca,'on')
  plot(hca,x_xline,-y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.5:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.5:numel(S(is).Y));%S(is).Y;
    plot(hca,sx,sy,'color',[0 0 0],'linewidth',2)
  end
elseif 1 % 3D
  plot3(hca,x_xline*0,x_xline,y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  hold(hca,'on')
  plot3(hca,x_xline*0,x_xline,-y_xline,'linewidth',1,'linestyle','--','color',[0,0,0])
  % Draw field lines
  AY = AY0 - dA*t(it);
  S = contourcs(x,y,AY,AYlev0);
  for is = 1:numel(S)
    sx = interp1(1:numel(S(is).X),S(is).X,1:0.5:numel(S(is).X));
    sy = interp1(1:numel(S(is).Y),S(is).Y,1:0.5:numel(S(is).Y));%S(is).Y;
    plot3(hca,sx*0,sx,sy,'color',[0 0 0],'linewidth',2)
  end
end

pause(0.1)
drawnow
hold(hca,'off')
hca.Visible = 'off';
hca.Position = [0 0 1 1];









%% PIC: map of pxy
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(21);
twpe = 16000;
%zlim_ds = [-0.25 0.25];
%xlim_ds = [100.7 104.3];
xpick = 100.8:0.4:104.2;
xpick = 100.0:0.4:105.0;
zpick = -0.3:0.1:0.3;

%ds = ds100.twpelim(twpe).xlim(xlim_ds).zlim(zlim_ds);
ds = ds100.twpelim(twpe).xfind(xpick).zfind(zpick);

%ispecies = [2 4 6];
ispecies = [6];
%h = ds.plot_map(ispecies,1,'off-diag','bline',no02m);
%h = ds.plot_map(ispecies,3,'bline',no02m);
h = ds.plot_map(ispecies,1,'ratio',[4 6],'bline',no02m);

compact_panels(h.ax,0.,0.)
%colormap(pic_colors('blue_red'));
c_eval('h.ax(?).XGrid = ''off''; h.ax(?).YGrid = ''off'';',1:numel(h.ax))




%
ht = findobj(gcf,'type','text'); ht = ht(end:-1:1);
delete(ht)
for ip = 1:numel(h.ax)
  axis(h.ax(ip),'square');
  h.ax(ip).XTick = -10:5:10;
  h.ax(ip).YTick = -10:5:10;
  h.ax(ip).XTickLabelRotation = 0;
  
  % sum of integrand
  if 0
    hi = findobj(h.ax(ip).Children,'type','image');
    sumf = sum(hi.CData(:));
    if sumf > 0
      color = [1 0 0];
    else
      color = [0 0 1];
    end
    irf_legend(h.ax(ip),sprintf('%.3f',sumf),[0.02 0.98],'color',color)
    h.ax(ip).Color = color;
  end
end

%% Plot location of boxes
figure(22);
pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);

h = pic.plot_map({'pexy'});
h.CLim = 0.012*[-1 1];
colormap(pic_colors('blue_red'));
hold(h(1),'on')
ds.plot_boxes(h(1));
hold(h(1),'off')

%% PIC: peij
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
zlim = 1.5*[-1 1];

varstrs = {'pexy','peyz'}';
cbarlabels = {'P^e_{LM}','P^e_{MN}'};
clims = {9.9e-3*[-1 1],9.9e-3*[-1 1]};
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',0,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

c_eval('h(?).FontSize = 14;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

%% PIC: peij, + gradients
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
zlim = 1.5*[-1 1];

varstrs = {'pexy','peyz';'dxPexy','dzPezy'}';
cbarlabels = {'P^e_{LM}','P^e_{MN}';'\partial_LP^e_{LM}','\partial_NP^e_{MN}'}';
%cbarlabels = {'P^e_{LM}','P^e_{MN}';'\partial P^e_{LM}/\partial L','\partial P^e_{MN}/\partial N'}';

clims = {9.9e-3*[-1 1],9.9e-3*[-1 1];0.049*[-1 1],0.049*[-1 1]}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',0,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

c_eval('axis(h(?),''equal'');',1:numel(h))
compact_panels(h,0.01,0.11)

hb = findobj(gcf,'type','colorbar'); ht = ht(end:-1:1);
c_eval('hb(?).FontSize = 16;',1:numel(hb))
c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

% hl = findobj(gcf,'type','contour');
% c_eval('hl(?).LineWidth = 1;',1:numel(hl))

%% PIC: peij, + gradients divided by density
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
%zlim = 1.5*[-1 1];
zlim = 0.99*[-1 1];

varstrs = {'pexy','peyz';'-dxPexy./ne','-dzPezy./ne'}';
cbarlabels = {'P^e_{LM}','P^e_{MN}';'-\partial_LP^e_{LM}/n','-\partial_NP^e_{MN}/n'}';
%cbarlabels = {'P^e_{LM}','P^e_{MN}';'\partial P^e_{LM}/\partial L','\partial P^e_{MN}/\partial N'}';

clims = {9.9e-3*[-1 1],9.9e-3*[-1 1];0.199*[-1 1],0.199*[-1 1]}';
cmaps = {pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red'),pic_colors('blue_red')};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',3,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

c_eval('axis(h(?),''equal'');',1:numel(h))
compact_panels(h,0.01,0.11)

hb = findobj(gcf,'type','colorbar'); ht = ht(end:-1:1);
c_eval('hb(?).FontSize = 16;',1:numel(hb))
c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

%hl = findobj(gcf,'type','contour');
%c_eval('hl(?).LineWidth = 1;',1:numel(hl))
%c_eval('h(?).LineWidth = 1;',1:numel(h))
%c_eval('h(?).Position(2) = h(?).Position(2) + 0.02;',1:numel(h))

%% PIC: electron ohms law
pic = no02m;

twpe = 16000;
xlim = mean(pic.xi) + 3*[-1 1];
zlim = 0.99*[-1 1];

varstrs = {'Ey','vexBy','-(1/100)*vdvey','-dxPexy./ne','-dzPezy./ne'}';
cbarlabels = {'E','vxB','-m v\cdot\nabla v','-\partial_LP^e_{LM}/ne','-\partial_NP^e_{MN}/ne','sum'}';

varstrs = {'Ey+vexBy','-(1/100)*vdvey','-dxPexy./ne','-dzPezy./ne'}';
cbarlabels = {'E+vxB','-m v\cdot\nabla v','-\partial_LP^e_{LM}/ne','-\partial_NP^e_{MN}/ne','sum'}';

%varstrs = {'Ey+vexBy','-(1/100)*vdvey-dxPexy./ne-dzPezy./ne'}';
%cbarlabels = {'E+vxB','-v\cdot\nabla v - \partial_LP^e_{LM}/n - \partial_NP^e_{MN}/n','sum'}';


clims = {0.199*[-1 1],0.199*[-1 1],0.199*[-1 1],0.199*[-1 1],0.199*[-1 1],0.199*[-1 1]}';
cmapbr = pic_colors('blue_red');
cmaps = {cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr,cmapbr};

h = pic.twpelim(twpe).xlim(xlim).zlim(zlim).plot_map(varstrs,'A',0.2,'sep','smooth',3,'clim',clims,'cbarlabels',cbarlabels,'cmap',cmaps);

%c_eval('axis(h(?),''equal'');',1:numel(h))
%compact_panels(h,0.01,0.11)

hb = findobj(gcf,'type','colorbar'); ht = ht(end:-1:1);
c_eval('hb(?).FontSize = 16;',1:numel(hb))
c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

c_eval('h(?).Position(3) = 0.7;',1:numel(h))

if 0

  %%
  hl = findobj(gcf,'type','contour');
  c_eval('hl(?).LineWidth = 0.5;',1:numel(hl))
  c_eval('h(?).LineWidth = 0.5;',1:numel(h))
end

%% PIC: electron anisotropy
%no02m = PIC('/Volumes/Fountain/cno062/data/PIC/no_hot_bg_n02_m100/data_h5/fields.h5');
%ds100 = PICDist('/Volumes/DataRaid/cno062/no_hot_bg_n02_m100/data_h5/dists_new.h5');
figure(22);
twpe = 16000;

pic = no02m.twpelim(twpe).xlim(mean(no02m.xi)+[-4 4]).zlim([-2 2]);

varstrs = {'log10(tepar./teperp)','tepar./teperp-1'}';
varstrs = {'log10(tepar./teperp)'}';
cbarlabels = {'log_{10}T^e_{||}/T^e_{\perp}'}';
h = pic.plot_map(varstrs,'A',0.2,'sep','cbarlabels',cbarlabels);
colormap(pic_colors('blue_red'));

c_eval('h(?).FontSize = 16;',1:numel(h))

c_eval('h(?).YLabel.String = ''N (d_i)'';',1:numel(h))
c_eval('h(?).XLabel.String = ''L (d_i)'';',1:numel(h))

h(1).CLim = 1.5*[-1 1];

if 1
  %%
  
  hca = h(1);
  xpicks = 102.4;
  zpicks = -0.4;
  ds = ds100.twpelim(16000).xfind(xpicks).zfind(zpicks);
  hold(hca,'on')
  %[hd,hdl] = ds100.twpelim(twpe).xfind(xpick).zfind(zpick).plot_boxes(hca);
  [hd,hdl] = ds.plot_boxes(hca);
  hold(hca,'off')
end

%% Compare pressure to stress tensor

h = irf_plot(8);

hca = irf_panel('ve');
c_eval('irf_plot(hca,{mvaVe?.x,mvaVe?.y,mvaVe?.z},''comp'');',ic)
hca.YLabel.String = 'v_e (km/s)';

hca = irf_panel('Te anis');
c_eval('irf_plot(hca,{Temprat?},''comp'');',ic)
hca.YLabel.String = 'T^e_{||}/T^e_{\perp}';

for comp = ["xx","yy","zz","xy","xz","yz"]
  hca = irf_panel(char(comp));
  c_eval('irf_plot(hca,{mvaPe?.(comp),mvaSe?.(comp)},''comp'')',ic)
  hca.YLabel.String = comp;
  irf_legend(hca,{'P','S'}',[1.02,0.95])
end
