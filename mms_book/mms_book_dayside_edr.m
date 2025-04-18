%% Load data
ic = [3];

tint = irf.tint('2015-10-16T13:06:30.00Z/2015-10-16T13:07:30.00Z');

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
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('gsmB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint);',ic);

%c_eval('gseB?scm = mms.get_data(''B_gse_scm_brst_l2'',tint,?);',ic)

% Electric field
disp('Loading electric field...') 
c_eval('gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint);',ic);
%c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
%c_eval('tic; gseE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_gse_brst_l2'',tint); toc',ic);
%c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);
%
% Load spacecraft position
disp('Loading spacecraft position...')
c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4);

% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('dcv? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint);',ic);

% Particle moments
% Skymap distributions
if 1
  %%
  ic_ = ic;
  %ic = 3;
disp('Loading skymaps...')
c_eval('tic; ePDist? = mms.get_data(''PDe_fpi_brst_l2'',tint,?); toc;',ic)
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
c_eval('agyro? = 2*(facPeqq?.yy-facPeqq?.zz)/(facPeqq?.yy+facPeqq?.zz); agyro? = agyro?.abs;',ic);

% Compute temperature ratio An
c_eval('Temprat? = facPepp?.xx/(facPepp?.yy);',ic);

c_eval('Te?par = facTe?.xx;',ic);
c_eval('Te?perp = (facTe?.yy+facTe?.yy)/2;',ic);

c_eval('vte? = sqrt(2*units.eV*(facTe?.trace)/3/units.me)*1e-3;',ic) % km/s
c_eval('wce? = units.me*gseB?.abs*1e-9/units.me;',ic) % rad/s
c_eval('rce? = vte?/wce;',ic) % km

c_eval('wpi? = (ne?*1e6*units.e^2/units.eps0/units.mp).^0.5;',ic); % rad/s
c_eval('wpe? = (ne?*1e6*units.e^2/units.eps0/units.me).^0.5;',ic); % rad/s
c_eval('di? = units.c/wpi?*1e-3;',ic); % km
c_eval('de? = units.c/wpe?*1e-3;',ic); % km

c_eval('UB? = gseB?.abs2*1e-18/2/units.mu0;',ic); % SI
c_eval('betae? = gsePe?.trace*1e-9/3/UB?.resample(gsePe?);',ic);
c_eval('betai? = gsePi?.trace*1e-9/3/UB?.resample(gsePi?);',ic);
c_eval('beta? = betae? + betai?.resample(betae?);',ic);



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

%% Rotate things into new coordinate system 
L = -[0.3665,-0.1201,-0.9226];
M = -[0.5694,-0.7553,-0.3245];
N = -[0.7358,0.6443,-0.2084];
lmn = [L;M;N];


%LMN = [0.3665,-0.1201,-0.9226; 0.5694,-0.7553,-0.3245; 0.7358,0.6443,-0.2084]';
%L = LMN(1,:);
%M = LMN(2,:);
%N = LMN(3,:);

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

%% Reduce distributions
ic_dist = 3;
ic = ic_dist;
disp('Preparing reduced distributions.')
vint = [-Inf Inf];
elim = [30 40000];
vg = -10000:200:10000;


c_eval('ef1Dx? = ePDist?.elim(elim).reduce(''1D'',L,''vint'',vint,''vg'',vg);',ic)
c_eval('ef1Dy? = ePDist?.elim(elim).reduce(''1D'',M,''vint'',vint,''vg'',vg);',ic)
c_eval('ef1Dz? = ePDist?.elim(elim).reduce(''1D'',N,''vint'',vint,''vg'',vg);',ic)

%% Figure: Overview 1
ic = 3;

tint_edr = irf.tint('2017-07-11T22:33:58.00Z/2017-07-11T22:34:09.00Z'); %20151112071854

npanels = 8;
h = irf_plot(npanels);
iisub = 0;
cmap = colormap(pic_colors('candy4'));
fontsize = 18;
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
  irf_legend(hca,{'L','M','N'}',[1.02,0.98],'fontsize',fontsize);
end 
if 1 % E LMN
  isub = isub + 1;
  zoomy = [zoomy isub];
  hca = irf_panel('E LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  fhigh = 10;
  %c_eval('irf_plot(hca,{mvaE?.x.filt(0,fhigh,[],3),mvaE?.y.filt(0,fhigh,[],3),mvaE?.z.filt(0,fhigh,[],3)},''comp'');',ic)
  c_eval('irf_plot(hca,{mvaE?.x,mvaE?.y,mvaE?.z},''comp'');',ic)
  hca.YLabel.String = {'E (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{comps(1),comps(2),comps(3)}',[1.02,0.98],'fontsize',fontsize);
  %irf_legend(hca,{sprintf('f<%g Hz',fhigh)},[0.98 0.98],'fontsize',fontsize,'color',[0 0 0]);
end 

if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  %c_eval('irf_plot(hca,{-1*vte?.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  
  hca.YLabel.String = {'u_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{['u_' comps(1)],['u_' comps(2)],['u_' comps(3)]}',[1.02,0.98],'fontsize',fontsize);
end
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  set(hca,'ColorOrder',mms_colors('1'))
  %c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  %c_eval('irf_plot(hca,{-1*vte?.tlim(tint).smooth(30)*1e-3},''comp'',''--'');',ic)  
  
  hca.YLabel.String = {'u_e (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{['u_' comps(1)],['u_' comps(2)],['u_' comps(3)]}',[1.02,0.98],'fontsize',fontsize);
end
if 1 % e DEF omni
  isub = isub + 1;
  hca = irf_panel('i DEF omni');  
  c_eval('[hout,hcb] = irf_spectrogram(hca,ePDist?.omni.deflux.specrec,''log'');',ic)  
  set(hca,'yscale','log');
  %set(hca,'ytick',[1e1 1e2 1e3 1e4]);  
  colormap(hca,cmap) 
  if 0 % E_ExB
    hold(hca,'on')
    c_eval('vexb = gseVExB?_srvy.abs.resample(iPDist?_fast);',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    hexb_p = irf_plot(hca,Eexb,'color',0.2*[1 1 1],'linestyle','-');
    hexb_o = irf_plot(hca,Eexb*16,'color',0.5*[1 1 1],'linestyle','-');
    hold(hca,'off')
    %irf_legend([hexb_p,hexb_o],{'m_pv_{ExB}^2/2','m_Ov_{ExB}^2/2'}',[0.98 0.98])
  end
  if 0 % vi_par
    hold(hca,'on')
    c_eval('vexb = gseVi?_fast_par.abs;',ic)
    Eexb = vexb.^2*1e6*units.mp/2/units.eV;
    irf_plot(hca,Eexb,'b')
    hold(hca,'off')
  end
  hca.YLabel.String = {'E_e','(eV)'};   
  hca.YLabel.Interpreter = 'tex';
end
if 1 % e psd x,y,z, 3 panels
  for comp = ['x','y','z']
    isub = isub + 1;
    hca = irf_panel(['f1D ' comp]);
    c_eval(sprintf('f1D = ef1D%s?;',comp),ic)
    irf_spectrogram(hca,f1D.specrec('velocity_1D'));  
    hca.YLim = f1D.depend{1}(1,[1 end]);  
    if 0 % % Vi, VExB
      hold(hca,'on')    
      c_eval('hl = irf_plot(hca,{gseVi?.(comp),gseVExB?.(comp).resample(gseVi?)},''comp'');',ic)      
      hold(hca,'off')
      irf_legend(hl,{'v_i','v_{ExB}'},[0.02 0.12])
      hl = findobj(hca,'type','line');
    end
    hca.YLabel.String = {sprintf('v_{e%s}',comp),'(km/s)'}; 
    hca.YLabel.Interpreter = 'tex';
  end
end

if 0 % beta
  hca = irf_panel('beta');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{betae?, betai?, beta?},''comp'');',ic)  
  hca.YScale = 'log';
  hca.YLabel.String = {'\beta'};
  hca.YTick = 10.^[-3:5];
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'\beta_e','\beta_e','\beta'}',[1.02,0.98],'fontsize',fontsize);
end
if 0 % lengths
  hca = irf_panel('length scales');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{di?, de?, di?*sqrt(betai?.resample(di?)), rce?},''comp'');',ic)  
  hca.YScale = 'log';
  hca.YLabel.String = {'L (km)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'d_i','d_e','d_i\beta^{1/2}','\rho_{e}'}',[1.02,0.98],'fontsize',fontsize);
  irf_legend(hca,{'Hall','e- inert.','e- pres.'},[0.05,0.98],'fontsize',fontsize);
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:npanels
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  h(ii).FontSize = fontsize;
end

irf_zoom(h,'x',gseB3.time)
irf_zoom(h,'y')
irf_plot_axis_align