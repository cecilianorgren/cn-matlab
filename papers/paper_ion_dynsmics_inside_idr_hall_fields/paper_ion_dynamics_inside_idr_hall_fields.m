% Load data
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
mms.db_init('local_file_db','/Volumes/mms');
%db_info = datastore('mms_db');

units = irf_units;
% Torbert event 
ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z');
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('iPDistErr? = mms.get_data(''PDERRi_fpi_brst_l2'',tint,?);',ic) % missing some ancillary data
c_eval('iPDist?_nobg = iPDist?; iPDist?_nobg.data(iPDist?_nobg.data < iPDistErr?.data*1.01) = 0;',ic)
c_eval('iPDist?_counts = iPDist?; iPDist?_counts.data = (iPDist?.data./iPDistErr?.data).^2;',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',1:4);
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',1:4);
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('gseVExB? = cross(gseE?.resample(gseB?.time),gseB?)/gseB?.abs/gseB?.abs*1e3; gseVExB?.units = '''';',ic) % km/s
c_eval('gseVixB? = cross(gseVi?*1e3,gseB?.resample(gseVi?.time)*1e-9)*1e3; gseVixB?.units = '''';',ic) % mV/m
c_eval('gseVexB? = cross(gseVe?*1e3,gseB?.resample(gseVe?.time)*1e-9)*1e3; gseVexB?.units = '''';',ic) % mV/m

c_eval('gseVexB? = cross(gseVe?*1e3,gseB?.resample(gseVe?.time)*1e-9)*1e3; gseVexB?.units = '''';',ic) % mV/m
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic);

c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)

c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',1:4);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',1:4);

c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);


%c_eval('[enflux_new, enflux_BG, idist_new, idist_BG, Ni_new, gseVi_new, gsePi_new,Ni_bg, EnergySpectr_bg, Pres_bg, EnergySpectr_bg_self]= mms.remove_ion_penetrating_radiation_bg(iPDist?);',ic)

c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic)

c_eval('B?inf = irf.ts_scalar(gseB?.time,sqrt(gseB?.abs2.data*1e-18 + 1e-9*gsePi?.resample(gseB?).trace.data/3*(2*units.mu0)))*1e9;',ic)


c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)

c_eval('gseJe? = -units.e*ne?*gseVe?*1e3*1e6*1e9; gseJe?.units = ''nA/m^2''; gseJe?.coordinateSystem = ''GSE'';',1:4);
c_eval('gseJi? = units.e*ne?*gseVi?.resample(ne?.time)*1e3*1e6*1e9; gseJi?.units = ''nA/m^2''; gseJi?.coordinateSystem = ''GSE'';',1:4);
c_eval('gseJ? = (gseJe?+gseJi?);',1:4);

c_eval('gseR?brsttime = gseR?.resample(gseB?);',1:4)
[Jcurl,divB,gseB,JxB,gseCurvB,gseDivPb] = c_4_j('gseR?brsttime','gseB?');

gseJcurl = irf.ts_vec_xyz(Jcurl.time,Jcurl.data); gseJcurl.coordinateSystem = 'GSE';
gseJcurl.data = gseJcurl.data*1e9; Jcurl.units = 'nAm^{-2}';
gseJcurl.time = EpochTT(gseJcurl.time); gseJcurl.name = '4sc current density';

gseBav = (gseB1.resample(gseB2) + gseB2.resample(gseB2) + gseB3.resample(gseB2) + gseB4.resample(gseB2))/4; gseBav.name = 'B1234'; 
neav = (ne1.resample(ne1) + ne2.resample(ne1) + ne3.resample(ne1) + ne4.resample(ne1))/4; % gseE1 not there?
gseJxB = gseJcurl.cross(gseBav.resample(gseJcurl)); gseJxB.name = 'JxB'; gseJxB.units = 'nAm^-2 nT';
gseJxBne_mVm = (gseJxB*1e-18)/(neav.resample(gseJxB)*1e6)/units.e*1e3; gseJxBne_mVm.name = 'JxB/ne';
%gseJxBne_mVm.data(abs(gseJxBne_mVm.data)>100) = NaN;


c_eval('gseJxB? = gseJ?.cross(gseB?.resample(gseJ?)); gseJxB?.name = ''JxB''; gseJxB?.units = ''nA/m^2 nT'';',1:4)
c_eval('gseJxBy? =    gseJ?.x*gseB?.y.resample(gseJ?); gseJxBy?.name = ''Jx*By''; gseJxBy?.units = ''nA/m^2 nT'';',1:4)
c_eval('gseJyBx? = -1*gseJ?.y*gseB?.x.resample(gseJ?); gseJyBx?.name = ''-Jy*Bx''; gseJyBx?.units = ''nA/m^2 nT'';',1:4)
c_eval('gseJyBz? =    gseJ?.y*gseB?.z.resample(gseJ?); gseJyBz?.name = ''Jy*Bz''; gseJyBz?.units = ''nA/m^2 nT'';',1:4)
c_eval('gseJzBy? = -1*gseJ?.z*gseB?.y.resample(gseJ?); gseJzBy?.name = ''-Jz*By''; gseJzBy?.units = ''nA/m^2 nT'';',1:4)
c_eval('gseJzBx? =    gseJ?.z*gseB?.x.resample(gseJ?); gseJzBx?.name = ''Jz*Bx''; gseJxBy?.units = ''nA/m^2 nT'';',1:4)
c_eval('gseJxBz? = -1*gseJ?.x*gseB?.z.resample(gseJ?); gseJxBz?.name = ''-Jx*Bz''; gseJxBz?.units = ''nA/m^2 nT'';',1:4)

c_eval('gseJxBne?_mVm = (gseJxB?*1e-18)/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJxBne?_mVm.name = ''JxB/ne'';',1:4)

% FPI
c_eval('gseJxBy?ne_mVm = gseJxBy?*1e-18/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJxBy?ne_mVm.units = ''mV/m'';',1:4)
c_eval('gseJyBx?ne_mVm = gseJyBx?*1e-18/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJyBx?ne_mVm.units = ''mV/m'';',1:4)
c_eval('gseJyBz?ne_mVm = gseJyBz?*1e-18/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJyBz?ne_mVm.units = ''mV/m'';',1:4)
c_eval('gseJzBy?ne_mVm = gseJzBy?*1e-18/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJzBy?ne_mVm.units = ''mV/m'';',1:4)
c_eval('gseJzBx?ne_mVm = gseJzBx?*1e-18/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJzBx?ne_mVm.units = ''mV/m'';',1:4)
c_eval('gseJxBz?ne_mVm = gseJxBz?*1e-18/(ne?.resample(gseJxB?)*1e6)/units.e*1e3; gseJxBz?ne_mVm.units = ''mV/m'';',1:4)

% Curl
gseJxBy_ne_mVm_curl =    gseJcurl.x*gseB.y*1e-18/(neav.resample(gseJcurl)*1e6)/units.e*1e3; gseJxBy_ne_mVm_curl.units = 'mV/m';
gseJyBx_ne_mVm_curl = -1*gseJcurl.y*gseB.x*1e-18/(neav.resample(gseJcurl)*1e6)/units.e*1e3; gseJxBy_ne_mVm_curl.units = 'mV/m';
gseJzBx_ne_mVm_curl =    gseJcurl.z*gseB.x*1e-18/(neav.resample(gseJcurl)*1e6)/units.e*1e3; gseJxBy_ne_mVm_curl.units = 'mV/m';
gseJxBz_ne_mVm_curl = -1*gseJcurl.x*gseB.z*1e-18/(neav.resample(gseJcurl)*1e6)/units.e*1e3; gseJxBy_ne_mVm_curl.units = 'mV/m';
gseJyBz_ne_mVm_curl =    gseJcurl.y*gseB.z*1e-18/(neav.resample(gseJcurl)*1e6)/units.e*1e3; gseJxBy_ne_mVm_curl.units = 'mV/m';
gseJzBy_ne_mVm_curl = -1*gseJcurl.z*gseB.y*1e-18/(neav.resample(gseJcurl)*1e6)/units.e*1e3; gseJxBy_ne_mVm_curl.units = 'mV/m';


% counts = (P/Perr)^2
% sum counts for all time steps,
% remove every bin which has counts<1.5

%% Local coordinate system
L = [0.9482,-0.255,-0.1893];
M = [0.1818,0.9245,-0.3350];
N = [0.2604,0.2832,0.9230];
lmn_edr = [L;M;N];

L_vi = -[-0.8906    0.4548    0.0045];
M_vi = [ 0.4539    0.8893   -0.0559];
N_vi = -[-0.0294   -0.0477   -0.9984];
lmn_vi = [L_vi; M_vi; N_vi];


L_gse = [1 0 0];
M_gse = [0 1 0];
N_gse = [0 0 1];
lmn_gse = [L_gse; M_gse; N_gse];

L_gse = [1 0 -.2]; L_gse = L_gse/norm(L_gse);
M_gse = [0 1 -0.2]; M_gse = cross(L_gse,cross(M_gse,L));
%N_gse = [0 0 1];
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

lmn = lmn_gse;
lmn = lmn_edr;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);


c_eval('mvaVixB? = gseVixB?*lmn''; mvaVixB?.name = ''Vi x B LMN'';',ic)
c_eval('mvaVexB? = gseVexB?*lmn''; mvaVexB?.name = ''Ve x B LMN'';',ic)
c_eval('mvaVExB? = gseVExB?*lmn''; mvaVExB?.name = ''E LMN'';',ic)
c_eval('mvaE? = gseE?*lmn''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaB? = gseB?*lmn''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaVi? = gseVi?*lmn''; mvaVi?.name = ''Vi LMN'';',ic)
c_eval('mvaVe? = gseVe?*lmn''; mvaVe?.name = ''Ve LMN'';',ic)
c_eval('mvaPi? = lmn*gsePi?*lmn''; mvaPi?.units = gsePi?.units;',ic)
c_eval('mvaVe?perp = gseVe?perp*lmn''; mvaVe?perp.name = ''Ve perp LMN'';',ic)


c_eval('mvaJ? = gseJ?*lmn'';',1:4)
c_eval('mvaJxBne?_mVm = gseJxBne?_mVm*lmn'';',1:4)
mvaJxBne_mVm = gseJxBne_mVm*lmn';

c_eval('mvaJxB? = mvaJ?.cross(mvaB?.resample(mvaJ?)); mvaJxB?.name = ''JxB''; mvaJxB?.units = ''nA/m^2 nT'';',ic)
c_eval('mvaJLBM? =    mvaJ?.x*mvaB?.y.resample(mvaJ?); mvaJLBM?.name = ''JL*BM''; mvaJLBM?.units = ''nA/m^2 nT'';',ic)
c_eval('mvaJMBL? = -1*mvaJ?.y*mvaB?.x.resample(mvaJ?); mvaJMBL?.name = ''-JM*BL''; mvaJMBL?.units = ''nA/m^2 nT'';',ic)
c_eval('mvaJMBN? =    mvaJ?.y*mvaB?.z.resample(mvaJ?); mvaJMBN?.name = ''JM*BN''; mvaJMBN?.units = ''nA/m^2 nT'';',ic)
c_eval('mvaJNBM? = -1*mvaJ?.z*mvaB?.y.resample(mvaJ?); mvaJNBM?.name = ''-JN*BM''; mvaJNBM?.units = ''nA/m^2 nT'';',ic)
c_eval('mvaJNBL? =    mvaJ?.z*mvaB?.x.resample(mvaJ?); mvaJNBL?.name = ''JN*BL''; mvaJNBL?.units = ''nA/m^2 nT'';',ic)
c_eval('mvaJLBN? = -1*mvaJ?.x*mvaB?.z.resample(mvaJ?); mvaJLBN?.name = ''-JL*BN''; mvaJLBN?.units = ''nA/m^2 nT'';',ic)


% FPI
c_eval('mvaJLBM?ne_mVm = mvaJLBM?*1e-18/(ne?.resample(mvaJxB?)*1e6)/units.e*1e3; mvaJLBM?ne_mVm.units = ''mV/m'';',ic)
c_eval('mvaJMBL?ne_mVm = mvaJMBL?*1e-18/(ne?.resample(mvaJxB?)*1e6)/units.e*1e3; mvaJMBL?ne_mVm.units = ''mV/m'';',ic)
c_eval('mvaJMBN?ne_mVm = mvaJMBN?*1e-18/(ne?.resample(mvaJxB?)*1e6)/units.e*1e3; mvaJMBN?ne_mVm.units = ''mV/m'';',ic)
c_eval('mvaJNBM?ne_mVm = mvaJNBM?*1e-18/(ne?.resample(mvaJxB?)*1e6)/units.e*1e3; mvaJNBM?ne_mVm.units = ''mV/m'';',ic)
c_eval('mvaJNBL?ne_mVm = mvaJNBL?*1e-18/(ne?.resample(mvaJxB?)*1e6)/units.e*1e3; mvaJNBL?ne_mVm.units = ''mV/m'';',ic)
c_eval('mvaJLBN?ne_mVm = mvaJLBN?*1e-18/(ne?.resample(mvaJxB?)*1e6)/units.e*1e3; mvaJLBN?ne_mVm.units = ''mV/m'';',ic)


c_eval('mvaVeMBL? = -1*mvaVe?.y*1e3*mvaB?.x.resample(mvaVe?.time)*1e-9*1e3; mvaVeMBL?.name = ''-VeM x BL'';',ic) % mV/m
c_eval('mvaVeLBM? = +1*mvaVe?.x*1e3*mvaB?.y.resample(mvaVe?.time)*1e-9*1e3; mvaVeLBM?.name = ''VeL x BM'';',ic) % mV/m

c_eval('mvaEVexB? = mvaE?.resample(mvaVexB?) + mvaVexB?; mvaEVexB?.name = ''E + ve x B'';',ic) % mV/m

c_eval('tsLgse? = irf.ts_vec_xyz(iPDist?.time,repmat(L,iPDist?.length,1));',ic)
c_eval('tsMgse? = irf.ts_vec_xyz(iPDist?.time,repmat(M,iPDist?.length,1));',ic)
c_eval('tsNgse? = irf.ts_vec_xyz(iPDist?.time,repmat(N,iPDist?.length,1));',ic)

c_eval('tsLdsl? = mms_dsl2gse(tsLgse?,defatt?,-1);',ic)
c_eval('tsMdsl? = mms_dsl2gse(tsMgse?,defatt?,-1);',ic)
c_eval('tsNdsl? = mms_dsl2gse(tsNgse?,defatt?,-1);',ic)

%c_eval('tsLdsl? = tsLgse?;',ic)
%c_eval('tsMdsl? = tsMgse?;',ic)
%c_eval('tsNdsl? = tsNgse?;',ic)

%% Reduce distributions
nMovMean = 5;
elim = [200 Inf];
c_eval('fi?_L = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L);',ic)

nMovMean = 5;
elim = [00 Inf];
c_eval('fi?_L_000_5 = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L);',ic)

nMovMean = 5;
elim = [100 Inf];
c_eval('fi?_L_100 = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L);',ic)

nMovMean = 2;
elim = [100 Inf];
c_eval('fi?_L_100_mm2 = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L);',ic)


elim = [100 Inf];open 
c_eval('fi?_L_100_only = iPDist?.elim(elim).reduce(''1D'',L);',ic)

elim = [200 Inf];
c_eval('fi?_L_200_only = iPDist?.elim(elim).reduce(''1D'',L);',ic)

elim = [00 Inf];
c_eval('fi?_L_000_only = iPDist?.elim(elim).reduce(''1D'',L);',ic)


%%
nMovMean = 2;
elim = [200 Inf];
c_eval('fi?_L = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L);',ic)
c_eval('fi?_M = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',M);',ic)
c_eval('fi?_N = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',N);',ic)


%% Reduce distribution in direction of maximum pressure
[Prot,Trot] = maximum_shear_direction(mvaPi3);
e1 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,1,:)));
e2 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,2,:)));
e3 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,3,:)));


tsTrot = irf.ts_vec_xyz(gsePi3.time,[gsePi3.data(:,1,1), gsePi3.data(:,2,2), gsePi3.data(:,3,3)*0]);
e1 = tsTrot.norm;

c_eval('fi?_e1 = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',e1);',ic)


%pdist_rot2 = pdist.shift(squeeze(vel), 10, R2, 'mms');

%% Define times, etc... things that are common for the entire study
tint_figure = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:35:00.00Z');
time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');
time_xline_ion = irf_time('2017-07-11T22:34:02.00Z','utc>EpochTT');
v_xline = -170;
time_vdf = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
tint_figure_zoom = irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:30.00Z');
tint_figure_edr = irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:12.00Z');
tint_figure_zoom_inner_idr = irf.tint('2017-07-11T22:33:50.00Z/2017-07-11T22:34:20.00Z');
nMovMean = 5;

%% Figure: Overview

h = irf_plot(10);

fontsize = 12;

if 1 % B gse
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
end 
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % Pi
  hca = irf_panel('Pi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPi?.xx.tlim(tint),mvaPi?.yy.tlim(tint),mvaPi?.zz.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'P_i (nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % Pi
  hca = irf_panel('Pi off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPi?.xy.tlim(tint),mvaPi?.xz.tlim(tint),mvaPi?.yz.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'P_i (nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'},[0.98,0.98],'fontsize',fontsize);
end

if 0 % dEFlux ion
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
end
if 0 % dEFlux ion
  hca = irf_panel('ion dEF omni movmean rem onecounts');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.1],'fontsize',fontsize,'color','k');
end

if 1 % fi red L
  hca = irf_panel('fi L');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_L.specrec;',ic)
  %specrec.f = specrec.f - v_xline;
  c_eval('irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
    
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.x,''k-'')',ic)
  %c_eval('irf_plot(hca,irf.ts_scalar(mvaVi?.time,repmat(v_xline,[mvaVi?.length,1])),''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iL} (km/s)'};

  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 1 % fi red M
  hca = irf_panel('fi M');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_M.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.y,''k-'')',ic)
  hca.NextPlot = "replace";
  
  hca.YLabel.String = {'v_{iM} (km/s)'};
end
if 1 % fi red L
  hca = irf_panel('fi N');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_N.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iN} (km/s)'};
end
if 1 % max shear direction
  hca = irf_panel('fi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_e1.specrec,''donotfitcolorbarlabel'',''lin'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{i,Pmax} (km/s)'};
end
if 1 % e1
  hca = irf_panel('Vi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{e1.norm.x,e1.norm.y,e1.norm.z},''comp'');',ic)  
  
  hca.YLabel.String = {'e_1 (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end

irf_plot_axis_align
%irf_zoom(h,'x',tint_figure)
irf_zoom(h,'x',tint_figure_zoom)
irf_zoom(h(1:3),'y')
irf_pl_mark(h,time_xline,'black','linestyle',':')
%irf_pl_mark(h,time_xline_ion,'red','linestyle',':')
%colormap(irf_colormap('thermal'))
colormap([pic_colors('candy_gray'); 1 1 1])
irf_plot_axis_align
h(end).XTickLabelRotation = 0;
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
%hlinks1 = linkprop(h(3:4),{'CLim'});
%hlinks2 = linkprop(h(5:7),{'CLim'});
%hlinks2 = linkprop(h(4:6),{'CLim'});
h(1).Title.String = sprintf('N_{movmean} = %g',nMovMean);
c_eval('h(?).FontSize = 14;',1:numel(h))

if 0
  %%
  h(3).YLim = 1499*[-1 1];
  h(4).YLim = 1499*[-1 1];
  h(5).YLim = 1499*[-1 1];
  h(5).CLim = [-3 -1];
  colormap([0.9 0.9 0.9; pic_colors('candy4')])
end

%% Figure: Hall fields
ic = 3;
h = irf_plot(4);

fontsize = 12;
fhighE = 10;

if 1 % B lmn
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end 
if 1 % E lmn
  hca = irf_panel('E LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaE?.x.filt(0,fhighE,[],5),mvaE?.y.filt(0,10,[],5),mvaE?.z.filt(0,10,[],5)},''comp'');',ic)
  hca.YLabel.String = {'E (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic),sprintf('E < %g Hz',fhighE)},[0.08 0.98],'fontsize',fontsize,'color','k');
end 
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint)*1e-3,mvaVe?.y.tlim(tint)*1e-3,mvaVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end
if 0 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % JxB
  hca = irf_panel('JxB LMN');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{gseJxB1.z,gseJxB2.z,gseJxB3.z,gseJxB4.z,gseJxB.z},'comp');
  
  hca.YLabel.String = {'JxB_N (km/s)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'1','2','3','4','curl'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % JxB/ne mV/m
  hca = irf_panel('E Hall LMN');
  set(hca,'ColorOrder',mms_colors('1243b'))
  c_eval('irf_plot(hca,{mvaE?.z.filt(0,fhighE,[],5),mvaJxBne?_mVm.z,mvaJLBM?ne_mVm,mvaJMBL?ne_mVm},''comp'');',ic)
  
  hca.YLabel.String = {'E_N (mV/m)'};
  set(hca,'ColorOrder',mms_colors('1243b'))
  irf_legend(hca,{,sprintf('E (< %g Hz)',fhighE),'JxB/ne','J_LB_M/ne','-J_MB_L/ne'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end
if 0 % JxB/ne mV/m
  hca = irf_panel('JxB/ne LMN');
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_plot(hca,{mvaJxBne1_mVm.z,mvaJxBne2_mVm.z,mvaJxBne3_mVm.z,mvaJxBne4_mVm.z,mvaJxBne_mVm.z},'comp');
  
  hca.YLabel.String = {'JxB_N/ne (mV/m)'};
  set(hca,'ColorOrder',mms_colors('1234b'))
  irf_legend(hca,{'1','2','3','4','curl'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % E + vexB LMN
  hca = irf_panel('E + VexB LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaEVexB?.x.tlim(tint),mvaEVexB?.y.tlim(tint),mvaEVexB?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'E + v_e x B (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end
if 0 % E + vexB N
  hca = irf_panel('E + VexB N');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_patch(hca,{mvaEVexB?.z.tlim(tint),0});',ic)  
  
  hca.YLabel.String = {'E + v_e x B (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end
if 0 % E + vexB
  hca = irf_panel('ve and vExB LMN');
  colors = [mms_colors('xyz'); (mms_colors('xyz')*1).^0.5];
  set(hca,'ColorOrder',colors)
  c_eval('irf_plot(hca,{mvaVe?perp.x.tlim(tint),mvaVe?perp.y.tlim(tint),mvaVe?perp.z.tlim(tint)},''comp'');',ic)  
  hold(hca,'on')
  set(hca,'ColorOrder',colors(4:6,:))
  c_eval('htmp = irf_plot(hca,{mvaVExB?.x.resample(mvaVe?),mvaVExB?.y.resample(mvaVe?),mvaVExB?.z.resample(mvaVe?)},''comp'');',ic)  
  c_eval('htmp.Children(?).LineStyle = ''--'';',1:3)  
  hold(hca,'off')
  
  hca.YLabel.String = {'v_L (km/s)'};
  set(hca,'ColorOrder',colors)
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end
if 0 % E + vexB L
  hca = irf_panel('ve and vExB L');
  colors = [mms_colors('123'); (mms_colors('xyz')*1).^0.5];
  set(hca,'ColorOrder',colors)
  c_eval('irf_plot(hca,{mvaVExB?.x.resample(mvaVe?),mvaVe?perp.x.tlim(tint),mvaVe?},''comp'');',ic)    
  hca.YLabel.String = {'v_L (km/s)'};
  set(hca,'ColorOrder',colors)
  irf_legend(hca,{'v_{ExB}','v_{e\perp}','v_e'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end
if 0 % E + vexB M
  hca = irf_panel('ve and vExB M');
  colors = [mms_colors('123'); (mms_colors('xyz')*1).^0.5];
  set(hca,'ColorOrder',colors)
  c_eval('irf_plot(hca,{mvaVExB?.y.resample(mvaVe?),mvaVe?perp.y.tlim(tint),mvaVe?.y.tlim(tint)},''comp'');',ic)    
  hca.YLabel.String = {'v_M (km/s)'};
  set(hca,'ColorOrder',colors)
  irf_legend(hca,{'v_{ExB}','v_{e\perp}','v_e'},[0.98,0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('MMS %g',ic)},[0.08 0.98],'fontsize',fontsize,'color','k');
end

%irf_zoom(h,'x',tint_figure)
irf_zoom(h,'x',tint_figure_edr)
irf_zoom(h,'y')
irf_pl_mark(h,time_xline,'black','linestyle',':')
%irf_pl_mark(h,time_xline_ion,'red','linestyle',':')
%colormap(irf_colormap('thermal'))
colormap([pic_colors('candy_gray'); 1 1 1])
irf_plot_axis_align
h(end).XTickLabelRotation = 0;
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))

c_eval('h(?).FontSize = 16;',1:numel(h))

legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;

for ii = 1:numel(h)
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0])
  nInd = nInd + 1;
  h(ii).FontSize = 12;
end

if 0
  %%
  h(3).YLim = 1499*[-1 1];
  h(4).YLim = 1499*[-1 1];
  h(5).YLim = 1499*[-1 1];
  h(5).CLim = [-3 -1];
  colormap([0.9 0.9 0.9; pic_colors('candy4')])
end


%% Figure: Reduced distributions, 2D


%% Fit of reduced
tint_fit =  irf.tint('2017-07-11T22:33:20.00Z/2017-07-11T22:34:30.00Z');
tint_fit =  irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:20.00Z');
vdf = fi1_N.tlim(tint_fit);
%vdf = if1DN1_nobg.tlim(tint_fit);
vdf = vdf(1:1:vdf.length);
%vdf = vdf(fix(vdf.length/2));
%vdf = vdf([1 vdf.length]);
nPop = 3; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  0,      1000e3];


nPop = 2; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3];

% nPop = 2; % Two populations looks ok, but beams are warmer and slower to cover region around v=0
% X0 = [0.01e6, -1000e3, 500e3,...
%       0.01e6,  1000e3, 500e3];

if 0
nPop = 4; % Three populations give colder and faster counterstreaming beams
X0 = [0.01e6, -1000e3, 500e3,...
      0.01e6,  1000e3, 500e3,...
      0.01e6,  -1000e3,      1000e3,...
      0.01e6,  1000e3,      1000e3];
end
tic; [fitdata,ts] = funFitVDF(vdf,'nPop',nPop,'plot',0,'guessprevious',0,'X0',X0,'weight',repmat([0 1 1],1,nPop)); toc;


%% Plot fit
h = irf_plot(7);
hca = irf_panel('vdf_obs');
irf_spectrogram(hca,vdf.specrec('velocity'),'lin')
hca = irf_panel('vdf_fit');
irf_spectrogram(hca,ts.f.specrec('velocity'),'lin')
hold(hca,'on'); irf_plot(hca,ts.vd,':'); hold(hca,'off'); 
hca = irf_panel('n_fit');
irf_plot(hca,ts.n)
hca = irf_panel('vd_fit');
irf_plot(hca,ts.vd)
hca = irf_panel('T_fit');
irf_plot(hca,ts.T)
hca = irf_panel('T*n_fit');
irf_plot(hca,ts.T.*ts.n)
hca = irf_panel('cost function');
irf_plot(hca,ts.cf)

irf_plot_axis_align(h)
irf_zoom(h,'x',vdf.time([1 vdf.length]))
hlinks = linkprop([irf_panel('vdf_obs'),irf_panel('vdf_fit')],{'CLim','YLim'});
colormap(pic_colors('candy4'))



%% Ion trajectory from +vM to -vM, integrate in Harris current sheet, yz
units = irf_units;

% Define magnetic field
a = 1.5*1e3;
b = 150*1e3;
xvec = a*linspace(-60,60,500);
zvec = 1e3*linspace(-1000,1000,400);
yvec = linspace(-0,0,1);

[X,Y,Z] = ndgrid(xvec,yvec,zvec);

% Integrate orbits
B0 = 20e-9;
E0 = 3e-3;
m = units.mp;
q = units.e;
lz = .5*b;

Bx = @(x,y,z) B0*2*z/b^2*0;
Bx = @(x,y,z) B0*tanh(z/b);
By = @(x,y,z) B0*x*0;
Bz = @(x,y,z) B0*2*x/a^2;
Ex = @(x,y,z) E0*x*0;
Ey = @(x,y,z) 1*(E0*x*0 + E0*z*0 + E0);
Ez = @(x,y,z) 1*-70e-3*(z/(2*lz)).*exp(-(z/(2*lz)).^2);

Ay = @(x,y,z) (x/a).^2 - (z/b).^2;
Ay =  @(x,y,z) -b*B0*log(cosh(z/b)) + B0*(x/a).^2 + 0*y; %.*exp(-z.^2/b.^2)
AY0 = Ay(X,Y,Z);

%Bx = @(x,y,z) 2*z/b^2;
%Bz = @(x,y,z) 2*x/a^2;
%contour(X,Z,Ay(X,Z))
%hold on
%BX = Bx(X,Z);
%BZ = Bz(X,Z);
%quiver(X,Z,BX,BZ)

options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,xvec([1 end])),...
                 'AbsTol',1e-6);
%options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,1,0.1*a*[-1 1]),...
%                 'RelTol',1e-6,'AbsTol',1e-6);
options = odeset('Events', @(t,xyz) eom_box_edge(t,xyz,1,0.5*xvec([1 end])),...
                 'RelTol',1e-12,'AbsTol',1e-12);
EoM = @(t,xyz) eom(t,xyz,m,q,Ex,Ey,Ez,Bx,By,Bz); 
tstart = 0;
tstop = 8;

vz0 = E0/B0;
x_init_all = [0 0 0 0 -500e3 700e3;...
              ];

x_init_all = [0 0 -1000e3 0 10e3 E0/B0;...
              ];



clear p
for ip = 1:size(x_init_all,1)
  x_init = x_init_all(ip,:);
  [t,x_sol] = ode45(EoM,[tstart tstop],x_init,options); % 
  p(ip).t = t;
  p(ip).x = x_sol(:,1);
  p(ip).y = x_sol(:,2);
  p(ip).z = x_sol(:,3);
  p(ip).vx = x_sol(:,4);
  p(ip).vy = x_sol(:,5);
  p(ip).vz = x_sol(:,6);
end

h = setup_subplots(4,4,'vertical');
isub = 1;
dotsize = 10;



if 1 % Bx, Ez
  hca = h(isub); isub = isub + 1;
  plot(hca,zvec*1e-3,Bx(0,0,zvec)*1e9,zvec*1e-3,Ez(0,0,zvec)*1e3,zvec*1e-3,Ey(0,0,zvec)*1e3,'linewidth',1)  
  irf_legend(hca,{'B_L','E_N','E_M'},[0.02 0.98])
  hca.XLabel.String = 'N (km)';
end
if 1 % potential from Ez
  hca = h(isub); isub = isub + 1;
  plot(hca,zvec*1e-3,-cumtrapz(zvec,Ez(0,0,zvec)*1e3)/1000,'linewidth',1)     
  hold(hca,'on')
  plot(hca,zvec*1e-3,zvec*0 - 5000,'linewidth',1)
  hold(hca,'off')

  hca.YLabel.String = '\phi = -\int E_N dN (V)';
  hca.XLabel.String = 'N (km)';
end
if 1 % Jy from Bx
  hca = h(isub); isub = isub + 1;
  dz = zvec(2)-zvec(1); % m
  Jy = units.mu0^(-1)*gradient(Bx(0,0,zvec))/(dz); % Jy = dBx/dz
  plot(hca,zvec*1e-3,Jy*1e9,'linewidth',1)      
  hca.YLabel.String = 'J_M (nA/m^2)';
  hca.XLabel.String = 'N (km)';
end
if 1 % vey from Jy from Bx
  hca = h(isub); isub = isub + 1;
  dz = zvec(2)-zvec(1);
  Jy = units.mu0^(-1)*gradient(Bx(0,0,zvec))/(dz);
  vey = -Jy/units.e/0.04e6; % J = nev
  vey = vey*1e-6; % 10^3 km/s
  plot(hca,zvec*1e-3,vey,'linewidth',1)      
  hold(hca,'on')
  plot(hca,zvec*1e-3,zvec*0 - 18,'linewidth',1)
  hold(hca,'off')
  hca.YLabel.String = 'v_{eM} = -J_y/ne (10^3 km/s)';
  hca.XLabel.String = 'N (km)';
end
if 0 % particle
  hca = h(isub); isub = isub + 1;
  hca.NextPlot = "add";
  for ip = 1:numel(p)
    plot(hca,p(ip).y*1e-3,p(ip).z*1e-3,p(ip).y(1)*1e-3,p(ip).z(1)*1e-3,'go')
  end
  hca.NextPlot = "replaceall";
  hca.XLabel.String = 'L (km)';
  hca.YLabel.String = 'N (km)';
end
if 1 % particle, vy
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)
    scatter(hca,p(ip).y*1e-3,p(ip).z*1e-3,dotsize,p(ip).vy*1e-3)
  end
  hcb = colorbar(hca);
  hcb.YLabel.String = 'v_M (km/s)';
  colormap(hca,pic_colors('bluegrayred'))  
  hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'M (km)';
  hca.YLabel.String = 'N (km)';
end
if 1 % particle, |vz|
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)
    scatter(hca,p(ip).y*1e-3,p(ip).z*1e-3,dotsize,abs(p(ip).vz)*1e-3)
  end
  hcb = colorbar(hca);
  hcb.YLabel.String = '|v_N| (km/s)';
  colormap(hca,pic_colors('candy_gray'))
  hca.XLabel.String = 'M (km)';
  hca.YLabel.String = 'N (km)';
end
if 1 % particle, sqrt(vy^2+vz^2)
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)
    scatter(hca,p(ip).y*1e-3,p(ip).z*1e-3,dotsize,sqrt(p(ip).vy.^2 + p(ip).vz.^2)*1e-3)    
  end
  hcb = colorbar(hca);
  hcb.YLabel.String = '(v_M^2+v_N^2)^{1/2} (km/s)';
  colormap(hca,pic_colors('candy_gray'))
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'M (km)';
  hca.YLabel.String = 'N (km)';
end
if 1 % particle, eV
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)
    eV = units.mp*(p(ip).vy.^2 + p(ip).vz.^2)/2/units.eV;
    scatter(hca,p(ip).y*1e-3,p(ip).z*1e-3,dotsize,eV)    
  end
  hcb = colorbar(hca);
  hcb.YLabel.String = 'Energy (eV)';
  colormap(hca,pic_colors('candy_gray'))
  %hca.CLim = max(abs(hca.CLim))*[-1 1];
  hca.XLabel.String = 'M (km)';
  hca.YLabel.String = 'N (km)';
end
if 1 % particle, time, vm
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)    
    plot(hca,p(ip).t,p(ip).vy*1e-3)    
  end    
  hca.YLabel.String = 'v_M (km/s)';
  hca.XLabel.String = 'time (s)';
end
if 1 % particle, time, vn
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p) 
    plot(hca,p(ip).t,(p(ip).vz*1e-3))    
  end    
  hca.YLabel.String = 'v_N (km/s)';
  hca.XLabel.String = 'time (s)';
end
if 1 % particle, time, vtot
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)    
    plot(hca,p(ip).t,sqrt(p(ip).vy.^2+p(ip).vz.^2)*1e-3)    
  end    
  hca.YLabel.String = '(v_M^2+v_N^2)^{1/2} (km/s)';
  hca.XLabel.String = 'time (s)';
end
if 1 % particle, time
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)
    eV = units.mp*(p(ip).vy.^2 + p(ip).vz.^2)/2/units.eV;
    plot(hca,p(ip).t,eV)    
  end    
  hca.YLabel.String = 'Energy (eV)';
  hca.XLabel.String = 'time (s)';
end
if 1 % particle, color, phase space
  hca = h(isub); isub = isub + 1;
  for ip = 1:numel(p)    
    scatter(hca,p(ip).vy*1e-3,p(ip).vz*1e-3,dotsize,p(ip).z*1e-3)
  end    
  axis(hca,'equal')
  axis(hca,'square')
  hca.YLabel.String = 'v_M';
  hca.XLabel.String = 'v_N';
  hca.XLim = 1500*[-1 1];
  hca.YLim = 1500*[-1 1];
  hcb = colorbar(hca);
  colormap(hca,pic_colors('blue_red'))
  hca.CLim = 2*b*1e-3*[-1 1];
end

c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).Box = ''on'';',1:numel(h))
%c_eval('axis(h(?),''equal'');',2:4)
hlinks = linkprop(h(5:8),{'XLim','YLim'});

%% Single (or a few) 1D reduced ion population

fontsize = 14;
nMovMean = 7;
elim = [500 Inf];
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
dt = 0;
time = time + dt;
pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
c_eval('pdist_with_noise = iPDist?.movmean(nMovMean).tlim(time+0.5*0.15*[-1 1]);',ic)

vlim = [-1500 1500];

h = setup_subplots(1,1);

t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
Ldsl = Tdsl(1,:);
Mdsl = Tdsl(2,:);
Ndsl = Tdsl(3,:);

hca = h(1);

vdf = pdist.reduce('1D',Mdsl);
vdf_with_noise = pdist_with_noise.reduce('1D',Mdsl);
plot(hca,vdf.depend{1},vdf.data,vdf_with_noise.depend{1},vdf_with_noise.data)
hca.XLim = vlim;
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'f_i (s/m^4)';
hl = findobj(gcf,'type','line');
hl.LineWidth = 2;
hca.XGrid = 'on';
hca.YGrid = 'on';
t1t2 = pdist.time + 0.5*nMovMean*0.150*[-1 1];
irf_legend(hca,sprintf('%s-%s',t1t2(1).utc('HH:MM:SS.mmm'),t1t2(2).utc('HH:MM:SS.mmm')),[0.02,0.98],'fontsize',fontsize,'color','k');

%% Ion distributions, One dist at the time
%h = setup_subplots(2,2);

% 
elim = [00 Inf];
pdist_all = iPDist1.movmean(nMovMean,'RemoveOneCounts',iPDist1_counts).elim(elim);
time = time_vdf;

for dt = 15%-2:1:25
[h1,h] = initialize_combined_plot('topbottom',2,1,4,0.4,'vertical');

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
irf_plot(hca,{mvaVi1})
hca = h1(isub); isub = isub + 1;
irf_plot(hca,{mvaVExB1.resample(pdist)})
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];

isub = 1;


time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
time = time + dt;
pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
%pdist_nobg = iPDist3_nobg.tlim(time + 7*0.5*0.150*[-1 1]).elim([300 Inf]);

c_eval('B = mvaB?.resample(pdist);',ic)  
c_eval('E = mvaE?.resample(pdist);',ic)  
c_eval('vExB = mvaVExB?.resample(pdist);',ic)  
 

c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  
nSmooth = 1;
hca = h(isub); isub = isub + 1;
vdf = pdist.reduce('2D',M,L);
%vdf = pdist_nobg.reduce('2D',M,L);
vdf.plot_plane(hca','smooth',nSmooth,'vectors',{dmpaB1.resample(pdist),'B'})
axis(hca,'square')
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'v_L (km/s)';
if 1 % plot B direction
  %%
  xlim = hca.XLim;
  ylim = hca.YLim;
  hold(hca,'on')
  %dt_distx = 0.030;
  %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
  B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
  B_ = mean(B__.data,1);
  B_std = std(B__.data,1);
  b = B_/norm(B_);
  B_std_inplane = std(B__.data(:,1:2),1);
  B_inplane = sqrt(sum(B_(1:2).^2));
  b = b*2000;
  quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
    hold(hca,'off')
    hca.XLim = xlim;
    hca.YLim = ylim;
    B_ = B_(1:2);
end
if 0 % plot B direction
  %%
  xlim = hca.XLim;
  ylim = hca.YLim;
  hold(hca,'on')
  %dt_distx = 0.030;
  %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
  B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
  B_ = mean(B__.data,1);
  B_std = std(B__.data,1);
  b = B_/norm(B_);
  B_std_inplane = std(B__.data(:,1:2),1);
  B_inplane = sqrt(sum(B_(1:2).^2));
  if B_inplane > 2*norm(B_std_inplane)
    k = b(1)/b(2);
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
%    irf_legend(hca,sprintf('B_{%s} = %.1f nT', comps,B_inplane),[0.98 0.98],'color','k','fontsize',fontsize_B_amp)
  
end
if 1 % plot ExB
  %%
  hold(hca,'on')
  hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.x.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
  hold(hca,'off')    
end
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];



hca = h(isub); isub = isub + 1;
vdf = pdist.reduce('2D',[M_vi],[N_vi]);
%vdf = pdist_nobg.reduce('2D',[M],[N]);
vdf.plot_plane(hca','smooth',nSmooth)
axis(hca,'square')
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'v_N (km/s)';
%xlim = hca.XLim;
%ylim = hca.YLim;
if 1 % plot ExB
  %%
  hold(hca,'on')
  hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
  hold(hca,'off')    
end
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];


hca = h(isub); isub = isub + 1;
vdf = pdist.reduce('2D',[L_vi],[N_vi]);
%vdf = pdist_nobg.reduce('2D',[L],[N]);
vdf.plot_plane(hca','smooth',nSmooth)
axis(hca,'square')
hca.XLabel.String = 'v_L (km/s)';
hca.YLabel.String = 'v_N (km/s)';
if 1 % plot ExB
  %%
  hold(hca,'on')
  hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'+k')
  %plot(hca,mean(ve.y.data,1)*1e-3,mean(ve.z.data,1)*1e-3,'ow')
  hold(hca,'off')    
end
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];



times_exact{1} = vdf.time;

if 0 % with bg
  hca = h(isub); isub = isub + 1;
  vdf = pdist.reduce('2D',L,M);
  vdf.plot_plane(hca','smooth',2);
end

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill');
%hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
%hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
axis(hca,'square')
vlim = 2000;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];
hca.ZLim = vlim*[-1 1];
%camlight(gca,90,-45)
%view(hca,[-1 -1 0.5])
%view(hca,[1 0.2 0.2])
view(hca,[0 1 0.2])
camlight(gca,0,0)

%hca = h(isub); isub = isub + 1;
%hs = pdist.plot_isosurface(hca,'rotate',lmn);

c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
%cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end




%% Several at the same time.
%h = setup_subplots(2,2);
dt_all = [-15 -10 0 10 15];
dt_all = [-15:5:15]+0;
dt_all = [-15:5:15]+-1.5;
dt_all = [-6:2:6]+.5;
dt_all = [-6:2:6]+.5;
dt_all = [-6:3:6]+.5;

dt_all = [-6:3:6]+.5;
dt_all = [-30 -3 1 5 15]-.5;
dt_all = [-30 -3 1 5 15]-.5;
%dt_all = [-6:2:6]+0;
%dt_all = [-6:2:6]+25;
%dt_all = [-6:2:6]-00;
times_utc = ['2017-07-11T22:33:50.582Z';...
             '2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:07.940Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

[h1,h] = initialize_combined_plot('topbottom',2,3,numel(dt_all),0.2,'vertical');

vL_Xline = 1*-170;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{mvaVi?.x-vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf(['v_L-(%g) (km/s)'],vL_Xline),'M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{mvaVExB?.resample(iPDist?)})',ic)
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';


isub = 1;
% nSmooth = 1; % specified further doen


elim = [0 Inf];
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)
time = time_vdf;



vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];
elim = [200 Inf];

% time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
% time = time + dt;
times = EpochTT(times_utc);
for it = 1:times.length%(1)
  time = times(it)  ;
  
  pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
  tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
  
  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))

  if 0
    %%
    edges = -0.5:1:9;
    centers = edges(2:end)-0.5*(edges(2)-edges(1));
    N = histcounts(counts.data(:),edges);
    Nsum = histcounts(count_sum(:),edges);
    %pdist_1crem.data()
    hca = subplot(1,1,1);
    bar(centers,[N; Nsum]',2)
    legend(hca,{'Not summed','Summed'})
    hca.YScale = 'log';
  end

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)    
  
  nSmooth = 0;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    if 1 % plot B direction
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
      B_ = mean(B__.data,1);
      B_std = std(B__.data,1);
      b = B_/norm(B_);
      B_std_inplane = std(B__.data(:,1:2),1);
      B_inplane = sqrt(sum(B_(1:2).^2));
      b = b*2000;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end

  if 1 % f(L,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl],'vint',vint_M);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  

  nSmooth = 3;
  iso_values = 1*10.^[-27:-15];
  iso_values = [6.5e-28];
  iso_values = 20e-28;
  vlim = 3000;
  if 0 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',lmn);
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')
    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 0 1])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.ZLabel.String = 'v_N (km/s)';
    hca.Title = [];
  end
  if 0 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill');
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 1 0])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.ZLabel.String = 'v_N (km/s)';
    hca.Title = [];
  end
  if 0 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill');
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[1 0 0])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.ZLabel.String = 'v_N (km/s)';
    hca.Title = [];
  end
  

  %times_exact{1} = vdf.time;
  times_exact{1} = pdist.time;

  c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
  %cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end


c_eval('h(?).FontSize = 8;',1:numel(h))


%colormap(pic_colors('candy_gray'))
colormap(pic_colors('thermal'))

hlinks_LM = linkprop(h(1:numel(dt_all):end),{'CLim'});
hlinks_LN = linkprop(h(2:numel(dt_all):end),{'CLim'});
hlinks_MN = linkprop(h(3:numel(dt_all):end),{'CLim'});

%hlinks_LM = linkprop(h(1:3:end),{'View'});
%hlinks_LN = linkprop(h(2:3:end),{'View'});
%hlinks_MN = linkprop(h(3:3:end),{'View'});

%compact_panels(h,0.04,0.01)
drawnow
hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
%delete(hb(1:end-3))
delete(hb)
hcb = colorbar(h(numel(h)-2),'location','north');

compact_panels(h(1:3:end),0,0)
compact_panels(h(2:3:end),0,0)
compact_panels(h(3:3:end),0,0)

hlinks = linkprop(h,{'XLim','YLim','CLim'});


c_eval('h(?).YTickLabel = [];',4:numel(h))
c_eval('h(?).YLabel = [];',4:numel(h))

%c_eval('h(?).YTickLabel = [];',4:numel(h))
%c_eval('h(?).YLabel = [];',4:numel(h))
c_eval('h(?).XTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).YTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))


c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).GridLineWidth = 1;',1:numel(h))

c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

%c_eval('h1(?).Position(3) = 0.3;',1:numel(h1))
c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

h1(1).Position = [0.1700    0.8530    0.3000    0.0848];
h1(2).Position = [0.1700    0.7672    0.3000    0.0848];

%hca.CLim = [-10 -7];

%% Several at the same time. forces on VDF
%h = setup_subplots(2,2);
dt_all = [-15  -10 0 10 15];
dt_all = [-15:5:15]+0;
dt_all = [-15:5:15]+-1.5;
%dt_all = [-6:2:6]+0;
%dt_all = [-6:2:6]+25;
%dt_all = [-6:2:6]-00;

[h1,h] = initialize_combined_plot('topbottom',2,3,numel(dt_all),0.2,'vertical');

vL_Xline = -0;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVi3.x-vL_Xline,mvaVi3.y,mvaVi3.z},'comp')
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf(['v_L-(%g) (km/s)'],vL_Xline),'M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVExB3.resample(iPDist3)})
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';


isub = 1;
%nSmooth = 1;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];
elim = [200 Inf];

for dt = dt_all%(1)

  time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
  time = time + dt;
  tint_dist = time + 10*0.5*0.150*[-1 1];
  %pdist = iPDist3.tlim(tint_dist).elim([600 Inf]);
  %pdist = iPDist3_nobg.tlim(time + 2*0.5*0.150*[-1 1]).elim([200 Inf]);
  
  counts = iPDist3_counts.tlim(tint_dist);
  counts.data(isnan(counts.data)) = 0;  
  count_sum = sum(counts.data,1);

  pdist_1crem = iPDist3.tlim(tint_dist);
  pdist_1crem.data(:,count_sum<1.5) = 0;
  pdist_1crem = pdist_1crem.elim(elim);
  pdist = pdist_1crem;

  %pdist = iPDist3_nobg.tlim(tint_dist).elim(elim);
  %pdist = iPDist3.tlim(tint_dist).elim(elim);
  if 0
    %%
    edges = -0.5:1:9;
    centers = edges(2:end)-0.5*(edges(2)-edges(1));
    N = histcounts(counts.data(:),edges);
    Nsum = histcounts(count_sum(:),edges);
    %pdist_1crem.data()
    hca = subplot(1,1,1);
    bar(centers,[N; Nsum]',2)
    legend(hca,{'Not summed','Summed'})
    hca.YScale = 'log';
  end

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)    
  
  nSmooth = 1;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    %vdf.depend{1} = vdf.depend{1} - vL_Xline;
    %vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf_data = squeeze(mean(vdf.data,1));
    
    vL_data = squeeze(vdf(1).depend{1});
    vM_data = squeeze(vdf(1).depend{2});
    fEN = E.resample(pdist).z.data*1e-3;
    fBM = B.resample(pdist).y.data*1e-9;
    fBL = B.resample(pdist).x.data*1e-9;
    [VL,VM] = ndgrid(vL_data,vM_data);
    %[id,VL,VM] = ndgrid(1:pdist.length,vL_data,vM_data);
    force_data = mean(fEN) + VL*1e3*mean(fBM) - VM*1e3*mean(fBL);
    force_data(vdf_data==0) = 0;
    % f = EN + vL*BM - vM*BL;

    %vdf.plot_plane(hca,'smooth',nSmooth)
    pcolor(hca,vL_data,vM_data,force_data')
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';

    colormap(hca,pic_colors('blue_red'))
    shading(hca,'flat');

    hold(hca,'on')
    contour(hca,vL_data,vM_data,log10(vdf_data)','color','k')
    hold(hca,'off')

    if 1 % plot B direction      
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
      B_ = mean(B__.data,1);
      B_std = std(B__.data,1);
      b = B_/norm(B_);
      B_std_inplane = std(B__.data(:,1:2),1);
      B_inplane = sqrt(sum(B_(1:2).^2));
      b = b*2000;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(hca,-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end

  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    if 1 % plot B direction
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
      B_ = mean(B__.data,1);
      B_std = std(B__.data,1);
      b = B_/norm(B_);
      B_std_inplane = std(B__.data(:,1:2),1);
      B_inplane = sqrt(sum(B_(1:2).^2));
      b = b*2000;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end

  if 1 % f(L,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl],'vint',vint_M);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 0 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  

  nSmooth = 3;
  iso_values = 1*10.^[-27:-15];
  iso_values = [6.5e-28];
  iso_values = 20e-28;
  vlim = 3000;
  if 0 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill','rotate',lmn);
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')
    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 0 1])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.ZLabel.String = 'v_N (km/s)';
    hca.Title = [];
  end
  if 0 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill');
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 1 0])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.ZLabel.String = 'v_N (km/s)';
    hca.Title = [];
  end
  if 0 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'val',iso_values,'smooth',nSmooth,'fill');
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[1 0 0])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    hca.ZLabel.String = 'v_N (km/s)';
    hca.Title = [];
  end
  

  times_exact{1} = vdf.time;

  c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
  %cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end


c_eval('h(?).FontSize = 8;',1:numel(h))


colormap(pic_colors('candy_gray'))

%hlinks_LM = linkprop(h(1:numel(dt_all):end),{'CLim'});
%hlinks_LN = linkprop(h(2:numel(dt_all):end),{'CLim'});
%hlinks_MN = linkprop(h(3:numel(dt_all):end),{'CLim'});

hlinks_LM = linkprop(h(1:3:end),{'View'});
hlinks_LN = linkprop(h(2:3:end),{'View'});
hlinks_MN = linkprop(h(3:3:end),{'View'});

%compact_panels(h,0.04,0.01)
drawnow
hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
delete(hb(1:end-3))

%c_eval('h(?).YTickLabel = [];',4:numel(h))
%c_eval('h(?).YLabel = [];',4:numel(h))
c_eval('h(?).XTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).YTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))


c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).GridLineWidth = 1;',1:numel(h))

c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

%c_eval('h1(?).Position(3) = 0.3;',1:numel(h1))
c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

h1(1).Position = [0.1700    0.8530    0.3000    0.0848];
h1(2).Position = [0.1700    0.7672    0.3000    0.0848];





%% One time but with more overview panels
dt_all = [-30:2:30];

%dt_all = 0;

[h1,h] = initialize_combined_plot('leftright',3,2,2,0.4,'vertical');
%[h1,h] = initialize_combined_plot('topbottom',2,3,numel(dt_all),0.2,'vertical');

isub = 1;

tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVi3})
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))},[0.01 1.02]);
hca.YLabel.Interpreter = 'tex';

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaVExB3.resample(iPDist3)})
irf_zoom(h1,'x',tint_zoom)
%h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';


hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
irf_plot(hca,{mvaE3.resample(iPDist3)})
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';
 
for dt = dt_all%(1)
%%
  delete(findobj(gcf,'type','patch'))
  isub = 1;
  nSmooth = 1;
  elim = [600 Inf];
  time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
  time = time + dt;
  tint_dist = time + 10*0.5*0.150*[-1 1];
  pdist = iPDist3.tlim(tint_dist).elim([600 Inf]);
  %pdist = iPDist3_nobg.tlim(time + 2*0.5*0.150*[-1 1]).elim([200 Inf]);
  
  counts = iPDist3_counts.tlim(tint_dist);
  counts.data(isnan(counts.data)) = 0;  
  count_sum = sum(counts.data,1);

  pdist_1crem = iPDist3.tlim(tint_dist);
  pdist_1crem.data(:,count_sum<1.5) = 0;
  pdist_1crem = pdist_1crem.elim(elim);
  pdist = pdist_1crem;

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)    

  isub = 1;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',[-inf inf]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    if 1 % plot B direction
      %%
      xlim = hca.XLim;
      ylim = hca.YLim;
      hold(hca,'on')
      %dt_distx = 0.030;
      %B_ = B.tlim(dist.time([1 end]) + 0.5*dt_dist*[-1 1]);
      B__ = B.tlim(pdist.time([1 end]) + 0.5*0.03*[-1 1]);
      B_ = mean(B__.data,1);
      B_std = std(B__.data,1);
      b = B_/norm(B_);
      B_std_inplane = std(B__.data(:,1:2),1);
      B_inplane = sqrt(sum(B_(1:2).^2));
      b = b*2000;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',1)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  isub = 2;
  if 1 % f(L,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  isub = 4;
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Mdsl],[Ndsl]);
    vdf.plot_plane(hca','smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      %%
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  times_exact{1} = vdf.time;
  

  nSmooth = 3;
  iso_values = 1*10.^[-27:-15];
  iso_values = [6.5e-28];
  iso_values = 7e-28;
  vlim = 2500;

  isub = 3;
  if 1 % isuorface
    hca = h(isub); isub = isub + 1;
    hca.ColorOrder = pic_colors('matlab');
    hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill','rotate',lmn);
    %hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
    %hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
    c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
    axis(hca,'square')
    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.ZLim = vlim*[-1 1];
    %camlight(gca,90,-45)
    %view(hca,[-1 -1 0.5])
    view(hca,[0 0 1])
    view(hca,[2 -1 0.2])
    camlight(gca,0,0)
    
    %h(isub-4).Title = hca.Title;
    %h(isub-4).Title.FontSize = 8;
    hca.Title = [];
  end
  colormap(pic_colors('candy_gray'))
  c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
  drawnow
  cn.print(sprintf('torbert_fi_ref_dt_%02.0f',dt-dt_all(1)))
end



%% 3D isosurface, to illustrate 3D tilt



elim = [2000 Inf];
nMovMean = 7;
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts);',ic)
pdist_all.data(:,1:20,:,:) = 0;
%pdist_all = pdist_all.elim(elim);

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
dt = -3;
time = time + dt;

h = setup_subplots(1,1);
isub = 1;

pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
%iso_values = [6.5e-27];
hs = pdist.plot_isosurface(hca,'vals',iso_values,'smooth',nSmooth,'fill');
%hs = pdist_nobg.plot_isosurface(hca,'smooth',nSmooth);
%hs = pdist.plot_isosurface(hca,'smooth',nSmooth,'fill','rotate',lmn);
c_eval('hs.Patch(?).FaceAlpha = 1;',1:numel(hs.Patch))
axis(hca,'square')
vlim = 2500;
hca.XLim = vlim*[-1 1];
hca.YLim = vlim*[-1 1];
hca.ZLim = vlim*[-1 1];
%camlight(gca,90,-45)
%view(hca,[-1 -1 0.5])
%view(hca,[1 0.2 0.2])
view(hca,[0 1 0.2])
camlight(gca,0,0)

if 1 % Add max energy
  %%
  Emax = max(pdist_all.depend{1}(:));
  %Emax = 32000;
  vmax = sqrt(2*units.eV*Emax/units.mp)*1e-3;
  
  [X,Y,Z] = sphere(20);
  hold(hca,'on')
  hsurf = surf(hca,vmax*X,vmax*Y,vmax*Z,'facealpha',0.1,'facecolor',[0 0 0],'EdgeColor','none');
  hold(hca,'off')
end


%% Figure: EGU

h = irf_plot(4);

fontsize = 12;

if 0 % B gse
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
end 
if 0 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  

  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,-170*[1 1]),'k--')
  hca.XGrid = 'off'; hca.YGrid = 'off';
  hold(hca,'off')

  if 0
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,0*[1 1]),'color',[0.5 0.5 0.5])
  hca.XGrid = 'off'; hca.YGrid = 'off';
  hold(hca,'off')
  end
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);


end
if 0 % Pi
  hca = irf_panel('Pi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPi?.xx.tlim(tint),mvaPi?.yy.tlim(tint),mvaPi?.zz.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'P_i (nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 0 % Pi
  hca = irf_panel('Pi off LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaPi?.xy.tlim(tint),mvaPi?.xz.tlim(tint),mvaPi?.yz.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'P_i (nPa)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'LM','LN','MN'},[0.98,0.98],'fontsize',fontsize);
end

if 0 % dEFlux ion
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
end
if 0 % dEFlux ion
  hca = irf_panel('ion dEF omni movmean rem onecounts');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.1],'fontsize',fontsize,'color','k');
end

if 1 % fi red L
  hca = irf_panel('fi L');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_L.specrec;',ic)
  specrec.p(specrec.p==0) = NaN;
  %specrec.f = specrec.f - v_xline;
  c_eval('[hs, hcb] = irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  hcb.YLabel.String = 'f_i(v_L) (s/m^4)';
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
    
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.x,''k-'')',ic)
  %c_eval('irf_plot(hca,irf.ts_scalar(mvaVi?.time,repmat(v_xline,[mvaVi?.length,1])),''k-'')',ic)
  hca.NextPlot = "replace";


  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,-170*[1 1]),'k--')
  hca.XGrid = 'off'; hca.YGrid = 'off';
  hold(hca,'off')


  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,0*[1 1]),'color',[0.5 0.5 0.5])
  hca.XGrid = 'off'; hca.YGrid = 'off';
  hold(hca,'off')

  hca.YLabel.String = {'v_{iL} (km/s)'};
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 1 % fi red M
  hca = irf_panel('fi M');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_M.specrec;',ic)
  specrec.p(specrec.p==0) = NaN;
  c_eval('[hs, hcb] = irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  hcb.YLabel.String = 'f_i(v_M) (s/m^4)';
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.y,''k-'')',ic)
  hca.NextPlot = "replace";
  

  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,0*[1 1]),'color',[0.5 0.5 0.5])
  hca.XGrid = 'off'; hca.YGrid = 'off';
  hold(hca,'off')

  hca.YLabel.String = {'v_{iM} (km/s)'};
end
if 1 % fi red L
  hca = irf_panel('fi N');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_N.specrec;',ic)
  specrec.p(specrec.p==0) = NaN;
  c_eval('[hs, hcb] = irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  hcb.YLabel.String = 'f_i(v_N) (s/m^4)';
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";


  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,0*[1 1]),'color',[0.5 0.5 0.5])
  hca.XGrid = 'off'; hca.YGrid = 'off';
  hold(hca,'off')

  hca.YLabel.String = {'v_{iN} (km/s)'};
end
if 0 % max shear direction
  hca = irf_panel('fi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_e1.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iN} (km/s)'};
end
if 0 % e1
  hca = irf_panel('Vi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{e1.x,e1.y,e1.z},''comp'');',ic)  
  
  hca.YLabel.String = {'e_1 (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end




%irf_zoom(h,'x',tint_figure)
irf_zoom(h,'x',tint_figure_zoom)
tint_figure_zoom_inner_idr = irf.tint('2017-07-11T22:33:50.00Z/2017-07-11T22:34:15.00Z');
%irf_zoom(h,'x',tint_figure_zoom_inner_idr)



irf_zoom(h(1),'y')
h(1).YLim = [-600 600]-170;

irf_pl_mark(h,time_xline,'black','linestyle','--')
%irf_pl_mark(h,time_xline_ion,'red','linestyle',':')
%colormap(irf_colormap('thermal'))
irf_plot_axis_align
h(end).XTickLabelRotation = 0;
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
hlinks1 = linkprop(h(2:4),{'CLim'});
h(2).CLim = [-5 -1];
h(2).CLim = log10(prctile((specrec.p(:)),[3 99]));
colormap([pic_colors('candy4')])

hlinks2 = linkprop(h(2:4),{'YLim'});
h(3).YLim = 2200*[-1 1];
%hlinks2.Targets(1).YLim = [-1500 1500];

%hlinks1 = linkprop(h(3:4),{'CLim'});
%hlinks2 = linkprop(h(5:7),{'CLim'});
%hlinks2 = linkprop(h(4:6),{'CLim'});
h(1).Title.String = sprintf('N_{movmean} = %g',nMovMean);
c_eval('h(?).FontSize = 14;',1:numel(h))
colormap(pic_colors('candy4'))

%hb = findobj(gcf,'type','colorbar'); %hb = hb(end:-1:1);
%hb(3).YLabel.String = 'f_i(v_L) (s/m^4)';
%hb(2).YLabel.String = 'f_i(v_M) (s/m^4)';
%hb(1).YLabel.String = 'f_i(v_N) (s/m^4)';

%%