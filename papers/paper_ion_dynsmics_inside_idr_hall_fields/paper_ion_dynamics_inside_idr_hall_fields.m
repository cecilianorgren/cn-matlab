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
%%
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

% intertial length
c_eval('wpi? = (ne?*1e6*units.e^2/units.mp/units.eps0).^0.5;',ic)
c_eval('di? = 1e-3*(units.c/wpi?);',ic)

c_eval('vA? = (gseB?.resample(ne?).abs2*1e-18/units.mp/units.mu0/(ne?*1e6)).^0.5*1e-3;',ic)



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
M_gse = [0 1 -0.2]; M_gse = cross(L_gse,cross(M_gse,L_gse)); M_gse = M_gse/norm(M_gse);
%N_gse = [0 0 1];
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

lmn = lmn_vi;
lmn = lmn_gse;
%lmn = lmn_edr;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);

%% Rotate into LMN
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


c_eval('tsSCaxis?_dsl = irf.ts_vec_xyz(iPDist?.time,repmat([0 0 1],iPDist?.length,1));',ic)
c_eval('tsSCaxis?_gse = mms_dsl2gse(tsSCaxis?_dsl,defatt?,1);',ic)
c_eval('tsSCaxis?_lmn = tsSCaxis?_gse*lmn'';',ic)

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


elim = [100 Inf]; 
c_eval('fi?_L_100_only = iPDist?.elim(elim).reduce(''1D'',L);',ic)

elim = [200 Inf];
c_eval('fi?_L_200_only = iPDist?.elim(elim).reduce(''1D'',L);',ic)

elim = [00 Inf];
c_eval('fi?_L_000_only = iPDist?.elim(elim).reduce(''1D'',L);',ic)


%%
nMovMean = 5;
elim = [200 Inf];
%c_eval('fi?_L = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L);',ic)
%c_eval('fi?_M = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',M);',ic)
%c_eval('fi?_N = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',N);',ic)
data_tmp = load('/Users/cecilia/Data/Matlab/fi_reduced_N=5_elow=200.mat');

%%
%v_cut_L = -170;1
%c_eval('fi?_L_pos_neg = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',L,''vg_edges'',v_cut_L+[-3000, 0, 3000]);',ic)
load('/Users/cecilia/Data/Matlab/fi_reduced_posneg.mat')

% find_time_dependent_elow_for_iPDist.m

load('/Users/cecilia/Data/Matlab/fi_reduced_elow.mat')
%c_eval('fi?_L_elow = iPDist?.reduce(''1D'',L,''lowerelim'',tsElow);',ic)
%c_eval('fi?_M_elow = iPDist?.reduce(''1D'',M,''lowerelim'',tsElow);',ic)
%c_eval('fi?_N_elow = iPDist?.reduce(''1D'',N,''lowerelim'',tsElow);',ic)

%c_eval('fi?_M_elow_vnpos = iPDist?.reduce(''1D'',M,''lowerelim'',tsElow,''vint'',[-Inf Inf],''vyint'',[500 Inf]);',ic)
%%

PD_orig = iPDist3;
nMean = [3,3,3,3]; nThresh = 3;
PD_clean = PD_orig.remove_noise(nMean,nThresh,iPDist3_counts);
PD_diff = PD_orig+-PD_clean;

tsElow = PD_orig.find_noise_energy_limit(5).movmean(30);
emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit

PD_orig_notmasked = PD_orig;
PD_clean_notmasked = PD_clean;
PD_diff_notmasked = PD_diff;

PD_orig = PD_orig.mask({emask_mat});
PD_clean = PD_clean.mask({emask_mat});
PD_diff = PD_diff.mask({emask_mat});

c_eval('fi?_L = PD_clean.reduce(''1D'',L);',ic)
c_eval('fi?_M = PD_clean.reduce(''1D'',M);',ic)
c_eval('fi?_N = PD_clean.reduce(''1D'',N);',ic)

v_cut_L = -170;
c_eval('fi?_L_pos_neg = PD_clean.reduce(''1D'',L,''vg_edges'',v_cut_L+[-3000, 0, 3000]);',ic)


%% Reduce distribution in direction of maximum pressure
[Prot,Trot] = maximum_shear_direction(mvaPi3);
e1 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,1,:)));
e2 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,2,:)));
e3 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,3,:)));



tsTrot = irf.ts_vec_xyz(gsePi3.time,[gsePi3.data(:,1,1), gsePi3.data(:,2,2), gsePi3.data(:,3,3)*0]);
e1 = tsTrot.norm;

c_eval('fi?_e1 = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim).reduce(''1D'',e1);',ic)


%pdist_rot2 = pdist.shift(squeeze(vel), 10, R2, 'mms');

%% Reduce distribution in direction of maximum eigenvalue
nMovMean = 7;
elim = [500 Inf];
c_eval('PD_cleaned = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)
c_eval('PD = iPDist?.movmean(2).elim(elim);',ic)
p_c = mms.psd_moments(PD_cleaned,scPot3,'energyrange',elim);

tsPmoms_dsl = p_c.P_psd;
%tsPmoms_gse = mms_dsl2gse(tsPmoms_dsl,defatt3,1);
%mvaPi3_moms = lmn*tsPmoms_gse*lmn';
mvaPi3_c = lmn*tsPmoms_dsl*lmn';

[tsEig_val_c, tsEig_v1_c, tsEig_v2_c] = mvaPi3_c.eig([1 2]);

v1_rot = tsEig_v1_c.data;
v1_rot(v1_rot(:,2)<0,:) = -v1_rot(v1_rot(:,2)<0,:);

tsEig_v1_c_3D = irf.ts_vec_xyz(tsEig_v1_c.time,[v1_rot v1_rot(:,1)*0]);

%e1 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,1,:)));
%e2 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,2,:)));
%e3 = irf.ts_vec_xyz(gsePi3.time,squeeze(Trot(:,3,:)));

%tsTrot = irf.ts_vec_xyz(gsePi3.time,[gsePi3.data(:,1,1), gsePi3.data(:,2,2), gsePi3.data(:,3,3)*0]);
%e1 = tsTrot.norm;

c_eval('fi?_e1 = PD_cleaned.elim(elim).reduce(''1D'',tsEig_v1_c_3D);',ic)
c_eval('fi?_M = PD_cleaned.elim(elim).reduce(''1D'',M);',ic)
c_eval('fi?_e1_ = PD.elim(elim).reduce(''1D'',M);',ic)
c_eval('fi?_M_ = PD.elim(elim).reduce(''1D'',M);',ic)

%% Randomly remove single counts based on an average number of counts per energy level at low energies
tint = time_xline + [-10 10];
elim = [0 500];
pdist_counts_elim = iPDist3_counts.elim(elim).tlim(tint);
pdist_count = iPDist3_counts.elim([0 Inf]).tlim(tint);

C_TE_elim = nansum(pdist_counts_elim.data,[3 4]); % sum over angles
C_TE = nansum(pdist_counts.data,[3 4]); % sum over angles

h = setup_subplots(2,1);
isub = 1;

hca = h(isub); isub = isub + 1;
plot(hca,movmean(C_TE,10,1))


hca = h(isub); isub = isub + 1;
plot(hca,movmean(C_TE_elim,10,1))

hlinks = linkprop(h,{'YLim'});

%% Define times, etc... things that are common for the entire study
tint_figure = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:35:00.00Z');
%time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');
time_xline = irf_time('2017-07-11T22:34:02.60Z','utc>EpochTT');
time_xline_ion = irf_time('2017-07-11T22:34:02.00Z','utc>EpochTT');
v_xline = -170;
time_vdf = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
tint_figure_zoom = irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:30.00Z');
tint_figure_zoom_incl_sep = irf.tint('2017-07-11T22:33:23.00Z/2017-07-11T22:34:30.00Z');
tint_figure_edr = irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:12.00Z');
tint_figure_zoom_inner_idr = irf.tint('2017-07-11T22:33:50.00Z/2017-07-11T22:34:20.00Z');
nMovMean = 5;

%% Calculate new shifted PDist (to be able to make proper pitchangles and azimuthala angles in NL plane)
aa = pdist_cleaned(1).shift([-170 0 0],200,[1 0 0; 0 1 0; 0 0 1],ic);

%% Calculate moments based on pdist with removed one counts. Then get eigenvectors from that and reduce in that direction.
nMovMean = 7;
c_eval('pdist_cleaned = iPDist?.movmean(nMovMean,''removeonecounts'',iPDist?_counts).tlim(tint);',ic)
c_eval('pdist_movmean = iPDist?.movmean(nMovMean).tlim(tint);',ic)
c_eval('pdist_diff = pdist_cleaned; pdist_diff.data = pdist_movmean.data - pdist_cleaned.data;',ic)
p_c = mms.psd_moments(pdist_cleaned,scPot3,'energyrange',[300 Inf]);
pmoms_m = mms.psd_moments(pdist_movmean,scPot3,'energyrange',[300 Inf]);
p_d = mms.psd_moments(pdist_diff,scPot3,'energyrange',[300 Inf]);

tsPmoms_dsl = p_c.P_psd;
tsPmoms_gse = mms_dsl2gse(tsPmoms_dsl,defatt3,1);

mvaPi3_moms = lmn*tsPmoms_gse*lmn';
mvaPi3_moms = lmn*tsPmoms_dsl*lmn';

mvaPi3_moms.data = (mvaPi3_moms.data + permute(mvaPi3_moms.data,[1 3 2]))/2;

%% Calculate pitchangle of ions with respect to current sheet normal
lmn1 = lmn_vi;
lmn2 = lmn_gse;
lmn3 = lmn_edr;

iPitch3_c_1 = pdist_cleaned.pitchangles(irf.ts_vec_xyz(pdist_cleaned.time,repmat(lmn1(3,:),pdist_cleaned.length,1)),16);
iPitch3_c_2 = pdist_cleaned.pitchangles(irf.ts_vec_xyz(pdist_cleaned.time,repmat(lmn2(3,:),pdist_cleaned.length,1)),16);
iPitch3_c_3 = pdist_cleaned.pitchangles(irf.ts_vec_xyz(pdist_cleaned.time,repmat(lmn3(3,:),pdist_cleaned.length,1)),16);

pdist_cleaned_shifted = pdist_cleaned.shift([-170 0 0],200,[1 0 0; 0 1 0; 0 0 1],ic);
iPitch3_c_1_s = pdist_cleaned_shifted.pitchangles(irf.ts_vec_xyz(pdist_cleaned.time,repmat(lmn1(3,:),pdist_cleaned.length,1)),16);
iPitch3_c_2_s = pdist_cleaned_shifted.pitchangles(irf.ts_vec_xyz(pdist_cleaned.time,repmat(lmn2(3,:),pdist_cleaned.length,1)),16);
iPitch3_c_3_s = pdist_cleaned_shifted.pitchangles(irf.ts_vec_xyz(pdist_cleaned.time,repmat(lmn3(3,:),pdist_cleaned.length,1)),16);

%% Calculate % of ions moving towards or away from the X line
tint_left = [fi3_L.time(1) time_xline_ion];
tint_right = [time_xline_ion fi3_L.time(fi3_L.length)];

%[f_left,ind] = fi3_L.tlim([fi3_L.time(1) time_xline_ion]);
fi3_L_neg = fi3_L.data(:,1:size(fi3_L.data,2)/2);
fi3_L_pos = fi3_L.data(:,size(fi3_L.data,2)/2+1:end);
[ind_left,f_left] = fi3_L.time.tlim([fi3_L.time(1) time_xline_ion]);
[ind_right,f_right] = fi3_L.time.tlim([time_xline_ion fi3_L.time(fi3_L.length)]);

if 0
sum_all_left = sum(fi3_L_elow.data(ind_left,:),2);
sum_all_right = sum(fi3_L_elow.data(ind_right,:),2);
sum_left = sum(fi3_L_neg(ind_left,:),2);
sum_right = sum(fi3_L_neg(ind_right,:),2); 
clear ratio_towa ratio_away
ratio_towa(ind_left) = sum_left./sum_all_left;
ratio_towa(ind_right) = sum_right./sum_all_right;
ratio_away(ind_left) = sum_left./sum_all_left;
ratio_away(ind_right) = sum_right./sum_all_right;
end

if 0
fi3_L_neg = fi3_L_elow.data(:,1:47);
fi3_L_pos = fi3_L_elow.data(:,48:end);


sum_all = sum(fi3_L_elow.data,2);
sum_neg = sum(fi3_L_neg(:,:),2);
sum_pos = sum(fi3_L_pos(:,:),2); 
ratio_neg = sum_neg./sum_all;
ratio_pos = sum_pos./sum_all;
tsFneg = irf.ts_scalar(fi3_L.time,ratio_neg);
tsFpos = irf.ts_scalar(fi3_L.time,ratio_pos);

tsFaway = tsFneg.tlim(tint_left).combine(tsFpos.tlim(tint_right)); tsFaway.name = 'Away';
tsFtowa = tsFpos.tlim(tint_left).combine(tsFneg.tlim(tint_right)); tsFtowa.name = 'Towards';

%pos = fi3_L_pos_neg.tlim(tint_left).data(:,1)./

tsFLratio = irf.ts_scalar(fi3_L_elow.time,f_L_ratio);
end


sum_neg = fi3_L_pos_neg.data(:,1);
sum_pos = fi3_L_pos_neg.data(:,2);
sum_all = sum_neg + sum_pos;
ratio_neg = sum_neg./sum_all;
ratio_pos = sum_pos./sum_all;
tsFneg = irf.ts_scalar(fi3_L.time,ratio_neg);
tsFpos = irf.ts_scalar(fi3_L.time,ratio_pos);

tsFaway = tsFneg.tlim(tint_left).combine(tsFpos.tlim(tint_right)); tsFaway.name = 'Away';
tsFtowa = tsFpos.tlim(tint_left).combine(tsFneg.tlim(tint_right)); tsFtowa.name = 'Towards';

%% Calculate eigenvalues and vector for P
PP = gsePi3;
%PP = mvaPi3_moms;
%PP = pmoms_c.P_psd;
%PP = mvaP_c;
PP = mvaPi3;

nt = PP.length;
all_eig_vals = zeros(nt,3);
all_eig_vec1 = zeros(nt,3);
all_eig_vec2 = zeros(nt,3);
all_eig_vec3 = zeros(nt,3);
nMovMean = 1;
data_P = movmean(PP.data,nMovMean,1);
data_P = PP.data;
for it = 1:nt
  [V,D] = eig(squeeze(data_P(it,:,:)));
  all_eig_vals(it,:) = diag(D);
  all_eig_vec1(it,:) = V(:,1);
  all_eig_vec2(it,:) = V(:,2);
  all_eig_vec3(it,:) = V(:,3);
end
tsP = irf.ts_tensor_xyz(PP.time,data_P); tsP.name = 'P (smoothed)';
tsP_eig_vals = irf.ts_scalar(PP.time,all_eig_vals);  tsP_eig_vals.name = 'Eigenvalues';
tsP_eig_vec1 = irf.ts_vec_xyz(PP.time,all_eig_vec1); tsP_eig_vec1.name = 'V_1';
tsP_eig_vec2 = irf.ts_vec_xyz(PP.time,all_eig_vec2); tsP_eig_vec2.name = 'V_2';
tsP_eig_vec3 = irf.ts_vec_xyz(PP.time,all_eig_vec3); tsP_eig_vec3.name = 'V_3';

nt = PP.length;
all_eig2d_vals = zeros(nt,2);
all_eig2d_vec1 = zeros(nt,2);
all_eig2d_vec2 = zeros(nt,2);
all_eig2d_vec3 = zeros(nt,2);
nMovMean = 7;
data_P = movmean(PP.data,nMovMean,1);
for it = 1:nt
  [V,D] = eig(squeeze(data_P(it,1:2,1:2)));
  [Dsort, idsort] = sort(diag(D),'descend');
  all_eig2d_vals(it,:) = Dsort;
  all_eig2d_vec1(it,:) = V(:,idsort(1));
  all_eig2d_vec2(it,:) = V(:,idsort(2));
end

tsP_2d = irf.ts_tensor_xyz(PP.time,data_P(:,1:3,1:3));  tsP_2d.name = 'P (smoothed)';
tsP_eig2d_vals = irf.ts_scalar(PP.time,all_eig2d_vals); tsP_eig2d_vals.name = 'Eigenvalues';
tsP_eig2d_vec1 = irf.ts_vec_xy(PP.time,all_eig2d_vec1); tsP_eig2d_vec1.name = 'V_1';
tsP_eig2d_vec2 = irf.ts_vec_xy(PP.time,all_eig2d_vec2); tsP_eig2d_vec2.name = 'V_2';



[tsEig_val_3d, tsEig_v1_3d, tsEig_v2_3d, tsEig_v3_3d] = mvaPi3.eig([1 2 3]);
[tsEig_val_2d, tsEig_v1_2d, tsEig_v2_2d] = mvaPi3.eig([1 2]);

%% Calculation of boune times etc
f_wb = @(B,L,v) sqrt(units.e*B*v/units.mp/L);
f_fb = @(B,L,v) f_wb(B,L,v)/2/pi;
f_Tb = @(B,L,v) f_fb(B,L,v).^(-1);

B = 5e-9;
L = 600e3;
L = 200e3;
L = 1170e3;
L = 1170e3/4;
L = 400e3;
L = 800e3;
v = 1000e3;
wb = f_wb(B0,L,v);
fb = f_fb(B0,L,v);
Tb = f_Tb(B0,L,v);

disp(sprintf('B = %g nT, L = %g km, v = %g km/s: wb = %.2f rad/s, fb = %.2f Hz, Tb = %.2f s, Tb/2 = %.2f s',B0*1e9,L*1e-3,v*1e-3,wb,fb,Tb,Tb/2))


%% Noise-handling: Figure 1
h = irf_plot(7);

tsElow = iPDist3.find_noise_energy_limit(5).movmean(30);

PD1 = iPDist3;
PD2 = iPDist3.mask({[tsElow.data tsElow.data*0+Inf]});
PD3 = iPDist3.mask({[tsElow.data*0 tsElow.data]});
PD5 = iPDist3.mask({[tsElow.data*0 tsElow.data*0+500]});
PD4 = iPDist3.mask({[tsElow.data*0+500 tsElow.data*0+Inf]});

nMovMean = 21;
ne = ne3.movmean(nMovMean);
ni = ni3.movmean(nMovMean);
n1 = PD1.n.movmean(nMovMean);
n2 = PD2.n.movmean(nMovMean);
n3 = PD3.n.movmean(nMovMean);
n4 = PD4.n.movmean(nMovMean);
n5 = PD5.n.movmean(nMovMean);

colors = pic_colors('matlab');

if 1 % dEFlux ion PD1
  hca = irf_panel('ion dEF omni PD1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD1.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  irf_legend(hca,{'PD1'}',[0.01 0.1],'color','k','backgroundcolor',colors(1,:))
end
if 1 % dEFlux ion PD2
  hca = irf_panel('ion dEF omni PD2');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD2.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  irf_legend(hca,{'PD2'}',[0.01 0.1],'color','k','backgroundcolor',colors(2,:))
end
if 1 % dEFlux ion PD3
  hca = irf_panel('ion dEF omni PD3');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD3.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  irf_legend(hca,{'PD3'}',[0.01 0.1],'color','k','backgroundcolor',colors(3,:))
end
if 1 % dEFlux ion PD4
  hca = irf_panel('ion dEF omni PD4');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD4.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  irf_legend(hca,{'PD4'}',[0.01 0.1],'color','k','backgroundcolor',colors(4,:))
end
if 1 % dEFlux ion PD5
  hca = irf_panel('ion dEF omni PD5');
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_spectrogram(hca,PD5.deflux.omni.specrec,'donotfitcolorbarlabel');
  hca.YScale = 'log'; 
  irf_legend(hca,{'PD5'}',[0.01 0.1],'color','k','backgroundcolor',colors(5,:))
end

colors_density = circshift(colors,2);
if 1 % density for different energy ranges
  hca = irf_panel('n A');
  hca.ColorOrder = colors_density;
  irf_plot(hca,{ne,ni,n1,n2,n3,n4,n5},'comp')
  hca.YLabel.String = 'n (cm^{-3})';
  hca.ColorOrder = colors_density;
  irf_legend(hca,{'n_e (FPI)','n_i (FPI)','n_1','n_2','n_3','n_4','n_5'}',[1.01 0.98])
  irf_legend(hca,{sprintf('%g-point moving average',nMovMean)}',[0.01 0.05],'color','k')
  hl = findobj(hca,'type','line');
  c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
end
if 1 % density for different energy ranges
  hca = irf_panel('n A _');
  hca.ColorOrder = colors_density([1 2 5 7],:);
  irf_plot(hca,{ne,ni,n3,n5},'comp')
  hca.YLabel.String = 'n (cm^{-3})';
  hca.ColorOrder = colors_density([1 2 5 7],:);
  irf_legend(hca,{'n_e (FPI)','n_i (FPI)','n_3','n_5'}',[1.05 0.98])
  irf_legend(hca,{sprintf('%g-point moving average',nMovMean)}',[0.01 0.05],'color','k')
  hl = findobj(hca,'type','line');
  c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
end

irf_plot_axis_align
irf_zoom(h,'x',[ne3.time.start ne3.time.stop])

c_eval('h(?).FontSize = 16;',1:numel(h))
compact_panels(h,0.002)
h(end).XTickLabelRotation = 0;
c_eval('h(?).Layer = ''top'';',1:numel(h))
hlinks = linkprop(h(1:5),{'CLim','YLim'});


%% Noise-handling: Illustrate distribution of counts in the angular plane

PD = iPDist3;

PD2 = iPDist3.mask({[tsElow.data tsElow.data*0+Inf]});
PD3 = iPDist3.mask({[tsElow.data*0 tsElow.data]});

tint1 = time_xline + [-5 0];

tint2 = time_xline + [-30 30];

h = setup_subplots(2,2,'vertical');
isub = 1;

if 1 % low E  
  hca = h(isub); isub = isub + 1;
  %data = mean(PD3.data,[1 2]);
  mms.plot_skymap(hca,PD3.tlim(tint1),'flat');
  hca.Title.String = sprintf('high E: %s - %s',tint1(1).utc('HH:MM:SS'),tint1(2).utc('HH:MM:SS'));
end
if 1 % high E  
  hca = h(isub); isub = isub + 1;
  %data = mean(PD2.data,[1 2]);
  mms.plot_skymap(hca,PD2.tlim(tint1),'flat');
  hca.Title.String = sprintf('Low E: %s - %s',tint1(1).utc('HH:MM:SS'),tint1(2).utc('HH:MM:SS'));
end
if 1 % low E  
  hca = h(isub); isub = isub + 1;
  %data = mean(PD3.data,[1 2]);
  mms.plot_skymap(hca,PD3.tlim(tint2),'flat');
  hca.Title.String = sprintf('High E: %s - %s',tint1(1).utc('HH:MM:SS'),tint1(2).utc('HH:MM:SS'));
end
if 1 % high E  
  hca = h(isub); isub = isub + 1;
  %data = mean(PD2.data,[1 2]);
  mms.plot_skymap(hca,PD2.tlim(tint2),'flat');
  hca.Title.String = sprintf('Low E: %s - %s',tint1(1).utc('HH:MM:SS'),tint1(2).utc('HH:MM:SS'));
end

%hlinks = linkprop(h,{'CLim','YLim','XLim'});
c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).FontSize = 16;',1:numel(h))


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

if 1 % dEFlux ion
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
end
if 1 % dEFlux ion
  hca = irf_panel('ion dEF omni movmean rem onecounts');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.1],'fontsize',fontsize,'color','k');
end

if 0 % fi red L
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
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 0 % fi red M
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
if 0 % max shear direction
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
if 0 % e1
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
irf_zoom(h,'x',tint_figure_zoom_incl_sep)

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

%% Figure: Overview, paper

h = irf_plot(8);
c_eval('h(?).Position(2) = h(?).Position(2)-0.02;',1:numel(h))
fontsize = 14;
color_nan = [0.9 0.9 0.9]-0.1;
color_nan = [0.9 0.9 0.9]+0.1;

if 1 % B gse
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B','(nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
end 
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint),mvaVe?.y.tlim(tint),mvaVe?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint),mvaVi?.y.tlim(tint),mvaVi?.z.tlim(tint)},''comp'');',ic)  
  
  hca.YLabel.String = {'v_i','(km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end

if 0 % dEFlux ion
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  hca.YLabel.String = {'E_i','(eV)'};
end
if 0 % dEFlux ion
  hca = irf_panel('ion dEF omni movmean rem onecounts');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).deflux.omni.specrec,''donotfitcolorbarlabel'');',ic)  
  hca.YScale = 'log'; 
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.1],'fontsize',fontsize,'color','k');
  irf_legend(hca,{sprintf('one-counts removed')},[0.1,0.1],'fontsize',fontsize,'color','k');
  %irf_legend(hca,{sprintf('lower energy limit','for reduced VDFs')},[0.98,0.5],'fontsize',fontsize,'color','k');
  irf_legend(hca,{'lower energy limit','for reduced VDFs'}',[0.98,0.5],'fontsize',fontsize,'color','k');
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint,elim(1)*[1 1]),'k--')
  hold(hca,'off')
  hca.YLabel.String = {'E_i','(eV)'};
end

if 1 % dEFlux ion
  hca = irf_panel('ion dEF omni');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = PD_clean_notmasked.deflux.omni.specrec;',ic)
  specrec.p(specrec.p==0) = NaN;  
  irf_spectrogram(hca,specrec,'donotfitcolorbarlabel')
  hca.YScale = 'log'; 

  hold(hca,'on')
  irf_plot(hca,tsElow,'k','linewidth',2)
  irf_legend(hca,{'E_{low}'},[0.98 0.5],'color','k')
  hold(hca,'off')
  %hca.YLabel.String = 'E_{low} (eV)';
  hca.YLabel.String = {'E_i','(eV)'};
  hca.YLabel.Interpreter = 'tex';
  hca.Color = color_nan;  
end

if 0 % movmean counts
  hca = irf_panel('omni counts movmean 5');
  specrec = omni_counts2_movmean.specrec;
  [~,hcb] = irf_spectrogram(hca,specrec,'lin','donotfitcolorbarlabel');
  hca.YScale = 'log';
  hcb.YLabel.String = 'Total counts';
  
  irf_legend(hca,{sprintf('%g-point averaged counts',nMovMean),sprintf('%g-point averaged E_{low}',nSmoothElow)}',[0.1 0.1],'color','k','fontsize',fontsize,'fontweight','light','backgroundcolor','none')

  hca.CLim(1) = 0;
  hold(hca,'on')
  irf_plot(hca,tsElow,'k','linewidth',2)
  hold(hca,'off')
  hca.YLabel.String = 'E_{low} (eV)';
  hca.YLabel.Interpreter = 'tex';
end
if 1 % fi red L elow
  hca = irf_panel('fi L');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_L.specrec;',ic)
  specrec.p(specrec.p==0) = NaN;  
  %specrec.f = specrec.f - v_xline;
  c_eval('irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
    
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.x,''k-'')',ic)
  %c_eval('irf_plot(hca,irf.ts_scalar(mvaVi?.time,repmat(v_xline,[mvaVi?.length,1])),''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iL}','(km/s)'};
  hca.Color = color_nan;

  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 1 % fi red M elow
  hca = irf_panel('fi M');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_M.specrec; ',ic)
  specrec.p(specrec.p==0) = NaN;
  c_eval('irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)   
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.y,''k-'')',ic)
  hca.NextPlot = "replace";
  
  hca.YLabel.String = {'v_{iM}','(km/s)'};
  hca.Color = color_nan;
end
if 1 % fi red L elow
  hca = irf_panel('fi N');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('specrec = fi?_N.specrec; ',ic)
  specrec.p(specrec.p==0) = NaN;
  c_eval('irf_spectrogram(hca,specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iN}','(km/s)'};
  hca.Color = color_nan;
end

if 0 % fi red L
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

  hca.YLabel.String = {'v_{iL}','(km/s)'};

  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 0 % fi red M
  hca = irf_panel('fi M');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_M.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.y,''k-'')',ic)
  hca.NextPlot = "replace";
  
  hca.YLabel.String = {'v_{iM}','(km/s)'};
end
if 0 % fi red L
  hca = irf_panel('fi N');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_spectrogram(hca,fi?_N.specrec,''donotfitcolorbarlabel'');',ic)  
  
  set(hca,'ColorOrder',mms_colors('xyza'))
  %irf_legend(hca,{'x','y','z'},[0.98,0.3],'fontsize',fontsize);
  
  
  hca.NextPlot = "add";
  %c_eval('irf_plot(hca,mvaVi?.z,''k-'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iN}','(km/s)'};
end

if 1
  hca = irf_panel('fraction away from x line');
  set(hca,'ColorOrder',mms_colors('24'))
  ts1 = irf.ts_scalar(fi3_L_pos_neg.time,fi3_L_pos_neg.data(:,1)./sum(fi3_L_pos_neg.data(:,:),2));
  ts2 = irf.ts_scalar(fi3_L_pos_neg.time,fi3_L_pos_neg.data(:,2)./sum(fi3_L_pos_neg.data(:,:),2));
  data_away = zeros(ts1.length,1);
  idxt = find(ts1.time>time_xline);
  data_away = ts1.data;
  data_away(idxt) = ts2.data(idxt);
  ts_away = irf.ts_scalar(fi3_L_pos_neg.time,data_away);

  
  irf_plot(hca,{ts1,ts2},'comp')
  hold(hca,'on')
  hl = irf_patch(hca,{ts_away,0},'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','none');
  hl.EdgeColor = 'none';
  hl.FaceColor = [0.5 0.5 0.5];
  hl.FaceAlpha = 0.2;
  hold(hca,'off')
  
  hca.YLabel.String = {'Fraction','of ions'};
  
  set(hca,'ColorOrder',mms_colors('241'))  
  %irf_legend(hca,{'moving tailward','moving Earthward','moving away from X line'}',[0.2 0.7])
  irf_legend(hca,{'moving tailward','    moving Earthward','    moving away from X line'},[0.1 0.98],'fontsize',fontsize-2)
  hca.YTick = [0:0.25:1];
end
if 0 % max shear direction
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
if 0 % e1
  hca = irf_panel('Vi e1');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{e1.norm.x,e1.norm.y,e1.norm.z},''comp'');',ic)  
  
  hca.YLabel.String = {'e_1 (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end

irf_plot_axis_align
%irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint_figure_zoom)
irf_zoom(h,'x',tint_figure_zoom_incl_sep)

irf_zoom(h(1:3),'y')
irf_pl_mark(h,time_xline,'black','linestyle',':')
%irf_pl_mark(h,time_xline_ion,'red','linestyle',':')
%colormap(irf_colormap('thermal'))
%colormap([pic_colors('candy_gray'); 1 1 1])
irf_plot_axis_align
h(end).XTickLabelRotation = 0;
c_eval('h(?).YLabel.Interpreter = ''tex'';',1:numel(h))
%hlinks1 = linkprop(h(3:4),{'CLim'});
%hlinks2 = linkprop(h(5:7),{'CLim'});
%hlinks2 = linkprop(h(4:6),{'CLim'});
%h(1).Title.String = sprintf('N_{movmean} = %g',nMovMean);
c_eval('h(?).FontSize = fontsize;',1:numel(h))

hlinks = linkprop(h([5 6 7]),{'CLim'});
h(6).CLim = [-4   -1];
%colormap(irf_colormap('waterfall'))
%colormap(pic_colors('candy_gray'))
colormap(pic_colors('candy6'))
%hlinks = linkprop(h([4 5]),{'CLim'});
h(4).CLim = [ 3   6.0];
%hb = findobj(gcf,'type','colorbar');
%delete(hb(end))
%hb(4).Position(4) = hb(4).Position(4)*2.1;

% Add length on top
if 1
  %ax2 = axes('position',h(1).position);
  ax2 = axes('Position',h(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none','TickDir','out');
  %userdata = get(gcf,'userdata');
  %%
  di = mean(di3.tlim(time_xline+5*[-1 1]).data,1);
  dx = time_xline - tint_figure_zoom_incl_sep(1);

  xlim_s = h(1).XLim-dx; % time
  
  xlim_km = xlim_s*170;
  xlim_di = xlim_km/di;
  
  %xdata = (xlim_km(1):1000:xlim_km(end))/170;
  %xdata = (-10000:1000:10000)/170;
  
  xdata_di = -100:1:100;
  
  ax2.XTick = xdata_di;

  ax2.XLim = xlim_di;

  ax2.XLabel.String = 'd_i';
  ax2.XAxisLocation = 'top';
  ax2.YAxisLocation = 'right';
  ax2.YTick = [];
  ax2.XTickLabelRotation = 0;
  ax2.FontSize = fontsize;

%%
end


legends = {'a)','b)','c)','d)','e)','f)','g)','h)','i)','j)','k)','l)','m)'};
nInd = 1;
for ii = 1:numel(h)
  irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize+1)
  nInd = nInd + 1;
  %h(ii).FontSize = 16;
end


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
units = irf_units;
colors = pic_colors('matlab');

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
times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:07.940Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:10.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

if 0
times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:51.082Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:57.062Z';...
             '2017-07-11T22:34:02.000Z';...
             '2017-07-11T22:34:10.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];
end

nMovMean = 7;

if 1 % Calculate eigenvectors
  %%
  fT = 90;
  tint = time_xline_ion + 0.5*fT*[-1 1];
  
  c_eval('pdist_cleaned = iPDist?.movmean(nMovMean,''removeonecounts'',iPDist?_counts).tlim(tint);',ic)
  %pdist_cleaned = pdist_cleaned.tlim(tint);
  c_eval('pdist_original = iPDist?.movmean(nMovMean).tlim(tint);',ic)
  c_eval('pdist_diff = pdist_cleaned; pdist_diff.data = pdist_original.data - pdist_cleaned.data;',ic)
  pdist_diff.data(pdist_diff.data<0) = 0;
  
  p_c = mms.psd_moments(pdist_cleaned,scPot3,'energyrange',[1000 Inf]);
  p_o = mms.psd_moments(pdist_original,scPot3,'energyrange',[1000 Inf]);
  p_d = mms.psd_moments(pdist_diff,scPot3,'energyrange',[1000 Inf]);
  %

  dslP_c = p_c.P_psd;
  [tsEig_val_c_dsl, tsEig_v1_c_dsl, tsEig_v2_c_dsl] = dslP_c.eig([1 2]);
  
  mvaP_c = lmn*p_c.P_psd*lmn';
  mvaP_o = lmn*p_o.P_psd*lmn';
  mvaP_d = lmn*p_d.P_psd*lmn';
  
  
  [tsEig_val_c, tsEig_v1_c, tsEig_v2_c] = mvaP_c.eig([1 2]);
  [tsEig_val_o, tsEig_v1_o, tsEig_v2_o] = mvaP_o.eig([1 2]);
  [tsEig_val_d, tsEig_v1_d, tsEig_v2_d] = mvaP_d.eig([1 2]);
  %tsEig_v1_c = tsEig_v1_c;
  ss = sign(tsEig_v1_c.data(:,2)); % rotate into preferred direction (sign ambiguity)

  tsEig_v1_c.data = tsEig_v1_c.data.*[ss ss];
  tsEig_v2_c.data = tsEig_v2_c.data.*[ss ss];

end

nRows = 4;
nCols = size(times_utc,1);
[h1,h] = initialize_combined_plot('topbottom',3,nRows,nCols,0.3,'vertical');

vL_Xline = 0*-170;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
tint_zoom = irf.tint('2017-07-11T22:33:34.00Z/2017-07-11T22:34:30.00Z'); %20151112071854

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
%irf_legend(hca,{'L','M','N'},[0.98 0.98]);
%irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
irf_legend(hca,{sprintf(['v_L'],vL_Xline),'v_M','v_N'},[0.01 0.98]);
%irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';

if 0
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{mvaVExB?.resample(iPDist?)})',ic)
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';
end

if 1 % tsEig_v1_c
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  irf_plot(hca,{tsEig_v1_c.tlim(tint)})
  hold(hca,'on')

  hold(hca,'off')
  
  hca.YLabel.String = 'v_1 ';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'L','M'},[0.98 0.98]);
  hca.YLabel.Interpreter = 'tex';
end
if 1 % tsEig_val ratio
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  eig_ratio = irf.ts_scalar(tsEig_val_c.time,tsEig_val_c.data(:,1)./tsEig_val_c.data(:,2));
  irf_plot(hca,{eig_ratio.tlim(tint)})
  
  hca.YLabel.String = '\lambda_1/\lambda_2';
  %hca.ColorOrder = mms_colors('xyz');
  %irf_legend(hca,{'L','M'},[0.98 0.98]);
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1 3];
end

isub = 1;
% nSmooth = 1; % specified further doen


irf_zoom(h1,'x',tint_figure_zoom_incl_sep+[+10 -10])

elim = [000 Inf];
%elim = [200 Inf];
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)
%c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim);',ic)
pdist_all = pdist_cleaned;
time = time_vdf;

fontsize_leg = 9;
fontsize = 10;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];

% time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
% time = time + dt;
times = EpochTT(times_utc);
for it = 1:times.length%(1)
  time = times(it);
  time = time+-12;
  
  
  pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
  tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
  elow = max(tsElow.tlim(tint_dist).data);  
  %elow = 2000;
  pdist = pdist.elim([elow Inf]);
  %pdist.data(:,:,:,1) = 0;
  %pdist.data(:,:,:,end) = 0;

  % Calculate eigenvectors and values
  pmoms_tmp = mms.psd_moments(pdist,scPot3);
  mvaP = irf.ts_tensor_xyz(pmoms_tmp.P_psd.time,lmn*squeeze(pmoms_tmp.P_psd.data)*lmn');
  [tsEig_val_2d_tmp, tsEig_v1_2d_tmp, tsEig_v2_2d_tmp] = mvaP.eig([1 2]);


  ss = tsEig_v1_2d_tmp.data(2);
  
  %tsEig_v1_2d_tmp.data = tsEig_v1_2d_tmp*ss;
  %tsEig_v2_2d_tmp.data = tsEig_v2_2d_tmp*ss;

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
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  %scaxis
  

  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    if 0 % plot B direction
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
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    if 1 % eigenvectors, based on actual pdist chosen in the loop
      hold(hca,'on')     

      %EV1 = tsP_eig2d_vec1.resample(pdist.time).data;
      %EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      EV1 = tsEig_v1_2d_tmp.data;
      EV2 = tsEig_v2_2d_tmp.data;
      
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,+2*EV1(1)*scaxis_scale,+2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',1.5)
      quiver(hca,+EV1(1)*scaxis_scale,+EV1(2)*scaxis_scale,-2*EV1(1)*scaxis_scale,-2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',1.5)

      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,+2*EV2(1)*scaxis_scale,+2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',1.5)
      quiver(hca,+EV2(1)*scaxis_scale,+EV2(2)*scaxis_scale,-2*EV2(1)*scaxis_scale,-2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',1.5)
      hold(hca,'off') 
      irf_legend(hca,sprintf('l_1/l_2=%.2f',tsEig_val_2d_tmp.data(1,1)./tsEig_val_2d_tmp.data(1,2)),[0.98 0.98],'fontsize',10,'color','k')
    end
    if 0 % eigenvectors, calculated as a timeseries before .... they are the same now
      hold(hca,'on')     

      %EV1 = tsP_eig2d_vec1.resample(pdist.time).data;
      %EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      EV1 = tsEig_v1_c.resample(pdist.time).data;
      EV2 = tsEig_v2_c.resample(pdist.time).data;
      
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,+2*EV1(1)*scaxis_scale,+2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,+EV1(1)*scaxis_scale,+EV1(2)*scaxis_scale,-2*EV1(1)*scaxis_scale,-2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)

      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,+2*EV2(1)*scaxis_scale,+2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      quiver(hca,+EV2(1)*scaxis_scale,+EV2(2)*scaxis_scale,-2*EV2(1)*scaxis_scale,-2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
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
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
      
  
  if 0 % 1D at a certain vn range    
  %%
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [500 inf];
    vint2 = [-inf -500];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
  end
  if 1 % 1D along maximum eigenvector
  %%
    hca = h(isub); isub = isub + 1;
    EV1 = tsEig_v1_c_dsl.resample(pdist.time).data;
      %EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      
    v1 = [EV1 0];
    v3 = [EV2 0];
    v2 = cross(v3,v1);
    vdf = pdist.reduce('2D',v1,v2);
    %vdf = vdf_MN;

    vint1 = [00 inf];
    vint2 = [-inf -00];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
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



for ip = 3:nRows:numel(h)
  patch(h(ip),[-3000 3000 3000 -3000],[-3000 -3000 -500 -500],colors(2,:),...
    'facealpha',0.1,'edgecolor',colors(2,:));
  %h(ip).Children = circshift(h(ip).Children,1);

  patch(h(ip),[-3000 3000 3000 -3000],[3000 3000 500 500],colors(1,:),...
    'facealpha',0.1,'edgecolor',colors(1,:));
  h(ip).Children = circshift(h(ip).Children,1);
end

%colormap(pic_colors('candy_gray'))
%colormap(pic_colors('thermal'))
%colormap(pic_colors('candy6'))
colormap(pic_colors('candy4'))

hlinks_LM = linkprop(h(1:nRows:end),{'CLim'});
hlinks_LN = linkprop(h(2:nRows:end),{'CLim'});
hlinks_MN = linkprop(h(3:nRows:end),{'CLim'});
hlinks_M = linkprop(h(4:nRows:end),{'YLim','YTick','XTick'});

%hlinks_LM = linkprop(h(1:3:end),{'View'});
%hlinks_LN = linkprop(h(2:3:end),{'View'});
%hlinks_MN = linkprop(h(3:3:end),{'View'});

h(4).YTickLabelMode = 'auto';

%compact_panels(h,0.04,0.01)
drawnow

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
%delete(hb(1:end-3))
delete(hb)
hcb = colorbar(h(numel(h)-nRows+1),'location','northoutside');
hcb.Label.String = 'PSD (s^2m^{-4})';
hcb.FontSize = 11;

%
compact_panels(h(1:nRows:end),0,0.0)
compact_panels(h(2:nRows:end),0,0)
compact_panels(h(3:nRows:end),0,0)
compact_panels(h(4:nRows:end),0,0)
%
c_eval('h(?).Position(2) = h(?).Position(2)+0.08;',1:nRows:numel(h))
c_eval('h(?).Position(2) = h(?).Position(2)+0.04;',2:nRows:numel(h))
c_eval('h(?).Position(2) = h(?).Position(2)-0.04;',4:nRows:numel(h))

i2d = sort([1:nRows:numel(h), 2:nRows:numel(h), 3:nRows:numel(h)]);
hlinks = linkprop(h(i2d),{'XLim','YLim','CLim'});


c_eval('h(?).YTickLabel = [];',4:numel(h))
c_eval('h(?).YLabel = [];',4:numel(h))

%c_eval('h(?).YTickLabel = [];',4:numel(h))
%c_eval('h(?).YLabel = [];',4:numel(h))
c_eval('h(?).XTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).YTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))

c_eval('h(?).FontSize = 12;',1:numel(h))

c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).GridLineWidth = 1;',1:numel(h))

c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

%c_eval('h1(?).Position(3) = 0.3;',1:numel(h1))
c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

dy = 0.0848;
dy = 0.05;
h1(1).Position = [0.300    0.8730    0.4000    dy];
%h1(2).Position = [0.1700    0.7672    0.3000    0.0848];
h1(2).Position = [0.300    0.8730-dy    0.4000    dy];
h1(3).Position = [0.300    0.8730-dy*2    0.4000    dy];

%hca.CLim = [-9.5 -7.5];
h(1).CLim = [-10 -7.5];

h(4).YTickMode = 'auto';
h(4).YLabel.String = 'f_i(v_M) (s/m^4)';
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).FontSize = 8;',1:numel(h))


%% Figure: Reduced distributions, 2D, using mask and new clean
units = irf_units;
colors = pic_colors('matlab');

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
times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:07.940Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:10.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

if 0
times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:51.082Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:57.062Z';...
             '2017-07-11T22:34:02.000Z';...
             '2017-07-11T22:34:10.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];
end

PD_orig = iPDist3;
nMean = [7,3,3,3]; nThresh = 4;
PD_clean = PD_orig.remove_noise(nMean,nThresh,iPDist3_counts);
%PD_clean = PD_orig.remove_noise([4 3 3 3],4,iPDist3_counts);
%PD_clean = iPDist3.remove_noise([7 1 1 1],2,iPDist3_counts);
%PD_clean = PD_orig.remove_noise([3 3 3 3],2,iPDist3_counts);
PD_diff = PD_orig+-PD_clean;

tsElow = PD_orig.find_noise_energy_limit(5).movmean(30);

emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit

PD_orig = PD_orig.mask({emask_mat});
PD_clean = PD_clean.mask({emask_mat});
PD_diff = PD_diff.mask({emask_mat});

if 1 % Calculate eigenvectors
  %%
  fT = 90;
  tint = time_xline_ion + 0.5*fT*[-1 1];
    
  p_c = PD_clean.p;
  p_o = PD_orig.p;
  p_d = PD_diff.p;
  %
  
  [tsEig_val_c_dsl, tsEig_v1_c_dsl, tsEig_v2_c_dsl] = p_c.eig([1 2]);
  
  mvaP_c = lmn*p_c*lmn';
  mvaP_o = lmn*p_o*lmn';
  mvaP_d = lmn*p_d*lmn';
  
  
  [tsEig_val_c, tsEig_v1_c, tsEig_v2_c] = mvaP_c.eig([1 2]);
  [tsEig_val_o, tsEig_v1_o, tsEig_v2_o] = mvaP_o.eig([1 2]);
  [tsEig_val_d, tsEig_v1_d, tsEig_v2_d] = mvaP_d.eig([1 2]);
  %tsEig_v1_c = tsEig_v1_c;
  ss = sign(tsEig_v1_c.data(:,2)); % rotate into preferred direction (sign ambiguity)

  tsEig_v1_c.data = tsEig_v1_c.data.*[ss ss];
  tsEig_v2_c.data = tsEig_v2_c.data.*[ss ss];

end

%% Plot
nRows = 4;
nCols = size(times_utc,1);
[h1,h] = initialize_combined_plot('topbottom',1,nRows,nCols,0.3,'vertical');

vL_Xline = 1*-170;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
tint_zoom = irf.tint('2017-07-11T22:33:34.00Z/2017-07-11T22:34:30.00Z'); %20151112071854

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
%irf_legend(hca,{'L','M','N'},[0.98 0.98]);
%irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
irf_legend(hca,{sprintf(['v_L'],vL_Xline),'v_M','v_N'},[0.01 0.98]);
%irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';

if 0
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
c_eval('irf_plot(hca,{mvaVExB?.resample(iPDist?)})',ic)
irf_zoom(h1,'x',tint_zoom)
h1(2).YLim = [-2000 2000];
hca.YLabel.String = 'v_e (km/s)';
hca.ColorOrder = mms_colors('xyz');
irf_legend(hca,{'L','M','N'},[0.98 0.98]);
hca.YLabel.Interpreter = 'tex';
end

if 0 % tsEig_v1_c
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  irf_plot(hca,{tsEig_v1_c.tlim(tint)})
  hold(hca,'on')

  hold(hca,'off')
  
  hca.YLabel.String = 'v_1 ';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{'L','M'},[0.98 0.98]);
  hca.YLabel.Interpreter = 'tex';
end
if 0 % tsEig_val ratio
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  eig_ratio = irf.ts_scalar(tsEig_val_c.time,tsEig_val_c.data(:,1)./tsEig_val_c.data(:,2));
  irf_plot(hca,{eig_ratio.tlim(tint)})
  
  hca.YLabel.String = '\lambda_1/\lambda_2';
  %hca.ColorOrder = mms_colors('xyz');
  %irf_legend(hca,{'L','M'},[0.98 0.98]);
  hca.YLabel.Interpreter = 'tex';
  hca.YLim = [1 3];
end

isub = 1;
% nSmooth = 1; % specified further doen


irf_zoom(h1,'x',tint_figure_zoom_incl_sep+[+10 -10])

pdist_all = PD_clean;
pdist_all.data(:,:,:,[ ]) = 0;
time = time_vdf;

fontsize_leg = 9;
fontsize = 10;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];

% time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
% time = time + dt;
times = EpochTT(times_utc);
for it = 1:times.length%(1)
  time = times(it);
  time = time+-12*0+-2;
  
  
  pdist = pdist_all.tlim(time+nMean(1)*0.5*0.151*[-1 1]);
  tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*1*[-1 1];
  
  % Calculate eigenvectors and values
  p_tmp = pdist.p;
  mvaP = lmn*p_tmp*lmn';
  mvaP = irf.ts_tensor_xyz(pdist.time.start,mean(mvaP.data,1));
  [tsEig_val_2d_tmp, tsEig_v1_2d_tmp, tsEig_v2_2d_tmp] = mvaP.eig([1 2]);
  

  ss = tsEig_v1_2d_tmp.data(2);
  
  %tsEig_v1_2d_tmp.data = tsEig_v1_2d_tmp*ss;
  %tsEig_v2_2d_tmp.data = tsEig_v2_2d_tmp*ss;

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
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  %scaxis
  

  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    if 0 % plot B direction
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
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    if 0 % eigenvectors, based on actual pdist chosen in the loop
      hold(hca,'on')     

      %EV1 = tsP_eig2d_vec1.resample(pdist.time).data;
      %EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      EV1 = mean(tsEig_v1_2d_tmp.data,1);
      EV2 = mean(tsEig_v2_2d_tmp.data,1);
      
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,+2*EV1(1)*scaxis_scale,+2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',1.5)
      quiver(hca,+EV1(1)*scaxis_scale,+EV1(2)*scaxis_scale,-2*EV1(1)*scaxis_scale,-2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',1.5)

      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,+2*EV2(1)*scaxis_scale,+2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',1.5)
      quiver(hca,+EV2(1)*scaxis_scale,+EV2(2)*scaxis_scale,-2*EV2(1)*scaxis_scale,-2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',1.5)
      hold(hca,'off') 
      irf_legend(hca,sprintf('l_1/l_2=%.2f',tsEig_val_2d_tmp.data(1,1)./tsEig_val_2d_tmp.data(1,2)),[0.98 0.98],'fontsize',10,'color','k')
    end
    if 0 % eigenvectors, calculated as a timeseries before .... they are the same now
      hold(hca,'on')     

      %EV1 = tsP_eig2d_vec1.resample(pdist.time).data;
      %EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      EV1 = tsEig_v1_c.resample(pdist.time).data;
      EV2 = tsEig_v2_c.resample(pdist.time).data;
      
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,+2*EV1(1)*scaxis_scale,+2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,+EV1(1)*scaxis_scale,+EV1(2)*scaxis_scale,-2*EV1(1)*scaxis_scale,-2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)

      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,+2*EV2(1)*scaxis_scale,+2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      quiver(hca,+EV2(1)*scaxis_scale,+EV2(2)*scaxis_scale,-2*EV2(1)*scaxis_scale,-2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
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
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
      
  
  if 1 % 1D at a certain vn range    
  %%
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vmin = 000;
    vint1 = [500 inf];
    vint2 = [-inf -500];
    
    vint1 = [vmin inf];
    vint2 = [-inf -vmin];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = mean(vdf.data,1);
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = mean(vdf.data,1);
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
  end
  if 0 % 1D along maximum eigenvector
  %%
    hca = h(isub); isub = isub + 1;
    EV1 = mean(tsEig_v1_c_dsl.resample(pdist.time).data,1);
      %EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      
    v1 = [EV1 0];
    v3 = [EV2 0];
    v2 = cross(v3,v1);
    vdf = pdist.reduce('2D',v1,v2);
    %vdf = vdf_MN;

    vint1 = [00 inf];
    vint2 = [-inf -00];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = mean(vdf.data,1);
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = mean(vdf.data,1);
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
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



for ip = 3:nRows:numel(h)
  vmax = -2000;
  patch(h(ip),[-3000 vmax vmax -3000],[-3000 -3000 -vmin -vmin],colors(2,:),...
    'facealpha',0.1,'edgecolor',colors(2,:));
  %h(ip).Children = circshift(h(ip).Children,1);

  patch(h(ip),[-3000 vmax vmax -3000],[3000 3000 vmin vmin],colors(1,:),...
    'facealpha',0.1,'edgecolor',colors(1,:));
  h(ip).Children = circshift(h(ip).Children,1);
end

%colormap(pic_colors('candy_gray'))
%colormap(pic_colors('thermal'))
%colormap(pic_colors('candy6'))
%colormap(pic_colors('candy4'))
colormap(irf_colormap('magma'))

hlinks_LM = linkprop(h(1:nRows:end),{'CLim'});
hlinks_LN = linkprop(h(2:nRows:end),{'CLim'});
hlinks_MN = linkprop(h(3:nRows:end),{'CLim'});
hlinks_M = linkprop(h(4:nRows:end),{'YLim','YTick','XTick'});

h(4).YLim = [0 0.04];
%hlinks_LM = linkprop(h(1:3:end),{'View'});
%hlinks_LN = linkprop(h(2:3:end),{'View'});
%hlinks_MN = linkprop(h(3:3:end),{'View'});

h(4).YTickLabelMode = 'auto';

%compact_panels(h,0.04,0.01)
drawnow

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
%delete(hb(1:end-3))
delete(hb)
hcb = colorbar(h(numel(h)-nRows+1),'location','northoutside');
hcb.Label.String = 'PSD (s^2m^{-4})';
hcb.FontSize = 11;

%
compact_panels(h(1:nRows:end),0,0.0)
compact_panels(h(2:nRows:end),0,0)
compact_panels(h(3:nRows:end),0,0)
compact_panels(h(4:nRows:end),0,0)
%
c_eval('h(?).Position(2) = h(?).Position(2)+0.08;',1:nRows:numel(h))
c_eval('h(?).Position(2) = h(?).Position(2)+0.04;',2:nRows:numel(h))
c_eval('h(?).Position(2) = h(?).Position(2)-0.04;',4:nRows:numel(h))

i2d = sort([1:nRows:numel(h), 2:nRows:numel(h), 3:nRows:numel(h)]);
hlinks = linkprop(h(i2d),{'XLim','YLim','CLim'});


c_eval('h(?).YTickLabel = [];',nRows+1:numel(h))
c_eval('h(?).YLabel = [];',nRows+1:numel(h))

%c_eval('h(?).YTickLabel = [];',4:numel(h))
%c_eval('h(?).YLabel = [];',4:numel(h))
c_eval('h(?).XTick = -2000:1000:2000;',1:numel(h))
c_eval('h(?).YTick = -2000:1000:2000;',setdiff(1:numel(h),4:4:numel(h)))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))

c_eval('h(?).FontSize = 12;',1:numel(h))

c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).GridLineWidth = 1;',1:numel(h))

c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

%c_eval('h1(?).Position(3) = 0.3;',1:numel(h1))
c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

dy = 0.0848;
%dy = 0.05;
h1(1).Position = [0.300    0.7730    0.4000    dy];
%h1(1).Position = [0.300    0.8730    0.4000    dy];
%h1(2).Position = [0.1700    0.7672    0.3000    0.0848];
%h1(2).Position = [0.300    0.8730-dy    0.4000    dy];
%h1(3).Position = [0.300    0.8730-dy*2    0.4000    dy];

%hca.CLim = [-9.5 -7.5];
h(1).CLim = [-10 -7.5];

h(4).YTickMode = 'auto';
h(4).YLabel.String = 'f_i(v_M) (s/m^4)';
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).FontSize = 8;',1:numel(h))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
hall = findobj(gcf,'type','axes'); hall = hall(end:-1:1);
c_eval('hall(?).LineWidth = 1.5;',1:numel(hall))


%% Make movie of 2d VDFs, to get a better feeling of the turning.

fT = 30;
tint = time_xline_ion + 0.5*fT*[-1 1];

nMovMean = 7;
%c_eval('pdist_all = iPDist?.movmean(nMovMean,''removeonecounts'',iPDist?_counts).tlim(tint);',ic)
c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim).tlim(tint);',ic)
c_eval('pdist_all = pdist_diff.tlim(tint);',ic)
elows_tmp = tsElow; elows_tmp.data = movmean(elows_tmp.data,nMovMean);
elows_all = elows_tmp.tlim(tint);

nRows = 2;
nCols = 2;
[h1,h] = initialize_combined_plot('topbottom',1,nRows,nCols,0.3,'vertical');

vL_Xline = 0*-170;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854

if 1 % vi
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'},[.98 0.05]);
  hca.YLabel.Interpreter = 'tex';
end

isub = 1;
% nSmooth = 1; % specified further doen


irf_zoom(h1,'x',tint_figure_zoom_incl_sep)

elim = [000 Inf];
time = time_vdf;

fontsize_leg = 9;
fontsize = 10;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];

% time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
% time = time + dt;
directory_ = strrep(printpath,'\','');
vidfile = VideoWriter([directory_ 'testmovie.mp4'],'MPEG-4');
open(vidfile);
     
clear F
times = pdist_all.time;
for it = 1:3:times.length %1:nMovMean:times.length
  time = times(it);
  
  pdist = pdist_all.tlim(time+0.5*0.151*[-1 1]);
  tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
  elow = max(tsElow.tlim(tint_dist).data);  
  pdist = pdist.elim([elow Inf]);
  
  if exist('hmark','var'); delete(hmark); end
  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))
  
  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  scaxis
  

  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;
  isub = 1;

  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    if 0 % plot B direction
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
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end

    if 1 % plot eigenvectors
      hold(hca,'on')     
      eigval = tsP_eig2d_vals.resample(pdist.time).data;
      EV1 = tsP_eig2d_vec1.resample(pdist.time).data;
      EV2 = tsP_eig2d_vec2.resample(pdist.time).data;
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,2*EV1(1)*scaxis_scale,2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,2*EV2(1)*scaxis_scale,2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      hold(hca,'off')
      irf_legend(hca,sprintf('eig1/eig2=%.2f',eigval(1)/eigval(2)),[0.98 0.98],'color','k')
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
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
    
  if 1 % 1D at a certain vn range    
  %%
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [500 inf];
    vint2 = [-inf -500];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('%g km/s < v_N',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
    hca.YLim = [0 0.025]
    hca.XGrid = 'on';
    hca.YGrid = 'on';
  end
  
  irf_legend(h(1),{sprintf('N = %g',nMovMean)},[0.02 1.01])
  hlinks = linkprop(h(1:3),{'CLim'});
  hlinks.Targets(1).CLim = [-10 -7.5];
  c_eval('h(?).FontSize = 12;',1:numel(h))
  colormap(pic_colors('candy6'))

  drawnow
  pause(0.1)
  F(it) = getframe(gcf); 
  writeVideo(vidfile,F(it));  
end
close(vidfile);




%% Make movie of 2d VDFs, with original, cleaned, and diff pdists

fT = 30;
tint = time_xline_ion + 0.5*fT*[-1 1];

nMovMean = 4;
%c_eval('pdist_all = iPDist?.tlim(tint);',ic)
%c_eval('pdist_cleaned = iPDist?.movmean(nMovMean,''removeonecounts'',iPDist?_counts).tlim(tint);',ic)
%pdist_cleaned = PD.tlim(tint);
%c_eval('pdist_movmean = iPDist?.movmean(nMovMean).tlim(tint);',ic)
%c_eval('pdist_diff = pdist_cleaned; pdist_diff.data = pdist_movmean.data - pdist_cleaned.data;',ic)
%pdist_diff.data(pdist_diff.data<0) = 0;
c_eval('pdist_all = iPDist?.tlim(tint);',ic)
nMean = [4 3 3 3];
nThresh = 4;
pdist_cleaned = iPDist3.remove_noise(nMean,nThresh,iPDist3_counts).tlim(tint);

pdist_diff = pdist_all+-pdist_cleaned;

tsElow = pdist_all.find_noise_energy_limit(5).movmean(30);

emask_mat = [tsElow.data*0 tsElow.data]; % setting all datapoints within these energy bounds to nan, effectively applying a lower energy limit

pdist_cleaned = pdist_cleaned.mask({emask_mat});
pdist_all = pdist_all.mask({emask_mat});
pdist_diff = pdist_diff.mask({emask_mat});

%%
if 1
  nMovMean = nMean(1);

%  pdist_original = PD_orig;
%  pdist_cleaned = PD_new;
%  pdist_diff = PD_diff;

  pdist_original = pdist_all;
  pdist_cleaned = pdist_cleaned;
  pdist_diff = pdist_diff;

  %pdist_original = PD2_orig.movMean(nMovMean);
  %pdist_cleaned = PD2_new;
  %pdist_diff = PD2_diff;
  
end

%%
p_c = mms.psd_moments(pdist_cleaned,scPot3,'energyrange',[300 Inf]);
p_o = mms.psd_moments(pdist_original,scPot3,'energyrange',[300 Inf]);
p_d = mms.psd_moments(pdist_diff,scPot3,'energyrange',[300 Inf]);
%%
mvaP_c = lmn*p_c.P_psd*lmn';
mvaP_o = lmn*p_o.P_psd*lmn';
mvaP_d = lmn*p_d.P_psd*lmn';


[tsEig_val_c, tsEig_v1_c, tsEig_v2_c] = mvaP_c.eig([1 2]);
[tsEig_val_o, tsEig_v1_o, tsEig_v2_o] = mvaP_o.eig([1 2]);
[tsEig_val_d, tsEig_v1_d, tsEig_v2_d] = mvaP_d.eig([1 2]);

if 0
tsEig_val_c = tsEig_val_c.movmean(nMean(1));
tsEig_val_o = tsEig_val_o.movmean(nMean(1));
tsEig_val_d = tsEig_val_d.movmean(nMean(1));
tsEig_v1_c = tsEig_v1_c.movmean(nMean(1));
tsEig_v1_o = tsEig_v1_o.movmean(nMean(1));
tsEig_v1_d = tsEig_v1_d.movmean(nMean(1));
tsEig_v2_c = tsEig_v2_c.movmean(nMean(1));
tsEig_v2_o = tsEig_v2_o.movmean(nMean(1));
tsEig_v2_d = tsEig_v2_d.movmean(nMean(1));
end
%%
units = irf_units;
elows_tmp = tsElow; elows_tmp.data = movmean(elows_tmp.data,nMovMean);
elows_all = elows_tmp.tlim(tint);

nRows = 3;
nCols = 3;
[h1,h] = initialize_combined_plot('topbottom',1,nRows,nCols,0.15,'horizontal');
h1.Position(2) = h1.Position(2) - 0.05;
h1.Position(4) = h1.Position(4)*6;

vL_Xline = 0*-170;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854

if 1 % vi
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'},[.98 0.05]);
  hca.YLabel.Interpreter = 'tex';
end

isub = 1;
% nSmooth = 1; % specified further down


irf_zoom(h1,'x',tint_figure_zoom_incl_sep)

elim = [000 Inf];
time = time_vdf;

fontsize_leg = 9;
fontsize = 10;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];

% time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
% time = time + dt;
directory_ = strrep(printpath,'\','');
vidfile = VideoWriter([directory_ 'testmovie.mp4'],'MPEG-4');
open(vidfile);
     
clear F
nDistToMean = nMean(1);
times = pdist_all.time;
for it = 1:7:times.length %1:nMovMean:times.length
  time = times(it);
  
  pdist_c = pdist_cleaned.tlim(time+nDistToMean*0.5*0.151*[-1 1]);
  pdist_o = pdist_original.tlim(time+nDistToMean*0.5*0.151*[-1 1]);
  pdist_d = pdist_diff.tlim(time+nDistToMean*0.5*0.151*[-1 1]);
  pdist = pdist_o;
  tint_dist = [pdist.time.start pdist.time.stop] + nDistToMean*0.5*0.150*[-1 1];
  elow = max(tsElow.tlim(tint_dist).data);  
  
  pdist_c = pdist_c.elim([elow Inf]);
  pdist_o = pdist_o.elim([elow Inf]);
  pdist_d = pdist_d.elim([elow Inf]);
  
  if exist('hmark','var'); delete(hmark); end
  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))
  
  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  scaxis
  

  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;
  isub = 1;
  
  pdist = pdist_o;
  tsEig = tsEig_val_o;
  tsEigV1 = tsEig_v1_o;
  tsEigV2 = tsEig_v2_o;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end

    if 1 % plot eigenvectors
      hold(hca,'on')     
      eigval = mean(tsEig.resample(pdist.time).data,1);
      EV1 = mean(tsEigV1.resample(pdist.time).data,1);
      EV2 = mean(tsEigV2.resample(pdist.time).data,1);
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,2*EV1(1)*scaxis_scale,2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,2*EV2(1)*scaxis_scale,2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      hold(hca,'off')
      irf_legend(hca,sprintf('eig1/eig2=%.2f',eigval(1)/eigval(2)),[0.98 0.98],'color','k','fontsize',fontsize_leg)
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
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  pdist = pdist_c;
  tsEig = tsEig_val_c;
  tsEigV1 = tsEig_v1_c;
  tsEigV2 = tsEig_v2_c;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end

    if 1 % plot eigenvectors
      hold(hca,'on')     
      eigval = mean(tsEig.resample(pdist.time).data,1);
      EV1 = mean(tsEigV1.resample(pdist.time).data,1);
      EV2 = mean(tsEigV2.resample(pdist.time).data,1);
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,2*EV1(1)*scaxis_scale,2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,2*EV2(1)*scaxis_scale,2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      hold(hca,'off')
      irf_legend(hca,sprintf('eig1/eig2=%.2f',eigval(1)/eigval(2)),[0.98 0.98],'color','k','fontsize',fontsize_leg)
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
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  
  pdist = pdist_d;
  tsEig = tsEig_val_d;
  tsEigV1 = tsEig_v1_d;
  tsEigV2 = tsEig_v2_d;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    if 1 % plot eigenvectors
      hold(hca,'on')     
      eigval = mean(tsEig.resample(pdist.time).data,1);
      EV1 = mean(tsEigV1.resample(pdist.time).data,1);
      EV2 = mean(tsEigV2.resample(pdist.time).data,1);
      quiver(hca,-EV1(1)*scaxis_scale,-EV1(2)*scaxis_scale,2*EV1(1)*scaxis_scale,2*EV1(2)*scaxis_scale,0,'color',colors(1,:),'linewidth',2)
      quiver(hca,-EV2(1)*scaxis_scale,-EV2(2)*scaxis_scale,2*EV2(1)*scaxis_scale,2*EV2(2)*scaxis_scale,0,'color',colors(2,:),'linewidth',2)
      hold(hca,'off')
      irf_legend(hca,sprintf('eig1/eig2=%.2f',eigval(1)/eigval(2)),[0.98 0.98],'color','k','fontsize',fontsize_leg)
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
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    %hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  
  irf_legend(h(1),{sprintf('Window = [%g,%g,%g,%g]',nMean(1),nMean(2),nMean(3),nMean(4)),sprintf('N < %g',nThresh)},[0.02 1.05],'color','k')
  %irf_legend(h(1),{sprintf('N = %g',nMovMean)},[0.02 1.01])
  hlinks = linkprop(h,{'CLim'});
  hlinks.Targets(1).CLim = [-10 -7.5];
  c_eval('h(?).FontSize = 10;',1:numel(h))
  colormap(pic_colors('candy6'))
  hb = findobj(gcf,'type','colorbar');
  delete(hb(2:end))

  drawnow
  pause(0.1)
  F(it) = getframe(gcf); 
  writeVideo(vidfile,F(it));  
end
close(vidfile);




%
%% Distributions to go along wth sketch of trajectory in NM plane

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:07.940Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

time_pdist = EpochTT('2017-07-11T22:34:03.000Z');

elim = [000 Inf];
%elim = [200 Inf];
nMovMean = 5;
c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim);',ic)

fontsize_leg = 9;
fontsize = 10;

time = time_pdist;


pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
elow = max(tsElow.tlim(tint_dist).data);  

pdist = pdist.elim([elow Inf]);

t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
Ldsl = Tdsl(1,:);
Mdsl = Tdsl(2,:);
Ndsl = Tdsl(3,:);

c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  
 
colors = pic_colors('matlab');
nSmooth = 0;
nContours = 0;

h = setup_subplots(1,2);
isub = 1;

  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl]);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    hca.YLabel.String = 'v_N (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    if 1
      hold(hca,'on')
      patch(hca,[-3000 3000 3000 -3000],[-3000 -3000 -500 -500],colors(2,:),...
        'facealpha',0.1,'edgecolor',colors(2,:));
      %h(ip).Children = circshift(h(ip).Children,1);
    
      patch(hca,[-3000 3000 3000 -3000],[3000 3000 500 500],colors(1,:),...
        'facealpha',0.1,'edgecolor',colors(1,:));
      hca.Children = circshift(hca.Children,1);
      hold(hca,'off')
    end
    hca.CLim = [-10 -7.5];
    colormap(pic_colors('candy6'))
  end
    
  if 1 % 1D at a certain vn range    
  %
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [500 inf];
    vint2 = [-inf -500];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    
    vdf1 = PDist(time,data1,'1Dcart',vdf.depend{1});
    vdf2 = PDist(time,data2,'1Dcart',vdf.depend{1});


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2,'linewidth',2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1499;
    hca.XLim = vlim*[-1 1];  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = [-2000:500:2000];
    hca.XTickLabelRotation = 0;
    if 0 % dashed line
      %%
      hold(hca,'on')
      fmax = @(v,n,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
      
      ff = fmax(v_center*1e3,0.005*1e6,200e3,500e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',colors(3,:))
      ff = fmax(v_center*1e3,0.009*1e6,-500e3,300e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',[0 0.7 0.7])
      
      ff = fmax(v_center*1e3,0.005*1e6,100e3,500e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',colors(5,:))
      ff = fmax(v_center*1e3,0.009*1e6,-400e3,500e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',[0.5 0.0 0.9])
      
      hold(hca,'off')
    end
    if 0 % patch
      %%
      hold(hca,'on')
      fmax = @(v,n,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
      
      ff1 = fmax(v_center*1e3,0.004*1e6,200e3,400e3);
      aa = [v_center, v_center(end:-1:1)];
      bb = [ff1 v_center*0];
      patch(hca,aa,bb,colors(3,:),...
        'facealpha',0.2,'edgecolor','none')

      ff2 = fmax(v_center*1e3,0.006*1e6,-500e3,200e3);
      aa = [v_center, v_center(end:-1:1)];
      bb = [ff2 v_center*0];
      patch(hca,aa,bb,[0 0.7 0.7],...
        'facealpha',0.3,'edgecolor','none')

      if 0
      bb = [ff1+ff2 v_center*0];
      patch(hca,aa,bb,[0 0.7 0.7],...
        'facealpha',0.3,'edgecolor',[0 0 0])      
      end

      %ff = fmax(v_center*1e3,0.005*1e6,100e3,500e3);
      %plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',colors(5,:))
      %ff = fmax(v_center*1e3,0.009*1e6,-400e3,500e3);
      %plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',[0.5 0.0 0.9])
      
      hold(hca,'off')
    end
  end
  
 %
for ip = 1:numel(h)
  hca = h(ip);
  %axis(hca,'square')
  hca.Position(3) = 0.25;
end
h(1).Title.String = sprintf('%s - %s',tint_dist(1).utc('HH:MM:SS.mmm'),tint_dist(2).utc('SS.mmm'));
%compact_panels(h,0.01,0.2)

%% Distributions to illustrste gyroturning, divide f(v_L,v_M<>0)

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:07.940Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

time_pdist = EpochTT('2017-07-11T22:34:03.000Z')+0;

elim = [000 Inf];
%elim = [200 Inf];
nMovMean = 5;
c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim);',ic)

fontsize_leg = 9;
fontsize = 10;



h = setup_subplots(2,3,'vertical');
isub = 1;

for dt = [-10 0 7]
time = time_pdist+dt;


pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
elow = max(tsElow.tlim(tint_dist).data);  

pdist = pdist.elim([elow Inf]);

t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
Ldsl = Tdsl(1,:);
Mdsl = Tdsl(2,:);
Ndsl = Tdsl(3,:);

c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  
 
colors = pic_colors('matlab');
nSmooth = 0;
nContours = 0;


  vint1 = [0 inf];
  vint2 = [-inf -0];
  vlim_1d = 1499;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Ldsl],[Mdsl]);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours,'log10',0)
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.YLabel.String = 'v_M (km/s)';
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)*1e0,mean(vExB.z.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    vlim = 2500;
    hca.XTick = [-2000:1000:2000];
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    if 1
      hold(hca,'on')
      patch(hca,vlim_1d*[-1 1 1 -1],[3000 3000 vint1(1) vint1(1)],colors(1,:),...
        'facealpha',0.1,'edgecolor',colors(1,:));
      %h(ip).Children = circshift(h(ip).Children,1);
    
      patch(hca,vlim_1d*[-1 1 1 -1],[-3000 -3000 vint2(2) vint2(2)],colors(2,:),...
        'facealpha',0.1,'edgecolor',colors(2,:));
      hca.Children = circshift(hca.Children,1);
      hold(hca,'off')
    end
    %hca.CLim = [-10 -7.5];
    colormap(hca,pic_colors('candy6'))
    hca.Title.String = sprintf('%s - %s',tint_dist(1).utc('HH:MM:SS.mmm'),tint_dist(2).utc('SS.mmm'));
  end
    
  if 1 % 1D at a certain vn range    
  %
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2));
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2));
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;
    
    vdf1 = PDist(time,data1,'1Dcart',vdf.depend{1});
    vdf2 = PDist(time,data2,'1Dcart',vdf.depend{1});


    hca.ColorOrder = colors;
    plot(hca,v_center,data1,v_center,data2,'linewidth',2)    
    hca.XLabel.String = 'v_L (km/s)';
    
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_L) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,...
      {sprintf('v_M > %g km/s',vint1(1)),sprintf('v_M < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    
    hca.XLim = vlim_1d*[-1 1];  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = [-2000:500:2000];
    hca.XTickLabelRotation = 0;
    if 0 % dashed line
      %%
      hold(hca,'on')
      fmax = @(v,n,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
      
      ff = fmax(v_center*1e3,0.005*1e6,200e3,500e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',colors(3,:))
      ff = fmax(v_center*1e3,0.009*1e6,-500e3,300e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',[0 0.7 0.7])
      
      ff = fmax(v_center*1e3,0.005*1e6,100e3,500e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',colors(5,:))
      ff = fmax(v_center*1e3,0.009*1e6,-400e3,500e3);
      plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',[0.5 0.0 0.9])
      
      hold(hca,'off')
    end
    if 0 % patch
      %%
      hold(hca,'on')
      fmax = @(v,n,vd,vt) n.*(1/pi./vt.^2)^(1/2)*exp(-(v-vd).^2./vt.^2);
      
      ff1 = fmax(v_center*1e3,0.004*1e6,200e3,400e3);
      aa = [v_center, v_center(end:-1:1)];
      bb = [ff1 v_center*0];
      patch(hca,aa,bb,colors(3,:),...
        'facealpha',0.2,'edgecolor','none')

      ff2 = fmax(v_center*1e3,0.006*1e6,-500e3,200e3);
      aa = [v_center, v_center(end:-1:1)];
      bb = [ff2 v_center*0];
      patch(hca,aa,bb,[0 0.7 0.7],...
        'facealpha',0.3,'edgecolor','none')

      if 0
      bb = [ff1+ff2 v_center*0];
      patch(hca,aa,bb,[0 0.7 0.7],...
        'facealpha',0.3,'edgecolor',[0 0 0])      
      end

      %ff = fmax(v_center*1e3,0.005*1e6,100e3,500e3);
      %plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',colors(5,:))
      %ff = fmax(v_center*1e3,0.009*1e6,-400e3,500e3);
      %plot(hca,v_center,ff,'linewidth',2,'linestyle','--','color',[0.5 0.0 0.9])
      
      hold(hca,'off')
    end
  end
end
c_eval('h(?).YLabel = []; h(?).YTickLabel = [];',3:numel(h)) 
hb = findobj(gcf,'type','colorbar');
delete(hb(2:3))
%%
compact_panels(h(1:2:end),0.01,0.01)
compact_panels(h(2:2:end),0.01,0.01)
compact_panels(h(1:2:end),0.01,0.01)
c_eval('h(?).Position(1) = h(?-1).Position(1);',2:2:numel(h))
 %%
for ip = 1:numel(h)
  hca = h(ip);
  %axis(hca,'square')
  %hca.Position(3) = 0.25;
  hca.XLim = 2500*[-1 1];
end

%compact_panels(h,0.01,0.2)


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
nMovMean = 20;
elim = [500 Inf];
c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).elim(elim);',ic)

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
dt = 0;
time = time + dt;
pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
c_eval('pdist_with_noise = iPDist?.movmean(nMovMean).tlim(time+0.5*0.15*[-1 1]);',ic)
c_eval('pdist_with_noise_elim = iPDist?.elim(elim).movmean(nMovMean).tlim(time+0.5*0.15*[-1 1]);',ic)

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
vdf_with_noise_elim = pdist_with_noise_elim.reduce('1D',Mdsl);
plot(hca,vdf.depend{1},vdf.data*1.,...
  vdf_with_noise.depend{1},vdf_with_noise.data,...
  vdf_with_noise_elim.depend{1},vdf_with_noise_elim.data)
hca.XLim = vlim;
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'f_i (s/m^4)';
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))
hca.XGrid = 'on';
hca.YGrid = 'on';
t1t2 = pdist.time + 0.5*nMovMean*0.150*[-1 1];
irf_legend(hca,{sprintf('one-counts removed, nMovMean=%g',nMovMean),'original',sprintf('E > %g eV',elim(1))}',[0.98,0.98],'fontsize',fontsize);

irf_legend(hca,sprintf('%s-%s',t1t2(1).utc('HH:MM:SS.mmm'),t1t2(2).utc('HH:MM:SS.mmm')),[0.02,0.98],'fontsize',fontsize,'color','k');

%% Single (or a few) 1D reduced ion population

fontsize = 14;

nMovMean = [2 5 10 20];
elow = 100:100:2000;
c_eval('pdist_orig = iPDist?;', ic)
c_eval('pdist_mean! = iPDist?.movmean(!,''RemoveOneCounts'',iPDist?_counts);',ic, nMovMean)

%vint = [];
%%

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
dt = 0;
time = time + dt;

% Pick out the time
pdist_ = pdist_orig.tlim(time+0.5*0.15*[-1 1]);
c_eval('pdist_mean?_ = pdist_mean?.tlim(time+0.5*0.15*[-1 1]);',nMovMean)

t_dist_center = pdist_.time.start + (pdist_.time.stop - pdist_.time.start)/2;
c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
Ldsl = Tdsl(1,:);
Mdsl = Tdsl(2,:);
Ndsl = Tdsl(3,:);

vg_edges = -2000:20:2000;


h = setup_subplots(2,3,'horizontal');
isub = 1;

if 1 % Different elow, M
  hca = h(isub); isub = isub + 1;
  
  vdf_orig = pdist_.reduce('1D',Mdsl,'nMC',1000,'vg_edges',vg_edges);
  plot(hca,vdf_orig.depend{1},vdf_orig.data)
  hold(hca,'on')
  estep = 2;
  for iE = 1:estep:numel(elow)  
    vdf_elim = pdist_.elim([elow(iE) Inf]).reduce('1D',Mdsl);
    plot(hca,vdf_elim.depend{1},vdf_elim.data)  
  end
  legs = arrayfun(@(x)sprintf('E>%g eV',x),elow(1:estep:end),'UniformOutput',false);
  legend(hca,{'Original', legs{:}},'Box','off','fontsize',7,'location','northeast')
  %irf_legend(hca,{'Original', legs{:}}',[0.98 0.98])
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_M (km/s)';
  hca.YLabel.String = 'f(v_M) (s/m^4)';
  hca.YLim = [0 0.1];
  hca.XLim = [-2000 2000];
end
if 1 % Different elow, N
  hca = h(isub); isub = isub + 1;
  
  vdf_orig = pdist_.reduce('1D',Ndsl,'nMC',1000,'vg_edges',vg_edges);
  plot(hca,vdf_orig.depend{1},vdf_orig.data)
  hold(hca,'on')
  estep = 2;
  for iE = 1:estep:numel(elow)  
    vdf_elim = pdist_.elim([elow(iE) Inf]).reduce('1D',Ndsl);
    plot(hca,vdf_elim.depend{1},vdf_elim.data)  
  end
  legs = arrayfun(@(x)sprintf('E>%g eV',x),elow(1:estep:end),'UniformOutput',false);
  %legend(hca,{'Original', legs{:}},'Box','off')
  legend(hca,{'Original', legs{:}},'Box','off','fontsize',7,'location','northeast')
  %irf_legend(hca,{'Original', legs{:}}',[0.98 0.98])
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_N (km/s)';
  hca.YLabel.String = 'f(v_N) (s/m^4)';
  hca.YLim = [0 0.1];
  hca.XLim = [-2000 2000];
end
if 1 % Different elow, L
  hca = h(isub); isub = isub + 1;
  
  vdf_orig = pdist_.reduce('1D',Ldsl,'nMC',1000,'vg_edges',vg_edges);
  plot(hca,vdf_orig.depend{1},vdf_orig.data)
  hold(hca,'on')
  estep = 2;
  for iE = 1:estep:numel(elow)  
    vdf_elim = pdist_.elim([elow(iE) Inf]).reduce('1D',Ldsl);
    plot(hca,vdf_elim.depend{1},vdf_elim.data)  
  end
  legs = arrayfun(@(x)sprintf('E>%g eV',x),elow(1:estep:end),'UniformOutput',false);
  %legend(hca,{'Original', legs{:}},'Box','off')
  legend(hca,{'Original', legs{:}},'Box','off','fontsize',7,'location','northeast')
  %irf_legend(hca,{'Original', legs{:}}',[0.98 0.98])
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'f(v_L) (s/m^4)';
  hca.YLim = [0 0.1];
  hca.XLim = [-2000 2000];
end
if 0 % Different nMovMean
  hca = h(isub); isub = isub + 1;
  
  vdf_orig = pdist_.reduce('1D',Mdsl,'nMC',1000,'vg_edges',vg_edges);
  plot(hca,vdf_orig.depend{1},vdf_orig.data)
  hold(hca,'on')
  for iN = 1:numel(nMovMean)  
    c_eval('vdf_mean = pdist_mean?_.reduce(''1D'',Mdsl);',nMovMean(iN))
    plot(hca,vdf_mean.depend{1},vdf_mean.data)  
  end
  legs = arrayfun(@(x)sprintf('nMovMean=%g',x),nMovMean,'UniformOutput',false);
  legend(hca,{'Original', legs{:}},'Box','off')
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_M (km/s)';
  hca.YLabel.String = 'f(v_M) (s/m^4)';
  hca.YLim = [0 0.1];
  hca.XLim = [-2000 2000];
end

if 1 % Different elow, 2d map m
  hca = h(isub); isub = isub + 1;
    
  vdf_orig = pdist_.reduce('1D',Mdsl,'nMC',1000,'vg_edges',vg_edges);
  all_vdfs = zeros(numel(elow)+1,size(vdf_orig.depend{1},2));
  all_vdfs(1,:) = vdf_orig.data;
  hold(hca,'on')
  for iE = 1:numel(elow)  
    vdf_elim = pdist_.elim([elow(iE) Inf]).reduce('1D',Mdsl,'nMC',1000,'vg_edges',vg_edges);
    %plot(hca,vdf_elim.depend{1},vdf_elim.data)  
    all_vdfs(1+iE,:) = vdf_elim.data;
  end
  e_center = [0 elow]; de = diff(e_center);
  e_edges = -50:100:(elow(end)+50);
  v_edges = vdf_orig.ancillary.v_edges*1e-3;
  %e_edges = [e_edges(1)-0.5*de(1) e_edges+0.5*de];
  [VE,EE] = ndgrid(v_edges,e_edges);
  surf(hca,VE,EE,VE*0,(all_vdfs)')
  view(hca,[0 0 1])
  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f(v_M) (s/m^4)';
  hca.CLim = [-2 -1];
  hca.CLim = [0 0.04];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  shading(hca,'flat');
  hca.Layer = 'top';
  hca.Box = 'on';

  %legs = arrayfun(@(x)sprintf('E>%g eV',x),elow,'UniformOutput',false);
  %legend(hca,{'Original', legs{:}},'Box','off')
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_M (km/s)';
  hca.YLabel.String = 'E_{low} (eV)';
  %hca.YLim = [0 0.1];
  hca.YLim = [0 elow(end)];
  hca.XLim = [-1000 1000];

  hold(hca,'on')
  v = sqrt(2*units.eV*e_center/units.mp)*1e-3;
  plot(hca,v,e_center,'w')
  plot(hca,-v,e_center,'w')
  hold(hca,'off')
end


if 1 % Different elow, n
  hca = h(isub); isub = isub + 1;
    
  vdf_orig = pdist_.reduce('1D',Ndsl,'nMC',1000,'vg_edges',vg_edges);
  all_vdfs = zeros(numel(elow)+1,size(vdf_orig.depend{1},2));
  all_vdfs(1,:) = vdf_orig.data;
  hold(hca,'on')
  for iE = 1:numel(elow)  
    vdf_elim = pdist_.elim([elow(iE) Inf]).reduce('1D',Ndsl,'nMC',1000,'vg_edges',vg_edges);
    %plot(hca,vdf_elim.depend{1},vdf_elim.data)  
    all_vdfs(1+iE,:) = vdf_elim.data;
  end
  e_center = [0 elow]; de = diff(e_center);
  e_edges = -50:100:(elow(end)+50);
  v_edges = vdf_orig.ancillary.v_edges*1e-3;
  %e_edges = [e_edges(1)-0.5*de(1) e_edges+0.5*de];
  [VE,EE] = ndgrid(v_edges,e_edges);
  surf(hca,VE,EE,VE*0,(all_vdfs)')
  view(hca,[0 0 1])
  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f(v_N) (s/m^4)';
  hca.CLim = [-2 -1];
  hca.CLim = [0 0.04];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  shading(hca,'flat');
  hca.Layer = 'top';
  hca.Box = 'on';

  %legs = arrayfun(@(x)sprintf('E>%g eV',x),elow,'UniformOutput',false);
  %legend(hca,{'Original', legs{:}},'Box','off')
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_N (km/s)';
  hca.YLabel.String = 'E_{low} (eV)';
  %hca.YLim = [0 0.1];
  hca.YLim = [0 elow(end)];
  hca.XLim = [-1000 1000];

  hold(hca,'on')
  v = sqrt(2*units.eV*e_center/units.mp)*1e-3;
  plot(hca,v,e_center,'w')
  plot(hca,-v,e_center,'w')
  hold(hca,'off')
end


if 1 % Different elow, L
  hca = h(isub); isub = isub + 1;
    
  vdf_orig = pdist_.reduce('1D',Ldsl,'nMC',1000,'vg_edges',vg_edges);
  all_vdfs = zeros(numel(elow)+1,size(vdf_orig.depend{1},2));
  all_vdfs(1,:) = vdf_orig.data;
  hold(hca,'on')
  for iE = 1:numel(elow)  
    vdf_elim = pdist_.elim([elow(iE) Inf]).reduce('1D',Ldsl,'nMC',1000,'vg_edges',vg_edges);
    %plot(hca,vdf_elim.depend{1},vdf_elim.data)  
    all_vdfs(1+iE,:) = vdf_elim.data;
  end
  e_center = [0 elow]; de = diff(e_center);
  e_edges = -50:100:(elow(end)+50);
  v_edges = vdf_orig.ancillary.v_edges*1e-3;
  %e_edges = [e_edges(1)-0.5*de(1) e_edges+0.5*de];
  [VE,EE] = ndgrid(v_edges,e_edges);
  surf(hca,VE,EE,VE*0,(all_vdfs)')
  view(hca,[0 0 1])
  
  hcb = colorbar('peer',hca);
  hcb.YLabel.String = 'f(v_L) (s/m^4)';
  hca.CLim = [-2 -1];
  hca.CLim = [0 0.04];
  hca.XGrid = 'on';
  hca.YGrid = 'on';
  shading(hca,'flat');
  hca.Layer = 'top';
  hca.Box = 'on';

  %legs = arrayfun(@(x)sprintf('E>%g eV',x),elow,'UniformOutput',false);
  %legend(hca,{'Original', legs{:}},'Box','off')
  %irf_legend(hca,sprintf('%g<v<%g',vint(1),vint(2)),[0.02 0.98]);
  hold(hca,'off')
  hca.XLabel.String = 'v_L (km/s)';
  hca.YLabel.String = 'E_{low} (eV)';
  %hca.YLim = [0 0.1];
  hca.YLim = [0 elow(end)];
  hca.XLim = [-1000 1000];

  hold(hca,'on')
  v = sqrt(2*units.eV*e_center/units.mp)*1e-3;
  plot(hca,v,e_center,'w')
  plot(hca,-v,e_center,'w')
  hold(hca,'off')
end

if 0
  %%
  irf_plot_axis_align(h([1 4]))
  irf_plot_axis_align(h([2 5]))
  irf_plot_axis_align(h([3 6]))
end
 

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
c_eval('h(?).LineWidth = 1.5;',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))
%%
time = time_xline + 15;

nMovMean = 20;

c_eval('pdist_orig = iPDist?.movmean(nMovMean).tlim(time+0.5*0.15*[-1 1]);',ic)
c_eval('pdist_mm = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts).tlim(time+0.5*0.15*[-1 1]);',ic)





omni = pdist_orig.omni;
omni_mm = pdist_mm.omni;

% make fit
%x_data_mm = log10(omni_mm.depend{1});
%y_data_mm = log10(omni_mm.data);
%FO_mm = fit(x_data_mm', y_data_mm', 'linearinterp')

elim_fit = [0 1000];
x_data = log10(omni.elim(elim_fit).depend{1});
y_data = log10(omni.elim(elim_fit).data);
y_data(isnan(y_data)) = 0;
FO = fit(x_data', y_data', 'poly1');
FO_fit = FO.p1*x_data + FO.p2;
FO_fit_all = FO.p1*log10(omni.depend{1}) + FO.p2;

%%
h = setup_subplots(1,1);
isub = 1;

if 0
hca = h(isub); isub = isub + 1;
loglog(hca,omni.depend{1},omni.data,...
           omni_mm.depend{1},omni_mm.data)
end

hca = h(isub); isub = isub + 1;
plot(hca,log10(omni.depend{1}),log10(omni.data),...
           x_data,y_data,'*',...
           x_data,FO_fit,...
           log10(omni_mm.depend{1}),log10(omni_mm.data),...
           log10(omni_mm.depend{1}),log10(omni_mm.data+10.^FO_fit_all),'--',...
           'linewidth',2)
legend(hca,{'Original','Fitted data',sprintf('Fit: y = %g*x + %g',FO.p1,FO.p2),'moving mean removal','moving mean removal + fit'},'box','off')
hca.XLabel.String = 'log_{10} E (eV)';
hca.YLabel.String = 'log_{10} PSD (s^3/cm^6)';
irf_legend(hca,sprintf('N = %g',nMovMean),[0.98 0.4],'k')

h(1).Title.String = omni.time.utc('yyyy-mm-ddTHH:MM:SS.mmmZ');
%%
hca.XLim = vlim;
hca.XLabel.String = 'v_M (km/s)';
hca.YLabel.String = 'f_i (s/m^4)';
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 2;',1:numel(hl))
hca.XGrid = 'on';
hca.YGrid = 'on';
t1t2 = pdist.time + 0.5*nMovMean*0.150*[-1 1];
irf_legend(hca,{sprintf('one-counts removed, nMovMean=%g',nMovMean),'original',sprintf('E > %g eV',elim(1))}',[0.98,0.98],'fontsize',fontsize);

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
nMovMean = 5;

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
dt_all = [-10:5:10]+-1.5;
%dt_all = [-6:2:6]+0;
%dt_all = [-6:2:6]+25;
%dt_all = [-6:2:6]-00;

[h1,h] = initialize_combined_plot('topbottom',2,2,numel(dt_all),0.2,'vertical');

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
  if 0 % f(L,M)
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

  if 0 % f(L,M)
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
  

  if 1 % f(L,N), force EL + vN*BM
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Ndsl);
    %vdf.depend{1} = vdf.depend{1} - vL_Xline;
    %vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    vdf_data = squeeze(mean(vdf.data,1));
    
    vL_data = squeeze(vdf(1).depend{1});
    vN_data = squeeze(vdf(1).depend{2});
    fEL = E.resample(pdist).x.data*1e-3;
    fBM = B.resample(pdist).y.data*1e-9;
    fBL = B.resample(pdist).x.data*1e-9;
    [VL,VN] = ndgrid(vL_data,vM_data);
    %[id,VL,VM] = ndgrid(1:pdist.length,vL_data,vM_data);
    force_data = mean(fEL) - VN*1e3*mean(fBM);% - VM*1e3*mean(fBL);
    force_data(vdf_data==0) = 0;
    % f = EN + vL*BM - vM*BL;

    %vdf.plot_plane(hca,'smooth',nSmooth)
    pcolor(hca,vL_data,vN_data,force_data')
    shading(hca,'flat')
    colormap(hca,pic_colors('blue_gray_red'))
    hca.CLim = max(abs(hca.CLim))*[-1 1];

    %vdf.plot_plane(hca,'smooth',nSmooth)
    axis(hca,'square')
    hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_N (km/s)';
    hcb = colorbar('peer',hca);
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


%colormap(pic_colors('candy_gray'))

%hlinks_LM = linkprop(h(1:numel(dt_all):end),{'CLim'});
%hlinks_LN = linkprop(h(2:numel(dt_all):end),{'CLim'});
%hlinks_MN = linkprop(h(3:numel(dt_all):end),{'CLim'});

%hlinks_LM = linkprop(h(1:2:end),{'View'});
%hlinks_LN = linkprop(h(2:2:end),{'View'});
%hlinks_MN = linkprop(h(3:2:end),{'View'});

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



elim = [3000 Inf];
nMovMean = 7;
%c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts);',ic)
c_eval('pdist_all = iPDist?.movmean(nMovMean);',ic)
%pdist_all.data(:,1:20,:,:) = 0;
pdist_all = pdist_all.elim(elim);

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
dt = -3;
dt = 0;
time = time + dt;

h = setup_subplots(1,1);
isub = 1;

pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
iso_values = [1e-27];
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

%% Compare single time omni spectra for moving averages and not

it  = 3;
time = times(it);
  
c_eval('pdist_with_noise = iPDist?.movmean(nMovMean).tlim(time+0.5*0.15*[-1 1]);',ic)
c_eval('pdist_with_noise_elim = iPDist?.elim(elim).movmean(nMovMean).tlim(time+0.5*0.15*[-1 1]);',ic)


pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);
tint_dist = pdist.time + 0.5*0.150*nMovMean*[-1 1];
elow = max(tsElow.tlim(tint_dist).data);  
%elow = 2000;
pdist = pdist.elim([elow Inf]);

%% Reduced f(vM), how to best illustrate the presence of two populations

fT = 11;
tint = time_xline_ion + 0.5*fT*[-1 1];

nMovMean = 5;
c_eval('pdist_all = iPDist?.movmean(nMovMean).tlim(tint);',ic)
%c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim);',ic)
%c_eval('elows_tmp = elows_all = tsElow.movmean(nMovMean).tlim(tint);',ic)
elows_tmp = tsElow; elows_tmp.data = movmean(elows_tmp.data,nMovMean);
elows_all = elows_tmp.tlim(tint);

%%
n_dists = pdist_all.length;

t_per_panel = 5;
npanels = ceil(n_dists/t_per_panel);
nrows = 3;
ncols = ceil(npanels/nrows);

[h1,h] = initialize_combined_plot('topbottom',2,nrows,ncols,0.3,'vertical');


isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); 

if 1 % Vi
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  %irf_legend(hca,{'L','M','N'},[0.98 0.98]);
  irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
  %irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
  hca.YLabel.Interpreter = 'tex';
end
if 2 % Vi
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{tsElow},''comp'')',ic)
  hca.YLabel.String = 'E_{low} (eV)';
  hca.ColorOrder = mms_colors('1');
  %irf_legend(hca,{'L','M','N'},[0.98 0.98]);
  %irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
  %irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
  hca.YLabel.Interpreter = 'tex';
end

colors = pic_colors('matlab');
isub = 1;
times_all = pdist_all.time;
for ip = 1:npanels
  its = t_per_panel*(ip-1)+(1:t_per_panel);
  its(its>n_dists) = []; 
  times = times_all(its);    
  pdist = pdist_all(its);
  %pdist.data(:,:,:,[1 end]) = 0;
  %elows1 = tsElow.resample(pdist);
  elows2 = elows_all(its);
  %elows1.data
  elows2.data;
  elows = elows2*1+00;
  
  %elow = 2000;  
  %pdist.data(:,:,:,1) = 0;
  %pdist.data(:,:,:,end) = 0;
  
  tint_dist = [times.start times.stop];
  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))

  

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  
  hca = h(isub); isub = isub + 1;

  if 1 % 2D at a certain vn range
    
    vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows,'nMC',500);

    vint = [500 inf];
    vint = [-inf -500];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    vidx = find(all([v_center>vint(1); v_center<vint(2)]',2));
    dv = v_center(2)-v_center(1);
    data = vdf.data;
    data = data(:,:,vidx);
    data = sum(data,3)*dv*1e3;
    


    hca.ColorOrder = colors;
    plot(hca,v_center,data)
    %axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    
    
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    irf_legend(hca,sprintf('%g < v_N < %g',vint(1),vint(2)),[0.2 0.98],'fontsize',10)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
  end
  if 0 % 1D
    vdf = pdist.reduce('1D',Mdsl,'lowerelim',elows);
      
    hca.ColorOrder = colors;
    plot(hca,vdf.depend{1}(1,:),vdf.data)
    %axis(hca,'square')
    hca.XLabel.String = 'v_M (km/s)';
    
    
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
  end

end

drawnow
compact_panels(h,0.01,0.005)

hlinks = linkprop(h,{'YLim','XLim'});
h(1).YLim = [0 0.03];


c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).XMinorGrid = ''on'';',1:numel(h))


delete(h((npanels+1):end))


irf_zoom(h1,'x',tint_zoom)






%% Reduced f(vL), how to best illustrate the gyroturning

fT = 2;
tint = time_xline_ion + 0.5*fT*[-1 1] +0;

nMovMean = 5;
c_eval('pdist_all = iPDist?.movmean(nMovMean).tlim(tint);',ic)
%c_eval('pdist_all = iPDist?.movmean(nMovMean).elim(elim);',ic)
%c_eval('elows_tmp = elows_all = tsElow.movmean(nMovMean).tlim(tint);',ic)
elows_tmp = tsElow; elows_tmp.data = movmean(elows_tmp.data,nMovMean);
elows_all = elows_tmp.tlim(tint);

%%
n_dists = pdist_all.length;

t_per_panel = 1;
npanels = ceil(n_dists/t_per_panel);
nrows = 3;
ncols = ceil(npanels/nrows);

[h1,h] = initialize_combined_plot('topbottom',1,nrows,ncols,0.3,'vertical');


isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); 

if 1 % Vi
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  %irf_legend(hca,{'L','M','N'},[0.98 0.98]);
  irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
  %irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
  hca.YLabel.Interpreter = 'tex';
end
if 0 % Elow
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{tsElow},''comp'')',ic)
  hca.YLabel.String = 'E_{low} (eV)';
  hca.ColorOrder = mms_colors('1');
  %irf_legend(hca,{'L','M','N'},[0.98 0.98]);
  %irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
  %irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
  hca.YLabel.Interpreter = 'tex';
end

colors = pic_colors('matlab');
isub = 1;
times_all = pdist_all.time;
for ip = 1:npanels
  its = t_per_panel*(ip-1)+(1:t_per_panel);
  its(its>n_dists) = []; 
  times = times_all(its);    
  pdist = pdist_all(its);
  pdist.data(:,:,:,[1 end]) = 0;
  %elows1 = tsElow.resample(pdist);
  elows2 = elows_all(its);
  %elows1.data
  elows2.data;
  elows = elows2*1+00;
  
  %elow = 2000;  
  %pdist.data(:,:,:,1) = 0;
  %pdist.data(:,:,:,end) = 0;
  
  tint_dist = [times.start times.stop];
  c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))

  

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  
  hca = h(isub); isub = isub + 1;

  if 1 % 1D at a certain vn range
    v_edges = -1000:1000:1000;
    nv = numel(v_edges) - 1;

    vdf = pdist.reduce('2D',Ldsl,Mdsl,'lowerelim',elows);        
    
    data = zeros(nv,size(vdf.data,3));
    legs_edges = cell(nv,1);
    v_center = vdf.depend{1}(1,:);   
    dv = v_center(2)-v_center(1);
    for iv = 1:nv
      vint = [v_edges(iv) v_edges(iv+1)];
      vidx = find(all([v_center>vint(1); v_center<vint(2)]',2));      
      data_tmp = vdf.data;
      data_tmp = data_tmp(:,:,vidx);
      data(iv,:) = sum(data_tmp,3)*dv*1e3;
      legs_edges{iv} = sprintf('%g<v_M<%g',vint(1),vint(2));
    end
    
    
       
    

    hca.ColorOrder = [colors; colors];
    plot(hca,v_center,data)    
    hca.XLabel.String = 'v_M (km/s)';      
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',7)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = [colors; colors];
    irf_legend(hca,legs_edges,[0.98 0.98],'fontsize',10)
    %irf_legend(hca,sprintf('%g < v_N < %g',vint(1),vint(2)),[0.2 0.98],'fontsize',10)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
  end

end

drawnow
compact_panels(h,0.01,0.005)

hlinks = linkprop(h,{'YLim','XLim'});
h(1).YLim = [0 0.03];

hl = findobj(gcf,'type','line'); 
c_eval('hl(?).LineWidth = 1;',1:numel(hl))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).XMinorGrid = ''on'';',1:numel(h))


delete(h((npanels+1):end))


irf_zoom(h1,'x',tint_zoom)





