%% Set up database 
%mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cecilianorgren/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
mms.db_init('local_file_db','/Volumes/mms');
%db_info = datastore('mms_db');

units = irf_units;
% Torbert event 
ic = 3;
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z');

%% Load data
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


c_eval('EdotB? = (gseE?.resample(gseB?).x*gseB?.x + gseE?.resample(gseB?).y*gseB?.y + gseE?.resample(gseB?).z*gseB?.z)*1e-3*1e-9;',ic)
%% Clean distributions, do it once hear so it will be the same for all
nMean = [5,3,3,3]; nThresh = 3;

PD_orig = iPDist3;
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

%% Define times, etc... things that are common for the entire study
tint_figure = irf.tint('2017-07-11T22:33:00.00Z/2017-07-11T22:35:00.00Z');
%time_xline = irf_time('2017-07-11T22:34:03.00Z','utc>EpochTT');
time_xline = irf_time('2017-07-11T22:34:02.60Z','utc>EpochTT');
time_xline_ion = irf_time('2017-07-11T22:34:02.00Z','utc>EpochTT');
time_vdf = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
tint_figure_zoom = irf.tint('2017-07-11T22:33:30.00Z/2017-07-11T22:34:30.00Z');
tint_figure_zoom_incl_sep = irf.tint('2017-07-11T22:33:23.00Z/2017-07-11T22:34:30.00Z');
tint_figure_edr = irf.tint('2017-07-11T22:33:55.00Z/2017-07-11T22:34:12.00Z');
tint_figure_zoom_inner_idr = irf.tint('2017-07-11T22:33:50.00Z/2017-07-11T22:34:20.00Z');

vL_xline = -170;

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

% This one is defined by eye so that the up/down going double M-populations
% align in M
L_gse = [1 0 -.2]; L_gse = L_gse/norm(L_gse);
M_gse = [0 1 -0.2]; M_gse = cross(L_gse,cross(M_gse,L_gse)); M_gse = M_gse/norm(M_gse);
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

% This one makes the counterstreaming (vN) ions summetric about vL-vL_xline
% (they're both at vL-vL_xline = 0).
L_gse = [1 0 -.1]; L_gse = L_gse/norm(L_gse);
M_gse = [0 1 -0.2]; M_gse = cross(L_gse,cross(M_gse,L_gse)); M_gse = M_gse/norm(M_gse);
N_gse = cross(L_gse,M_gse);
lmn_gse = [L_gse; M_gse; N_gse];

% Based on eigenvalues of p
dbcsP = PD_clean.p;

[tsEigVal,tsEigV1,tsEigV2,tsEigV3] = dbcsP.eig;
tint_N = time_xline+0.1*[-1 1];
tint_N = EpochTT('2017-07-11T22:34:03.000Z')+0.5*nMean(1)*0.151*[-1 1];
c_eval('gseEigV1 = mms_dsl2gse(tsEigV1,defatt?,1);',ic)
c_eval('gseEigV2 = mms_dsl2gse(tsEigV2,defatt?,1);',ic)
c_eval('gseEigV3 = mms_dsl2gse(tsEigV3,defatt?,1);',ic)

N_gse = mean(gseEigV1.tlim(tint_N).data,1); N_gse = N_gse/norm(N_gse);
L_gse = cross(cross(N_gse,[1 0 0]),N_gse); L_gse = L_gse/norm(L_gse); % L mostly towards [1 0 0] GSE
M_gse = cross(N_gse,L_gse); M_gse = M_gse/norm(M_gse);
lmn_eig = [L_gse; M_gse; N_gse];

N_gse = mean(gseEigV1.tlim(tint_N).data,1); N_gse = N_gse/norm(N_gse);
M_gse = mean(gseEigV2.tlim(tint_N).data,1); M_gse = M_gse/norm(M_gse);
M_gse = cross(cross(N_gse,M_gse),N_gse); M_gse = M_gse/norm(M_gse);
L_gse = cross(M_gse,N_gse); L_gse = L_gse/norm(L_gse);
lmn_eig2 = [L_gse; M_gse; N_gse];

% Based on eigenvalues of p
dbcsP = PD_clean.movmean(nMean(1)).p;

[tsEigVal,tsEigV1,tsEigV2,tsEigV3] = dbcsP.eig;
tint_N = time_xline+0.1*[-1 1];
tint_N = EpochTT('2017-07-11T22:34:03.000Z')+0.5*1*0.151*[-1 1];
c_eval('gseEigV1 = mms_dsl2gse(tsEigV1,defatt?,1); gseEigV1 = gseEigV1.tlim(tint_N);',ic)
c_eval('gseEigV2 = mms_dsl2gse(tsEigV2,defatt?,1); gseEigV2 = gseEigV2.tlim(tint_N);',ic)
c_eval('gseEigV3 = mms_dsl2gse(tsEigV3,defatt?,1); gseEigV3 = gseEigV3.tlim(tint_N);',ic)
EigVal = tsEigVal.tlim(tint_N);

N_gse = mean(gseEigV1.tlim(tint_N).data,1); N_gse = N_gse/norm(N_gse);
L_gse = cross(cross(N_gse,[1 0 0]),N_gse); L_gse = L_gse/norm(L_gse); % L mostly towards [1 0 0] GSE
M_gse = cross(N_gse,L_gse); M_gse = M_gse/norm(M_gse);
lmn_eig = [L_gse; M_gse; N_gse];

% The eigenvectors directly
N_gse = mean(gseEigV1.tlim(tint_N).data,1); N_gse = N_gse/norm(N_gse);
M_gse = mean(gseEigV2.tlim(tint_N).data,1); M_gse = M_gse/norm(M_gse);
M_gse = cross(cross(N_gse,M_gse),N_gse); M_gse = M_gse/norm(M_gse);
L_gse = cross(M_gse,N_gse); L_gse = L_gse/norm(L_gse);
lmn_eig2 = [L_gse; M_gse; N_gse];

lmn = lmn_vi;
lmn = lmn_gse;
%lmn = lmn_edr;
lmn = lmn_eig;
lmn = lmn_eig2;
L = lmn(1,:);
M = lmn(2,:);
N = lmn(3,:);
lmn = [L; M; N];
%lmn


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
mvaCurvB = gseCurvB*lmn';
mvaJcurl = gseJcurl*lmn';

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


%% Reduce 1D distributions
c_eval('fi?_L = PD_clean.reduce(''1D'',L);',ic)
c_eval('fi?_M = PD_clean.reduce(''1D'',M);',ic)
c_eval('fi?_N = PD_clean.reduce(''1D'',N);',ic)

v_cut_L = -170;
c_eval('fi?_L_pos_neg = PD_clean.reduce(''1D'',L,''vg_edges'',v_cut_L+[-3000, 0, 3000]);',ic)



%% Calculate % of ions moving towards or away from the X line
tint_left = [fi3_L.time(1) time_xline_ion];
tint_right = [time_xline_ion fi3_L.time(fi3_L.length)];

%[f_left,ind] = fi3_L.tlim([fi3_L.time(1) time_xline_ion]);
fi3_L_neg = fi3_L.data(:,1:size(fi3_L.data,2)/2);
fi3_L_pos = fi3_L.data(:,size(fi3_L.data,2)/2+1:end);
[ind_left,f_left] = fi3_L.time.tlim([fi3_L.time(1) time_xline_ion]);
[ind_right,f_right] = fi3_L.time.tlim([time_xline_ion fi3_L.time(fi3_L.length)]);


sum_neg = fi3_L_pos_neg.data(:,1);
sum_pos = fi3_L_pos_neg.data(:,2);
sum_all = sum_neg + sum_pos;
ratio_neg = sum_neg./sum_all;
ratio_pos = sum_pos./sum_all;
tsFneg = irf.ts_scalar(fi3_L.time,ratio_neg);
tsFpos = irf.ts_scalar(fi3_L.time,ratio_pos);

tsFaway = tsFneg.tlim(tint_left).combine(tsFpos.tlim(tint_right)); tsFaway.name = 'Away';
tsFtowa = tsFpos.tlim(tint_left).combine(tsFneg.tlim(tint_right)); tsFtowa.name = 'Towards';

tsFracEarthward = irf.ts_scalar(fi3_L_pos_neg.time,fi3_L_pos_neg.data(:,1)./sum(fi3_L_pos_neg.data(:,:),2));
tsFracTailward = irf.ts_scalar(fi3_L_pos_neg.time,fi3_L_pos_neg.data(:,2)./sum(fi3_L_pos_neg.data(:,:),2));
data_away = zeros(tsFracEarthward.length,1);
idxt = find(tsFracEarthward.time>time_xline);
data_away = tsFracEarthward.data;
data_away(idxt) = tsFracTailward.data(idxt);
tsFracAway = irf.ts_scalar(fi3_L_pos_neg.time,data_away);

%% Calculation of bounce times etc
f_wb = @(B,L,v) sqrt(units.e*B*v/units.mp/L);
f_fb = @(B,L,v) f_wb(B,L,v)/2/pi;
f_Tb = @(B,L,v) f_fb(B,L,v).^(-1);

B0 = 5e-9;
L_ = 600e3;
L_ = 200e3;
L_ = 1170e3;
L_ = 1170e3/4;
L_ = 400e3;
L_ = 800e3;
v = 1000e3;
wb = f_wb(B0,L_,v);
fb = f_fb(B0,L_,v);
Tb = f_Tb(B0,L_,v);

disp(sprintf('B = %g nT, L = %g km, v = %g km/s: wb = %.2f rad/s, fb = %.2f Hz, Tb = %.2f s, Tb/2 = %.2f s',B0*1e9,L_*1e-3,v*1e-3,wb,fb,Tb,Tb/2))

%% Calculate azimuthal distribution in the LM plane (using macroparticles)
units = irf_units;
PD = PD_clean;

nt = PD.length;
Vsc = scPot3.resample(PD);
nMP = 5000;
MP = PD.macroparticles('ntot',nMP,'scpot',Vsc);
vL_shift = vL_xline*1;
%vL_shift = 0;
%vL_shift = -999;
E_edges = [PD.ancillary.energy(1,1) - PD.ancillary.delta_energy_minus(1,1), PD.ancillary.energy(1,:) + PD.ancillary.delta_energy_plus(1,:)];
E_minus = E_edges(1:end-1);
E_plus = E_edges(2:end);
v_minus = sqrt(2*units.e*E_minus/units.mp); % m/s
v_plus = sqrt(2*units.e*E_plus/units.mp); % m/s
% Cylindrical volume element v*dv*dtheta
% vol = (v2^2-v1^2)/2*Delta_theta
d_vel = (v_plus.^2 - v_minus.^2)/2; % (m/s)^2

nEnergy = numel(E_edges)-1;
dazim = 10;
azimuth_edges = -180:dazim:180;
nAzimuth = numel(azimuth_edges)-1;

% 2D volume, radial-azmiuthal bins
d_vel_mat = repmat(d_vel,nt,1,nAzimuth);
d2v = d_vel_mat.*dazim*pi/180; % (m/s)^2


dn_tot_all = zeros(nt,nEnergy,nAzimuth);
%df_tot_all = zeros(nt,nEnergy,nAzimuth);
%dv_tot_all = zeros(nt,nEnergy,nAzimuth);

for it = 1:nt
  Ltmp = tsLdsl3(it);
  Mtmp = tsMdsl3(it);
  Ntmp = tsNdsl3(it);

  dv = MP(it).dv; % same units as the original f (I think)
  df = MP(it).df;
  dn = df.*dv;
  vx = MP(it).vx;
  vy = MP(it).vy;
  vz = MP(it).vz;
  vL = [vx vy vz]*Ltmp.data' - 1*vL_shift;
  %vL = [vx vy vz]*Ltmp.data' - -1000;
  vM = [vx vy vz]*Mtmp.data';
  vN = [vx vy vz]*Ntmp.data';
  theta_NL = atan2d(vN,vL);
  v2 = (vL.^2 + vM.^2 + vN.^2)*1e6; % (cm/s)^2 -> (m/s)^2
  vLN = sqrt(vL.^2 + vN.^2)*1e3; % (m/s)^1 -- It's the energy corresponding to this one should use, because the end results will be a reduced cylinder
  E = units.mp*vLN.^2/2/units.eV; % eV
  % E(E<elow) = 0;
  % dn(vM<0) = 0;
  % dn(vLN<1500) = 0;  
  fun = @sum;
  [dn_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dn,'Fun',fun);
  %[dv_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dv,'Fun',fun);
  %[df_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',df,'Fun',fun);

  dn_tot_all(it,:,:) = dn_tot;
  %dv_tot_all(it,:,:) = dv_tot;
  %df_tot_all(it,:,:) = df_tot;
end

%f_tot = dn_tot_all./d2v; % [1/cm^3]/[m^2/s^2]

% [dn_tot_all] = cm^-3
% [d2v] = m^2/s^2
% cm^-3/(cm^2/s^2) = s^2/cm^5, reduced along one dimension
f_tot = dn_tot_all*1e6./d2v; 

iAzim = PDist(PD.time,f_tot,'azimuthangle',PD.depend{1},mid{2}); % scaling factor applied above, units in s^2/cm^5
iAzim.ancillary.v_shift = vL_shift;
to_SI = 1;%(1e2)^5; % s^2/cm^5 -> s^2/m^5
iAzim = iAzim*to_SI;
iAzim.units = 's^2/m^5';

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

h = irf_plot(8);
c_eval('h(?).Position(2) = h(?).Position(2)-0.02;',1:numel(h))
fontsize = 13;
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
  scale = 1e3;
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint)/scale,mvaVe?.y.tlim(tint)/scale,mvaVe?.z.tlim(tint)/scale},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e','(10^3 km/s)'};
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
  irf_plot(hca,tsElow,'color',[0 0 0],'linewidth',3)
  %irf_legend(hca,{'E_{low}'},[0.02 0.5],'color','k','fontweight','bold','fontsize',fontsize)
  %irf_legend(hca,{'E_{low}'},[0.02 0.5],'color','k','fontweight','bold','fontsize',fontsize,'backgroundcolor','w')
  hold(hca,'off')
  %hca.YLabel.String = 'E_{low} (eV)';
  hca.YLabel.String = {'E_i','(eV)'};
  hca.YLabel.Interpreter = 'tex';
  hca.Color = color_nan;  
  hca.XGrid = 'off';
  hca.YGrid = 'off';
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
if 1 % fi red L
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
  c_eval('irf_plot(hca,irf.ts_scalar(mvaVi?.time,repmat(vL_xline,[mvaVi?.length,1])),''k--'')',ic)
  hca.NextPlot = "replace";

  hca.YLabel.String = {'v_{iL}','(km/s)'};
  hca.Color = color_nan;

  %irf_legend(hca,sprintf('N_{mean} = %g',nMovMean),[0.02,0.1],'fontsize',fontsize);
  %irf_legend(hca,{sprintf('N_{mean} = %g',nMovMean),' one-counts removed'},[0.02,0.98],'fontsize',fontsize,'color','k');
end
if 1 % fi red M
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
if 1 % fi red N
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

if 1 % fraction of ions
  hca = irf_panel('fraction away from x line');
  set(hca,'ColorOrder',mms_colors('24'))
  

  
  irf_plot(hca,{tsFracEarthward,tsFracTailward},'comp')
  hold(hca,'on')
  hl = irf_patch(hca,{tsFracAway,0},'facecolor',[0 0 0],'facealpha',0.2,'edgecolor','none');
  hl.EdgeColor = 'none';
  hl.FaceColor = [0.5 0.5 0.5];
  hl.FaceAlpha = 0.2;
  hold(hca,'off')
  
  hca.YLabel.String = {'Fraction','of ions'};
  
  colors = [mms_colors('24'); 0.3 0.3 0.3];  
  set(hca,'ColorOrder',colors)  
  %irf_legend(hca,{'moving tailward','moving Earthward','moving away from X line'}',[0.2 0.7])
  irf_legend(hca,{'moving tailward','    moving Earthward','           moving away from X line'},[0.045 0.98],'fontsize',fontsize-3)
  %hca.YTick = [0:0.25:1];
end

irf_plot_axis_align
%irf_zoom(h,'x',tint_figure)
%irf_zoom(h,'x',tint_figure_zoom)
irf_zoom(h,'x',tint_figure_zoom_incl_sep)

irf_zoom(h(1:3),'y')
irf_pl_mark(h,time_xline+0,'black','linestyle','-.')
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

c_eval('h(?).XGrid = ''off''; h(?).YGrid = ''off'';',1:numel(h))

hlinks = linkprop(h([5 6 7]),{'CLim'});
h(6).CLim = [-4   -1];
%colormap(irf_colormap('waterfall'))
%colormap(pic_colors('candy_gray'))
colormap(pic_colors('candy6'))
colormap(irf_colormap('magma'))
%hlinks = linkprop(h([4 5]),{'CLim'});
h(4).CLim = [ 3   6.0];
%hb = findobj(gcf,'type','colorbar');
%delete(hb(end))
%hb(4).Position(4) = hb(4).Position(4)*2.1;

% Add length on top
if 1
  %ax2 = axes('position',h(1).position);
  delete(ax2)
  ax2 = axes('Position',h(1).Position,'XAxisLocation','top','YAxisLocation','right','color','none','TickDir','out');
  userdata = get(gcf,'userdata');
  %
  di = mean(di3.tlim(time_xline+5*[-1 1]).data,1);
  tstart = irf_time(userdata.t_start_epoch,'epoch>EpochTT');
  dt = time_xline - tint_figure_zoom_incl_sep(1);
  dt = time_xline - tstart;
%dt = 0;



  xlim_s = h(1).XLim-1*dt; % time
  
  xlim_km = xlim_s*170;
  xlim_di = xlim_km/di;
  
  %xdata = (xlim_km(1):1000:xlim_km(end))/170;
  %xdata = (-10000:1000:10000)/170;
  
  xdata_di = -100:1:100;
  
  ax2.XTick = xdata_di+0*dt/di;

  ax2.XLim = xlim_di+0*dt/di;

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
  %irf_legend(h(ii),legends{nInd},[0.01 0.98],'color',[0 0 0],'fontsize',fontsize+1)
  irf_legend(h(ii),legends{nInd},[-0.18 0.999],'color',[0 0 0],'fontsize',fontsize+1)
  nInd = nInd + 1;
  %h(ii).FontSize = 16;
  h(ii).YLabel.Position(1) = -0.10;
end


if 0
  %%
  h(3).YLim = 1499*[-1 1];
  h(4).YLim = 1499*[-1 1];
  h(5).YLim = 1499*[-1 1];
  h(5).CLim = [-3 -1];
  colormap([0.9 0.9 0.9; pic_colors('candy4')])
end

irf_legend(0,{sprintf('L=[%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3)),sprintf('nMean=[%g,%g,%g,%g], nThresh = %g',nMean(1),nMean(2),nMean(3),nMean(4),nThresh)}',[0.05 0.01],'fontsize',fontsize-5)

%% Figure: Hall fields
ic = 3;
h = irf_plot(4);

fontsize = 14;
fhighE = 10;

if 1 % B lmn
  hca = irf_panel('B LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaB?.x,mvaB?.y,mvaB?.z},''comp'');',ic)
  hca.YLabel.String = {'B (nT)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
end 
if 1 % E lmn
  hca = irf_panel('E LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))  
  c_eval('irf_plot(hca,{mvaE?.x.filt(0,fhighE,[],5),mvaE?.y.filt(0,10,[],5),mvaE?.z.filt(0,10,[],5)},''comp'');',ic)
  hca.YLabel.String = {'E (mV/m)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98 0.98],'fontsize',fontsize);
  irf_legend(hca,{sprintf('E < %g Hz',fhighE)},[0.08 0.98],'fontsize',fontsize,'color','k');
end 
if 1 % Ve
  hca = irf_panel('Ve LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVe?.x.tlim(tint)*1e-3,mvaVe?.y.tlim(tint)*1e-3,mvaVe?.z.tlim(tint)*1e-3},''comp'');',ic)  
  
  hca.YLabel.String = {'v_e (10^3 km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
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
irf_pl_mark(h,time_xline,'black','linestyle','-.')
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
  irf_legend(h(ii),legends{nInd},[0.02 0.98],'color',[0 0 0],'fontsize',fontsize)
  nInd = nInd + 1;
  %h(ii).FontSize = 12;
end

if 0
  %%
  h(3).YLim = 1499*[-1 1];
  h(4).YLim = 1499*[-1 1];
  h(5).YLim = 1499*[-1 1];
  h(5).CLim = [-3 -1];
  colormap([0.9 0.9 0.9; pic_colors('candy4')])
end

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
c_eval('h(?).LineWidth = 1.5;',1:numel(h))

irf_legend(0,{sprintf('L=[%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3)),sprintf('nMean=[%g,%g,%g,%g], nThresh = %g',nMean(1),nMean(2),nMean(3),nMean(4),nThresh)},[0.05 1])
irf_plot_axis_align

%% Figure: Azimuthal angle (LN plane) spectrogram
 
tint_plot = [iAzim_0170.time.start iAzim_0170.time.stop];
tint_plot = tint_figure_zoom_incl_sep;

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:50.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:10.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

times_utc = [...%'2017-07-11T22:33:45.000Z';...
             '2017-07-11T22:33:54.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             %'2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:00.502Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:05.500Z';...
             %'2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:10.000Z';...
             '2017-07-11T22:34:20.000Z'];

iVDF = 1;
time = EpochTT(times_utc(iVDF,:)); % Time of example distribution

h1 = irf_plot(2);

h1(1).Position = [0.5 0.75 0.44 0.15];
h1(2).Position = [0.5 0.15 0.44 0.6];

h2 = subplot(1,2,1);
h2.Position = [0.1 0.15 0.5 0.5];
%h2.Position = [0.1181    0.2759    0.2267    0.5810];
h2.Position = [0.1581    0.2759    0.2267    0.5810];

fontsize = 10;
elim = [2000 10000];
vlim_data = sqrt(2*elim*units.eV/(units.mp))*1e-3; % km/s

if 1 % Vi
  hca = irf_panel('Vi LMN');
  set(hca,'ColorOrder',mms_colors('xyza'))
  c_eval('irf_plot(hca,{mvaVi?.x.tlim(tint_plot),mvaVi?.y.tlim(tint_plot),mvaVi?.z.tlim(tint_plot)},''comp'');',ic)    
  hca.YLabel.String = {'v_i (km/s)'};
  set(hca,'ColorOrder',mms_colors('xyza'))
  irf_legend(hca,{'L','M','N'},[0.98,0.98],'fontsize',fontsize);
end
if 1 % iAzim_170
  hca = irf_panel('azim iAzim_170');
  ts = iAzim.mask({[tsElow.data*0 tsElow.data*1]});
  ts = ts.elim(elim);
  %ts = iAzim_4;
  specrec = ts.tlim(tint_plot).specrec('pitchangle');
  specrec.p = smooth2(specrec.p,0,0);
  %specrec.df = 1;
  [hca_, hcb] = irf_spectrogram(hca,specrec,'donotfitcolorbarlabel'); 
  hcb.YLabel.String = 'log_{10} f_i (s^3/m^6)';
  hca.Layer = 'top';
  hca.YTick = -180:30:180;
  hold(hca,'on')
  irf_plot(hca,irf.ts_scalar(tint_plot,90*[1 1]),'k--')
  irf_plot(hca,irf.ts_scalar(tint_plot,-90*[1 1]),'k--')
  hold(hca,'off')  
  hca.YLabel.String = '\theta_{LN} (deg)';
  hca.YLabel.Interpreter = 'tex';
  %irf_legend(hca,{sprintf('v_L -> v_L-(%g km/s)',iAzim_0170.ancillary.v_shift)},[0.98 0.98],'color','k','fontsize',12)
  hca.XGrid = 'on'; hca.YGrid = 'on';
  %hca.XGrid = 'off'; hca.YGrid = 'off';
  %hca.CLim = [-30 -26]+30;
  %hca.CLim = [1 4];
end

irf_plot_axis_align(h1)
h1(1).YLim(1) = [-799];

vL_xline_use = iAzim.ancillary.v_shift;


pdist = PD_clean.tlim(time+nMean(1)*0.15*0.5*[-1 1]);
tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*[-1 1];


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


scaxis_scale = 2000;
nSmooth = 0;
nContours = 0;

  isub = 1;
    vlim = vlim_data(2)*2.3;
  if 1 % f(L,N)
    hca = h2(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl],'vint',[-Inf Inf]);
    vdf.depend{1} = vdf.depend{1} - 1*vL_xline_use;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - 1*vL_xline_use;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    %hca.XLabel.String = 'v_L (km/s)';
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
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
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    hca.XTick = -4000:2000:4000;
    hca.YTick = -4000:2000:4000;
    hca.XMinorGrid = 'on';
    hca.YMinorGrid = 'on';
    %hca.XMinorTick = -4000:2000:4000;
    %hca.YMinorTick = -3000:1000:3000;
    if 1 % add elim circular markings and angular marking combined
      hold(hca,'on')      
      ang = linspace(0,360,361);
      color_grid = ([0.9 0.9 0.9]-0.3)*0;
      plot(hca,vlim_data(1)*cosd(ang), vlim_data(1)*sind(ang),'color',color_grid,'linestyle','-')
      plot(hca,vlim_data(2)*cosd(ang), vlim_data(2)*sind(ang),'color',color_grid,'linestyle','-')

      angs = -180:30:180;
      for ia = 1:numel(angs)
        xx = vlim_data*cosd(angs(ia));
        yy = vlim_data*sind(angs(ia));
        if abs(angs(ia)) == 90
          plot(hca,xx,yy,'color',color_grid,'linestyle','--','linewidth',2)
        else
          plot(hca,xx,yy,'color',color_grid,'linestyle','-')

        end
        text_margin = 1.9;
        if angs(ia) == 180
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,{sprintf('%4.f^o',angs(ia))}, ...
            'color',[0 0 0]+0.0,'VerticalAlignment','bottom','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        elseif angs(ia) == -180
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,{sprintf('%g^o',angs(ia))}, ...
            'color',[0 0 0]+0.0,'VerticalAlignment','top','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        else
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,sprintf('%4.f^o',angs(ia)), ...
            'color',[0 0 0]+0.0,'VerticalAlignment','middle','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        end

      end
      
      hold(hca,'off')
    end
    if 0 % add elim circular markings
      hold(hca,'on')      
      ang = linspace(0,360,361);
      color_grid = ([0.9 0.9 0.9]-0.3)*0;
      plot(hca,vlim_data(1)*cosd(ang), vlim_data(1)*sind(ang),'color',color_grid,'linestyle','--')
      plot(hca,vlim_data(2)*cosd(ang), vlim_data(2)*sind(ang),'color',color_grid,'linestyle','--')
      hold(hca,'off')
    end
    if 0 % add circular grid
      %%
      hold(hca,'on')      
      ang = linspace(0,360,361);
      r = 2500;
      color_grid = ([0.9 0.9 0.9]-0.3)*0;
      x0 = vL_xline-vL_xline_use;
      x0 = 0;
      y0 = 0;
      %plot(hca,x0 + r*cosd(ang), y0 + r*sind(ang),'color',color_grid)
      plot(hca,x0 + 2500*cosd(ang), y0 + 2500*sind(ang),'color',color_grid)
      %plot(hca,x0 + 2000*cosd(ang), y0 + 2000*sind(ang),'color',color_grid)
      plot(hca,x0,0,'kx',x0,0,'ko','linewidth',2)
      angs = -180:30:180;
      for ia = 1:numel(angs)
        xx = x0 + [0 r*cosd(angs(ia))];
        yy = y0 + [0 r*sind(angs(ia))];
        if abs(angs(ia)) == 90
          plot(hca,xx,yy,'color',color_grid,'linestyle','--','linewidth',2)
        else
          plot(hca,xx,yy,'color',color_grid,'linestyle','-')

        end
        text_margin = 1.15;
        if angs(ia) == 180
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,{sprintf('%4.f^o',angs(ia))}, ...
            'color','k','VerticalAlignment','bottom','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        elseif angs(ia) == -180
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,{sprintf('%g^o',angs(ia))}, ...
            'color','k','VerticalAlignment','top','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        else
          ht = text(hca,xx(2)*text_margin,yy(2)*text_margin,sprintf('%4.f^o',angs(ia)), ...
            'color','k','VerticalAlignment','middle','HorizontalAlignment','center', ...
            'backgroundcolor','none');
        end

      end
      
      hold(hca,'off')
    end
  
  end


h1(end).XTickLabelRotation = 0;
%colormap('parula')
cmap = irf_colormap('magma');
colormap(cmap)
%colormap(flipdim(irf_colormap('Spectral'),1))
%irf_colormap('vik')
%colormap(pic_colors('candy4'))
%colormap(flipdim(pic_colors('thermal'),1))


irf_zoom(h1,'x',tint_plot)
h1(1).YLim = [-799 399];
%h1(2).CLim = [1 4];

irf_legend(h2,'(a)',[-0.13 1.05],'color','k','fontsize',fontsize+2)
irf_legend(h1(1),'(b)',[-0.13 1.01],'color','k','fontsize',fontsize+2)
irf_legend(h1(2),'(c)',[-0.13 1.0],'color','k','fontsize',fontsize+2)

%h(1).Title.String = sprintf('L = [%.2f,%.2f,%.2f]; M = [%.2f,%.2f,%.2f]; N = [%.2f,%.2f,%.2f];',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3));
h = findobj(gcf,'type','axes'); h = h(end:-1:1);
c_eval('h(?).FontSize = fontsize;',1:numel(h))

%c_eval('hmark(?) = irf_pl_mark(h1(?),time_xline,''k'',''linestyle'',''-.'',''linewidth'',2); ;',1:numel(h1))
drawnow
c_eval('hmark(?) = irf_pl_mark(h1(?),time_xline,''k'',''linestyle'','':'',''linewidth'',4); ;',1:numel(h1))

h = findobj(gcf,'type','axes'); h = h(end:-1:1);
hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
c_eval('h(?).LineWidth = 1.5;',1:numel(h))

irf_legend(0,{sprintf('L=[%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3)),sprintf('nMean=[%g,%g,%g,%g], nThresh = %g',nMean(1),nMean(2),nMean(3),nMean(4),nThresh)},[0.05 1])

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

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:53.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];

times_utc = [...%'2017-07-11T22:33:25.000Z';...
             '2017-07-11T22:33:53.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             '2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:04.000Z';...
             '2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:15.940Z'];


times_utc = ['2017-07-11T22:33:45.000Z';...
             '2017-07-11T22:33:54.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             %'2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:00.502Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:05.500Z';...
             %'2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:10.000Z'];

times_utc = [...%'2017-07-11T22:33:45.000Z';...
             '2017-07-11T22:33:54.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             %'2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:00.502Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:05.500Z';...
             %'2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:10.000Z';...
             '2017-07-11T22:34:20.000Z'];
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

vdf_legs = {'I','II','III','IV','V','VI'};
vdf_legs = {'A','B','C','D','E','F'};
% Plot
nRows = 4;
nCols = size(times_utc,1);
[h1,h] = initialize_combined_plot('topbottom',2,nRows,nCols,0.3,'vertical');

vL_Xline = 1*vL_xline;
fontsize = 9;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854
tint_zoom = irf.tint('2017-07-11T22:33:34.00Z/2017-07-11T22:34:30.00Z'); %20151112071854

hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
%c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
hca.YLabel.String = 'v_i (km/s)';
hca.ColorOrder = mms_colors('xyz');
%irf_legend(hca,{'L','M','N'},[0.98 0.98]);
%irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
%irf_legend(hca,{sprintf(['v_L'],vL_Xline),'   v_M','   v_N'},[0.01 0.12],'fontsize',fontsize);
irf_legend(hca,{sprintf(['L'],vL_Xline),'  M','   N'},[0.08 0.10],'fontsize',fontsize);
%irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';

if 1
hca = h1(isub); isub = isub + 1;
hca.ColorOrder = mms_colors('xyz');
%c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
c_eval('irf_plot(hca,{mvaB?.x-0*vL_Xline,mvaB?.y,mvaB?.z},''comp'')',ic)
hca.YLabel.String = 'B (nT)';
hca.ColorOrder = mms_colors('xyz');
%irf_legend(hca,{'L','M','N'},[0.98 0.98]);
%irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'}',[1.01 0.98]);
%irf_legend(hca,{sprintf(['v_L'],vL_Xline),'   v_M','   v_N'},[0.01 0.12],'fontsize',fontsize);
irf_legend(hca,{sprintf(['L'],vL_Xline),'  M','   N'},[0.08 0.10],'fontsize',fontsize);
%irf_legend(hca,{sprintf('L=[%.2f,%.2f,%.2f]',L(1),L(2),L(3)),sprintf('M=[%.2f,%.2f,%.2f]',M(1),M(2),M(3)),sprintf('N=[%.2f,%.2f,%.2f]',N(1),N(2),N(3))}',[1.01 0.98]);
hca.YLabel.Interpreter = 'tex';
end

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


isub = 1;
% nSmooth = 1; % specified further doen


irf_zoom(h1,'x',tint_figure_zoom_incl_sep+[+8 -8])
irf_zoom(h1,'y')

pdist_all = PD_clean;
pdist_all.data(:,:,:,[ ]) = 0;
time = time_vdf;

fontsize_leg = fontsize;
%fontsize = 11;


vint_L = [-Inf -170];
vint_L = [-Inf Inf];
vint_M = [-Inf Inf];
vint_N = [-Inf Inf];

% time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
% time = time + dt;
times = EpochTT(times_utc);
for it = 1:times.length%(1)
  time = times(it);
  time = time+12*+-0;
  time = time+0;
  
  
  pdist = pdist_all.tlim(time+nMean(1)*0.5*0.151*[-1 1]);
  tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*1*[-1 1];
  elow = mean(tsElow.tlim(tint_dist).data,1);
 
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
  
  B_scale = 1/1000; % 2 nT per 1000 km/s
  scaxis_scale = 2000;
  vscale = 1e3;
  vlim = 2500/vscale;
  B_scale = B_scale*vscale;
  nSmooth = 0;
  nContours = 0;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours, '10^3 km/s')
    hca.Position = position;
    axis(hca,'square')
    
    if vL_Xline == 0
      hca.XLabel.String = 'v_L (10^3 km/s)';
    else
      hca.XLabel.String = sprintf(['v_L-(%g) (10^3 km/s)'],vL_Xline);
      hca.XLabel.String = 'v_L-v_{L}^{Xline} (10^3 km/s)';
    end
    hca.YLabel.String = 'v_M (10^3 km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
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
      b = b*B_inplane/B_scale;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',2)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)/vscale,mean(vExB.y.data,1)/vscale,'ok','MarkerFaceColor','w','markersize',5);
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

    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end

  %irf_legend(hca,sprintf('VDF %s',vdf_legs{it}),[0.02 0.98],'fontweight','bold','fontsize',fontsize-0,'color','k')
  hca.Title.String = sprintf('VDF %s',vdf_legs{it});
  if 1 % f(L,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[L_vi],[N_vi]);
    vdf = pdist.reduce('2D',[Ldsl],[Ndsl],'vint',vint_M);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours, '10^3 km/s')
    hca.Position = position;
    axis(hca,'square')

    hca.YLabel.String = 'v_N (10^3 km/s)';
    if vL_Xline == 0
      hca.XLabel.String = 'v_L (10^3 km/s)';
    else
      hca.XLabel.String = sprintf(['v_L-(%g) (10^3 km/s)'],vL_Xline);
      hca.XLabel.String = 'v_L-v_{L}^{Xline} (10^3 km/s)';
    end

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
      B_std_inplane = std(B__.data(:,[1 3]),1);
      B_inplane = sqrt(sum(B_([1 3]).^2));
      b = b*B_inplane/B_scale;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(3),2*b(1),2*b(3),0,'k','linewidth',2)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_([1 3]);
    end  
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)/vscale,mean(vExB.z.data,1)/vscale,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end
    
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end
  
  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl],'vint',vint_L);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours, '10^3 km/s')
    hca.Position = position;
    axis(hca,'square')
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
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
      B_std_inplane = std(B__.data(:,[2 3]),1);
      B_inplane = sqrt(sum(B_([2 3]).^2));
      b = b*B_inplane/B_scale;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(2),-b(3),2*b(2),2*b(3),0,'k','linewidth',2)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_([2 3]);
    end  
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.y.data,1)/vscale,mean(vExB.z.data,1)/vscale,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(2)*scaxis_scale,-scaxis(3)*scaxis_scale,2*scaxis(2)*scaxis_scale,2*scaxis(3)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end

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
    if it == 6
      data1 = data1/2;
      data2 = data2/2;
    end
    plot(hca,v_center/vscale,data1,v_center/vscale,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    if it == 6
      irf_legend(hca,'(1/2)f_i(v_M)',[0.98 0.98],'color','k','fontsize',fontsize_leg)
    end
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    if it == 1
    %irf_legend(hca,...
    %  {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))},...
    %  [0.98 0.98],'fontsize',fontsize_leg)
    irf_legend(hca,...
      {'v_N > 0','v_N < 0'}',...
      [0.98 0.98],'fontsize',fontsize_leg)
    end
    vlim = 1500/vscale;
    hca.XLim = vlim*[-1 1];  
  end

  %times_exact{1} = vdf.time;
  times_exact{1} = pdist.time;

  c_eval('hmark(?) = irf_pl_mark(h1(!),[times_exact{?}(1).epochUnix times_exact{?}(end).epochUnix] + 0.5*0.03*[-1 1],[0.5 0.5 0.5]+0.2);',1,1:numel(h1))
  %cn.print(sprintf('torbert_fi_ref_dt_%g',dt))
end


c_eval('h1(?).FontSize = fontsize;',1:numel(h1))
c_eval('h(?).FontSize = fontsize;',1:numel(h))

hq = findobj(gcf,'type','quiver'); hq = hq(end:-1:1);
c_eval('hq(?).MaxHeadSize=0.5;',1:numel(hq))

for ip = 3:nRows:numel(h)
  vmax = -2000;
  patch(h(ip),-[-3000 vmax vmax -3000]/vscale,[-3000 -3000 -vmin -vmin]/vscale,colors(2,:),...
    'facealpha',0.1,'edgecolor',colors(2,:));
  %h(ip).Children = circshift(h(ip).Children,1);

  patch(h(ip),-[-3000 vmax vmax -3000]/vscale,[3000 3000 vmin vmin]/vscale,colors(1,:),...
    'facealpha',0.1,'edgecolor',colors(1,:));
  h(ip).Children = circshift(h(ip).Children,4);
end

irf_legend(h1(1),'(a)',[0.02 0.98],'color','k','fontsize',fontsize)
irf_legend(h1(2),'(b)',[0.02 0.98],'color','k','fontsize',fontsize)
for it = 1:times.length  
    %irf_legend(h(1+(it-1)*4),sprintf('(c%g)',it),[0.02 0.98],'k')
    %irf_legend(h(2+(it-1)*4),sprintf('(d%g)',it),[0.02 0.98],'k')
    %irf_legend(h(3+(it-1)*4),sprintf('(e%g)',it),[0.02 0.98],'k')
    %irf_legend(h(4+(it-1)*4),sprintf('(f%g)',it),[0.02 0.98],'k')

    irf_legend(h(1+(it-1)*4),sprintf('(c%s)',vdf_legs{it}),[0.02 0.98],'color','k','fontsize',fontsize)
    irf_legend(h(2+(it-1)*4),sprintf('(d%s)',vdf_legs{it}),[0.02 0.98],'color','k','fontsize',fontsize)
    irf_legend(h(3+(it-1)*4),sprintf('(e%s)',vdf_legs{it}),[0.02 0.98],'color','k','fontsize',fontsize)
    irf_legend(h(4+(it-1)*4),sprintf('(f%s)',vdf_legs{it}),[0.02 0.98],'color','k','fontsize',fontsize)
end

if 1 % add B scale quiver
  %%
  hold(h(17),'on')   
  quiver(h(17),1,2,1,0,0,'k','linewidth',2,'MaxHeadSize',0.5)
  %irf_legend(h(17),sprintf('B=%g nT',1*B_scale),[0.5 0.98],'k')
  text(h(17),0.9,2,sprintf('%g nT',1*B_scale),'color','k','HorizontalAlignment','right','fontsize',fontsize_leg)
  hold(h(17),'on')
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

%h(4).YLim = [0 0.035];
h(4).YLim = [0 0.025];

h(4).YTickLabelMode = 'auto';

%compact_panels(h,0.04,0.01)
drawnow

hb = findobj(gcf,'type','colorbar'); hb = hb(end:-1:1);
%delete(hb(1:end-3))
delete(hb)
hcb = colorbar(h(numel(h)-nRows+1),'location','northoutside');
hcb.Label.String = 'PSD (s^2m^{-4})';
hcb.FontSize = fontsize-0;
%hcb.Position(2)=0.741;

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
c_eval('h(?).XTick = (-2000:1000:2000)/vscale;',1:numel(h))
c_eval('h(?).YTick = (-2000:1000:2000)/vscale;',setdiff(1:numel(h),4:4:numel(h)))
c_eval('h(?).XTickLabelRotation = 0;',1:numel(h))


c_eval('h(?).Layer = ''top'';',1:numel(h))
c_eval('h(?).GridLineWidth = 0.5;',1:numel(h))
c_eval('h(?).LineWidth = 0.5;',1:numel(h))
c_eval('h1(?).LineWidth = 0.5;',1:numel(h1))

c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

%c_eval('h1(?).Position(3) = 0.3;',1:numel(h1))
c_eval('h1(?).XTickLabelRotation = 0;',1:numel(h1))
c_eval('h1(?).XLabel = [];',1:(numel(h1)-1))

h1(2).YLim(2) = 4;

dy = 0.0848;
%dy = 0.05;

%for ih1 = 1:numel(h1)
%  h1(ih1).Position(2) = h1(ih1).Position(2) - 0.05;
%  h1(ih1).Position(4) = h1(ih1).Position(4)*6;
%end
h1(1).Position = [0.300    0.7730+dy    0.4000    dy];
h1(2).Position = [0.300    0.7730    0.4000    dy];
%h1(1).Position = [0.300    0.8730    0.4000    dy];
%h1(2).Position = [0.1700    0.7672    0.3000    0.0848];
%h1(2).Position = [0.300    0.8730-dy    0.4000    dy];
%h1(3).Position = [0.300    0.8730-dy*2    0.4000    dy];

%hca.CLim = [-9.5 -7.5];
h(1).CLim = [-10 -7.5];
h(1).CLim = [-10 -7.3];

h(4).YTickMode = 'auto';
h(4).YLabel.String = 'f_i(v_M) (s/m^4)';
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))

c_eval('h1(?).FontSize = fontsize;',1:numel(h1))

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.;',1:numel(hl))
hall = findobj(gcf,'type','axes'); hall = hall(end:-1:1);
c_eval('hall(?).LineWidth = 1.;',1:numel(hall))

irf_legend(0,{sprintf('L=[%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3)),sprintf('nMean=[%g,%g,%g,%g], nThresh = %g',nMean(1),nMean(2),nMean(3),nMean(4),nThresh)},[0.05 1])

%% Single time distribution f(vL,vM)
times_utc = [...%'2017-07-11T22:33:45.000Z';...
             '2017-07-11T22:33:54.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             %'2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:00.502Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:05.500Z';...
             %'2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:10.000Z';...
             '2017-07-11T22:34:20.000Z'];

h = subplot(1,1,1);
isub = 1;
it = 4;
 time = times(it);
  
  
  pdist = pdist_all.tlim(time+nMean(1)*0.5*0.151*[-1 1]);
  tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*1*[-1 1];
  elow = mean(tsElow.tlim(tint_dist).data,1);
 

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  
  B_scale = 1/1000; % 2 nT per 1000 km/s
  scaxis_scale = 2000;
  vscale = 1e3;
  vlim = 2500/vscale;
  B_scale = B_scale*vscale;
  nSmooth = 0;
  nContours = 0;
  if 1 % f(L,M)
    hca = h(isub); isub = isub + 1;
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours, '10^3 km/s')
    
    axis(hca,'square')
    
    if vL_Xline == 0
      hca.XLabel.String = 'v_L (10^3 km/s)';
    else
      hca.XLabel.String = sprintf(['v_L-(%g) (10^3 km/s)'],vL_Xline);
      hca.XLabel.String = 'v_L-v_{L}^{Xline} (10^3 km/s)';
    end
    hca.YLabel.String = 'v_M (10^3 km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
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
      b = b*B_inplane/B_scale;
      %quiver(-b(2),-b(1),2*b(2),2*b(1),0,'k','linewidth',1)
      quiver(-b(1),-b(2),2*b(1),2*b(2),0,'k','linewidth',2)
        hold(hca,'off')
        hca.XLim = xlim;
        hca.YLim = ylim;
        B_ = B_(1:2);
    end     
    if 1 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)/vscale,mean(vExB.y.data,1)/vscale,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end
    if 0 % sc spin axis
      hold(hca,'on')      
      quiver(-scaxis(1)*scaxis_scale,-scaxis(2)*scaxis_scale,2*scaxis(1)*scaxis_scale,2*scaxis(2)*scaxis_scale,0,'k','linewidth',1)
      hold(hca,'off') 
    end

    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
  end

  hca.FontSize = 12;
  hca.LineWidth = 1.5;
  h(1).CLim = [-10 -7.3];
  colormap(irf_colormap('magma'))
  %irf_legend(hca,sprintf('VDF %s',vdf_legs{it}),[0.02 0.98],'fontweight','bold','fontsize',fontsize-0,'color','k')
  hca.Title.String = sprintf('VDF %s',vdf_legs{it});


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
time_pdist = EpochTT(times_utc(3,:));

elim = [000 Inf];
%elim = [200 Inf];
nMovMean = 5;
c_eval('pdist_all = PD_clean;',ic)

fontsize_leg = 10;
fontsize = 10;

time = time_pdist;

pdist_tint = time + 0.5*0.150*nMean(1)*[-1 1];
pdist = pdist_all.tlim(pdist_tint);
tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*[-1 1]; % this just reshifts to the actual FPI timeline


pdist = pdist;
times = pdist.time;

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

%h = setup_subplots(2,2);
clear h;
h(1) = subplot(2,4,[1 2 5 6]);
h(2) = subplot(2,4,[3 4]);
h(3) = subplot(2,4,[7 8]);
isub = 1;

  if 1 % f(M,N)
    hca = h(isub); isub = isub + 1;
    %vdf = pdist_nobg.reduce('2D',[M_vi],[N_vi]);
    vdf_MN = pdist.reduce('2D',[Mdsl],[Ndsl]);
    vdf = vdf_MN;
    position = hca.Position;
    vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours, '10^3 km/s')
    hca.Position = position;
    hcb = findobj(gcf,'type','colorbar');
    hcb.Location = 'northoutside';
    axis(hca,'square')
    hca.XLabel.String = 'v_M (10^3 km/s)';
    hca.YLabel.String = 'v_N (10^3 km/s)';
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
    vlim = 2500*1e-3;
    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    if 0
      hold(hca,'on')
      patch(hca,[2000 3000 3000 2000],[-3000 -3000 -0 -0],colors(2,:),...
        'facealpha',0.1,'edgecolor',colors(2,:));
      %h(ip).Children = circshift(h(ip).Children,1);
    
      patch(hca,[2000 3000 3000 2000],[3000 3000 00 00],colors(1,:),...
        'facealpha',0.1,'edgecolor',colors(1,:));
      hca.Children = circshift(hca.Children,1);
      hold(hca,'off')
    end
    if 1
      hold(hca,'on')
    
      patch(hca,[-3000 0000 0000 -3000],[3000 3000 00 00],colors(1,:),...
        'facealpha',0.3,'edgecolor',colors(1,:)*0);

      patch(hca,[0000 3000 3000 0000],[3000 3000 00 00],colors(2,:),...
        'facealpha',0.3,'edgecolor',colors(2,:)*0);

      patch(hca,[-3000 0000 0000 -3000],[-3000 -3000 -0 -0],colors(3,:),...
        'facealpha',0.3,'edgecolor',colors(3,:)*0);
      patch(hca,[0000 3000 3000 0000],[-3000 -3000 -0 -0],colors(4,:),...
        'facealpha',0.3,'edgecolor',colors(4,:)*0);
    


      hca.Children = circshift(hca.Children,1);
      hold(hca,'off')
    end
    hca.CLim = [-10 -7.5];
    %colormap(pic_colors('candy6'))
    colormap(irf_colormap('magma'))
  end
      
  vscale = 1e3;
  if 1 % 1D at a certain vn range, patches
    hca = h(isub); isub = isub + 1;

    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [0 inf];
    vint2 = [-inf -0];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);


    data = vdf.data;    
    data = sum(data,3)*dv*1e3;

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2)); % vN range
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;
    

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2)); % vN range
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;    
    
   % vdf1 = PDist(times,data1,'1Dcart',vdf.depend{1});
   % vdf2 = PDist(times,data2,'1Dcart',vdf.depend{1});

    data1 = mean(data1,1); % vN > 0
    data2 = mean(data2,1); % vN < 0

    hca.ColorOrder = colors;

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data1_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data1_y = [data1(idx) data1(idx)*0 0];
    
    idx = find(v_center>0); % vM range
    data2_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data2_y = [data1(idx) data1(idx)*0 0];

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data3_x = [v_center(idx) v_center(idx(end:-1:1))];
    data3_y = [data2(idx) data2(idx)*0];
    
    idx = find(v_center>0);
    data4_x = [v_center(idx) v_center(idx(end:-1:1))];
    data4_y = [data2(idx) data2(idx)*0];
    
    plot(hca,v_center/vscale,mean(data1,1),'k')
    hold(hca,'on')
    
    patch(hca,data1_x/vscale,data1_y,colors(1,:),...
        'facealpha',0.3,'edgecolor',colors(1,:)*0);
    patch(hca,data2_x/vscale,data2_y,colors(2,:),...
        'facealpha',0.3,'edgecolor',colors(2,:)*0);
    
    %patch(hca,data3_x,data3_y,colors(3,:),...
    %    'facealpha',0.3,'edgecolor',colors(3,:)*0);
    %patch(hca,data4_x,data4_y,colors(4,:),...
    %    'facealpha',0.3,'edgecolor',colors(4,:)*0);
    hold(hca,'off')

    %plot(hca,v_center,data1,v_center,data2,'linewidth',2)    
    hca.XLabel.String = 'v_M (10^3 km/s)';
    
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
      {sprintf('v_N > %g km/s',0)},...
      [0.98 0.98],'fontsize',fontsize_leg,'color','k')
  
    vlim = 1499;
    hca.XLim = vlim*[-1 1]/vscale;  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = [-2000:500:2000]/vscale;
    hca.XTickLabelRotation = 0;
    hca.Box = 'on';
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
  
  if 1 % 1D at a certain vn range, patches 
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [0 inf];
    vint2 = [-inf -0];
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
    
   % vdf1 = PDist(times,data1,'1Dcart',vdf.depend{1});
   % vdf2 = PDist(times,data2,'1Dcart',vdf.depend{1});

    data1 = mean(data1,1);
    data2 = mean(data2,1);

    hca.ColorOrder = colors;

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data1_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data1_y = [data1(idx) data1(idx)*0 0];
    
    idx = find(v_center>0); % vM range
    data2_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data2_y = [data1(idx) data1(idx)*0 0];

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data3_x = [v_center(idx) v_center(idx(end:-1:1))];
    data3_y = [data2(idx) data2(idx)*0];
    
    idx = find(v_center>0);
    data4_x = [v_center(idx) v_center(idx(end:-1:1))];
    data4_y = [data2(idx) data2(idx)*0];
    

    plot(hca,v_center/vscale,mean(data2,1),'k')
    hold(hca,'on')

    patch(hca,data3_x/vscale,data3_y,colors(3,:),...
        'facealpha',0.3,'edgecolor',colors(3,:)*0);    
    patch(hca,data4_x/vscale,data4_y,colors(4,:),...
        'facealpha',0.3,'edgecolor',colors(4,:)*0);
    hold(hca,'off')
    
    %patch(hca,data3_x,data3_y,colors(3,:),...
    %    'facealpha',0.3,'edgecolor',colors(3,:)*0);
    %patch(hca,data4_x,data4_y,colors(4,:),...
    %    'facealpha',0.3,'edgecolor',colors(4,:)*0);
    hold(hca,'off')

    %plot(hca,v_center,data1,v_center,data2,'linewidth',2)    
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
      {sprintf('v_N < %g km/s',0)},...
      [0.98 0.98],'fontsize',fontsize_leg,'color','k')
  
    vlim = 1499;
    hca.XLim = vlim*[-1 1]/vscale;  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = [-2000:500:2000]/vscale;
    hca.XTickLabelRotation = 0;
    hca.Box = 'on';
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
  

  if 0 % 1D at a certain vn range    
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
    
   % vdf1 = PDist(times,data1,'1Dcart',vdf.depend{1});
   % vdf2 = PDist(times,data2,'1Dcart',vdf.depend{1});

    data1 = mean(data1,1);
    data2 = mean(data2,1);

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
      {sprintf('v_N > %g km/s',0),sprintf('v_N < %g km/s',0)}',...
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
  

h(1).Position(3) = 0.32;
h(1).Position(1) = 0.17;
for ip = 1:numel(h)
  hca = h(ip);
  %axis(hca,'square')
%  hca.Position(3) = 0.25;
end
c_eval('h(?).FontSize = fontsize;',1:numel(h))
c_eval('h(?).LineWidth = 1.5;',1:numel(h))
%h(1).Title.String = sprintf('%s - %s',tint_dist(1).utc('HH:MM:SS.mmm'),tint_dist(2).utc('SS.mmm'));
%compact_panels(h,0.01,0.2)

irf_legend(0,{sprintf('L=[%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3)),sprintf('nMean=[%g,%g,%g,%g], nThresh = %g',nMean(1),nMean(2),nMean(3),nMean(4),nThresh)},[0.05 1])

%%
%%
%%
%% Figure, overplot different timesas patches
fontsize = 12;

h = subplot(1,1,1);
isub = 1;

colors = mms_colors('matlab');
roman_numerals = {'I','II','III','IV','V','VI'};

hca = h(isub); isub = isub + 1;
times = EpochTT(times_utc);
leg_strs = {};
leg_colors = [];
iCount = 0;
for it = [ 3 4 5 6]%1:times.length%(1)
  iCount = iCount + 1;
  time = times(it);
  time = time+-12*0+-0;
  
  
  pdist = pdist_all.tlim(time+nMean(1)*0.5*0.151*[-1 1]);
  tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*1*[-1 1];
  
  leg_strs{iCount} = sprintf('VDF %s',roman_numerals{it});
  leg_colors(iCount,1:3) = colors(it,:);
  
 
  %c_eval('hmark = irf_pl_mark(h1,tint_dist,[0.5 0.5 0.5]);',1:numel(h1))

  t_dist_center = pdist.time.start + (pdist.time.stop - pdist.time.start)/2;
  c_eval('Tdsl = [tsLdsl?.resample(t_dist_center).data; tsMdsl?.resample(t_dist_center).data; tsNdsl?.resample(t_dist_center).data];',ic)
  Ldsl = Tdsl(1,:);
  Mdsl = Tdsl(2,:);
  Ndsl = Tdsl(3,:);
  
  c_eval('B = mvaB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)  
  c_eval('E = mvaE?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic) 
  c_eval('vExB = mvaVExB?.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]);',ic)   
  c_eval('scaxis = mean(tsSCaxis?_lmn.tlim(pdist.time([1 end]) + 0.5*0.15*[-1 1]).data,1);',ic) 
  

  holdon = 0;
  
  scaxis_scale = 2000;
  nSmooth = 0;
  nContours = 0;
  nContours = -9; vlim = 2500;
  nContours = -8.2; vlim = 1700;
  if numel(nContours) == 1; nContours = [nContours nContours]; end
  
  if 1 % f(L,M)    
    vdf = pdist.reduce('2D',Ldsl,Mdsl,'vint',vint_N,'nMC',1000);
    vdf.depend{1} = vdf.depend{1} - vL_Xline;
    vdf.ancillary.vx_edges = vdf.ancillary.vx_edges - vL_Xline;
    position = hca.Position;
    data = squeeze(mean(vdf.data,1));
    data = smooth2(data,nSmooth);
    %data = smooth2(data,nSmooth);
    data = log10(data);
    v1_ = squeeze(mean(vdf.depend{1},1));
    v2_ = squeeze(mean(vdf.depend{1},1));
    htmp = contour(hca,v1_,v2_,data',nContours,'edgecolor',colors(it,:),'LineWidth',2);
    %vdf.plot_plane(hca,'smooth',nSmooth,'contour',nContours)
    hca.Position = position;
    axis(hca,'square')
    
    if vL_Xline == 0
      hca.XLabel.String = 'v_L (km/s)';
    else
      hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
      hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    end
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    
    if 0 % plot ExB
      hold(hca,'on')
      hbulk = plot(hca,mean(vExB.x.data,1)*1e0,mean(vExB.y.data,1)*1e0,'ok','MarkerFaceColor','w','markersize',5);
      hold(hca,'off')    
    end

    hca.XLim = vlim*[-1 1];
    hca.YLim = vlim*[-1 1];
    if not(holdon); holdon = 1; hold(hca,'on'); end

    contourdata = htmp;
    diffhtmp = diff(contourdata,1,2);

    dl = sqrt(sum(diffhtmp(1,:).^2 + diffhtmp(2,:).^2,1));
    jumps = find(dl > 200);
    contourdata(:,jumps) = NaN;

    polys = splitContours(htmp);

    %
    minArea = 0e3; % threshold
    
    for i = 1:numel(polys)
        x = polys{i}.x;
        y = polys{i}.y;
    
        A = polyarea(x, y);
    
        if A > minArea   % keep only “big” areas
            %patch(x, y, 'b', 'FaceAlpha', 0.4);
            hp = patch(hca,x,y,x*0+1,'facealpha',0.2,'facecolor',colors(it,:),'edgecolor','none','LineWidth',2);
        end
    end

    %hp = patch(hca,contourdata(1,:),contourdata(2,:),contourdata(1,:)*0+1,'facealpha',0.2,'facecolor',colors(it,:),'edgecolor','none','LineWidth',2);
  end
end
hold(hca,'off');


hca.ColorOrder = leg_colors;
irf_legend(hca,leg_strs',[0.02 0.02],'fontsize',fontsize,'fontweight','bold')
irf_legend(hca,sprintf('f(v_L,v_M) = 10^{%g}',nContours(1)),[0.02 0.98],'fontsize',fontsize,'fontweight','bold','color',[0 0 0])

hl = findobj(gcf,'type','line');
c_eval('hl(?).LineWidth = 1.5;',1:numel(hl))
hall = findobj(gcf,'type','axes'); hall = hall(end:-1:1);
c_eval('hall(?).LineWidth = 1.5;',1:numel(hall))
c_eval('h(?).XGrid = ''on''; h(?).YGrid = ''on'';',1:numel(h))
c_eval('h(?).FontSize = fontsize;',1:numel(h))
  

%% Make movie of 2d VDFs, to get a better feeling of the turning.

fT = 40;
tint = time_xline_ion + 0.5*fT*[-1 1];

nMovMean = 7;
%c_eval('pdist_all = iPDist?.movmean(nMovMean,''removeonecounts'',iPDist?_counts).tlim(tint);',ic)
c_eval('pdist_all = PD_clean;',ic)

nRows = 2;
nCols = 2;
[h1,h] = initialize_combined_plot('topbottom',1,nRows,nCols,0.3,'vertical');

compact_panels(h,0.1,0.1)

vL_Xline = 1*-170;

isub = 1;
tint_zoom = irf.tint('2017-07-11T22:33:24.00Z/2017-07-11T22:34:40.00Z'); %20151112071854

if 1 % vi
  hca = h1(isub); isub = isub + 1;
  hca.ColorOrder = mms_colors('xyz');
  c_eval('irf_plot(hca,{mvaVi?.x-0*vL_Xline,mvaVi?.y,mvaVi?.z},''comp'')',ic)
  hca.YLabel.String = 'v_i (km/s)';
  hca.ColorOrder = mms_colors('xyz');
  %irf_legend(hca,{sprintf(['v_L-(%g km/s)'],vL_Xline),'v_M','v_N'},[.98 0.05]);
  irf_legend(hca,{'v_L','v_M','v_N'},[.98 0.05]);
  hca.YLabel.Interpreter = 'tex';
end

isub = 1;
% nSmooth = 1; % specified further doen


irf_zoom(h1,'x',tint_figure_zoom_incl_sep+[+8 -8])
%irf_zoom(h1,'x',tint_figure_zoom_incl_sep)
irf_zoom(h1,'x',tint + [-3 3])

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
vidfile.FrameRate = 20;
open(vidfile);
     
clear F
times = pdist_all.tlim(tint).time;
for it = 1:1:times.length %1:nMovMean:times.length
  time = times(it);
  
  pdist = pdist_all.tlim(time+nMean(1)*0.5*0.151*[-1 1]);
  tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*[-1 1];
  pdist = pdist;
  
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
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
    hca.YLabel.String = 'v_M (km/s)';
    vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
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

    if 0 % plot eigenvectors
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
    hca.XLabel.String = 'v_L-v_{L}^{Xline} (km/s)';
    %hca.XLabel.String = sprintf(['v_L-(%g) (km/s)'],vL_Xline);
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
      {sprintf('%g km/s < v_N',vint1(1)),sprintf('v_N < %g km/s',vint2(2))}',...
      [0.98 0.98],'fontsize',fontsize_leg)
  
    vlim = 1500;
    hca.XLim = vlim*[-1 1];  
    hca.YLim = [0 0.025]
    hca.XGrid = 'on';
    hca.YGrid = 'on';
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
    if 0%it == 6
      data1 = data1/2;
      data2 = data2/2;
    end
    plot(hca,v_center/vscale,data1,v_center/vscale,data2)    
    hca.XLabel.String = 'v_M (km/s)';
    
    if it == 6
      irf_legend(hca,'(1/2)f_i(v_M)',[0.98 0.98],'k')
    end
    axis(hca,'square')
    hca.YLabel.String = 'f_i(v_M) (s/m^4)';
    %vmin = sqrt(2*units.eV*min(elow)/units.mp)*1e-3;
    %irf_legend(hca,sprintf('v>%.0f km/s',vmin),[0.02 0.98],'color','k','fontsize',10)
    %irf_legend(hca,sprintf('%g',ip),[0.02 0.98],'color','k','fontsize',10)
  
    %E_legs = arrayfun(@(x) sprintf('%.0f eV',x),elows.data,'UniformOutput',false);
    %E_legs{1} = {['E > ' E_legs{1}]};
    hca.ColorOrder = colors;
    %irf_legend(hca,E_legs,[0.98 0.98],'fontsize',10)
    if it == 1
    %irf_legend(hca,...
    %  {sprintf('v_N > %g km/s',vint1(1)),sprintf('v_N < %g km/s',vint2(2))},...
    %  [0.98 0.98],'fontsize',fontsize_leg)
    irf_legend(hca,...
      {'v_N > 0','v_N < 0'}',...
      [0.98 0.98],'fontsize',fontsize_leg)
    end
    vlim = 1500/vscale;
    hca.XLim = vlim*[-1 1];  
  end


  %irf_legend(h(1),{sprintf('N = %g',nMovMean)},[0.02 1.01])
  hlinks = linkprop(h(1:3),{'CLim'});
  hlinks.Targets(1).CLim = [-10 -7.5];
  c_eval('h(?).FontSize = 13;',1:numel(h))
  %colormap(pic_colors('candy6'))
  colormap(irf_colormap('magma'))
  h(1).CLim = [-10 -7.3];

  set(gcf,'color','white');

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
[h1,h] = initialize_combined_plot('topbottom',2,nRows,nCols,0.15,'horizontal');

for ih1 = 1:numel(h1)
  h1(ih1).Position(2) = h1(ih1).Position(2) - 0.05;
  h1(ih1).Position(4) = h1(ih1).Position(4)*6;
end
%h1.Position(2) = h1.Position(2) - 0.05;
%h1.Position(4) = h1.Position(4)*6;

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
time_pdist = EpochTT(times_utc(3,:));

elim = [000 Inf];
%elim = [200 Inf];
nMovMean = 5;
c_eval('pdist_all = PD_clean;',ic)

fontsize_leg = 9;
fontsize = 10;

time = time_pdist;

pdist_tint = time + 0.5*0.150*nMean(1)*[-1 1];
pdist = pdist_all.tlim(pdist_tint);
tint_dist = [pdist.time.start pdist.time.stop] + 0.5*0.150*[-1 1]; % this just reshifts to the actual FPI timeline


pdist = pdist;
times = pdist.time;

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

h = setup_subplots(2,2);
h(1) = subplot(2,4,[1 2 5 6]);
h(2) = subplot(2,4,[3 4]);
h(3) = subplot(2,4,[7 8]);
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
    if 0
      hold(hca,'on')
      patch(hca,[2000 3000 3000 2000],[-3000 -3000 -0 -0],colors(2,:),...
        'facealpha',0.1,'edgecolor',colors(2,:));
      %h(ip).Children = circshift(h(ip).Children,1);
    
      patch(hca,[2000 3000 3000 2000],[3000 3000 00 00],colors(1,:),...
        'facealpha',0.1,'edgecolor',colors(1,:));
      hca.Children = circshift(hca.Children,1);
      hold(hca,'off')
    end
    if 1
      hold(hca,'on')
    
      patch(hca,[-3000 0000 0000 -3000],[3000 3000 00 00],colors(1,:),...
        'facealpha',0.3,'edgecolor',colors(1,:)*0);

      patch(hca,[0000 3000 3000 0000],[3000 3000 00 00],colors(2,:),...
        'facealpha',0.3,'edgecolor',colors(2,:)*0);

      patch(hca,[-3000 0000 0000 -3000],[-3000 -3000 -0 -0],colors(3,:),...
        'facealpha',0.3,'edgecolor',colors(3,:)*0);
      patch(hca,[0000 3000 3000 0000],[-3000 -3000 -0 -0],colors(4,:),...
        'facealpha',0.3,'edgecolor',colors(4,:)*0);
    


      hca.Children = circshift(hca.Children,1);
      hold(hca,'off')
    end
    hca.CLim = [-10 -7.5];
    %colormap(pic_colors('candy6'))
    colormap(irf_colormap('magma'))
  end
      
  if 1 % 1D at a certain vn range, patches
    hca = h(isub); isub = isub + 1;

    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [0 inf];
    vint2 = [-inf -0];
    %vint = [-500 500];
    v_center = vdf.depend{1}(1,:);
    dv = v_center(2)-v_center(1);


    data = vdf.data;    
    data = sum(data,3)*dv*1e3;

    vidx1 = find(all([v_center>vint1(1); v_center<vint1(2)]',2)); % vN range
    data1 = vdf.data;
    data1 = data1(:,:,vidx1);
    data1 = sum(data1,3)*dv*1e3;
    

    vidx2 = find(all([v_center>vint2(1); v_center<vint2(2)]',2)); % vN range
    data2 = vdf.data;
    data2 = data2(:,:,vidx2);
    data2 = sum(data2,3)*dv*1e3;    
    
   % vdf1 = PDist(times,data1,'1Dcart',vdf.depend{1});
   % vdf2 = PDist(times,data2,'1Dcart',vdf.depend{1});

    data1 = mean(data1,1); % vN > 0
    data2 = mean(data2,1); % vN < 0

    hca.ColorOrder = colors;

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data1_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data1_y = [data1(idx) data1(idx)*0 0];
    
    idx = find(v_center>0); % vM range
    data2_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data2_y = [data1(idx) data1(idx)*0 0];

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data3_x = [v_center(idx) v_center(idx(end:-1:1))];
    data3_y = [data2(idx) data2(idx)*0];
    
    idx = find(v_center>0);
    data4_x = [v_center(idx) v_center(idx(end:-1:1))];
    data4_y = [data2(idx) data2(idx)*0];
    
    plot(hca,v_center,mean(data1,1),'k')
    hold(hca,'on')
    
    patch(hca,data1_x,data1_y,colors(1,:),...
        'facealpha',0.3,'edgecolor',colors(1,:)*0);
    patch(hca,data2_x,data2_y,colors(2,:),...
        'facealpha',0.3,'edgecolor',colors(2,:)*0);
    
    %patch(hca,data3_x,data3_y,colors(3,:),...
    %    'facealpha',0.3,'edgecolor',colors(3,:)*0);
    %patch(hca,data4_x,data4_y,colors(4,:),...
    %    'facealpha',0.3,'edgecolor',colors(4,:)*0);
    hold(hca,'off')

    %plot(hca,v_center,data1,v_center,data2,'linewidth',2)    
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
      {sprintf('v_N > %g km/s',0)},...
      [0.98 0.98],'fontsize',fontsize_leg,'color','k')
  
    vlim = 1499;
    hca.XLim = vlim*[-1 1];  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = [-2000:500:2000];
    hca.XTickLabelRotation = 0;
    hca.Box = 'on';
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
  
  if 1 % 1D at a certain vn range, patches 
    hca = h(isub); isub = isub + 1;
    %vdf = pdist.reduce('2D',Mdsl,Ndsl,'lowerelim',elows);
    vdf = vdf_MN;

    vint1 = [0 inf];
    vint2 = [-inf -0];
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
    
   % vdf1 = PDist(times,data1,'1Dcart',vdf.depend{1});
   % vdf2 = PDist(times,data2,'1Dcart',vdf.depend{1});

    data1 = mean(data1,1);
    data2 = mean(data2,1);

    hca.ColorOrder = colors;

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data1_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data1_y = [data1(idx) data1(idx)*0 0];
    
    idx = find(v_center>0); % vM range
    data2_x = [v_center(idx) v_center(idx(end:-1:1)) v_center(idx(1))];
    data2_y = [data1(idx) data1(idx)*0 0];

    idx = find(v_center<0); idx = [idx, idx(end)+1]; % vM range
    data3_x = [v_center(idx) v_center(idx(end:-1:1))];
    data3_y = [data2(idx) data2(idx)*0];
    
    idx = find(v_center>0);
    data4_x = [v_center(idx) v_center(idx(end:-1:1))];
    data4_y = [data2(idx) data2(idx)*0];
    

    patch(hca,data3_x,data3_y,colors(3,:),...
        'facealpha',0.3,'edgecolor',colors(3,:)*0);
    hold(hca,'on')
    patch(hca,data4_x,data4_y,colors(4,:),...
        'facealpha',0.3,'edgecolor',colors(4,:)*0);
    
    %patch(hca,data3_x,data3_y,colors(3,:),...
    %    'facealpha',0.3,'edgecolor',colors(3,:)*0);
    %patch(hca,data4_x,data4_y,colors(4,:),...
    %    'facealpha',0.3,'edgecolor',colors(4,:)*0);
    hold(hca,'off')

    %plot(hca,v_center,data1,v_center,data2,'linewidth',2)    
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
      {sprintf('v_N < %g km/s',0)},...
      [0.98 0.98],'fontsize',fontsize_leg,'color','k')
  
    vlim = 1499;
    hca.XLim = vlim*[-1 1];  
    hca.XGrid = 'on';
    hca.YGrid = 'on';
    hca.XTick = [-2000:500:2000];
    hca.XTickLabelRotation = 0;
    hca.Box = 'on';
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
  

  if 0 % 1D at a certain vn range    
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
    
   % vdf1 = PDist(times,data1,'1Dcart',vdf.depend{1});
   % vdf2 = PDist(times,data2,'1Dcart',vdf.depend{1});

    data1 = mean(data1,1);
    data2 = mean(data2,1);

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
      {sprintf('v_N > %g km/s',0),sprintf('v_N < %g km/s',0)}',...
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
c_eval('h(?).LineWidth = 1.5;',1:numel(h))
h(1).Title.String = sprintf('%s - %s',tint_dist(1).utc('HH:MM:SS.mmm'),tint_dist(2).utc('SS.mmm'));
%compact_panels(h,0.01,0.2)

irf_legend(0,{sprintf('L=[%.2f,%.2f,%.2f], M = [%.2f,%.2f,%.2f], N = [%.2f,%.2f,%.2f]',L(1),L(2),L(3),M(1),M(2),M(3),N(1),N(2),N(3)),sprintf('nMean=[%g,%g,%g,%g], nThresh = %g',nMean(1),nMean(2),nMean(3),nMean(4),nThresh)},[0.05 1])

%% Distributions to illustrate gyroturning, divide f(v_L,v_M<>0)

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
nMovMean = 5;
%c_eval('pdist_all = iPDist?.movmean(nMovMean,''RemoveOneCounts'',iPDist?_counts);',ic)
c_eval('pdist_all = PD_clean.movmean(nMovMean);',ic)
%pdist_all.data(:,1:20,:,:) = 0;
%pdist_all = pdist_all.elim(elim);

times_utc = [...%'2017-07-11T22:33:45.000Z';...
             '2017-07-11T22:33:54.582Z';...
             %'2017-07-11T22:34:00.582Z';...
             %'2017-07-11T22:33:58.062Z';...
             '2017-07-11T22:34:00.502Z';...
             '2017-07-11T22:34:03.000Z';...
             '2017-07-11T22:34:05.500Z';...
             %'2017-07-11T22:34:08.540Z';...
             %'2017-07-11T22:34:20.940Z';...
             '2017-07-11T22:34:10.000Z';...
             '2017-07-11T22:34:20.000Z'];

time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');
time = irf_time('2017-07-11T22:34:02.000Z','utc>EpochTT');

dt = -3;
dt = 0;
time = time + dt;

h = setup_subplots(1,1);
isub = 1;

pdist = pdist_all.tlim(time+0.5*0.15*[-1 1]);

hca = h(isub); isub = isub + 1;
hca.ColorOrder = pic_colors('matlab');
hca.ColorOrder = circshift(hca.ColorOrder,1);
nSmooth = 3;
iso_values = 1*10.^[-27:-15];
iso_values = [6.5e-28];
%iso_values = [1e-27];
iso_values = [2e-27];
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
camlight(gca,0,1)
camlight(gca,180,0)

if 0 % Add max energy
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








%% Calculate azimuthal distribution in the LM plane (using macroparticles), non cleaned
units = irf_units;
PD = PD_clean;

nt = PD.length;
Vsc = scPot3.resample(PD);
nMP = 5000;
MP = PD.macroparticles('ntot',nMP,'scpot',Vsc);
vL_shift = vL_xline*1;
%vL_shift = -999;
E_edges = [PD.ancillary.energy(1,1) - PD.ancillary.delta_energy_minus(1,1), PD.ancillary.energy(1,:) + PD.ancillary.delta_energy_plus(1,:)];
E_minus = E_edges(1:end-1);
E_plus = E_edges(2:end);
v_minus = sqrt(2*units.e*E_minus/units.mp); % m/s
v_plus = sqrt(2*units.e*E_plus/units.mp); % m/s
% Cylindrical volume element v*dv*dtheta
% vol = (v2^2-v1^2)/2*Delta_theta
d_vel = (v_plus.^2 - v_minus.^2)/2; % (m/s)^2
d_vel_tot = (v_plus(end).^2 - v_minus(1).^2)/2; % (m/s)^2

%d3v = d_vel_mat.*d_azim.*d_polar_mat; % (m/s)^3

nEnergy = numel(E_edges)-1;
dazim = 5;
azimuth_edges = -180:dazim:180;
nAzimuth = numel(azimuth_edges)-1;


d_vel_mat = repmat(d_vel,nt,1,nAzimuth);
d2v = d_vel_mat.*dazim*pi/180; % (m/s)^2

d_vel_tot_mat = repmat(d_vel_tot,nt,nAzimuth);
d2v_intE = d_vel_tot_mat.*dazim; % (m/s)^3 ????

dn_tot_all_ = zeros(nt,nAzimuth);
dv_tot_all_ = zeros(nt,nAzimuth);
dn_tot_all = zeros(nt,nEnergy,nAzimuth);
df_tot_all = zeros(nt,nEnergy,nAzimuth);
dv_tot_all = zeros(nt,nEnergy,nAzimuth);

for it = 1:nt
  %elow = tsElow.resample(PD(it).time).data;
  dv = MP(it).dv; % same units as the original f (I think)
  df = MP(it).df;
  dn = df.*dv;
  vx = MP(it).vx;
  vy = MP(it).vy;
  vz = MP(it).vz;
  vL = [vx vy vz]*L' - vL_shift;
  vM = [vx vy vz]*M';
  vN = [vx vy vz]*N';
  theta_NL = atan2d(vN,vL);
  v2 = (vL.^2 + vM.^2 + vN.^2)*1e6; % (cm/s)^2 -> (m/s)^2
  vLN = sqrt(vL.^2 + vN.^2)*1e3; % (m/s)^1 -- It's the energy corresponding to this one should use, because the end results will be a reduced cylinder
  E = units.mp*vLN.^2/2/units.eV; % eV
  %E(E<elow) = 0;
  %dn(vM<0) = 0;
  %dn(vLN<1500) = 0;  
  fun = @sum;
  [dn_tot_ edges_ mid_ loc_] = histcn([theta_NL],azimuth_edges,'AccumData',dn,'Fun',fun);
  [dv_tot_ edges_ mid_ loc_] = histcn([theta_NL],azimuth_edges,'AccumData',dv,'Fun',fun);
  [dn_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dn,'Fun',fun);
  [dv_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',dv,'Fun',fun);
  [df_tot edges mid loc] = histcn([E, theta_NL],E_edges,azimuth_edges,'AccumData',df,'Fun',fun);
  % for iloc = 1:nAzimuth    
  %   dv_tot_all(it,iloc) = sum(MP(it).dv(loc==iloc));
  %   df_tot_all(it,iloc) = sum(MP(it).df(loc==iloc));
  % end
  dn_tot_all_(it,:) = dn_tot_;
  dv_tot_all_(it,:) = dv_tot_;


  dn_tot_all(it,:,:) = dn_tot;
  dv_tot_all(it,:,:) = dv_tot;
  df_tot_all(it,:,:) = df_tot;
  
end

% Not seperated by E
d3v_tot = sum(PD.d3v('mat'),[2:4]);
d3v_tot_az_bin_ = repmat(d3v_tot/nAzimuth,[1,nAzimuth]);
f_tot_2_ = dn_tot_all_./d3v_tot_az_bin_;

% Not seperated by E, 2
d3v_tot = sum(PD.d3v('mat'),[2:4]);
d3v_tot_az_bin_ = repmat(d3v_tot/nAzimuth,[1,nAzimuth]);
f_tot_2_ = dn_tot_all_./d3v_tot_az_bin_;

% Not seperated by E, 2
f_tot_3_ = dn_tot_all./d2v;


% iAzim__ = PDist(PD.time,f_tot,'azimuthangle',mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
% iAzim__.ancillary.v_shift = vL_shift;
% to_SI = (1e2*1e3)^6;
% iAzim__ = iAzim_*to_SI;
% %iAzim_.units = 's^3/m^6';
% eval(sprintf('iAzim__%04.0f = iAzim__;',abs(vL_shift)));

% Seperated by E
d3v_tot_per_E = sum(PD.d3v('mat'),[3:4]);
d3v_tot_az_bin = repmat(d3v_tot_per_E/nAzimuth,[1,1,nAzimuth]);
%f_tot_1 = dn_tot_all./dv_tot_all;
f_tot = dn_tot_all./d3v_tot_az_bin; 

% Seperated by E 2
%f_tot_1 = dn_tot_all./dv_tot_all;
f_tot_3 = dn_tot_all./dv_tot_all;

% Seperated by E 3
% [dn_tot_all] = cm^-3
% [d2v] = m^2/s^2
f_tot_4_ = dn_tot_all./(d2v*1e4); % cm^-3/(cm^2/s^2) = s^2/cm^5, reduced along one dimension



% NOT Seperated by E 4
%f_tot_1 = dn_tot_all./dv_tot_all;
f_tot_5 = dn_tot_all_./d2v_intE; % (1/cm^3)/(m^2/s^2) = (s^2)*(m^2/cm^3) = (s^2)*(m^2/cm^3)

iAzim_ = PDist(PD.time,f_tot,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
iAzim_.ancillary.v_shift = vL_shift;
to_SI = (1e2*1e3)^6;
iAzim_ = iAzim_*to_SI;
%iAzim_.units = 's^3/m^6';
eval(sprintf('iAzim_%04.0f = iAzim_;',abs(vL_shift)));
%iAzim_ = PDist(PD.time,f_tot_1,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6

iAzim_2 = PDist(PD.time,f_tot_3,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
iAzim_2.ancillary.v_shift = vL_shift;
to_SI = (1e2*1e3)^6;
iAzim_2 = iAzim_2*to_SI;
%iAzim_.units = 's^3/m^6';
eval(sprintf('iAzim_2_%04.0f = iAzim_2;',abs(vL_shift)));

iAzim_3 = PDist(PD.time,f_tot_4_,'azimuthangle',PD.depend{1},mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
iAzim_3.ancillary.v_shift = vL_shift;
to_SI = (1e2*1e3)^6;
iAzim_3 = iAzim_3*to_SI;
%iAzim_.units = 's^3/m^6';
eval(sprintf('iAzim_3_%04.0f = iAzim_3;',abs(vL_shift)));

f_tmp = reshape(f_tot_5,size(f_tot_5,1),1,size(f_tot_5,2));
f_tmp = repmat(f_tmp,1,2,1);
iAzim_4 = PDist(PD.time,f_tmp,'azimuthangle',PD.depend{1}(:,[10 11]),mid{2}); % scaling factor to go from 1/m^6 -> 1/cm^6
iAzim_4.ancillary.v_shift = vL_shift;
to_SI = (1e2*1e3)^6;
iAzim_4 = iAzim_4*to_SI;
%iAzim_.units = 's^3/m^6';
eval(sprintf('iAzim_4_%04.0f = iAzim_4;',abs(vL_shift)));