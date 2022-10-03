ic = 1;
tint_fast = irf.tint('2017-07-25T20:00:00.00Z/2017-07-25T24:00:00.00Z');

% Load datastore
localuser = datastore('local','user');
%mms.db_init('local_file_db',['/Users/' localuser '/data']);
mms.db_init('local_file_db',['/Users/' localuser '/Data/MMS']);
db_info = datastore('mms_db');   

%%

c_eval('gseB?_srvy = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',tint_fast);',ic);
c_eval('gsmB?_srvy = mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gsm_srvy_l2'',tint_fast);',ic);
c_eval('iPDist?_fast = mms.get_data(''PDi_fpi_fast_l2'',tint_fast,?);',ic) % missing some ancillary data
c_eval('gseVi?_fast = mms.get_data(''Vi_gse_fpi_fast_l2'',tint_fast,?);',ic);
c_eval('gseVe?_fast = mms.get_data(''Ve_gse_fpi_fast_l2'',tint_fast,?);',ic);
c_eval('gseE?_fast = mms.get_data(''E_gse_edp_fast_l2'',tint_fast,?);',ic);

c_eval('gseVExB?_srvy = gseE?_fast.resample(gseB?_srvy).cross(gseB?_srvy);',ic)
%c_eval('gseVExB?_srvy = cross(gseE?_fast.resample(gseB?_srvy.time),gseB?)/gseB?_srvy.abs2*1e3; gseVExB?_srvy.units = '''';',ic) % km/s

c_eval('ne?_fast = mms.get_data(''Ne_fpi_fast_l2'',tint_fast,?);',ic);
c_eval('ni?_fast = mms.get_data(''Ni_fpi_fast_l2'',tint_fast,?);',ic);
c_eval('Pi?_fast = mms.get_data(''Pi_gse_fpi_fast_l2'',tint_fast,?);',ic);
c_eval('Pe?_fast = mms.get_data(''Pe_gse_fpi_fast_l2'',tint_fast,?);',ic);

c_eval('nOp?_srvy = mms.get_data(''Noplus_hpca_srvy_l2'',tint_fast,?);',ic);
c_eval('nHp?_srvy = mms.get_data(''Nhplus_hpca_srvy_l2'',tint_fast,?);',ic);

%c_eval('pOp?_srvy = mms.get_data(''Poplus_hpca_srvy_l2'',tint_fast,?);',ic);
%c_eval('pHp?_srvy = mms.get_data(''Phplus_hpca_srvy_l2'',tint_fast,?);',ic);

%c_eval('pOps?_srvy = irf.ts_scalar(pOp?_srvy.time,mean(pOp?_srvy.data,2));',ic);
%c_eval('pHps?_srvy = irf.ts_scalar(pHp?_srvy.time,mean(pHp?_srvy.data,2));',ic);

c_eval('vOp?_srvy = mms.get_data(''Voplus_gsm_hpca_srvy_l2'',tint_fast,?);',ic);
c_eval('vHp?_srvy = mms.get_data(''Vhplus_gsm_hpca_srvy_l2'',tint_fast,?);',ic);


c_eval('iPDist?_Hp_fast = mms.get_data(''Omnifluxhplus_hpca_srvy_l2'',tint_fast,?);',ic) % missing some ancillary data 
c_eval('iPDist?_Op_fast = mms.get_data(''Omnifluxoplus_hpca_srvy_l2'',tint_fast,?);',ic) % missing some ancillary data 
%c_eval('iPDist?_Opp_fast = mms.get_data(''Omnifluxoplusplus_hpca_srvy_l2'',tint_fast,?);',ic) % missing some ancillary data 

c_eval('[gseVi?_fast_par,gseVi?_fast_perp] = irf_dec_parperp(gseB?_srvy.resample(gseVi?_fast),gseVi?_fast); gseVi?_fast_par.name = ''Vi par''; gseVi?_fast_perp.name = ''Vi perp'';',ic)
c_eval('[gseVe?_fast_par,gseVe?_fast_perp] = irf_dec_parperp(gseB?_srvy.resample(gseVe?_fast),gseVe?_fast); gseVe?_fast_par.name = ''Ve par''; gseVe?_fast_perp.name = ''Ve perp'';',ic)

% Convective electric fields
c_eval('gseVexB?_fast = gseVe?_fast.cross(gseB?_srvy.resample(gseVe?_fast))*1e-3; gseVexB?_fast.units = ''mV/m'';',ic)
c_eval('gseVixB?_fast = gseVi?_fast.cross(gseB?_srvy.resample(gseVi?_fast))*1e-3; gseVixB?_fast.units = ''mV/m'';',ic)
% Non-ideal electric field, E+VexB
c_eval('gseEVexB?_fast = gseE?_fast.resample(gseVexB?_fast.time)+gseVexB?_fast; gseEVexB?_fast.name = ''E+VexB'';',ic)
c_eval('gseEVixB?_fast = gseE?_fast.resample(gseVixB?_fast.time)+gseVixB?_fast; gseEVixB?_fast.name = ''E+VixB'';',ic)

% Magnetic field pressure
c_eval('PB?_srvy = gseB?_srvy.abs2/2/units.mu0*1e-9; PB?_srvy.name = ''Magnetic pressure''; PB?_srvy.units =''nPa'';',ic)
% Plasma beta
%c_eval('betae?_srvy = PB?.resample(gsePi?)*3/(gsePe?.trace.resample(gsePi?));',ic)
%c_eval('betai?_srvy = PB?.resample(gsePi?)*3/(gsePi?.trace);',ic)
c_eval('beta?_srvy_fpi = PB?_srvy.resample(Pi?_fast)*3/(Pi?_fast.trace+Pe?_fast.trace.resample(Pi?_fast)); beta?_srvy_fpi.name = ''beta ie'';',ic)
%c_eval('betai?_srvy_hpca = PB?_srvy.resample(pOp?_srvy)*3/(pOp?_srvy.trace + pHp?_srvy.trace); betai?_srvy_hpca.name = ''beta op hp'';',ic)

% Energetic particles
c_eval('feeps_ion_omni?_srvy = mms.get_data(''Omnifluxion_epd_feeps_srvy_l2'',tint_fast,?);',ic)
c_eval('feeps_pa?_srvy = mms.get_data(''Pitchanglefluxion_epd_feeps_srvy_l2'',tint_fast,?);',ic)
c_eval('feeps_ele_omni?_srvy = mms.get_data(''Omnifluxelectron_epd_feeps_srvy_l2'',tint_fast,?);',ic)
%% Rotated coordinates
tint_mva = irf.tint('2017-07-25T20:14:08.398745849Z/2017-07-25T21:57:25.093696533Z'); % early
%tint_mva = irf.tint('2017-07-25T21:36:16.246635253Z/2017-07-25T23:52:59.490427001Z'); % later

[out,l,v]=irf_minvar(gseB1_srvy.tlim(tint_mva));

Radjust = [-1 0 0; 0 -1 0;0 0 1]; % early
%Radjust = [-1 0 0; 0 1 0;0 0 -1]; % later
%Radjust = [1 0 0; 0 1 0; 0 0 1];
R = v*Radjust';


c_eval('lmnB?_srvy = gseB?_srvy*R'';',ic)
c_eval('lmnE?_fast = gseE?_fast*R'';',ic)
c_eval('lmnVi?_fast = gseVi?_fast*R'';',ic)
c_eval('lmnVe?_fast = gseVe?_fast*R'';',ic)
c_eval('lmnVExB?_srvy = gseVExB?_srvy*R'';',ic)

c_eval('lmnVi? = gseVi?*R'';',ic)