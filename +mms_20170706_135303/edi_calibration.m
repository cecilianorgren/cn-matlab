%edi_calibration
units = irf_units;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
tint = tint + [-5 5]; % using the above edges causes problem with new EDI files because they have different versions that adjoining file
tint_zoom = irf.tint('2017-07-06T13:54:00.000Z/2017-07-06T13:54:10.000Z');
tint_zoom = irf.tint('2017-07-06T13:53:55.000Z/2017-07-06T13:54:15.000Z');
mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
db_info = datastore('mms_db');   

% Load data
c_eval('tic; gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
R = mms.get_data('R_gse',tint);
c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
c_eval('tic; [iPDist?,iPDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; [ePDist?,ePDistErr?] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint+[5 -5],?);',4)
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)

% Set up edi and fpi pitch angle data
% EDI energy and corresponding velocity
E_edi = 500; % eV
% v_edi = sqrt(2*units.e*E_edi./units.me); % m/s
% dE_edi = 25; % eV
% 
% E_edi_plus = E_edi + dE_edi;
% E_edi_minus = E_edi - dE_edi;
% v_edi_plus = sqrt(2*units.e*E_edi_plus./units.me); % m/s
% v_edi_minus = sqrt(2*units.e*E_edi_minus./units.me); % m/s
% v_edi_plusminus = v_edi_plus-v_edi_minus;
% dv_edi_minus = v_edi_minus - v_edi;
% dv_edi_plus = v_edi_plus - v_edi;
% dv_edi = dv_edi_plus - dv_edi_minus; % m/s
% 
% v_edi_edges = [v_edi_minus,v_edi_plus]*1e-3;
% pa_edi_edges = -45:11.25:45;
% pa_edi_centers = [(-45+11.25):11.25:0 0:11.25:(45-11.25)];
% az_edges = 0:10:360;

tint_zoom = irf.tint('2017-07-06T13:54:05.50Z/2017-07-06T13:54:05.65Z');
tint_zoom = tint_zoom + [-1 1];

c_eval('ePitch?_fpi = ePDist?.pitchangles(dmpaB?,16);',1:4);
c_eval('ePitch?_flux_fpi = ePitch?_fpi.flux;',1:4);
c_eval('ePitch?_flux_fpi_bin_closest_to_edi = ePitch?_flux_fpi.elim(E_edi);',1:4);

% Get FPI data from moms, never mind this, it integrates over too many
% pitch angles
%dobj = dataobj('/Volumes/Fountain/Data/MMS/mms1/fpi/brst/l2/des-moms/2017/07/06/mms1_fpi_brst_l2_des-moms_20170706135303_v3.3.0.cdf');