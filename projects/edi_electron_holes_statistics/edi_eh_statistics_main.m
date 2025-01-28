%% Load data
ic = 1:4;
tint = irf.tint('2017-07-06T13:53:03.00Z/2017-07-06T13:55:33.00Z');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
db_info = datastore('mms_db');

c_eval('gseB?  = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseE?  = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('scPot? = mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint);',ic);
c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('ePitch?_flux_edi = mms.get_data(''Flux-amb-pm2_edi_brst_l2'',tint,?);',ic)


%% Cross-correlate EDI and EDP

% High-pass filter integrated E field, to avoid difting due to DC field
c_eval('A?_000 = ePitch?_flux_edi.palim([0 15]);',ic)
c_eval('A?_180 = ePitch?_flux_edi.palim([165 180]);',ic)

c_eval('intEdt? = irf_integrate(gseE?par);');    
flow = 50;
c_eval('B? = gseE?par.filt(flow,0,[],5);',ic)
c_eval('B?_int = intEdt?.filt(flow,0,[],5);',ic)

T = 0.1;
dt = 0.02;
maxLag = 0;
c_eval('tsC?_000 = edi_eh_running_average_correlation(A?_000,B?_int,T,dt,maxLag);',ic)
c_eval('tsC?_180 = edi_eh_running_average_correlation(A?_180,B?_int,T,dt,maxLag);',ic)

maxLagE = 10;
c_eval('[tsCE?!,tsLE?!] = edi_eh_running_average_correlation(gseE?par,gseE!par,T,dt,maxLagE);',ic,ic)

irf_plot({A1_000,A1_180,B1,B1_int,tsC1_000,tsC1_180})
