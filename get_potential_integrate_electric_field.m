% Integrate electroic field
tint = 1;

mms.db_init('local_file_db','/Volumes/mms');
%db_info = datastore('mms_db');

units = irf_units;
ic = 1:4;
tint = irf.tint('2017-07-06T15:28:00.00Z/2017-07-06T15:30:00.00Z');
%%
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',1:4);
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',1:4);
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',1:4);
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);

c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic);

c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)

c_eval('[gseE?par,gseE?perp] = irf_dec_parperp(gseB?,gseE?); gseE?par.name = ''E par''; gseE?perp.name = ''E perp'';',ic)
c_eval('[gseVe?par,gseVe?perp] = irf_dec_parperp(gseB?,gseVe?); gseVe?par.name = ''Ve par''; gseVe?perp.name = ''Ve perp'';',ic)

%%
tint = irf.tint('2017-07-06T15:29:19.099354003Z/2017-07-06T15:29:19.195225830Z');
v_ph_gse = 430 * [-0.42 -0.48 -0.77]; %509*[-0.47 -0.40 -0.79];
c_eval('b?_loc = mean(gseB?.tlim(tint).norm.data,1);',ic)
b_av = (b1_loc + b2_loc + b3_loc + b4_loc)/4;
vph_par = dot(v_ph_gse,b_av);
ffilt = 2; % Hz
c_eval('Etoint? = gseE?par;',ic)
c_eval('Etoint? = Etoint?.filt(ffilt,0,[],3);',ic);    
c_eval('intEdt? = irf_integrate(Etoint?);',ic);    
%minpeakdistance = 50;
%c_eval('[PKS?,LOCS?,W?] = findpeaks(intEdt?.data,''MinPeakDistance'',minpeakdistance);',ic)
c_eval('intEdt?_detrend = intEdt?;',ic)

c_eval('phi? = intEdt?.tlim(tint+0.0*[-1 1])*vph_par;',ic);    


