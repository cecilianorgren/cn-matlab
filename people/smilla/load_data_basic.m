%% Set up database
mms.db_init('local_file_db','/data/mms');
%db_info = datastore('mms_db');

%% Specify time to load
units = irf_units;
% Torbert event 
ic = 3; % spacecraft number
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z');

%% Load data
% c_eval replaces all the ? with the number specifed at the end, e.g. 
% ic = 1:4 above to load the data from all spacecraft
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);
c_eval('iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic);
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic);

c_eval('gseR? = mms.get_data(''R_gse'',tint,?);',1:4)

%% Plot 1
% The most way to plot a set of quantities
irf_plot({gseB3,gseVi3,gseVe3,gseE3,gsePi3,iPDist3.deflux.omni.specrec('energy')})

%% Plot 2
% Plot components separately
irf_plot({gseB3.x,gseB3.y,gseB4.z},'comp')

