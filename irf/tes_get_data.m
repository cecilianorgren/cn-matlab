% test_get_data
% mms1_fpi_fast_l2_des-moms_20170725160000_v3.3.0.cdf

%mms.db_init('local_file_db','/data/mms/');
mms.db_init('local_file_db','/Users/cno062/Data/MMS');

tint = irf.tint('2017-07-25T16:00:00.00Z/2017-07-25T18:00:00.00Z');

ic = 1;

c_eval('gseR = mms.get_data(''R_gse'',tint,?);',ic)
c_eval('ne = mms.get_data(''Ne_fpi_fast_l2'',tint,?);',ic);  
c_eval('gseTe = mms.get_data(''Te_dbcs_fpi_fast_l2'',tint,?);',ic); 
c_eval('gseTi = mms.get_data(''Ti_dbcs_fpi_fast_l2'',tint,?);',ic);
c_eval('gsePe = mms.get_data(''Pe_dbcs_fpi_fast_l2'',tint,?);',ic);
c_eval('gsePi = mms.get_data(''Pi_dbcs_fpi_fast_l2'',tint,?);',ic);
c_eval('gseVi = mms.get_data(''Vi_dbcs_fpi_fast_l2'',tint,?);',ic);