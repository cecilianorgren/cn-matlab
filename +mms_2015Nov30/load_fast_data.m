ic = 1;
tint = irf.tint('2015-11-30T00:00:00.00Z/2015-11-30T02:00:00.00Z'); %20151112071854
if 1
if 0
  %load('/Users/Cecilia/Data/MMS/20151112071854_2016-08-23.mat') % has this dobj thing
  load('/Users/Cecilia/Data/MMS/20151112071854_2016-08-24.mat')
  return
end
%% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
mms.db_init('local_file_db','/data/mms');
db_info = datastore('mms_db');   

%% Particle distributions: electrons and ions
disp('Loading particle distributions...')
%c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_fast_l2_des-dist'',tint+[20 0])); toc',ic)
%c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_fast_l2_dis-dist'',tint+[20 0])); toc',ic)
c_eval('tic; ePDist?fast = mms.get_data(''PDe_fpi_fast_l2'',tint,?); toc',ic)
c_eval('tic; iPDist?fast = mms.get_data(''PDi_fpi_fast_l2'',tint,?); toc',ic)

%% Setting tint from ePDist1
%tint = [ePDist1([1 ePDist1.length]).time];
end
%% Make event directory
fileName = iPDist1fast.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
eventPath = ['/Users/Cecilia/Research/Events/' dirName '/'];
mkdir(eventPath)

%% Load defatt, for coordinate tranformation
disp('Loading defatt...')
%load /Users/Cecilia/Data/MMS/2015Oct16/defatt.mat
c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic);
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
c_eval('defatt? = mms_removerepeatpnts(defatt?);',ic)
    
%% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB?srvy=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_dmpa_srvy_l2'',tint); toc;',ic);
c_eval('tic; gseB?srvy=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gse_srvy_l2'',tint); toc;',ic);
c_eval('tic; gsmB?srvy=mms.db_get_ts(''mms?_fgm_srvy_l2'',''mms?_fgm_b_gsm_srvy_l2'',tint); toc;',ic);

%% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?fast=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_gse_fast_l2'',tint); toc',ic);
c_eval('tic; dslE?fast=mms.db_get_ts(''mms?_edp_fast_l2_dce'',''mms?_edp_dce_dsl_fast_l2'',tint); toc',ic);

%% Load spacecraft position
disp('Loading spacecraft position...')
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

%% Spacecraft potential
%disp('Loading spacecraft potential...')
%c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);

%% Particle moments

% Skymap distributions
%c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
%c_eval('iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0]));',ic)

% Pressure and temperature
disp('Loading pressure and temperature...')
c_eval('tic; dbcsPe?fast = mms.get_data(''Pe_dbcs_fpi_fast_l2'',tint,?); toc;',ic) 
c_eval('tic; dbcsTe?fast = mms.get_data(''Te_dbcs_fpi_fast_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsPi?fast = mms.get_data(''Pi_dbcs_fpi_fast_l2'',tint,?); toc;',ic) 
c_eval('tic; dbcsTi?fast = mms.get_data(''Ti_dbcs_fpi_fast_l2'',tint,?); toc;',ic)

c_eval('gsePe?fast = mms.rotate_tensor(dbcsPe?fast,''gse'',?); gsePe?fast.units = ''nPa''; gsePe?fast.coordinateSystem = ''GSE'';',ic)
c_eval('gseTe?fast = mms.rotate_tensor(dbcsTe?fast,''gse'',?);',ic)
c_eval('gsePi?fast = mms.rotate_tensor(dbcsPi?fast,''gse'',?); gsePi?fast.units = ''nPa''; gsePe?fast.coordinateSystem = ''GSE'';',ic)
c_eval('gseTi?fast = mms.rotate_tensor(dbcsTi?fast,''gse'',?);',ic)

c_eval('facPe?fast = mms.rotate_tensor(gsePe?fast,''fac'',gseB?srvy); facPe?fast.units = ''nPa''; facPe?fast.coordinateSystem = ''FAC'';',ic)
c_eval('facTe?fast = mms.rotate_tensor(gseTe?fast,''fac'',gseB?srvy); facTi?fast.units = ''eV''; facTi?fast.coordinateSystem = ''FAC'';',ic)
c_eval('facPi?fast = mms.rotate_tensor(gsePi?fast,''fac'',gseB?srvy); facPi?fast.units = ''nPa''; facPi?fast.coordinateSystem = ''FAC'';',ic)
c_eval('facTi?fast = mms.rotate_tensor(gseTi?fast,''fac'',gseB?srvy); facTi?fast.units = ''eV''; facTi?fast.coordinateSystem = ''FAC'';',ic)

% Density
disp('Loading density...')
c_eval('tic; ne?fast = mms.db_get_ts(''mms?_fpi_fast_l2_des-moms'',''mms?_des_numberdensity_dbcs_fast'',tint); toc;',ic);
c_eval('tic; ni?fast = mms.db_get_ts(''mms?_fpi_fast_l2_dis-moms'',''mms?_dis_numberdensity_dbcs_fast'',tint); toc;',ic);

% Velocity
disp('Loading bulk velocities...')
c_eval('tic; dbcsVe?fast = mms.get_data(''Ve_dbcs_fpi_fast_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsVi?fast = mms.get_data(''Vi_dbcs_fpi_fast_l2'',tint,?); toc;',ic)

%c_eval('tic; dbcsVe?fast = mms.get_data(''Ve_dbcs_fpi_fast_l2'',tint,?); toc;',ic)
%c_eval('tic; dbcsVi?fast = mms.get_data(''Vi_dbcs_fpi_fast_l2'',tint,?); toc;',ic)

c_eval('gseVe?fast = mms_dsl2gse(dbcsVe?fast,defatt?); gseVe?fast.coordinateSystem = ''GSE'';',ic)
c_eval('gseVi?fast = mms_dsl2gse(dbcsVi?fast,defatt?); gseVe?fast.coordinateSystem = ''GSE'';',ic)
