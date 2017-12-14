ic = 1:4;
tint = irf.tint('2015-11-12T07:18:54.00Z/2015-11-12T07:19:45.00Z'); %20151112071854

if 0
  %load('/Users/Cecilia/Data/MMS/20151112071854_2016-08-23.mat') % has this dobj thing
  load('/Users/Cecilia/Data/MMS/20151112071854_2016-08-24.mat')
  return
end
%% Load datastore
mms.db_init('local_file_db','/Volumes/Nexus/data');
db_info = datastore('mms_db');   

%% Particle distributions: electrons and ions
disp('Loading particle distributions...')
c_eval('tic; [ePDist?,ePDist?error] = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)

%% Setting tint from ePDist1
tint = ePDist1([1 ePDist1.length]).time;

%% Make event directory
fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
fileNameSplit = strsplit(fileName{1},'_'); numName = fileNameSplit{6};
dirName = sprintf('%s-%s-%s_%s',numName(1:4),numName(5:6),numName(7:8),numName(9:14));
[~,computername]=system('hostname');
%if strfind(computername,'ift0227887')
  eventPath = ['/Users/cno062/Research/Events/' dirName '/'];
%else
%  eventPath = ['/Users/Cecilia/Research/Events/' dirName '/'];
%end
mkdir(eventPath)

%% Load defatt, for coordinate tranformation
disp('Loading defatt...')
%load /Users/Cecilia/Data/MMS/2015Oct16/defatt.mat
c_eval('defatt? = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic);
c_eval('defatt?.zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint).zdec;',ic);
c_eval('defatt? = mms_removerepeatpnts(defatt?);',ic)
    
%% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);
c_eval('tic; gsmB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gsm_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB?scm=mms.db_get_ts(''mms?_scm_brst_l2_scb'',''mms?_scm_acb_gse_scb_brst_l2'',tint); toc',ic);

%% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?par=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_par_epar_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?hmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_dsl_brst_l2'',tint); toc',ic);
c_eval('tic; E?parhmfe=mms.db_get_ts(''mms?_edp_brst_l2_hmfe'',''mms?_edp_hmfe_par_epar_brst_l2'',tint); toc',ic);

%% Load spacecraft position
disp('Loading spacecraft position...')
R = mms.get_data('R_gse',tint);
if size(R.gseR1,2) == 4
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?(:,2:4));',1:4); % dfg_srvy_l2pre
else
  c_eval('gseR? = irf.ts_vec_xyz(R.time,R.gseR?);',1:4); % mec
end

%% Spacecraft potential
disp('Loading spacecraft potential...')
c_eval('tic; scPot?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_scpot_brst_l2'',tint); toc;',ic);
c_eval('tic; dcv?=mms.db_get_ts(''mms?_edp_brst_l2_scpot'',''mms?_edp_dcv_brst_l2'',tint); toc;',ic);

%% Particle moments

% Skymap distributions
%c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
%c_eval('iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0]));',ic)

% Pressure and temperature
disp('Loading pressure and temperature...'); tic
c_eval('gsePe? = mms.get_data(''Pe_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTe? = mms.get_data(''Te_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gsePi? = mms.get_data(''Pi_gse_fpi_brst_l2'',tint,?);',ic) 
c_eval('gseTi? = mms.get_data(''Ti_gse_fpi_brst_l2'',tint,?);',ic); toc

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...'); tic;
c_eval('ne? = mms.get_data(''Ne_fpi_brst_l2'',tint,?);',ic);
c_eval('ni? = mms.get_data(''Ni_fpi_brst_l2'',tint,?);',ic); toc

% Velocity
disp('Loading bulk velocities...'); tic
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVi? = mms.get_data(''Vi_gse_fpi_brst_l2'',tint,?);',ic); toc
c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?);',ic); toc

%c_eval('tic; gseVe?fast = mms.get_data(''Ve_gse_fpi_fast_l2'',fastTint,?); toc;',ic)
%c_eval('tic; gseVi?fast = mms.get_data(''Vi_gse_fpi_fast_l2'',fastTint,?); toc;',ic)
