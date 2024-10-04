ic = 1:4;
tint = irf.tint('2015-11-05T05:18:05.00Z/2015-11-05T05:19:00.00Z');
tint = irf.tint('2015-11-05T04:47:05.00Z/2015-11-05T04:48:00.00Z');
tint = irf.tint('2015-11-03T12:55:04.00Z/2015-11-03T12:58:00.00Z'); %20151103125504
tint = irf.tint('2015-11-03T15:49:14.00Z/2015-11-03T15:52:00.00Z'); %20151103154914

tint = irf.tint('2015-10-16T11:27:54.00Z/2015-10-16T11:29:13.00Z'); %20151016112754
%20151016130334
%20151016130524
%20151017162454
tint = irf.tint('2015-09-15T11:17:44.00Z/2015-09-15T11:21:24.00Z'); %20150915111744
tint = irf.tint('2016-01-01T09:03:14.00Z/2016-01-01T09:04:15.00Z'); %20160101090314
tint = irf.tint('2015-10-16T10:34:14.00Z/2015-10-16T10:35:55.00Z'); %20151016103414
tint = irf.tint('2015-11-12T07:18:54.00Z/2015-11-12T07:19:45.00Z'); %20151112071854
tint = irf.tint('2015-11-17T14:15:54.00Z/2015-11-17T14:18:14.00Z'); %20151117141554
tint = irf.tint('2016-02-02T09:03:14.00Z/2016-02-02T09:04:14.00Z'); %20160101090314
tint = irf.tint('2017-08-20T01:53:00.00Z/2017-08-20T01:57:00.00Z'); %20160101090314
%% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/data/mms');
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');   



%% Particle distributions: electrons and ions
disp('Loading particle distributions...')
c_eval('tic; ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0])); toc',ic)
c_eval('tic; iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0])); toc',ic)

%c_eval('tic; iPDist? = mms.get_data(''PDi_fpi_brst_l2'',tint,?); toc',ic)
%% Setting tint from ePDist1
tint = [ePDist1([1 ePDist1.length]).time];

%% Make event directory
fileName = ePDist1.userData.GlobalAttributes.Logical_file_id;
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
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint); toc;',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint); toc;',ic);

%% Electric field
disp('Loading electric field...')
c_eval('tic; gseE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint); toc',ic);
c_eval('tic; dslE?=mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint); toc',ic);

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

%% Particle moments
% load dataobject first
% electrons
if 0
  c_eval('edobj?=dataobj(mms.get_filepath(''mms?_fpi_brst_l2_des-moms'',tint(1)));',ic)

  c_eval('ne? = mms.variable2ts(get_variable(edobj?,''mms?_des_numberdensity_dbcs_brst''));',ic);

  c_eval('vex? = mms.variable2ts(get_variable(edobj?,''mms?_des_bulkx_dbcs_brst''));',ic);
  c_eval('vey? = mms.variable2ts(get_variable(edobj?,''mms?_des_bulky_dbcs_brst''));',ic);
  c_eval('vez? = mms.variable2ts(get_variable(edobj?,''mms?_des_bulkz_dbcs_brst''));',ic);
  c_eval('dbcsVe?=irf.ts_vec_xyz(vex?.time,[vex?.data vey?.data vez?.data]);',ic)
  c_eval('dbcsVe1.name = ''ve''; dbcsVe?.userData = vex?.userData;',ic)

  c_eval('Pexx? = mms.variable2ts(get_variable(edobj?,''mms?_des_presxx_dbcs_brst''));',ic);
  c_eval('Pexy? = mms.variable2ts(get_variable(edobj?,''mms?_des_presxy_dbcs_brst''));',ic);
  c_eval('Pexz? = mms.variable2ts(get_variable(edobj?,''mms?_des_presxz_dbcs_brst''));',ic);  
  c_eval('Peyy? = mms.variable2ts(get_variable(edobj?,''mms?_des_presyy_dbcs_brst''));',ic);
  c_eval('Peyz? = mms.variable2ts(get_variable(edobj?,''mms?_des_presyz_dbcs_brst''));',ic);  
  c_eval('Pezz? = mms.variable2ts(get_variable(edobj?,''mms?_des_preszz_dbcs_brst''));',ic);   
  c_eval('DataPe = zeros(Pexx?.length,3,3); DataPe(:,1,1) = Pexx?.data; DataPe(:,1,2) = Pexy?.data; DataPe(:,1,3) = Pexz?.data; DataPe(:,2,2) = Peyy?.data; DataPe(:,2,3) = Peyz?.data; DataPe(:,3,3) = Pezz?.data;',ic);
  c_eval('dbcsPe?=irf.ts_tensor_xyz(Pexx?.time,DataPe);',ic)
  c_eval('dbcsPe?.name = ''Pe''; dbcsPe?.units = ''nPa''; dbcsPe?.userData = Pexx?.userData; dbcsPe?.coordinateSystem = ''DBCS'';',ic)

  c_eval('Texx? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempxx_dbcs_brst''));',ic);
  c_eval('Texy? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempxy_dbcs_brst''));',ic);
  c_eval('Texz? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempxz_dbcs_brst''));',ic);  
  c_eval('Teyy? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempyy_dbcs_brst''));',ic);
  c_eval('Teyz? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempyz_dbcs_brst''));',ic);  
  c_eval('Tezz? = mms.variable2ts(get_variable(edobj?,''mms?_des_tempzz_dbcs_brst''));',ic);   
  c_eval('DataTe = zeros(Texx?.length,3,3); DataTe(:,1,1) = Texx?.data; DataTe(:,1,2) = Texy?.data; DataTe(:,1,3) = Texz?.data; DataTe(:,2,2) = Teyy?.data; DataTe(:,2,3) = Teyz?.data; DataTe(:,3,3) = Tezz?.data;',ic);
  c_eval('dbcsTe?=irf.ts_tensor_xyz(Texx?.time,DataTe);',ic)
  c_eval('dbcsTe?.name = ''Te''; dbcsTe?.units = ''eV''; dbcsTe?.userData = Texx?.userData; dbcsTe?.coordinateSystem = ''DBCS'';',ic)

  % ions
  c_eval('idobj?=dataobj(mms.get_filepath(''mms?_fpi_brst_l2_dis-moms'',tint(1)));',ic)

  c_eval('ni? = mms.variable2ts(get_variable(idobj?,''mms?_dis_numberdensity_dbcs_brst''));',ic);

  c_eval('vix? = mms.variable2ts(get_variable(idobj?,''mms?_dis_bulkx_dbcs_brst''));',ic);
  c_eval('viy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_bulky_dbcs_brst''));',ic);
  c_eval('viz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_bulkz_dbcs_brst''));',ic);
  c_eval('dbcsVi?=irf.ts_vec_xyz(vix?.time,[vix?.data viy?.data viz?.data]);',ic)
  c_eval('dbcsVi?.name = ''vi''; dbcsVi?.userData = vix?.userData; dbcsVi?.coordinateSystem = ''DBCS'';',ic)

  c_eval('Pixx? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presxx_dbcs_brst''));',ic);
  c_eval('Pixy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presxy_dbcs_brst''));',ic);
  c_eval('Pixz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presxz_dbcs_brst''));',ic);  
  c_eval('Piyy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presyy_dbcs_brst''));',ic);
  c_eval('Piyz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_presyz_dbcs_brst''));',ic);  
  c_eval('Pizz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_preszz_dbcs_brst''));',ic);   
  c_eval('DataPi = zeros(Pixx?.length,3,3); DataPi(:,1,1) = Pixx?.data; DataPi(:,1,2) = Pixy?.data; DataPi(:,1,3) = Pixz?.data; DataPi(:,2,2) = Piyy?.data; DataPi(:,2,3) = Piyz?.data; DataPi(:,3,3) = Pizz?.data;',ic);
  c_eval('dbcsPi?=irf.ts_tensor_xyz(Pixx?.time,DataPi);',ic)
  c_eval('dbcsPi?.name = ''Pi''; dbcsPi?.units = ''nPa''; dbcsPi?.userData = Pixx?.userData; dbcsPi?.coordinateSystem = ''DBCS'';',ic)

  c_eval('Tixx? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempxx_dbcs_brst''));',ic);
  c_eval('Tixy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempxy_dbcs_brst''));',ic);
  c_eval('Tixz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempxz_dbcs_brst''));',ic);  
  c_eval('Tiyy? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempyy_dbcs_brst''));',ic);
  c_eval('Tiyz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempyz_dbcs_brst''));',ic);  
  c_eval('Tizz? = mms.variable2ts(get_variable(idobj?,''mms?_dis_tempzz_dbcs_brst''));',ic);   
  c_eval('DataTi = zeros(Tixx?.length,3,3); DataTi(:,1,1) = Tixx?.data; DataTi(:,1,2) = Tixy?.data; DataTi(:,1,3) = Tixz?.data; DataTi(:,2,2) = Tiyy?.data; DataTi(:,2,3) = Tiyz?.data; DataTi(:,3,3) = Tizz?.data;',ic);
  c_eval('dbcsTi?=irf.ts_tensor_xyz(Tixx?.time,DataTi);',ic)
  c_eval('dbcsTi?.name = ''Ti''; dbcsTi?.units = ''eV''; dbcsTi?.userData = Tixx?.userData; dbcsTi?.coordinateSystem = ''DBCS'';',ic)
  
  % Rotate into GSE coordinates
  disp('Rotate into GSE...'); tic
  c_eval('gsePe? = mms.rotate_tensor(dbcsPe?,''gse'',?); gsePe?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
  c_eval('gsePi? = mms.rotate_tensor(dbcsPi?,''gse'',?); gsePi?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
  c_eval('gseTe? = mms.rotate_tensor(dbcsTe?,''gse'',?);',ic)
  c_eval('gseTi? = mms.rotate_tensor(dbcsTi?,''gse'',?);',ic)
  c_eval('gseVe? = mms_dsl2gse(dbcsVe?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic)
  c_eval('gseVi? = mms_dsl2gse(dbcsVi?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic); toc

  % Rotate into FAC coordinates
  disp('Rotate into FAC...'); tic
  c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?); facTe?.coordinateSystem = ''FAC'';',ic)
  c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.coordinateSystem = ''FAC'';',ic)
  c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?); facTi?.coordinateSystem = ''FAC'';',ic)
  c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.coordinateSystem = ''FAC'';',ic); toc  
end

% Skymap distributions
%c_eval('ePDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_des-dist'',tint+[20 0]));',ic)
%c_eval('iPDist? = mms.make_pdist(mms.get_filepath(''mms?_fpi_brst_l2_dis-dist'',tint+[20 0]));',ic)

% Pressure and temperature
disp('Loading pressure and temperature...')
c_eval('tic; dbcsPe? = mms.get_data(''Pe_dbcs_fpi_brst_l2'',tint,?); toc;',ic) 
c_eval('tic; dbcsTe? = mms.get_data(''Te_dbcs_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsPi? = mms.get_data(''Pi_dbcs_fpi_brst_l2'',tint,?); toc;',ic) 
c_eval('tic; dbcsTi? = mms.get_data(''Ti_dbcs_fpi_brst_l2'',tint,?); toc;',ic)

c_eval('gsePe? = mms.rotate_tensor(dbcsPe?,''gse'',?); gsePe?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseTe? = mms.rotate_tensor(dbcsTe?,''gse'',?);',ic)
c_eval('gsePi? = mms.rotate_tensor(dbcsPi?,''gse'',?); gsePi?.units = ''nPa''; gsePe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseTi? = mms.rotate_tensor(dbcsTi?,''gse'',?);',ic)

c_eval('facPe? = mms.rotate_tensor(gsePe?,''fac'',gseB?); facPe?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTe? = mms.rotate_tensor(gseTe?,''fac'',gseB?);',ic)
c_eval('facPi? = mms.rotate_tensor(gsePi?,''fac'',gseB?); facPi?.units = ''nPa''; facPe?.coordinateSystem = ''FAC'';',ic)
c_eval('facTi? = mms.rotate_tensor(gseTi?,''fac'',gseB?);',ic)

% Density
disp('Loading density...')
c_eval('tic; ne? = mms.db_get_ts(''mms?_fpi_brst_l2_des-moms'',''mms?_des_numberdensity_dbcs_brst'',tint); toc;',ic);
c_eval('tic; ni? = mms.db_get_ts(''mms?_fpi_brst_l2_dis-moms'',''mms?_dis_numberdensity_dbcs_brst'',tint); toc;',ic);

% Velocity
disp('Loading bulk velocities...')
c_eval('tic; dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?); toc;',ic)
c_eval('tic; dbcsVi? = mms.get_data(''Vi_dbcs_fpi_brst_l2'',tint,?); toc;',ic)

%c_eval('tic; dbcsVe?fast = mms.get_data(''Ve_dbcs_fpi_fast_l2'',fastTint,?); toc;',ic)
%c_eval('tic; dbcsVi?fast = mms.get_data(''Vi_dbcs_fpi_fast_l2'',fastTint,?); toc;',ic)

c_eval('gseVe? = mms_dsl2gse(dbcsVe?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic)
c_eval('gseVi? = mms_dsl2gse(dbcsVi?,defatt?); gseVe?.coordinateSystem = ''GSE'';',ic)
