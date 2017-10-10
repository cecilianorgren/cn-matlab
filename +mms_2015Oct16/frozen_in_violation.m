%% Load data
tint = irf.tint('2015-10-16T10:32:30.00Z/2015-10-16T10:34:10.00Z'); % magnetosphere-magnetosheath-magnetosphere
ic = 1:4;

% Load defatt, for coordinate tranformation
disp('Loading defatt...')
load /Users/Cecilia/Data/MMS/2015Oct16/defatt.mat

% Magnetic field
disp('Loading magnetic field...')
c_eval('tic; dmpaB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_dmpa'',tint); toc',ic);
c_eval('tic; gseB?=mms.db_get_ts(''mms?_dfg_brst_l2pre'',''mms?_dfg_brst_l2pre_gse'',tint); toc',ic);

% Electric field
disp('Loading electric field...')
Eoffs1 = 2.05; Eoffs2 = 2.58; Eoffs3 = 2.53; Eoffs4 = 1.40;
c_eval('tic; dslE?dce=mms.db_get_ts(''mms?_edp_brst_ql_dce'',''mms?_edp_dce_xyz_dsl'',tint); toc',ic);
c_eval('dslE?dce.data(:,1)=dslE?dce.data(:,1)-Eoffs?;',ic)
c_eval('gseE?dce = mms_dsl2gse(dslE?dce,defatt?);',ic)
c_eval('gseE? = gseE?dce;',ic)

% Spacecraft position
disp('Loading spacecraft position...')
R  = mms.get_data('R_gse',tint);
c_eval('gseR? = TSeries(R.time,R.gseR?,''to'',1);',ic);

% Electron velocities 

%% Check rotation of electron frame electric field, following Schudder2015
c_eval('gseVexB? = gseVe?.cross(gseB?).resample(gseVe?.time)*1e-3; gseVexB?.units = ''mV/m'';',ic)
c_eval('gseRe? = gseE?+gseVexB?.resample(gseVe.time);',ic)
c_eval('gseR?brsttime = gseR?.resample(gseVe?);',1:4)

[curlRe,divRe,Re,~,~,~] = c_4_j('gseR?brsttime','gseRe?');
curlRe = irf.ts_vec_xyz(curlRe.time,curlRe.data);
curlRe.data = curlRe.data; curlRe.units = 'mV/m/km';
%Jcurl.time = EpochTT(Jcurl.time); Jcurl.name = '4sc current density';


