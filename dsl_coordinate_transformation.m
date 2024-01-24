ic = [3];
tint = irf.tint('2017-07-11T22:31:00.00Z/2017-07-11T22:37:20.00Z'); %20151112071854

% Load datastore
%mms.db_init('local_file_db','/Volumes/Nexus/data');
%mms.db_init('local_file_db','/Volumes/Fountain/Data/MMS');
mms.db_init('local_file_db','/Users/cecilia/Data/MMS');
%mms.db_init('local_file_db','/Users/cno062/Data/MMS');
mms.db_init('local_file_db','/Volumes/mms');
db_info = datastore('mms_db');

%% Load data
c_eval('dmpaB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_dmpa_brst_l2'',tint);',ic);
c_eval('gseB? = mms.db_get_ts(''mms?_fgm_brst_l2'',''mms?_fgm_b_gse_brst_l2'',tint);',ic);

c_eval('dslE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_dsl_brst_l2'',tint);',ic);
c_eval('gseE? = mms.db_get_ts(''mms?_edp_brst_l2_dce'',''mms?_edp_dce_gse_brst_l2'',tint);',ic);

c_eval('dbcsVe? = mms.get_data(''Ve_dbcs_fpi_brst_l2'',tint,?);',ic)
c_eval('gseVe? = mms.get_data(''Ve_gse_fpi_brst_l2'',tint,?);',ic)

c_eval('zra = mms.db_get_variable(''mms?_ancillary_defatt'',''zra'',tint);',ic)
c_eval('zdec = mms.db_get_variable(''mms?_ancillary_defatt'',''zdec'',tint);',ic)
defatt = zra; defatt.zdec = zdec.zdec;

%% Rotate things into new coordinate system 
c_eval('x?_dsl = irf.ts_vec_xyz(gseE?.time,repmat([1 0 0],gseE?.length,1));',ic)
c_eval('y?_dsl = irf.ts_vec_xyz(gseE?.time,repmat([0 1 0],gseE?.length,1));',ic)
c_eval('z?_dsl = irf.ts_vec_xyz(gseE?.time,repmat([0 0 1],gseE?.length,1));',ic) % == spacecraft spin axis
c_eval('x_dsl?_gse = mms_dsl2gse(x?_dsl,defatt,1);',ic)
c_eval('y_dsl?_gse = mms_dsl2gse(y?_dsl,defatt,1);',ic)
c_eval('z_dsl?_gse = mms_dsl2gse(z?_dsl,defatt,1);',ic) % spacecraft spin axis in gse

if 0
  %%
  hca = subplot(1,1,1);
  for it = 1:2000:z_dsl3_gse.length
    quiver3(hca,0,0,0,x_dsl3_gse.data(it,1),x_dsl3_gse.data(it,2),x_dsl3_gse.data(it,3),0)
    hold(hca,'on')
    quiver3(hca,0,0,0,y_dsl3_gse.data(it,1),y_dsl3_gse.data(it,2),y_dsl3_gse.data(it,3),0)
    quiver3(hca,0,0,0,z_dsl3_gse.data(it,1),z_dsl3_gse.data(it,2),z_dsl3_gse.data(it,3),0)
    
    hca.XLim = [-1 1];
    hca.YLim = [-1 1];
    hca.ZLim = [-1 1];
    
    hold(hca,'off')
    pause(0.01)
  end
  
end

e1 = z3_dsl.norm.data; e2 = z_dsl3_gse.norm.data;
e1e2 = sum(e1.*e2,2);
ang = acosd(e1e2);
tsAng_spin_axis = irf.ts_scalar(gseE3.time,ang);

L_gse = [0.9482,-0.255,-0.1893]; % LMN in GSE
M_gse = [0.1818,0.9245,-0.3350]; % LMN in GSE
N_gse = [0.2604,0.2832,0.9230];  % LMN in GSE
lmn_gse = [L_gse;M_gse;N_gse];   % LMN in GSE

%tt = irf_time('2017-07-11T22:34:02.00Z','utc>EpochTT'); %20151112071854
%tint_defatt = tt + [-2 2];

tsL_gse = irf.ts_vec_xyz(tt,L_gse); % LMN in GSE
tsM_gse = irf.ts_vec_xyz(tt,M_gse); % LMN in GSE
tsN_gse = irf.ts_vec_xyz(tt,N_gse); % LMN in GSE

L_dsl = mms_dsl2gse(tsL_gse, defatt, -1); L_dsl = L_dsl.data; % LMN in DSL
M_dsl = mms_dsl2gse(tsM_gse, defatt, -1); M_dsl = M_dsl.data; % LMN in DSL
N_dsl = mms_dsl2gse(tsN_gse, defatt, -1); N_dsl = N_dsl.data; % LMN in DSL
lmn_dsl = [L_dsl; M_dsl; N_dsl];                              % LMN in DSL

% DSL from GSE
c_eval('gseE?_from_dsl = mms_dsl2gse(dslE?,defatt,1); gseE?_from_dsl.name = ''E gse from dsl'';',ic)

% LMN in GSE
c_eval('mvaB?_gse = gseB?*lmn_gse''; mvaB?.name = ''B LMN'';',ic)
c_eval('mvaE?_gse = gseE?*lmn_gse''; mvaE?.name = ''E LMN'';',ic)
c_eval('mvaVe?_gse = gseVe?*lmn_gse''; mvaVe?.name = ''Ve LMN'';',ic)

 % LMN in DSL, is this wrong?
c_eval('mvaVe?_dsl = mms_dsl2gse(mvaVe?_gse,defatt,-1); mvaVe?_dsl.name = ''Ve LMN dsl'';',ic)
c_eval('mvaE?_dsl  = mms_dsl2gse(mvaE?_gse,defatt,-1); mvaE?_dsl.name = ''E LMN dsl'';',ic)
c_eval('mvaB?_dsl  = mms_dsl2gse(mvaB?_gse,defatt,-1); mvaB?_dsl.name = ''B LMN dsl'';',ic)

% LMN in DSL, is this right?
c_eval('mvaE?_dsl2 = dslE?*lmn_dsl''; mvaE?_dsl2.name = ''E LMN DSL2'';',ic)
c_eval('mvaB?_dsl2 = dmpaB?*lmn_dsl''; mvaB?_dsl2.name = ''B LMN DSL2'';',ic)
c_eval('mvaVe?_dsl2 = dbcsVe?*lmn_dsl''; mvaVe?_dsl2.name = ''Ve LMN DSL2'';',ic)

%%
e1 = gseE3.norm.data; e2 = dslE3.norm.data;
e1e2 = sum(e1.*e2,2);
ang = acosd(e1e2);
tsAng_gseEdslE = irf.ts_scalar(gseE3.time,ang);

% Sanity check: comparing gseE_from_data_file and gseE = mms_dsl2gse(dslE_from_data_file,defatt,1).
e1 = gseE3_from_dsl.norm.data; e2 = gseE3.norm.data;
e1e2 = sum(e1.*e2,2);
ang = acosd(e1e2);
tsAng_from_cdf = irf.ts_scalar(gseE3.time,ang);

% Comparing (lmn) mvaE
e1 = mvaE3_dsl.norm.data; e2 = mvaE3_dsl2.norm.data;
e1e2 = sum(e1.*e2,2);
ang = acosd(e1e2);
tsAng_in_lmn = irf.ts_scalar(gseE3.time,ang);

h = irf_plot(4);
hca = irf_panel('E LMN from DSL x');
irf_plot(hca,{mvaE3_dsl.x,mvaE3_dsl2.x,mvaE3_dsl.x-1*mvaE3_dsl2.x},'comp')
hca = irf_panel('E LMN from DSL y');
irf_plot(hca,{mvaE3_dsl.y,mvaE3_dsl2.y,mvaE3_dsl.y-1*mvaE3_dsl2.y},'comp')
hca = irf_panel('E LMN from DSL z');
irf_plot(hca,{mvaE3_dsl.z,mvaE3_dsl2.z,mvaE3_dsl.z-1*mvaE3_dsl2.z},'comp')
hca = irf_panel('angle between Es');
irf_plot(hca,{tsAng,tsAng_in_lmn},'comp')


%hca = irf_panel('E diff');
%irf_plot({mvaE3_dsl-1*mvaE3_dsl2},'comp')

