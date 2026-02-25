% event20030903. From Zhou2009, to compare wavelengths
cd /Users/Cecilia/Data/Cluster/20030903/
tint = toepoch([2003 09 03 23 25 00; 2003 09 03 23 35 00])';

%% Download data
for kk=[7]%[1 2 4 5]
    dataproducts = {'C?_CP_FGM_FULL_ISR2',...
                    'C?_CP_EFW_L2_E3D_INERT',...
                    'C?_CP_PEA_PITCH_SPIN_DEFLUX',...
                    'C?_CP_STA_CWF_HBR_ISR2',...
                    'C?_CP_STA_CWF_ISR2',...
                    'C?_CP_STA_CWF_GSE',...
                    'C?_CP_PEA_MOMENTS',...
                    'C?_CP_CIS_HIA_ONBOARD_MOMENTS'};
    caa_download(tint,dataproducts{kk})
end

%% Load data
sclist=1:4;
c_eval('diE?inert = c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
c_eval('diE? = c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'',''mat'');',sclist);
c_eval('diB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_ISR2'',''mat'');',sclist);
c_eval('diB?fgm = c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);

c_eval('gseE?inert = c_caa_var_get(''E_Vec_xyz_GSE__C?_CP_EFW_L2_E3D_GSE'',''mat'');',sclist);
%c_eval('diE? = c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'',''mat'');',sclist);
c_eval('gseB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_GSE'',''mat'');',sclist);
c_eval('gseB?fgm = c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);
%c_eval('gseB?fgm=irf_abs(gseB?fgm);',sclist)

c_eval('gseB?=c_fgm_staff_combine(gseB?fgm(:,1:4),gseB?sta);',sclist)
c_eval('gseB?=irf_abs(gseB?);',sclist)
c_eval('gsmB?=c_coord_trans(''GSE'',''GSM'',gseB?,''cl_id'',?);',sclist);

c_eval('gseR? = c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);

c_eval('[diE?,ang?]=irf_edb(diE?,diB?fgm);',sclist);
c_eval('gseE?=c_coord_trans(''ISR2'',''GSE'',diE?,''cl_id'',?);',sclist);
   