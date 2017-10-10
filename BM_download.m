%% from tin to epoch
interval=13;
%% FGM
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_FGM_FULL'');',3:4);
%% EFW
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_EFW_L2_E3D_INERT'');',3:4);
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_EFW_L2_P'');',3:4);
%% PEACE
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_PEA_PITCH_SPIN_DPFlux'');',3:4);
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_PEA_PITCH_SPIN_DEFlux'');',3:4);
%% CIS n?t fel med dessa, finns aldrig...
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_CIS-HIA_ONBOARD_MOMENTS'');',3:4);
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_CIS-CODIF_HS_H1_MOMENTS'');',3:4);
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_CIS-HIA_HS_1D_PEF'');',3:4);
c_eval('caa_download([tint(interval,1) tint(interval,2)],''C?_CP_CIS-CODIF_H1_1D_PEF'');',3:4);