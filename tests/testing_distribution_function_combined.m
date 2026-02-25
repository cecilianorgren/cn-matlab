% testing cn_construct_distribution_data_combined.m and 
% cn_construct_distribution_data_combined.m
% See also cn_construct_distribution_data_combined.m, 
%          cn_plot_distribution_function_combined.m
%          f_sphere.m - to be implemented
sclist=4;
%% Construct data
%c_eval('e_psd?=c_caa_distribution_data(''C?_CP_CIS_HIA_HS_MAG_IONS_PSD'');',sclist);
%c_eval('hia_hs_psd?=c_caa_distribution_data(''C?_CP_CIS_HIA_HS_MAG_IONS_PSD'');',3);
%c_eval('rap_pad?=c_caa_distribution_data(''C?_CP_RAP_PAD_E3DD'');',sclist);
%% Plot data for defined time interval
t1=[2007 08 31 10 19 0.00]; 
t2=[2007 08 31 10 19 020.00]; % lhdw
tint=toepoch([t1; t2]);
%% Plot many in subplots with dt defined
for k=1;
    h=subplot(1,2,1);
    c_eval(['h',num2str(1),'=c_caa_plot_distribution_function(h,''tint'',tint+k,''cross-section'',e_psd?);'],sclist)
    h=subplot(1,2,2);
    c_eval(['h',num2str(2),'=c_caa_plot_distribution_function(h,''tint'',tint+k,''polar'',e_psd?);'],sclist)
end
%% compare with c_caa_distribution_function.m
% from help
rap_pad4=c_caa_distribution_data('C4_CP_RAP_PAD_E3DD');
hia_hs_psd3=c_caa_distribution_data('C3_CP_CIS_HIA_HS_MAG_IONS_PSD');;
pea_3dxph_psd4=c_caa_distribution_data('C4_CP_PEA_3DXPH_PSD');
pea_3dxh_psd4=c_caa_distribution_data('C4_CP_PEA_PITCH_3DXH_PSD');
pea_pitch_full_heea_psd3=c_caa_distribution_data('C3_CP_PEA_PITCH_FULL_PSD','HEEA');
pea_pitch_full_leea_psd3=c_caa_distribution_data('C3_CP_PEA_PITCH_FULL_PSD','LEEA');
pea_3drl_psd3=c_caa_distribution_data('C3_CP_PEA_PITCH_3DRL_PSD');
codif_hs_psd3=c_caa_distribution_data('C3_CP_CIS_CODIF_HS_H1_PSD');
%% C4_CP_PEA_PITCH_3DXH_PSD, o polar, x cross
close
%pea_3dxh_psd4=c_caa_distribution_data('C4_CP_PEA_PITCH_3DXH_PSD');
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',pea_3dxh_psd4,pea_3dxph_psd4);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',pea_3dxh_psd4);
c_caa_distribution_function(h(3),'C4_CP_PEA_PITCH_3DXH_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C4_CP_PEA_PITCH_3DXH_PSD','tint',tint,'cross-section')
%% C3_CP_PEA_PITCH_3DRL_PSD, o polar, x cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',pea_3drl_psd3);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',pea_3drl_psd3);
c_caa_distribution_function(h(3),'C3_CP_PEA_PITCH_3DRL_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C3_CP_PEA_PITCH_3DRL_PSD','tint',tint,'cross-section')
%% C3_CP_PEA_PITCH_FULL_PSD, HEEA, x polar, x cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',pea_pitch_full_heea_psd3);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',pea_pitch_full_heea_psd3);
%c_caa_distribution_function(h(3),'C3_CP_PEA_PITCH_FULL_PSD','tint',tint,'polar')
%c_caa_distribution_function(h(4),'C3_CP_PEA_PITCH_FULL_PSD','tint',tint,'cross-section')
%% C3_CP_PEA_PITCH_FULL_PSD, LEEA, x polar, x cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',pea_pitch_full_leea_psd3);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',pea_pitch_full_leea_psd3);
c_caa_distribution_function(h(3),'C3_CP_PEA_PITCH_FULL_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C3_CP_PEA_PITCH_FULL_PSD','tint',tint,'cross-section')
%% C3_CP_CIS_CODIF_HS_H1_PSD, o polar, o cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',codif_hs_psd3);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',codif_hs_psd3);
c_caa_distribution_function(h(3),'C3_CP_CIS_CODIF_HS_H1_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C3_CP_CIS_CODIF_HS_H1_PSD','tint',tint,'cross-section')
%% C3_CP_PEA_PITCH_3DRL_PSD, x polar, x cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',pea_3drl_psd3,codif_hs_psd3);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',pea_3drl_psd3);
c_caa_distribution_function(h(3),'C3_CP_PEA_PITCH_3DRL_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C3_CP_PEA_PITCH_3DRL_PSD','tint',tint,'cross-section')
%% C4_CP_RAP_PAD_E3DD, o polar, o cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',rap_pad4,rap_pad4,rap_pad4);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',rap_pad4,'pitchangle',[0 10 20 90]);
%c_caa_plot_distribution_function(h(4),'tint',tint,'cross-section',rap_pad4);
c_caa_distribution_function(h(3),'C4_CP_RAP_PAD_E3DD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C4_CP_RAP_PAD_E3DD','tint',tint,'cross-section')
%% C4_CP_PEA_3DXPH_PSD, o polar, o cross
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',pea_3dxph_psd4);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',pea_3dxph_psd4);
c_caa_distribution_function(h(3),'C4_CP_PEA_3DXPH_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C4_CP_PEA_3DXPH_PSD','tint',tint,'cross-section')
%% C3_CP_CIS_HIA_HS_MAG_IONS_PSD, o polar, o cross-section
close
for k=1:4; h(k)=subplot(2,2,k); end
c_caa_plot_distribution_function(h(1),'tint',tint,'polar',hia_hs_psd3);
c_caa_plot_distribution_function(h(2),'tint',tint,'cross-section',hia_hs_psd3);
c_caa_distribution_function(h(3),'C3_CP_CIS_HIA_HS_MAG_IONS_PSD','tint',tint,'polar')
c_caa_distribution_function(h(4),'C3_CP_CIS_HIA_HS_MAG_IONS_PSD','tint',tint,'cross-section')
%% general checks for both
% ticks f?r polar plots funkar inte... vet inte vad jag ?ndrat
% legend: med vanlig legend t'cker den f?r kurvan om l?dan ?r f?r liten.
% med irf_legend s? ?ker texten utanf?r om boxen ?r f?r liten
%% products to check
% polar
% o C4_CP_PEA_3DXPH_PSD
% o C4_CP_CIS_HIA_HS_MAG_IONS_PSD
% o C4_CP_RAP_PAD_E3DD
% x C3_CP_PEA_PITCH_3DXH_PSD
% x C3_CP_PEA_PITCH_FULL_PSD
% o C3_CP_PEA_PITCH_3DRL_PSD
% o C3_CP_CIS_CODIF_HS_H1_PSD
% x 
% x 
% x 
% x 
% x 
% x 


% cross-section
% o C4_CP_PEA_3DXPH_PSD
% o C4_CP_CIS_HIA_HS_MAG_IONS_PSD
% o C4_CP_RAP_PAD_E3DD
% x C3_CP_PEA_PITCH_3DXH_PSD
% x C3_CP_PEA_PITCH_FULL_PSD
% o C3_CP_PEA_PITCH_3DRL_PSD
% o C3_CP_CIS_CODIF_HS_H1_PSD
% x 
% x 
% x 
% x 
% x 
% x 
% x 