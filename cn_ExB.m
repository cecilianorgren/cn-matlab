%% Load data
cd /home/cecilia/data/20070617/
caa_load

%% Get vectors (from C2)
[caagseExB2,~,gseExB2]=c_caa_var_get('v_drift_GSE__C2_CP_EFW_L2_V3D_GSE');
[caadiE2,~,diE2]=c_caa_var_get('E_Vec_xyz_ISR2__C2_CP_EFW_L2_E3D_INERT');
[caagseB2,~,gseB2]=c_caa_var_get('B_vec_xyz_gse__C2_CP_FGM_FULL');

%% Convert diE to gseE 
gseE2=c_coord_trans('dsi','gse',diE2,'CL_ID',2);

%% Calculate B^2
B2abs=irf_abs(gseB2,1);
B2pow2=B2abs.^2;
B2pow2=[ones(size(gseB2,1),1) repmat(B2pow2,[1 3])*1000]; %*1000 is for unities

%% Calculate ExB-drift
ExBmine=irf_cross(gseE2,gseB2./B2pow2);

%% Plot E, B, ExB vx ExB
figure(1)
h=irf_plot(4);

irf_plot(h(1),gseE2);
irf_plot(h(2),gseB2);
irf_plot(h(3),gseExB2);
irf_plot(h(4),ExBmine);

irf_legend(h(1),{'E_x','E_y','E_z'},[0.02 0.07]);
irf_legend(h(2),{'B_x','B_y','B_z'},[0.02 0.07]);
irf_legend(h(3),{'v_x','v_y','v_z'},[0.02 0.07]);
irf_legend(h(4),{'v_x','v_y','v_z'},[0.02 0.07]);

ylabel(h(1),'E GSE [mVm^{-1}]');
ylabel(h(2),'B GSE [nT]');
ylabel(h(3),'ExB GSE [kms^{-1}] from caa');
ylabel(h(4),'ExB GSE [kms^{-1}] my own');

%%
%
%figure(2)
%h=irf_plot(1);
%
%irf_plot({gseExB2 ExB},'comp');