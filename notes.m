%% all
set(gcf,'PaperPositionMode','auto')
print -dpng delme.png
%% web address
%https://launchpad.net/irfu-matlab
%http://caa.estec.esa.int/caa/home.xml
%http://www.cluster.rl.ac.uk/csdsweb-cgi/csdsweb_pick

%% plott various
tint=[toepoch([2007 6 18 13 19 0]) toepoch([2007 6 18 13 21 0])];

h=irf_plot(3); 
irf_plot(h(1),'B_vec_xyz_gse__C4_CP_FGM_FULL')
irf_legend(h(1),'C4',[0.02 0.02]);
irf_plot(h(2),'E_Vec_xy_ISR2__C4_CP_EFW_L2_E')
irf_plot(h(3),'Spacecraft_potential__C4_CP_EFW_L2_P')

irf_zoom(tint,'x',h);

%% Loading E and B
cd /home/cecilia/data/20070617
caa_load
%%
[caaB3,~,B3]=c_caa_var_get('B_vec_xyz_gse__C3_CP_FGM_FULL');
[caaB4,~,B4]=c_caa_var_get('B_vec_xyz_gse__C4_CP_FGM_FULL');
[caaE3,~,E3]=c_caa_var_get('E_Vec_xy_ISR2__C3_CP_EFW_L2_E');
[caaE4,~,E4]=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_E');

%% compare electric field C3/C4 xvsx yvsy zvsz
figure(1)
h=irf_plot(3);

irf_plot(h(1),E3(:,[1 2]),'g'); hold(h(1),'on');
irf_plot(h(1),E4(:,[1 3])); hold(h(1),'on');

irf_plot(h(2),E3(:,[1 2]),'g'); hold(h(2),'on');
irf_plot(h(2),E4(:,[1 3])); hold(h(2),'on');

irf_plot(h(3),E3(:,[1 2]),'g'); hold(h(3),'on');
irf_plot(h(3),E4(:,[1 3])); hold(h(3),'on');

irf_legend(h(1),{'C3','C4'},[0.02 0.05]);
irf_legend(h(2),{'C3','C4'},[0.02 0.05]);
irf_legend(h(2),{'C3','C4'},[0.02 0.05]);

ylabel(h(1),'Ex [mV/m] ISR2');
ylabel(h(2),'Ex [mV/m] ISR2');
ylabel(h(3),'Ex [mV/m] ISR2');

irf_zoom(tint,'x',h);

%% E3D is empty
figure(2)
[caaE3D,~,E3D]=c_caa_var_get('C1_CP_EFW_L2_E3D_INERT');
h=irf_plot(3);
irf_plot(h(1),E3D)

%% plot E for all four spacecrafts
figure
[caaE1,~,E1]=c_caa_var_get('E_Vec_xy_ISR2__C1_CP_EFW_L2_E');
[caaE2,~,E2]=c_caa_var_get('E_Vec_xy_ISR2__C2_CP_EFW_L2_E');
[caaE3,~,E3]=c_caa_var_get('E_Vec_xy_ISR2__C3_CP_EFW_L2_E');
[caaE4,~,E4]=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_E');
h=irf_plot(4);
irf_plot(h(1),E1(:,[1 2 3]));
irf_plot(h(2),E2(:,[1 2 3]));
irf_plot(h(3),E3(:,[1 2 3]));
irf_plot(h(4),E4(:,[1 2 3]));
ylabel(h(1),'C1 [mV/m]');
ylabel(h(2),'C2 [mV/m]');
ylabel(h(3),'C3 [mV/m]');
ylabel(h(4),'C4 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),{'x','y'},[0.02 0.05]);
irf_legend(h(3),{'x','y'},[0.02 0.05]);
irf_legend(h(4),{'x','y'},[0.02 0.05]);

%%
[caaR1,~,R1]=c_caa_var_get('sc_r_xyz_gse__C1_CP_AUX_POSGSE_1M');
[caaR2,~,R2]=c_caa_var_get('sc_r_xyz_gse__C2_CP_AUX_POSGSE_1M');
[caaR3,~,R3]=c_caa_var_get('sc_r_xyz_gse__C3_CP_AUX_POSGSE_1M');
[caaR4,~,R4]=c_caa_var_get('sc_r_xyz_gse__C4_CP_AUX_POSGSE_1M');


% figure
% c_pl_sc_conf_xyz([2007 6 18 13 20 00])

%% plot B for sp1 sp2 sp3 sp4 spacecrafts
figure
[caaB1,~,B1]=c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_FULL');
[caaB2,~,B2]=c_caa_var_get('B_vec_xyz_gse__C2_CP_FGM_FULL');
[caaB3,~,B3]=c_caa_var_get('B_vec_xyz_gse__C3_CP_FGM_FULL');
[caaB4,~,B4]=c_caa_var_get('B_vec_xyz_gse__C4_CP_FGM_FULL');
h=irf_plot(4);
irf_plot(h(1),B1(:,[1 2 3 4]));
irf_plot(h(2),B2(:,[1 2 3 4]));
irf_plot(h(3),B3(:,[1 2 3 4]));
irf_plot(h(4),B4(:,[1 2 3 4]));
ylabel(h(1),'C1 [nT]');
ylabel(h(2),'C2 [nT]');
ylabel(h(3),'C3 [nT]');
ylabel(h(4),'C4 [nT]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),{'x','y'},[0.02 0.05]);
irf_legend(h(3),{'x','y'},[0.02 0.05]);
irf_legend(h(4),{'x','y'},[0.02 0.05]);

%% compare magnetic field C3/C4
figure(2)
[caaB3,~,B3]=c_caa_var_get('B_vec_xyz_gse__C3_CP_FGM_FULL');
[caaB4,~,B4]=c_caa_var_get('B_vec_xyz_gse__C4_CP_FGM_FULL');
[caaE3,~,E3]=c_caa_var_get('E_Vec_xy_ISR2__C3_CP_EFW_L2_E');
[caaE4,~,E4]=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_E');
h=irf_plot(1);
irf_plot(h(1),{E3 E4},'comp')

%%
irf_plot(h(2),E3(:,[1 2]));
hold(h(2),'on');
irf_plot(h(2),E4(:,[1 2]),'g');
irf_legend(h(2),{'C3','C4'},[0.02 0.05]);
ylabel(h(2),'Ex [mV/m] ISR2');
irf_plot(h(3),'Spacecraft_potential__C4_CP_EFW_L2_P')
%%
h=irf_plot(1);
irf_plot(h(1),'B_vec_xyz_gse__C4_CP_FGM_FULL')

%%
figure(3)
h=irf_plot(3);
irf_plot(h(1),'B_vec_xyz_gse__C4_CP_FGM_FULL')
irf_plot(h(2),E3(:,[1 2]));
hold(h(2),'on');
irf_plot(h(2),E4(:,[1 2]),'g');
irf_legend(h(2),{'C3','C4'},[0.02 0.05]);
ylabel(h(2),'Ex [mV/m] ISR2');
irf_plot(h(3),'Spacecraft_potential__C4_CP_EFW_L2_P')
