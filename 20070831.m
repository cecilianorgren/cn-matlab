%% 2007 08 31 10:15-10:25

%% Specifying time interval used for downloading data
tint=[toepoch([2007 8 31 10 15 0]) toepoch([2007 8 31 10 25 0])];
%%
% caa_download(tint,'C3_CP_FGM_FULL');
% caa_download(tint,'C4_CP_FGM_FULL');
% caa_download(tint,'C3_CP_EFW_L2_E3D_INERT');
% caa_download(tint,'C4_CP_EFW_L2_E3D_INERT');

%% Loading data
cd /home/cecilia/data/20070831/
caa_load


%% Loading E and B
[caaB3,~,B3]=c_caa_var_get('B_vec_xyz_gse__C3_CP_FGM_FULL');
[caaB4,~,B4]=c_caa_var_get('B_vec_xyz_gse__C4_CP_FGM_FULL');

[caaE3,~,E3]=c_caa_var_get('E_Vec_xyz_ISR2__C3_CP_EFW_L2_E3D_INERT');
[caaE4,~,E4]=c_caa_var_get('E_Vec_xyz_ISR2__C4_CP_EFW_L2_E3D_INERT');

%% Loading ion and electron data
[caaEi2,~,Ei2]=c_caa_var_get('energy_table__C2_CP_CIS_CODIF_H1_1D_PEF');
[caaEi3,~,Ei3]=c_caa_var_get('energy_table__C3c_CP_CIS_CODIF_H1_1D_PEF');
[caaEi4,~,Ei4]=c_caa_var_get('energy_table__C4_CP_CIS_CODIF_H1_1D_PEF');

[caaEe2,~,Ee2]=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_DEFlux');
[caaEe3,~,Ee3]=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_DEFlux');
[caaEe4,~,Ee4]=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_DEFlux');

%% Plotting E
figure(1)
h=irf_plot(2);
irf_plot(h(1),E3(:,[1 2 3]));
irf_plot(h(2),E4(:,[1 2 3]));
ylabel(h(1),'C3 [mV/m]');
ylabel(h(2),'C4 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),{'x','y'},[0.02 0.05]);
irf_zoom(tint,'x',h);
add_timeaxis(h);

%% Sampling frequency fs
dt=E3(2,1)-E3(1,1);
fs=1/dt;

%% Getting frequency spectrum with no overlap
overlap=0;
%E2fft=caa_powerfft(E2,1024,450,0);
E3fft=caa_powerfft(E3,512,fs,overlap);
E4fft=caa_powerfft(E4,512,fs,overlap);

%% Getting frequency spectrum with overlap
overlap=50;
%E2fft_olap=caa_powerfft(E2,1024,450,overlap);
E3fft_olap=caa_powerfft(E3,512,fs,overlap);
E4fft_olap=caa_powerfft(E4,512,fs,overlap);

% Det blir skillnad i tidsskalan, med overlap blir antalet frekvenser f?rre
% och f?rskjutna

%% Plotting E with spectrogram and printing to file
figure(2)
h=irf_plot(5);
irf_plot(h(1),E3(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E3fft);
caa_spectrogram([h(4) h(5)],E3fft_olap);
ylabel(h(1),'E [mV/m]');
irf_legend(h(1),{'E_x' 'E_y'},[0.02 0.05]);
irf_legend(h(2),'E3_x 0% overlap',[0.02 0.05]);
irf_legend(h(3),'E3_y 0% overlap',[0.02 0.05]);
irf_legend(h(4),'E3_x 50% overlap',[0.02 0.05]);
irf_legend(h(5),'E3_y 50% overlap',[0.02 0.05]);
irf_zoom(tint,'x',h);
add_timeaxis(h);
print -dpng 20070831_E3_olaps.png

%% Plotting E with spectrogram and printing to file
figure(3)
h=irf_plot(3);
irf_plot(h(1),E4(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E4fft);
%caa_spectrogram([h(4) h(5)],E3fft_olap);
ylabel(h(1),'C1 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),'E4_x',[0.02 0.05]);
irf_legend(h(3),'E4_y',[0.02 0.05]);
%irf_legend(h(4),'E3_x 50% overlap',[0.02 0.05]);
%irf_legend(h(5),'E3_y 50% overlap',[0.02 0.05]);
irf_zoom(tint,'x',h);
add_timeaxis(h);
print -dpng 20070831_E4ps.png


%% Plotting E2 & E4 with spectrogram and printing to file
figure(4)

h=irf_plot(6);
irf_plot(h(1),E3(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E3fft);
irf_plot(h(4),E4(:,[1 2 3]));
caa_spectrogram([h(5) h(6)],E4fft);

ylabel(h(1),'C3 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),'E3_x',[0.02 0.05]);
irf_legend(h(3),'E3_y',[0.02 0.05]);
ylabel(h(4),'C4 [mV/m]');
irf_legend(h(4),{'x','y'},[0.02 0.05]);
irf_legend(h(5),'E4_x',[0.02 0.05]);
irf_legend(h(6),'E4_y',[0.02 0.05]);
irf_zoom(tint,'x',h);
add_timeaxis(h);
print -dpng 20070831_E34ps_0olap.png

%% Getting spectrogram with wavelets

dumwaveplot(E3(:,1),E3(:,[2 3]));
%% 
dumwaveplot(E3(:,1),E3(:,2)); %  x-component
%% Wavelet plots of E, B and S, needs B_0 - ambient magnetic field
irf_pl_ebs(E3,B3,B3)

%% Ion and electron spectrum
% There is no ion data for this time

[caaEl3,~,El3]=c_caa_var_get('Data__C3_CP_PEA_PITCH_SPIN_DEFlux');
El3data=El3.data;
%%
[caaEle3,~,Ele3]=c_caa_var_get('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_DEFlux');
[caaElen3,~,Elen3]=('Sweep_Energy__C3_CP_PEA_PITCH_SPIN_DEFlux');
%%
h=irf_plot(1);
irf_plot(h(1),El3.t,El3.data(:,1,1);
%%
pcolor(El3.data(1,:,:));
%%
El3.data(:,:,3);
pcolor(El3.data(:,:,3));
%%






