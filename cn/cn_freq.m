%% Load data
cd /home/cecilia/data/20070617/
caa_load

%% Loading variables
[caaB1,~,B1]=c_caa_var_get('B_vec_xyz_gse__C1_CP_FGM_FULL');
[caaB2,~,B2]=c_caa_var_get('B_vec_xyz_gse__C2_CP_FGM_FULL');
[caaB3,~,B3]=c_caa_var_get('B_vec_xyz_gse__C3_CP_FGM_FULL');
[caaB4,~,B4]=c_caa_var_get('B_vec_xyz_gse__C4_CP_FGM_FULL');

[caaE1,~,E1]=c_caa_var_get('E_Vec_xy_ISR2__C1_CP_EFW_L2_E');
[caaE2,~,E2]=c_caa_var_get('E_Vec_xy_ISR2__C2_CP_EFW_L2_E');
[caaE3,~,E3]=c_caa_var_get('E_Vec_xy_ISR2__C3_CP_EFW_L2_E');
[caaE4,~,E4]=c_caa_var_get('E_Vec_xy_ISR2__C4_CP_EFW_L2_E');

%% Plotting E
figure(1)
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

%% Getting frequency spectrum with no overlap
E1fft=caa_powerfft(E1,1024,450,0);
E2fft=caa_powerfft(E2,1024,450,0);
E3fft=caa_powerfft(E3,1024,450,0);
E4fft=caa_powerfft(E4,1024,450,0);

%% Getting frequency spectrum with 50% overlap
E1fft_olap=caa_powerfft(E1,1024,450,20);
E2fft_olap=caa_powerfft(E2,1024,450,20);
E3fft_olap=caa_powerfft(E3,1024,450,20);
E4fft_olap=caa_powerfft(E4,1024,450,20);

%% Plotting E with spectrogram and printing to file
figure(2)
h=irf_plot(5);
irf_plot(h(1),E1(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E1fft);
caa_spectrogram([h(4) h(5)],E1fft_olap);
ylabel(h(1),'C1 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),'E1_x',[0.02 0.05]);
irf_legend(h(3),'E1_y',[0.02 0.05]);
irf_legend(h(4),'E1_x 50% overlap',[0.02 0.05]);
irf_legend(h(5),'E1_y 50% overlap',[0.02 0.05]);
print -dpng E1ps.png

%%
figure(3)
h=irf_plot(3);
irf_plot(h(1),E2(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E2fft);
ylabel(h(1),'C1 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),'E2_x',[0.02 0.05]);
irf_legend(h(3),'E2_y',[0.02 0.05]);
print -dpng E2ps.png

%%
figure(4)
h=irf_plot(3);
irf_plot(h(1),E3(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E3fft);
ylabel(h(1),'C1 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),'E3_x',[0.02 0.05]);
irf_legend(h(3),'E3_y',[0.02 0.05]);
print -dpng E3ps.png

%%
figure(5)
h=irf_plot(3);
irf_plot(h(1),E4(:,[1 2 3]));
caa_spectrogram([h(2) h(3)],E4fft);
ylabel(h(1),'C1 [mV/m]');
irf_legend(h(1),{'x','y'},[0.02 0.05]);
irf_legend(h(2),'E4_x',[0.02 0.05]);
irf_legend(h(3),'E4_y',[0.02 0.05]);
set(gcf,'PaperPositionMode','auto')
print -dpng E4ps.png

%%





%%
a=power{1,1};
b=power{1,2};
fftmatrix0=zeros(size(a,2),size(a,1));
for k=1:size(t,1)  
    fftmatrix0(:,k)=a(k,:);
end
fftmatrix=flipdim(fftmatrix0,1);
figure(2)
imagesc(fftmatrix,[0 0.5])
title('Power spectrum')
ts=datestr(datenum(fromepoch(fix((max(t)-min(t))/2))));
xlabel([])
ylabel('f')
colorbar

%%
figure
plot(f,b(1,:))

%%
h=irf_plot(1);
irf_plot(h(1),[t power])

%%
c=[t power];
size(c)

%%
fftE1=fft(E1(:,[1 2 3]));
h=irf_plot(2);
irf_plot(h(1),fftE1);