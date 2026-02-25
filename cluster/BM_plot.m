%% Going to the right directory
% 2007****
% (0712,0715,0724,0731,0803,0809,0814,0821,0831,0902,0909,0926,1008,1017)
date= '20070926';
cd /home/cecilia/data/BM/20070926/
%%
caa_load
%% ge
tint_event=tint(12,:);
%% Loading into variables
c_eval('[caaB?,~,gseB?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');',3:4);
c_eval('[caaE?,~,diE?]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');',3:4);
c_eval('[caaP?,~,P?]=c_caa_var_get(''Spacecraft_potential__C3_CP_EFW_L2_P'');',3:4);

%% Convert diE to gseE 
%c_eval('gseE?=c_coord_trans(''dsi'',''gse'',diE?,''CL_ID'',?);',3:4);

%% Sampling frequency fs
dt=diE3(2,1)-diE3(1,1);
fs=1/dt;

%% Getting frequency spectrum with no overlap
if 0
overlap=0;
c_eval('diE?fft=caa_powerfft(diE?,1024,fs,overlap)',3:4)
end
%% Plot
h=irf_plot(6);

irf_plot(h(1),gseB3);
ylabel(h(1),'M_{GSE} [nT]');
irf_legend(h(1),{'x','y','z'},[0.02 0.05]);

irf_plot(h(2),diE3);
ylabel(h(2),'E_{ISR2} [mV/m]');
irf_legend(h(2),{'x','y','z'},[0.02 0.05]);

irf_plot(h(3),diE4);
ylabel(h(3),'E_{ISR2} [mV/m]');
irf_legend(h(3),{'x','y','z'},[0.02 0.05]);

irf_plot(h(4),P3);
ylabel(h(4),'V_{SC} [V]');
irf_legend(h(4),'C3',[0.02 0.05]);

irf_plot(h(5),'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel',...
             'log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
               ylabel(h(5),'E [eV]');
               
irf_legend(h(5),{['C' num2str(3)]},[0.98 0.05],'color','c')   
set(h(5),'yscale','log');
set(h(5),'ytick',[1 1e1 1e2 1e3 1e4 1e5]);
         % Insert a colorbar and moves it to the right of the plot. 'peer' is
% mandatory, h1(4) is the plot handle, 'EastOutside' specifies the position
% relative to the plot.

irf_plot(h(6),'Data__C3_CP_PEA_PITCH_SPIN_DEFlux','sum_dim1','colorbarlabel',...
             'log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
         irf_legend(h(6),{['C' num2str(3)]},[0.98 0.05],'color','c')  
         set(h(6),'yscale','log','ylim',[100 3e4]);
         set(h(6),'ytick',[1 1e1 1e2 1e3 1e4 1e5]);
         % Insert a colorbar and moves it to the right of the plot. 'peer' is
% mandatory, h1(4) is the plot handle, 'EastOutside' specifies the position relative to the plot.
% colb2 = colorbar('peer',h(6),'EastOutside');
% cpos2 = get(colb2,'Position');
% % Add a small step to the right.
% cpos2(1) = cpos2(1) +.092;
% set(colb2,'Position',cpos2);
         
 
irf_zoom(h,'x',tint_event);
irf_plot_axis_align
add_timeaxis(h(6),'usefig');

%% Printing to file
print -dpng overview_20070926.png

%% What we want to plot
% Spacecraft potential, B and abs(B), E s/c 3 s/c 4