%% cis plot.

h=irf_plot(5);

irf_plot(h(1),'velocity_isr2__C3_CP_CIS_HIA_ONBOARD_MOMENTS');
irf_plot(h(2),'flux__C3_CP_CIS_HIA_HS_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
irf_plot(h(3),diE3);
irf_plot(h(4),'velocity__C3_CP_CIS_CODIF_HS_H1_MOMENTS');
irf_plot(h(5),'flux__C3_CP_CIS_CODIF_H1_1D_PEF','colorbarlabel','log_{10} dEF\newline keV/cm^2 s sr keV','fitcolorbarlabel');
irf_plot_axis_align;

irf_zoom(h,'x',toepoch([2007 08 31 9 45 00;2007 08 31 12 15 00])')

%% check if it seems reasonable with polar plots
product = 'C4_CP_CIS-CODIF_HS_H1_PSD';

codif_pef = c_caa_distribution_data(product)