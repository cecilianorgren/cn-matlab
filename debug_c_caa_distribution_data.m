cd /Users/Cecilia/Data/Cluster/20070831;
t1 = [2007 08 31 10 10 00]; t2=[2007 08 31 10 20 00];
tint = toepoch([t1;t2])';
%caa_download(tint,'C3_CP_PEA_3DXPH_PSD')
%caa_download(tint,'C3_CP_FGM_FULL')
%caa_download(tint,'CL_SP_AUX')
%%
data_structure = c_caa_distribution_data('C3_CP_PEA_3DXPH_PSD');
%%
h = c_caa_plot_distribution_function('tint',tint,'polar',data_structure);

%%
h = c_caa_plot_distribution_function('tint',tint,'cross-section',data_structure);

%%

cd /Users/Cecilia/Data/Cluster/20080422;

tint = toepoch([2008 04 22 17 55 00.0;2008 04 22 18 10 00.0])';
%C3_CP_CIS-HIA_HS_1D_PEF

data_structure = c_caa_distribution_data('C3_CP_PEA_PITCH_SPIN_PSD');
%%
tint = toepoch([2008 04 22 17 55 00.0;2008 04 22 17 55 50.0])';
h = c_caa_plot_distribution_function('tint',tint,'cross-section',data_structure);

%%
h = c_caa_plot_distribution_function('tint',tint,'cross-section',data_structure);

%% OOOOOOOOOOBBBBBBBSSSS
% Slighty change version concerning HIA ions
% c_caa_distribution_data_debugged.
% TODO : test and commit
