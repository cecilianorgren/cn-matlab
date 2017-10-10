cd /Users/Cecilia/data/BM/20070831
addpath /Users/Cecilia/cn-matlab/

load matlabE;
load matlabB;
load matlabBeta;
load matlabCIS;
load matlabBfft;
load matlabEfft;

%% Spacecraft potential
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',3:4);
%% Electron densities
c_eval('scpNe?=c_efw_scp2ne(P?);',3:4);
c_eval('scpNe?=irf_resamp(scpNe?,diE?);',3:4);
c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
c_eval('peaNe?hr=irf_resamp(peaNe?,diE?);',3:4);
%% HIA velocity for C3
c_eval('hiaVi?=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',3);
c_eval('gsmhiaVi?=c_coord_trans(''gse'',''gsm'',hiaVi?,''cl_id'',?);',3)
%% input gsmE4, output gsmE4lowres2
cn_resamp
%%
overlap=50;
nfft=512;
fs=450;
c_eval('gsmB?fft=irf_powerfft(gsmB?,nfft,fs,overlap);',3:4);
c_eval('diE?fft=irf_powerfft(diE?,nfft,fs,overlap);',4);
%% plot and see whats missing
cn_overview
irf_zoom(h,'x',tint2)
irf_plot_axis_align(h)
userdata=get(gcf,'userdata');
pos1=get(userdata.subplot_handles(1),'Position');
pos4=get(userdata.subplot_handles(4),'Position');
pos=pos1;pos(2)=pos4(2);
set(h2,'position',pos);
set(hcb,'xticklabel',[]);
set(h2,'xticklabel',[]);
