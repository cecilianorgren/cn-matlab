% Get hia data for 2007 08 31
cd /Users/Cecilia/Dropbox' (IRFU)'/Exjobb2011/20070831_1000-1030/

%% Magnetic field
c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
c_eval('gsmB?fgm=c_coord_trans(''gse'',''GSM'',gseB?fgm,''cl_id'',?);',3:4);

%% Spacecraft potential
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',3:4);

%% load ion velocitites hia
sc=[1 3];
c_eval('caa_load(''C?_CP_CIS-HIA_ONBOARD_MOMENTS'');',sc);
c_eval('gsehiaVi?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sc);
c_eval('gsmhiaVi?=c_coord_trans(''gse'',''gsm'',gsehiaVi?,''cl_id'',?);',sc);

c_eval('gsehiaVi?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sc);
c_eval('hiaTi?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''temperature__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sc);
c_eval('hiaNi?=getmat(C?_CP_CIS_HIA_ONBOARD_MOMENTS,''density__C?_CP_CIS_HIA_ONBOARD_MOMENTS'');',sc);

%% Electron moments
% K 10^6
c_eval('parTe?=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
c_eval('perTe?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_'',''mat'');',3:4);
c_eval('Te?=[parTe?(:,1) (parTe?(:,2)+2*perTe?(:,2))/3];',3:4);
MK2eV = 8.621738*1e-5*1e6;
c_eval('eVTe?=[Te?(:,1) Te?(:,2)*MK2eV];',3:4)

c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',3:4);

c_eval('peaVe?gse=c_caa_var_get(''Data_Velocity_GSE__C?_CP_PEA_MOMENTS'',''mat'');',3:4);

%% Electric field
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',3:4);
c_eval('gseE?=c_coord_trans(''dsi'',''gse'',diE?,''CL_ID'',?);',3:4);
c_eval('gsmE?=c_coord_trans(''gse'',''gsm'',gseE?,''CL_ID'',?);',3:4);
c_eval('gseE?slow=irf_resamp(gseE?,gsehiaVi?);',3);
c_eval('gseE?slow=irf_resamp(gseE?,gsehiaVi!);',4,3);

%% ExB drift
c_eval('gseB?slow=irf_resamp(gseB?fgm,gsehiaVi!);',3:4,3);
c_eval('gseVExB? = irf_cross(gseE?slow,gseB?slow);',3:4)
%c_eval('gseVExB?_ = irf_cross(gseE?slow,irf_resampgseB?fgm);',3:4)
c_eval('aa=irf_abs(gseB?slow); gseVExB?(:,2:4) = gseVExB?(:,2:4)./[aa(:,5).^2 aa(:,5).^2 aa(:,5).^2]*1e3;',3:4)

%% compare velocities

irf_plot({gseB3fgm,gsehiaVi3,codVi3,gseVExB3})
irf_plot({gsehiaVi3,codVi3,gseVExB3},'comp')
irf_plot({gsehiaVi3,gseVExB3},'comp')

%%

h = irf_plot(9);
isub = 1;

if 1 % Magnetic field
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,gseB3fgm)
  hca.YLabel.String = {'B_{GSE}','(nT)'};
  irf_legend(hca,{'x','y','z'},[0.95,0.95])
end
if 1 % Density
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{hiaNi3,peaNe3},'comp')
  hca.YLabel.String = {'n','(cc)'};
  irf_legend(hca,{'n_i','n_e'},[0.95,0.95])
end
if 1 % Spacecraft potential
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,P3)
  hca.YLabel.String = {'-V_{sc}','(V)'};
end
if 1 % Electric field
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,gseE3)
  hca.YLabel.String = {'E_{GSE}','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.95,0.95])
end
if 1 % Electric field, slow
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,gseE3slow)
  hca.YLabel.String = {'E_{GSE}','(mV/m)'};
  irf_legend(hca,{'x','y','z'},[0.95,0.95])
end
if 0 % Ion velocities comparison, x 
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{gsehiaVi3(:,[1 2]),codVi3(:,[1 2]),gseVExB3(:,[1 2])},'comp')
  hca.YLabel.String = {'V_{i,x,GSE}','(km/s)'};
  irf_legend(hca,{'hia','codif','ExB'},[0.95,0.95])
end
if 0 % Ion velocities comparison, y 
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{gsehiaVi3(:,[1 3]),codVi3(:,[1 3]),gseVExB3(:,[1 3])},'comp')
  hca.YLabel.String = {'V_{i,y,GSE}','(km/s)'};
  irf_legend(hca,{'hia','codif','ExB'},[0.95,0.95])
end
if 0 % Ion velocities comparison, z 
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{gsehiaVi3(:,[1 4]),codVi3(:,[1 4]),gseVExB3(:,[1 4])},'comp')
  hca.YLabel.String = {'V_{i,z,GSE}','(km/s)'};
  irf_legend(hca,{'hia','codif','ExB'},[0.95,0.95])
end
if 1 % Ion velocities comparison, x, no codif
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{gsehiaVi3(:,[1 2]),gseVExB3(:,[1 2])},'comp')
  hca.YLabel.String = {'V_{i,x,GSE}','(km/s)'};
  irf_legend(hca,{'hia','ExB'},[0.95,0.95])
end
if 1 % Ion velocities comparison, y, no codif
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{gsehiaVi3(:,[1 3]),gseVExB3(:,[1 3])},'comp')
  hca.YLabel.String = {'V_{i,y,GSE}','(km/s)'};
  irf_legend(hca,{'hia','ExB'},[0.95,0.95])
end
if 1 % Ion velocities comparison, z, no codif
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,{gsehiaVi3(:,[1 4]),gseVExB3(:,[1 4])},'comp')
  hca.YLabel.String = {'V_{i,z,GSE}','(km/s)'};
  irf_legend(hca,{'hia','ExB'},[0.95,0.95])
end
if 1 % Electron velocity
  hca = h(isub); isub = isub + 1;
  irf_plot(hca,peaVe3gse)
  hca.YLabel.String = {'v_{e,GSE}','(km/s)'};
  irf_legend(hca,{'x','y','z'},[0.95,0.95])
end