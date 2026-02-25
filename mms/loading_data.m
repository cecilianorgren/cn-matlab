tic
if 1 % Loading data
%% Magnetic field
c_eval('gseB?fgm=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);
c_eval('gsmB?fgm=c_coord_trans(''gse'',''GSM'',gseB?fgm,''cl_id'',?);',3:4);
c_eval('gsmB?fgm=irf_abs(gsmB?fgm);',3:4)
c_eval('gseB?fgm=irf_abs(gseB?fgm);',3:4)
c_eval('diB?fgm=c_coord_trans(''GSM'',''DSI'',gsmB?fgm,''cl_id'',?);',3:4);
% Calculate and apply offset
offset=cn_offset(3,gsmB3fgm,4,gsmB4fgm,'gsm');  % whole BM interval
%offset2=cn_offset(...
%    3,irf_tlim(gsmB3fgm,toepoch([2007 08 31 10 10 0]),toepoch([2007 08 31 10 13 0])),...
%    4,irf_tlim(gsmB4fgm,toepoch([2007 08 31 10 10 0]),toepoch([2007 08 31 10 13 0])),...
%    'gsm');                                     % out in lobes
% gsmB3fgm=cn_apply_offset(3,gsmB3fgm,'gsm',offset2);

if 1, % read STAFF data form all sc, needs event defined
    switch event
        case '1'
            load mBS
            tint=[dBS3(1,1) dBS3(end,1)];
        case '2a'
            load mBS_20070902_1430-1440
            tint=[toepoch([2007 09 02 14 30 00]) toepoch([2007 09 02 14 40 00])];
        case '2b'
            load mBS_20070902_1545-1550
            tint=[toepoch([2007 09 02 15 45 00]) toepoch([2007 09 02 15 50 00])];
        case '3a'
            load mBS_20070926_0945-1000
            tint=[toepoch([2007 09 26 09 45 00]) toepoch([2007 09 26 10 00 00])];   
        case '3b'
            load mBS_20070926_1013-1030
            tint=[toepoch([2007 09 26 10 13 00]) toepoch([2007 09 26 10 30 00])];
        case '3c'
            load mBS_20070926_1045-1055
            tint=[toepoch([2007 09 26 10 45 00]) toepoch([2007 09 26 10 55 00])];
    end
    c_eval('B?staff=c_coord_trans(''DSC'',''gsm'',dBS?,''cl_id'',?);',3:4);
end

c_eval('gsmB?=c_fgm_staff_combine(gsmB?fgm(:,1:4),B?staff);',3:4)
c_eval('gseB?=c_coord_trans(''GSM'',''GSE'',gsmB?,''cl_id'',?);',3:4);
c_eval('diB?=c_coord_trans(''GSM'',''DSI'',gsmB?,''cl_id'',?);',3:4);
c_eval('gsmB?=irf_abs(gsmB?);',3:4);
%% Electric field
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',3:4);
c_eval('gseE?=c_coord_trans(''dsi'',''gse'',diE?,''CL_ID'',?);',3:4);
c_eval('gsmE?=c_coord_trans(''gse'',''gsm'',gseE?,''CL_ID'',?);',3:4);

%c_eval('[caaE?lowres,~,diE?lowres]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L3_E3D_INERT'');',3:4);

%% Magnetic field again
%c_eval('gsmB?fgm=irf_resamp(gsmB?fgm,diE?);',3:4)
%% Spacecraft potential
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',3:4);
%% Electron densities
c_eval('scpNe?=c_efw_scp2ne(P?);',4);
c_eval('scpNe?=irf_resamp(scpNe?,diE?);',4);
c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
%c_eval('peaNe?hr=irf_resamp(peaNe?,diE?);',3:4);
%% Ion densities
c_eval('hiaNi?=c_caa_var_get(''density__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',3);
c_eval('codNi?=c_caa_var_get(''density__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''mat'');',3:4);
%% Ion and electron (ExB) velocities GSE
c_eval('peaVe?=c_caa_var_get(''Data_Velocity_GSE__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
c_eval('diExB?=c_caa_var_get(''v_drift_ISR2__C?_CP_EFW_L2_V3D_INERT'',''mat'');',3:4);
c_eval('gseExB?=c_coord_trans(''dsi'',''gse'',diExB?,''CL_ID'',?);',3:4);
c_eval('gsmExB?=c_coord_trans(''dsi'',''gsm'',diExB?,''CL_ID'',?);',3:4);
c_eval('hiaVi?=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',3);
c_eval('hiaVi?=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',3);
c_eval('codVi?=c_caa_var_get(''velocity__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''mat'');',3:4 );
c_eval('hiaVi?=irf_resamp(hiaVi?,diE3);',3);
c_eval('codVi?=irf_resamp(codVi?,diE3);',3:4);
%% Spacecraft positions
c_eval('gsePos?=c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',3:4);         
c_eval('diPos?=c_coord_trans(''gse'',''isr2'',gsePos?,''cl_id'',?);',3:4);         
%% Electron temperature
% K 10^6
c_eval('parTe?=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
c_eval('perTe?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_'',''mat'');',3:4);
c_eval('Te?spinres=[parTe?(:,1) (parTe?(:,2)+2*perTe?(:,2))/3];',3:4);
c_eval('parTe?=irf_resamp(parTe?,diE3);',3:4)
c_eval('perTe?=irf_resamp(perTe?,diE3);',3:4)
% Adding to one temperature only, using C3
c_eval('Te?=[parTe?(:,1) (parTe?(:,2)+2*perTe?(:,2))/3];',3:4);


% eV
K2eV = 8.621738*1e-5;
MK2eV = 8.621738*1e-5*1e6;
c_eval('eVTe?=[Te?(:,1) Te?(:,2)*MK2eV];',3:4)
%c_eval('eVTe?=[Te?(:,1) Te?(:,2)*8.61734*10];',3:4)
c_eval('eVTe?=irf_resamp(eVTe?,diE3);',3:4)

%Te4=[parTe4(:,1) (parTe4(:,2)+2*perTe4(:,2))/3];
%eVTe4=[Te4(:,1) Te4(:,2)*8.61734*10];
%eVTe4=irf_resamp(eVTe4,diE3);
%Tetot4=[Tepar4(:,1) (Tepar4(:,2)+2*Teper4(:,2))/3];
%Teper=cn_average(Teper3,Teper4);
%% Ion temperatures
% K 10^6
c_eval('hiaTi?=c_caa_var_get(''temperature__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',3);
c_eval('codTi?=c_caa_var_get(''T__C?_CP_CIS_CODIF_HS_H1_MOMENTS'',''mat'');',3:4)
% eV
c_eval('hiaTi?eV=[hiaTi?(:,1) hiaTi?(:,2)*8.61734e-2*1000];',3)
c_eval('codTi?eV=[codTi?(:,1) codTi?(:,2)*8.61734e-2*1000];',3:4)

c_eval('codTi?eVhr=irf_resamp(codTi?eV,diE?);',3:4)
c_eval('hiaTi?eVhr=irf_resamp(hiaTi?eV,diE?);',3)
%eVTi3=irf_resamp(eVTi3,diE3);
eVTi3=codTi4eVhr;
end
% 0831: CODIF finns p? b?da, HIA p? C3
% 0902: samma
toc