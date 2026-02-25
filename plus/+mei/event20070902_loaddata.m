% Making figures for Maria Elena Innocenti's paper
% This script loads the data, either from CAS .cfd files or
% data20070902.mat
cd /Users/Cecilia/Data/Cluster/20070902/

doLoad = 0; % 1: load from CAA/, 2: load from data20070902.mat

if doLoad % Load from CAA/
sclist=3:4;
% Electric field
c_eval('diE? = c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
c_eval('diE?(:,1) = diE?(:,1)-1/450;',sclist) % Adjust for time error in CAA data.

%% Magnetic field
c_eval('gseB?fgm = c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);
c_eval('gseB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_GSE'',''mat'');',sclist);
c_eval('diB?sta=c_coord_trans(''GSE'',''ISR2'',gseB?sta,''cl_id'',?);',sclist);
c_eval('gsmB?sta=c_coord_trans(''GSE'',''GSM'',gseB?sta,''cl_id'',?);',sclist);
c_eval('gseB?=c_fgm_staff_combine(gseB?fgm(:,1:4),gseB?sta,''cl_id'',?);',sclist)
c_eval('absB? = irf_abs(gseB?fgm); absB? = absB?(:,[1 5]);',sclist)    

%% Electron moments
c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('parTe?MK=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('parTe?=[parTe?MK(:,1) irf_convert(parTe?MK(:,2),''MK2eV'')];',sclist);
%c_eval('perTe?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('peaNe?hf = irf_resamp(peaNe?,diE?);',sclist);
c_eval('e_deflux?=c_caa_var_get(''Data__C?_CP_PEA_PITCH_SPIN_DEFlux'',''mat'');',sclist);

%% Ion moments
units=irf_units;
sclist = 3;
%c_eval('dobj = caa_load(''C?_CP_CIS-HIA_ONBOARD_MOMENTS'');',sclist)
c_eval('gseVi?=c_caa_var_get(''velocity_gse__C?_CP_CIS-HIA_ONBOARD_MOMENTS'',''mat'');',sclist);
c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?,''cl_id'',?);',sclist);
c_eval('Ti?MK=c_caa_var_get(''temperature__C?_CP_CIS-HIA_ONBOARD_MOMENTS'',''mat'');',sclist);
c_eval('Ti?=[Ti?MK(:,1) irf_convert(Ti?MK(:,2),''MK2eV'')];',sclist);
c_eval('vti? = [Ti?(:,1) sqrt(2*units.e*Ti?(:,2)/units.mp)*1e-3];',sclist) % km/s

%% ExB drift
sclist = 3:4;
c_eval('gseExB? = c_caa_var_get(''v_drift_GSE__C?_CP_EFW_L2_V3D_GSE'',''mat'');',sclist);
c_eval('gsmExB? = c_coord_trans(''GSE'',''GSM'',gseExB?,''cl_id'',?);',[sclist]);

% Position
c_eval('gseR? = c_caa_var_get(''sc_pos_xyz_gse__C?_CP_FGM_FULL'',''mat'');',sclist);
c_eval('gsmR? = c_coord_trans(''gse'',''GSM'',gseR?,''cl_id'',?);',sclist);

% Make 3D electric field data
c_eval('diB? = c_coord_trans(''gse'',''isr2'',gseB?,''cl_id'',?);',sclist);
c_eval('[diE?,ang?]=irf_edb(diE?,diB?);',sclist);

% Transform fields to GSM coordinates
c_eval('gsmE? = c_coord_trans(''ISR2'',''GSM'',diE?,''cl_id'',?);',sclist);
c_eval('gsmB? = c_coord_trans(''ISR2'',''GSM'',diB?,''cl_id'',?);',sclist);

%% Spacecraft potential
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',sclist);
% Electron densities
c_eval('scpNe?=c_efw_scp2ne(P?);',sclist);
%c_eval('scpNe?=irf_resamp(scpNe?,diE?);',3:4);
%c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
%c_eval('peaNe?hr=irf_resamp(peaNe?,diE?);',3:4);

% Whisper
c_eval('whiPow?=c_caa_var_get(''Electric_Wave_Form_Power_Density__C?_CP_WHI_WAVE_FORM_ENERGY'',''mat'');',sclist);
c_eval('whiNe?=c_caa_var_get(''Electron_Density__C?_CP_WHI_ELECTRON_DENSITY'',''mat'');',sclist);

% Derived quantities
units = irf_units;

c_eval('ePr?=irf_multiply(1,peaNe?,1,parTe?MK,1); ePr? = [ePr?(:,1) ePr?(:,2)*units.kB*1e6*1e6];',sclist)
c_eval('iPr?=irf_multiply(1,peaNe?,1,Ti?MK,1); iPr? = [iPr?(:,1) iPr?(:,2)*units.kB*1e6*1e6];',3)
c_eval('magPr? = [absB?(:,1) absB?(:,2).^2*1e-18/2/units.mu0];',sclist)
c_eval('partPr? = irf_multiply(1,irf_add(1,ePr?,1,iPr?),1,magPr?,-1);',3)

%%
c_eval('Vtp? = irf_plasma_calc(gsmB?,peaNe?,0,parTe?,Ti?,''Vtp'');',3)
c_eval('Va? = irf_plasma_calc(gsmB?,peaNe?,0,parTe?,Ti?,''Va'');',3)
c_eval('beta?=irf_multiply(1,Vtp?,2,Va?,-2);',3)

% Electron thermal velocity
c_eval('vte? = [parTe?(:,1) sqrt(2*units.e*parTe?(:,2)/units.me)*1e-3];',sclist) % km/s

% Electron cyclotron frequency
c_eval('fce? = [absB?(:,1) absB?(:,2)*1e-9*units.e/units.me/2/pi];',sclist)

% Ion cyclotron frequency
c_eval('fci? = [absB?(:,1) absB?(:,2)*1e-9*units.e/units.mp/2/pi];',sclist)

% Lower hybrid frequency
c_eval('flh? = [absB?(:,1) absB?(:,2)*1e-9*units.e/sqrt(units.me*units.mp)/2/pi];',sclist)
    
% Electron gyroradius
c_eval('re? = irf_multiply(1/2/pi/sqrt(2),vte?,1,fce?,-1);',sclist) % km
% Ion gyroradius
c_eval('ri? = irf_multiply(1/2/pi/sqrt(2)*1e-3,Vtp?,1,fci?,-1);',3) % km

% Integrate E
c_eval('intdiE? = irf_integrate(diE?);',sclist) % mVs/m

if 0 % Make power spectrograms
nfft = 1024;
sfreq = 450;
overlap = 0.50;
c_eval('fftE? = irf_powerfft(diE?,nfft,sfreq,[overlap]);?',sclist)
c_eval('fftB? = irf_powerfft(diB?,nfft,sfreq,[overlap]);?',sclist)
%c_eval('absE?=irf_abs(diE?); wavE? = irf_wavelet(absE?(:,[1 5]));?',sclist)
%c_eval('absB?=irf_abs(diB?); wavE? = irf_wavelet(absB?(:,[1 5]));?',sclist)
end
else 
    load data20070902
    % Ion cyclotron frequency
    c_eval('fci? = [absB?(:,1) absB?(:,2)*1e-9*units.e/units.mp/2/pi];',sclist)
    % Ion gyroradius
    c_eval('ri? = irf_multiply(1/2/pi/sqrt(2)*1e-3,Vtp?,1,fci?,-1);',3) % km

end