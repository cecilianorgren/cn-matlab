sc = 3;
detector = 'HEEA';
c_eval(['psd?Data=c_caa_var_get(''Data_' detector '__C?_CP_PEA_PITCH_FULL_PSD'',''mat'');'],sc);
c_eval(['psd?Energy=c_caa_var_get(''Sweep_Energy_' detector '__C3_CP_PEA_PITCH_FULL_PSD'',''mat'');'],sc);
c_eval(['psd?Pitchangle=c_caa_var_get(''Sweep_PitchAngle_' detector '__C3_CP_PEA_PITCH_FULL_PSD'',''mat'');'],sc);
c_eval(['psd?Azimuth=c_caa_var_get(''Sweep_Azimuth_' detector '__C3_CP_PEA_PITCH_FULL_PSD'',''mat'');'],sc);
