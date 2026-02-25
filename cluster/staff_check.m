%% download staff data
tst=irf_time([2009 09 16 07 25 00]);dt=1.5*60;
c_get_batch(tst,dt,'vars','bsc','noproc')
%% load fgm and staff data
sc=1:4;
load mBSCR
c_eval('dB?staff=wBSC?;',sc);
c_eval('diB?staff=c_coord_trans(''SR2'',''ISR2'',dB?staff,''cl_id'',?);',sc);
c_eval('diB?fgm=c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sc);
%% filter data to see if staff is loaded correctly
c_eval('diB?s=irf_filt(diB?staff,1,10);',sc);
c_eval('diB?f=irf_filt(diB?fgm,1,10);',sc);
%% plot data to see if staff is loaded correctly
irf_plot({diB1s,diB1f},'comp');
%% combine staff and fgm
c_eval('diB?=c_fgm_staff_combine(diB?fgm,diB?staff);',sc);
%% load E to see if they match
c_eval('diE?=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sc);
c_eval('diE?ib=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_EB'',''mat'');',sc);
%%
% look if staff is read in correctly
% filter both in range 1-10 Hz, they should in most cases be identical

% low pass filter for FGM and STAFF
    Fs_fgm = 66.7;  % Sampling Frequency
    Fs_staff = 450;  % Sampling Frequency
    
    Fpass = 0.5;          % Passband Frequency
    Fstop = 10;          % Stopband Frequency
    Apass = 1;           % Passband Ripple (dB)
    Astop = 20;          % Stopband Attenuation (dB) (80 for staff)
    match = 'stopband';  % Band to match exactly
    
    % Construct an FDESIGN object and call its BUTTER method.
    ff_fgm  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs_fgm);
    d1_fgm = design(ff_fgm, 'butter', 'MatchExactly', match);
    
    ff_staff = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs_staff);
    d1_staff = design(ff_staff, 'butter', 'MatchExactly', match);
    
% higpass filter for STAFF
%     Fs = 450;  % Sampling Frequency
%     
%     Fstop = 1;           % Stopband Frequency
%     Fpass = 3;           % Passband Frequency
%     Astop = 80;          % Stopband Attenuation (dB)
%     Apass = 1;           % Passband Ripple (dB)
%     match = 'stopband';  % Band to match exactly
%     
%     % Construct an FDESIGN object and call its BUTTER method.
%     ff  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
%     d1_staff = design(ff, 'butter', 'MatchExactly', match);
    
b_fgm=diB3fgm;
b_fgm_filt=diB3fgm;
b_staff=diB3staff;
b_staff_filt=diB3staff;
db_staff=dB3staff;
db_staff_filt=dB3staff;

for j=2:4,
    b_fgm_filt(:,j)=filtfilt(d1_fgm.sosMatrix,d1_fgm.ScaleValues,b_fgm(:,j));
    b_staff_filt(:,j)=filtfilt(d1_staff.sosMatrix,d1_staff.ScaleValues,b_staff(:,j));
    db_staff_filt(:,j)=filtfilt(d1_staff.sosMatrix,d1_staff.ScaleValues,db_staff(:,j));
    c_eval('b_fgm_filt_2=irf_filt(b_fgm,Fpass,Fstop,Fs_fgm,5);',1:4)
    c_eval('b_staff_filt_2=irf_filt(b_staff,Fpass,Fstop,Fs_staff,5);',1:4)
    c_eval('db_staff_filt_2=irf_filt(db_staff,Fpass,Fstop,Fs_staff,5);',1:4)
    %c_eval('b_fgm_filt_3=butter(b_fgm,[Fpass Fstop]/Fs_fgm,''stop'');',1:4)
    %c_eval('b_staff_filt_3=butter(b_staff,[Fpass Fstop]/Fs_staff,''stop'');',1:4)
    %c_eval('db_staff_filt_3=butter(db_staff,[Fpass Fstop]/Fs_staff,''stop'');',1:4)
end

h=irf_plot({b_staff_filt_2,b_fgm_filt},'comp');
legend('irf filt','butter','andris')
%h=irf_plot(3);
%irf_plot(h(1),b_staff_filt(:,[1 4]))
%irf_plot(h(2),db_staff_filt(:,[1 4]))
%irf_plot(h(3),b_fgm_filt(:,[1 4]))
%legend('fgm','staff','staff2')
%irf_zoom(h,'x',[toepoch([2009 09 16 07 25 35]) toepoch([2009 09 16 07 25 45])])