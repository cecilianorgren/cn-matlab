%% Compare C3/C4 at better than spin resolution
if 1, % initialize figure
    fn=figure(62);
    h=irf_plot(6);
    set(fn,'defaultLineLineWidth',1);
end
sclist=3:4;
if 0, % read FGM data form all sc
    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');');
    %    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_5VPS'');',sclist);
    % remove offset
    c_eval(['offs?=mean(B?(:,4))-mean(B' num2str(sclist(1)) '(:,4));'],sclist(2:end));
    c_eval('dB?=c_coord_trans(''GSE'',''DSI'',B?,''cl_id'',?);',sclist(2:end));
    c_eval('dB?(:,4)=dB?(:,4)-offs?;',sclist(2:end));
    c_eval('B?=c_coord_trans(''DSI'',''GSE'',dB?,''cl_id'',?);',sclist(2:end));
    c_eval('B?=irf_abs(B?);',sclist);
    c_eval('gsmB?=irf_gse2gsm(B?);',sclist);
end
if 0, % read STAFF data form all sc
    load mBS
    c_eval('Bsc?=c_coord_trans(''DSC'',''DSI'',dBS?,''cl_id'',?);',sclist);
end
if 0, % lowpass filter B field data
   
    % All frequency values are in Hz.
    Fs = 66.7;  % Sampling Frequency
   
    N     = 100;  % Order
    Fpass = 2;    % Passband Frequency
    Fstop = 4;    % Stopband Frequency
    Wpass = 1;    % Passband Weight
    Wstop = 1;    % Stopband Weight
    dens  = 20;   % Density Factor
   
    % Calculate the coefficients using the FIRPM function.
    b  = firpm(N, [0 Fpass Fstop Fs/2]/(Fs/2), [1 1 0 0], [Wpass Wstop], ...
        {dens});
    d1 = dfilt.dffir(b);
   
   
    c_eval('gsmBfilt?=gsmB?(:,1:4);',sclist);
    c_eval('gsmBfilt?(:,2)=filtfilt(d1.Numerator,1,gsmB?(:,2));',sclist);
    c_eval('gsmBfilt?(:,3)=filtfilt(d1.Numerator,1,gsmB?(:,3));',sclist);
    c_eval('gsmBfilt?(:,4)=filtfilt(d1.Numerator,1,gsmB?(:,4));',sclist);
    c_eval('gsmBfilt?=irf_abs(gsmBfilt?);',sclist);
end
if 0, % highpass filter B field data
    % All frequency values are in Hz.
    Fs = 450;  % Sampling Frequency
   
    N     = 1000;  % Order
    Fstop = 1.5;   % Stopband Frequency
    Fpass = 2.5;   % Passband Frequency
    Wstop = 1;     % Stopband Weight
    Wpass = 0.1;   % Passband Weight
    dens  = 20;    % Density Factor
   
    % Calculate the coefficients using the FIRPM function.
    b  = firpm(N, [0 Fstop Fpass Fs/2]/(Fs/2), [0 0 1 1], [Wstop Wpass], ...
        {dens});
    d1 = dfilt.dffir(b);
    c_eval('Bsc_filt?=Bsc?(:,1:4);',sclist);
    c_eval('Bsc_filt?(:,2)=filtfilt(d1.Numerator,1,Bsc?(:,2));',sclist);
    c_eval('Bsc_filt?(:,3)=filtfilt(d1.Numerator,1,Bsc?(:,3));',sclist);
    c_eval('Bsc_filt?(:,4)=filtfilt(d1.Numerator,1,Bsc?(:,4));',sclist);
    c_eval('gsmBsc_filt?=irf_gse2gsm(Bsc_filt?);',sclist);
end
if 0, % (Newer matlab version) lowpass filter B field data
    Fs = 66.7;  % Sampling Frequency
   
    Fpass = 2;           % Passband Frequency
    Fstop = 4;           % Stopband Frequency
    Apass = 1;           % Passband Ripple (dB)
    Astop = 20;          % Stopband Attenuation (dB)
    match = 'stopband';  % Band to match exactly
   
    % Construct an FDESIGN object and call its BUTTER method.
    ff  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
    d1 = design(ff, 'butter', 'MatchExactly', match);
    c_eval('gsmBfilt?=gsmB?(:,1:4);',sclist);
    c_eval('gsmBfilt?(:,2)=filtfilt(d1.sosMatrix,d1.ScaleValues,gsmB?(:,2));',sclist);
    c_eval('gsmBfilt?(:,3)=filtfilt(d1.sosMatrix,d1.ScaleValues,gsmB?(:,3));',sclist);
    c_eval('gsmBfilt?(:,4)=filtfilt(d1.sosMatrix,d1.ScaleValues,gsmB?(:,4));',sclist);
    c_eval('gsmBfilt?=irf_abs(gsmBfilt?);',sclist);
end
if 0, % (Newer matlab version) highpass filter B field data
    % All frequency values are in Hz.
    Fs = 450;  % Sampling Frequency
   
    Fstop = 1;           % Stopband Frequency
    Fpass = 3;           % Passband Frequency
    Astop = 80;          % Stopband Attenuation (dB)
    Apass = 1;           % Passband Ripple (dB)
    match = 'stopband';  % Band to match exactly
   
    % Construct an FDESIGN object and call its BUTTER method.
    ff  = fdesign.highpass(Fstop, Fpass, Astop, Apass, Fs);
    d1 = design(ff, 'butter', 'MatchExactly', match);
    c_eval('Bsc_filt?=Bsc?(:,1:4);',sclist);
    c_eval('Bsc_filt?(:,2)=filtfilt(d1.sosMatrix,d1.ScaleValues,Bsc?(:,2));',sclist);
    c_eval('Bsc_filt?(:,3)=filtfilt(d1.sosMatrix,d1.ScaleValues,Bsc?(:,3));',sclist);
    c_eval('Bsc_filt?(:,4)=filtfilt(d1.sosMatrix,d1.ScaleValues,Bsc?(:,4));',sclist);
    c_eval('gsmBsc_filt?=irf_gse2gsm(Bsc_filt?);',sclist);
end
if 1, % construct combined magnetic field
    c_eval('gsmBfull?=irf_add(1,gsmBsc_filt?,1,gsmBfilt?(:,1:4));',sclist);
    c_eval('gsmBfull?=irf_abs(gsmBfull?);',sclist);
    c_eval('gsmB?=gsmBfull?;',sclist);
end
if 1, % read EFW data
    %     c_eval('[caaE?,~,diE?,diE?units]=c_caa_var_get(''E_Vec_xy_ISR2__C?_CP_EFW_L2_E'');',sclist);
    %     c_eval('[caaE?,~,diE?,diE?units]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L3_E3D_INERT'');',sclist);
    c_eval('[caaE?,~,diE?,diE?units]=c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'');',sclist);
end
if 1,   % PANEL: FGM B single s/c
    hca=irf_panel('B');
    irf_plot(hca,gsmB3);
    ylabel(hca,'B [nT] GSM');
    irf_zoom(hca,'y',[0 30]);
    irf_legend(hca,{'B_X','B_Y','B_Z','B'},[0.02 0.3])
    irf_legend(hca,'C3',[0.98 0.98],'color','k')
end
if 1,   % PANEL: FGM Bx two s/c
    hca=irf_panel('B 2sc');
    irf_plot(hca,gsmB3(:,[1 2]),'g');
    ylabel(hca,'B_X [nT] GSM');
    irf_zoom(hca,'y');
    hold(hca,'on');
    irf_plot(hca,gsmB4(:,[1 2]),'b');
    irf_legend(hca,{'','','C3','C4'},[0.98, 0.98],'color','cluster')
    hold(hca,'off');
end
if 1,   % PANEL: Delta Bx two s/c
    hca=irf_panel('dB 2sc');
    c_eval('irf_plot(hca,irf_add(1,gsmB3(:,[1 :4]),-1,gsmB4(:,[1 :4])));');
    ylabel(hca,'\Delta B [nT] GSM');
    irf_zoom(hca,'y');
end
if 1,   % PANEL: Staff By two s/c
    hca=irf_panel('STAFF By');
    irf_plot(hca,gsmBsc_filt3(:,[1 3]),'g');
    ylabel(hca,'B_Y [nT] GSM');
    irf_zoom(hca,'y');
    hold(hca,'on');
    irf_plot(hca,gsmBsc_filt4(:,[1 3]),'b');
    irf_legend(hca,{'','','C3','C4'},[0.98, 0.98],'color','cluster')
    hold(hca,'off');
end
if 1,   % PANEL: Staff Bz two s/c
    hca=irf_panel('STAFF Bz');
    irf_plot(hca,gsmBsc_filt3(:,[1 4]),'g');
    ylabel(hca,'B_Z [nT] GSM');
    irf_zoom(hca,'y');
    hold(hca,'on');
    irf_plot(hca,gsmBsc_filt4(:,[1 4]),'b');
    irf_legend(hca,{'','','C3','C4'},[0.98, 0.98],'color','cluster')
    hold(hca,'off');
end
if 1,   % PANEL: EFW Ex two s/c
    hca=irf_panel('Ex');
    irf_plot(hca,diE3(:,[1 2]),'g');
    ylabel(hca,'E_X [mV/m] ISR2');
    irf_zoom(hca,'y');
    hold(hca,'on');
    irf_plot(hca,diE4(:,[1 2]),'b');
    irf_legend(hca,{'','','C3','C4'},[0.98, 0.98],'color','cluster')
    hold(hca,'off');
end
if 0,   % PANEL: EFW E field in ISR2 reference frame single s/c
    hca=irf_panel('E ISR2');
    c_eval('irf_plot(hca,diE?)',ic);
    ylabel(hca,'E [mV/m] ISR2');
    irf_zoom(hca,'ylim','smart');
    irf_legend(hca,{'E_X','E_Y'},[0.02 0.49])
    irf_legend(hca,{['C' num2str(ic)]},[0.98 0.98],'color','k')
end
irf_plot_axis_align
irf_zoom(h,'x',tint);
irf_legend(0,'Figure C3/C4, electron scale E field',[0.02 0.02],'fontsize',8,'color',[0.5 0.5 0.5]);    irf_plot_axis_align
irf_pl_number_subplots(h,[0.02,0.97],'fontsize',14);
irf_timeaxis(h);