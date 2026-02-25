function [gsmB3 gsmB4]=av_combine(event,off)
%% Andris script
sclist=3:4;
if strcmp(off,'offs')
    if 1, % read FGM data form all sc
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
else
    c_eval('[caaB?,~,B?]=c_caa_var_get(''B_vec_xyz_gse__C?_CP_FGM_FULL'');');
    c_eval('gsmB?=c_coord_trans(''GSE'',''GSM'',B?,''cl_id'',?);',sclist);
end
if 1, % read STAFF data form all sc
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
    c_eval('Bsc?=c_coord_trans(''DSC'',''DSI'',dBS?,''cl_id'',?);',sclist);
end
if 1, % lowpass filter B field data
    Fs = 66.7;  % Sampling Frequency
   
    Fpass = 2;           % Passband Frequency
    Fstop = 4;           % Stopband Frequency
    Apass = 1;           % Passband Ripple (dB)
    Astop = 20;          % Stopband Attenuation (dB)
    match = 'stopband';  % Band to match exactly
   
    % Construct an FDESIGN object and call its BUTTER method.
    ff  = fdesign.lowpass(Fpass, Fstop, Apass, Astop, Fs);
    d1 = design(ff, 'butter', 'MatchExactly', match);    
    [B,A]= sos2tf(d1.sosMatrix,d1.ScaleValues);
    c_eval('gsmBfilt?=gsmB?(:,1:4);',sclist);
    c_eval('gsmBfilt?(:,2)=filtfilt(B,A,gsmB?(:,2));',sclist);
    c_eval('gsmBfilt?(:,3)=filtfilt(B,A,gsmB?(:,3));',sclist);
    c_eval('gsmBfilt?(:,4)=filtfilt(B,A,gsmB?(:,4));',sclist);
    c_eval('gsmBfilt?=irf_abs(gsmBfilt?);',sclist);
end
if 1, % highpass filter B field data
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
    [B,A]= sos2tf(d1.sosMatrix,d1.ScaleValues);
    c_eval('Bsc_filt?=Bsc?(:,1:4);',sclist);
    c_eval('Bsc_filt?(:,2)=filtfilt(B,A,Bsc?(:,2));',sclist);
    c_eval('Bsc_filt?(:,3)=filtfilt(B,A,Bsc?(:,3));',sclist);
    c_eval('Bsc_filt?(:,4)=filtfilt(B,A,Bsc?(:,4));',sclist);
    c_eval('gsmBsc_filt?=irf_gse2gsm(Bsc_filt?);',sclist);
end
if 1, % construct combined magnetic field
    c_eval('gsmBfull?=irf_add(1,gsmBsc_filt?,1,gsmBfilt?(:,1:4));',sclist);
    c_eval('gsmBfull?=irf_abs(gsmBfull?);',sclist);
    c_eval('gsmB?=gsmBfull?;',sclist);
end