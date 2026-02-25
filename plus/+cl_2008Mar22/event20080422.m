cd /Users/Cecilia/Data/Cluster/20080422/
tint = toepoch([2008 04 22 17 00 00;2008 04 22 18 30 00])';
savePath = '/Users/Cecilia/Research/LH2/20080422/';
sclist=1:4;
doLoad = 1;
if doLoad

% Electric field
c_eval('diE? = c_caa_var_get(''E_Vec_xyz_ISR2__C?_CP_EFW_L2_E3D_INERT'',''mat'');',sclist);
c_eval('diE?(:,1) = diE?(:,1)-1/450;') % Adjust for time error in CAA data.

% Magnetic field
c_eval('diB?fgm = c_caa_var_get(''B_vec_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);
c_eval('gseB?sta = c_caa_var_get(''B_vec_xyz_Instrument__C?_CP_STA_CWF_GSE'',''mat'');',sclist);
c_eval('diB?sta=c_coord_trans(''GSE'',''ISR2'',gseB?sta,''cl_id'',?);',sclist);
c_eval('gsmB?sta=c_coord_trans(''GSE'',''GSM'',gseB?sta,''cl_id'',?);',sclist);
c_eval('diB?=c_fgm_staff_combine(diB?fgm(:,1:4),diB?sta,''cl_id'',?);',sclist)
c_eval('absB? = irf_abs(diB?fgm); absB? = absB?(:,[1 5]);',sclist)    

% Electron moments
c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('parTe?MK=c_caa_var_get(''Data_Temperature_ComponentParallelToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('parTe?=[parTe?MK(:,1) irf_convert(parTe?MK(:,2),''MK2eV'')];',sclist);
%c_eval('perTe?=c_caa_var_get(''Data_Temperature_ComponentPerpendicularToMagField__C?_CP_PEA_MOMENTS'',''mat'');',sclist);
c_eval('peaNe?hf = irf_resamp(peaNe?,diE?);',sclist);
c_eval('e_deflux?=c_caa_var_get(''Data__C?_CP_PEA_PITCH_SPIN_DEFlux'',''mat'');',sclist);

% Ion moments
units=irf_units;
c_eval('gseVi?=c_caa_var_get(''velocity_gse__C?_CP_CIS_HIA_ONBOARD_MOMENTS'',''mat'');',[1 3]);
c_eval('gsmVi? = c_coord_trans(''GSE'',''GSM'',gseVi?,''cl_id'',?);',[1 3]);
c_eval('Ti?MK=c_caa_var_get(''temperature__C?_CP_CIS-HIA_ONBOARD_MOMENTS'',''mat'');',[1 3]);
c_eval('Ti?=[Ti?MK(:,1) irf_convert(Ti?MK(:,2),''MK2eV'')];',[1 3]);
c_eval('vti? = [Ti?(:,1) sqrt(2*units.e*Ti?(:,2)/units.mp)*1e-3];',[1 3]) % km/s

%% ExB drift
c_eval('gseExB? = c_caa_var_get(''v_drift_GSE__C?_CP_EFW_L2_V3D_GSE'',''mat'');',sclist);
c_eval('gsmExB? = c_coord_trans(''GSE'',''GSM'',gseExB?,''cl_id'',?);',[sclist]);

% Position
c_eval('diR? = c_caa_var_get(''sc_pos_xyz_isr2__C?_CP_FGM_FULL_ISR2'',''mat'');',sclist);
c_eval('gsmR? = c_coord_trans(''ISR2'',''GSM'',diR?,''cl_id'',?);',sclist);

% Make 3D electric field data
c_eval('[diE?,ang?]=irf_edb(diE?,diB?fgm);',sclist);

% Transform fields to GSM coordinates
c_eval('gsmE? = c_coord_trans(''ISR2'',''GSM'',diE?,''cl_id'',?);',sclist);
c_eval('gsmB? = c_coord_trans(''ISR2'',''GSM'',diB?,''cl_id'',?);',sclist);

% Spacecraft potential
c_eval('P?=c_caa_var_get(''Spacecraft_potential__C?_CP_EFW_L2_P'',''mat'');',sclist);
% Electron densities
c_eval('scpNe?=c_efw_scp2ne(P?);',sclist);
%c_eval('scpNe?=irf_resamp(scpNe?,diE?);',3:4);
%c_eval('peaNe?=c_caa_var_get(''Data_Density__C?_CP_PEA_MOMENTS'',''mat'');',3:4);
%c_eval('peaNe?hr=irf_resamp(peaNe?,diE?);',3:4);
c_eval('whiNe?=c_caa_var_get(''Electron_Density__C?_CP_WHI_ELECTRON_DENSITY'',''mat'');',sclist);

% Whisper
c_eval('whiPow?=c_caa_var_get(''Electric_Wave_Form_Power_Density__C?_CP_WHI_WAVE_FORM_ENERGY'',''mat'');',sclist);

% Derived quantities
units = irf_units;

c_eval('partPr?=irf_multiply(1,peaNe?,1,parTe?MK,1); partPr? = [partPr?(:,1) partPr?(:,2)*units.kB*1e6*1e6];',sclist)
c_eval('magPr? = [absB?(:,1) absB?(:,2).^2*1e-18/2/units.mu0];',sclist)
c_eval('betaPr? = irf_multiply(1,partPr?,1,magPr?,-1);',sclist)

c_eval('Vtp? = irf_plasma_calc(gsmB?,peaNe?,0,parTe?,Ti?,''Vtp'');',[1 3])
c_eval('Va? = irf_plasma_calc(gsmB?,peaNe?,0,parTe?,Ti?,''Va'');',[1 3])
c_eval('beta?=irf_multiply(1,Vtp?,2,Va?,-2);',[1 3])

% Electron thermal velocity
c_eval('vte? = [parTe?(:,1) sqrt(2*units.e*parTe?(:,2)/units.me)*1e-3];',sclist) % km/s

% Electron cyclotron frequency
c_eval('fce? = [absB?(:,1) absB?(:,2)*1e-9*units.e/units.me/2/pi];',sclist)
    
% Lower hybrid frequency
c_eval('flh? = [absB?(:,1) absB?(:,2)*1e-9*units.e/sqrt(units.me*units.mp)/2/pi];',sclist)
    
% Electron gyroradius
c_eval('re? = irf_multiply(1/2/pi/sqrt(2),vte?,1,fce?,-1);',sclist) % km

% Integrate E
c_eval('intdiE? = irf_integrate(diE?);',sclist) % mVs/m

% Make power spectrograms
nfft = 1024;
sfreq = 450;
overlap = 0.50;
c_eval('fftE? = irf_powerfft(diE?,nfft,sfreq,[overlap]);?',sclist)
c_eval('fftB? = irf_powerfft(diB?,nfft,sfreq,[overlap]);?',sclist)
%c_eval('absE?=irf_abs(diE?); wavE? = irf_wavelet(absE?(:,[1 5]));?',sclist)
%c_eval('absB?=irf_abs(diB?); wavE? = irf_wavelet(absB?(:,[1 5]));?',sclist)
else 
    load data2
end


%% Overview plot for single SC 
sc=1;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
h = irf_plot(8);
isub = 1;

if 1 % Electron differential energy flux
    hca = h(isub); isub=isub+1;
    c_eval('deflux = e_deflux?;',sc)
    
    specrec.t = deflux.t;
    specrec.f = deflux.dep_x{2}.data(1,:)'*1e3; % eV
    specrec.p = {squeeze(nansum(deflux.data,2))};
    irf_spectrogram(hca,specrec.t,specrec.p{1},specrec.f)
    hca.YScale = 'log';
    hca.YTick = [0.01 0.1 1 10 100 1000 10000];
    hca.YLabel.String = 'E_e [eV]';
    hc = colorbar('peer',hca);
    hc.Label.String = 'keV/cm^2 s sr keV';
    
    %irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',sc),'sum_dim1','colorbarlabel','log_{10} dEF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    %hca.YLabel.String = 'E [nT]';
end
if 1 % n (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('peaNe?',sc));
    hca.YLabel.String = 'n_e [cc]';    
end
if 1 % Te (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('parTe?',sc));
    hca.YLabel.String = 'T_{e,||} [eV]';    
end

if 0 % B ISR2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
try
if 1 % V GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmVi?',sc));
    hca.YLabel.String = 'Vi [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
end
if 0 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftE?.t,fftE?.p{1},fftE?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(nT)^2/Hz';
    hca.NextPlot = 'add';
    c_eval('irf_plot(hca,flh?)',sc);
    c_eval('irf_plot(hca,fce?)',sc);
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';
end
if 1 % B power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftB?.t,fftB?.p{1},fftB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(mV/m)^2/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end
if 0
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftEB?.t,fftEB?.p{1},fftEB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = 'mV/m/nT/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end


title(h(1),irf_ssub('C?',sc))
irf_zoom(h,'x',tint)
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
%eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])
%% Overview plot for 2-4 SC, fields (E,V,Vi)
sc=1:4;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
h = irf_plot(10);
isub = 1;

if 1 % B2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB1);
    hca.YLabel.String = 'B1 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB2);
    hca.YLabel.String = 'B2 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB3);
    hca.YLabel.String = 'B3 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % B4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmB4);
    hca.YLabel.String = 'B4 [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E1
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE1);
    hca.YLabel.String = 'E1 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE2);
    hca.YLabel.String = 'E2 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE3);
    hca.YLabel.String = 'E3 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E4
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmE4);
    hca.YLabel.String = 'E4 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % Vi1
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmVi1);
    hca.YLabel.String = 'Vi1 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % Vi3
    hca = h(isub); isub=isub+1;
    irf_plot(hca,gsmVi3);
    hca.YLabel.String = 'Vi3 [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end

title(h(1),'')
irf_zoom(h,'x',tint)
%irf_plot_axis_align
%set(gcf,'paperpositionmode','auto');
%strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
%eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])
%% Overview plot for single SC 
sc=3;
c_eval('tint=[diB?fgm(1,1) diB?fgm(end,1)];',sc)
h = irf_plot(8);
isub = 1;

if 1; try % Electron differential energy flux
    hca = irf_panel('Electron DEFlux');
    c_eval('deflux = e_deflux?;',sc)    
    specrec.t = deflux.t;
    specrec.f = deflux.dep_x{2}.data(1,:)'*1e3; % eV
    specrec.p = {squeeze(nansum(deflux.data,2))};
    irf_spectrogram(hca,specrec.t,specrec.p{1},specrec.f)
    hca.YScale = 'log';
    hca.YTick = [0.01 0.1 1 10 100 1000 10000];
    hca.YLabel.String = 'E_e [eV]';
    hc = colorbar('peer',hca);
    hc.Label.String = 'keV/cm^2 s sr keV';
    
    %irf_plot(hca,irf_ssub('Data__C?_CP_PEA_PITCH_SPIN_DEFlux',sc),'sum_dim1','colorbarlabel','log_{10} dEF\newline #/cm^2 s sr keV','fitcolorbarlabel');
    %hca.YLabel.String = 'E [nT]';
end; end

if 1 % n (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('peaNe?',sc));
    hca.YLabel.String = 'n_e [cc]';    
end
if 0 % Te (PEACE)
    hca = h(isub); isub=isub+1;    
    irf_plot(hca,irf_ssub('parTe?',sc));
    hca.YLabel.String = 'T_{e,||} [eV]';    
end

if 0 % B ISR2
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % B GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmB?',sc));
    hca.YLabel.String = 'B [nT]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
try
if 1 % V GSM
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmVi?',sc));
    hca.YLabel.String = 'Vi [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
end
if 0 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('diE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'ISR2'},[0.02, 0.95]);
end
if 1 % E
    hca = h(isub); isub=isub+1;
    irf_plot(hca,irf_ssub('gsmE?',sc));
    hca.YLabel.String = 'E [mV/m]';
    irf_legend(hca,{'x','y','z'},[0.98, 0.95]);
    irf_legend(hca,{'GSM'},[0.02, 0.95]);
end
if 1 % E power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftE?.t,fftE?.p{1},fftE?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(nT)^2/Hz';
    hca.NextPlot = 'add';
    c_eval('irf_plot(hca,flh?)',sc);
    c_eval('irf_plot(hca,fce?)',sc);
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';
end
if 1 % B power spectrogram
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftB?.t,fftB?.p{1},fftB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = '(mV/m)^2/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end
if 0
    hca = h(isub); isub=isub+1;
    c_eval('irf_spectrogram(hca,fftEB?.t,fftEB?.p{1},fftEB?.f)',sc);
    hca.YScale = 'log';
    hc = colorbar('peer',hca);
    hc.Label.String = 'mV/m/nT/Hz';
    hca.NextPlot = 'add';
    irf_plot(hca,irf_ssub('flh?',sc));
    irf_plot(hca,irf_ssub('fce?',sc));
    hca.YTick = [1 10 100 1000 10000];
    hca.YLabel.String = 'f [Hz]';    
end


title(h(1),irf_ssub('C?',sc))
irf_zoom(h,'x',tint)
irf_plot_axis_align
set(gcf,'paperpositionmode','auto');
strPrint = [datestr(irf_time(tint(1),'epoch>datenum'),'yyyy-mm-dd') '_Overview'];
%eval(['print -depsc /Users/Cecilia/Research/LH2/BurstModeOverview/' strPrint '.eps'])
%%
% Add B
h = irf_plot(2,'newfigure');
hca = h(1);
irf_plot(hca,diE1)

hca = h(2); hold(hca,'on')
for kk=1:singleTT.numel
    irf_plot(hca,singleTT.Description{kk}.phiB); 
end
hold(hca,'off')

%% Rotate everything into the reference system of the ion outflow jet
%irf_plot({gsmB3,gsmVi3,gsmE3})
%tint = get_time(2); tint = torow(tint);
Viloc = irf_tlim(gsmVi3,tint);
Viloc = mean(Viloc(:,2:4),1);

z=irf_norm(Viloc); % B/z direction, tries*3
y=irf_norm(irf_cross(irf_cross(z,[-0.93 0.32 -0.19]),z)); % perp1
x=irf_norm(irf_cross(y,z)); % perp2

sysB3 = [gsmB3(:,1) gsmB3(:,2:4)*[x' y' z']];
sysVi3 = [gsmVi3(:,1) gsmVi3(:,2:4)*[x' y' z']];
sysE3 = [gsmE3(:,1) gsmE3(:,2:4)*[x' y' z']];
   
h=irf_plot(3,'newfigure');
irf_plot({sysB3,sysVi3,sysE3})

